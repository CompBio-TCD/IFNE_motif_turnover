import os,csv
import matplotlib.pyplot as plt
import numpy as np
import json
import pandas as pd 

from Bio import AlignIO 
from Bio import Phylo
from matplotlib.gridspec import GridSpec
from collections import defaultdict
from matplotlib.colors import ListedColormap, BoundaryNorm

def plot_tree_and_multi_heatmaps(
    tree_newick,
    ordered_species,
    M_by_tf,              # dict: tf -> matrix (n_species x aln_len)
    out_png,
    id_to_common=None,
    tf_order=None,
    max_tfs=2,            # e.g. 2 for ESR1 + PGR
    width_tree=1.0,
    width_tf=1.6,         # smaller = skinnier panels
    show_tree_tip_labels=False,
    global_vmax=None      # set if you want shared colour scale
):
    """
    Plots: tree (left) + multiple TF heatmaps (right, side-by-side).
    M_by_tf[tf] should already be ordered by ordered_species rows.
    """
    tfs = list(M_by_tf.keys()) if tf_order is None else [t for t in tf_order if t in M_by_tf]
    tfs = tfs[:max_tfs]

    PLACENTAL_EXCLUDE = {
        "monDom5",   # opossum
        "macEug2",   # wallaby
        "sarHar1",   # Tasmanian devil
        "ornAna1",   # platypus
    }

    # Common labels
    if id_to_common:
        ylabels = [id_to_common.get(sp, sp) for sp in ordered_species]
    else:
        ylabels = ordered_species

    n = len(ordered_species)
    fig = plt.figure(figsize=(2 + width_tf * len(tfs), 0.32 * n), constrained_layout=True)
    gs = GridSpec(1, 1 + len(tfs), width_ratios=[width_tree] + [width_tf] * len(tfs), figure=fig)

    # ---- Tree panel ----
    ax_tree = fig.add_subplot(gs[0, 0])
    tree = Phylo.read(tree_newick, "newick")

    # prune tree to ordered_species
    keep = set(ordered_species)
    for leaf in list(tree.get_terminals()):
        if leaf.name not in keep:
            tree.prune(leaf)

    # optionally relabel leaves to common names (tree tip labels)
    if id_to_common:
        for leaf in tree.get_terminals():
            if leaf.name in id_to_common:
                leaf.name = id_to_common[leaf.name]

    if show_tree_tip_labels:
        Phylo.draw(tree, axes=ax_tree, do_show=False)
    else:
        Phylo.draw(tree, axes=ax_tree, do_show=False, label_func=lambda x: None)

    ax_tree.set_axis_off()
    for spine in ax_tree.spines.values():
        spine.set_visible(False)
    ax_tree.patch.set_visible(False)

    # ---- TF heatmap panels ----
    # optionally use a shared vmax for consistent colour scale
    if global_vmax is None:
        # compute from finite values across selected TFs
        vals = []
        for tf in tfs:
            A = np.ma.masked_invalid(M_by_tf[tf]).compressed()
            if A.size:
                vals.append(np.nanmax(A))
        global_vmax = max(vals) if vals else None

    ims = []
    for i, tf in enumerate(tfs, start=1):
        ax = fig.add_subplot(gs[0, i])
        masked = np.ma.masked_invalid(M_by_tf[tf])
        im = ax.imshow(masked, aspect="auto", interpolation="nearest",
                       vmin=0, vmax=global_vmax)
        ims.append((ax, im, tf))

        # y labels only on first heatmap
        ax.set_yticks(range(n))
        if i == 1:
            ax.set_yticklabels(ylabels)
        else:
            ax.set_yticklabels([])
        ax.set_xlabel("Alignment column")
        ax.set_title(tf)

    # One shared colourbar on the last axis (clean)
    cax = ims[-1][0]
    fig.colorbar(ims[-1][1], ax=[ax for ax, _, _ in ims], label="-log10(p)")

    plt.savefig(out_png, dpi=400)
    plt.close()
    print(f"[tree+multi] wrote {out_png}")

def prune_tree_to_mammals(tree, species_in_data):
    """
    Return a pruned copy of the tree containing only the mammalian clade.

    Uses MRCA of mammal 'anchor' species present in the tree.
    This avoids relying on internal node labels which are usually absent.
    """
    # candidates spanning eutherians + (optionally) marsupials/monotremes
    anchors = ["hg38", "mm10", "rn6", "bosTau", "bosTau9", "canFam3", "canFam", "monDom5", "ornAna2", "ornAna1"]
    tree_leaf_names = {leaf.name for leaf in tree.get_terminals()}

    # Keep anchors that exist in THIS tree; allow prefix matches like bosTau9 vs bosTau
    present_anchors = []
    for a in anchors:
        if a in tree_leaf_names:
            present_anchors.append(a)
        else:
            # try prefix match for e.g. bosTau -> bosTau9
            matches = [n for n in tree_leaf_names if n.startswith(a)]
            if matches:
                present_anchors.append(matches[0])

    # Must have at least 2 anchors to define an MRCA clade
    present_anchors = [a for a in present_anchors if a in species_in_data]
    if len(present_anchors) < 2:
        raise RuntimeError(f"Not enough mammal anchors found in tree/data overlap: {present_anchors}")

    mrca = tree.common_ancestor(present_anchors)

    mammal_leaves = {t.name for t in mrca.get_terminals()}

    # prune tree to mammal leaves
    for leaf in list(tree.get_terminals()):
        if leaf.name not in mammal_leaves:
            tree.prune(leaf)

    return tree

def get_species_from_maf(maf_path):
    """
    Parse a MAF file and return a sorted list of species IDs.

    Species ID is taken as: second field on 's' lines, up to the first '.'
    e.g. 'hg38.chr9' -> 'hg38'
    """
    species = set()
    with open(maf_path) as f:
        for line in f:
            if not line.startswith("s "):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            src = fields[1]          # e.g. 'hg38.chr9'
            sp = src.split(".")[0]   # -> 'hg38'
            species.add(sp)
    return sorted(species)

def load_alignment(maf_path, species_prefixes=None):
    """
    Load aligned sequences for the species of interest from a MAF.
    
    If species_prefixes is None, derive them from 's' lines in the MAF.
    Returns:
        aln_seqs: dict species -> aligned sequence (may contain '-')
    """
    if species_prefixes is None:
        species_prefixes = get_species_from_maf(maf_path)
    
    aln_seqs = {sp: [] for sp in species_prefixes}
    
    for block in AlignIO.parse(maf_path, "maf"):
        for rec in block:
            src = rec.id        # e.g. 'hg38.chr9'
            sp = src.split(".")[0]
            if sp in aln_seqs:
                aln_seqs[sp].append(str(rec.seq))
    
    # Concatenate all blocks for each species
    for sp in aln_seqs:
        aln_seqs[sp] = "".join(aln_seqs[sp])

    return aln_seqs
    
def build_maps(aligned_seq):
    """
    aligned_seq: string with A/C/G/T/- of length L.
    Returns:
        aln2seq: list of length L; aln2seq[col] = seq_pos (0-based) or None if gap
        seq2aln: list of length ungapped_len; seq2aln[seq_pos] = alignment column (0-based)
    """
    aln2seq = []
    seq2aln = []
    seq_pos = -1

    for col, base in enumerate(aligned_seq):
        if base == "-":
            aln2seq.append(None)
        else:
            seq_pos += 1
            aln2seq.append(seq_pos)
            seq2aln.append(col)

    return aln2seq, seq2aln

def read_fimo_combined(path):
    rows = []
    with open(path) as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

def map_hits_to_alignment(fimo_rows, maps, aln_len, cap_logp=20.0):
    """
    fimo_rows: list of dicts from fimo_combined.tsv
    maps: dict species -> {"aln2seq":..., "seq2aln":...}
    aln_len: length of the alignment (same for all species)
    Returns:
        scores: dict[species] -> np.array(aln_len) of max -log10(p) per column (NaN if no hit)
    """
    scores = {}
    for sp, m in maps.items():
        scores[sp] = np.full(aln_len, np.nan, dtype=float)

    for r in fimo_rows:
        sp = r["species"]
        if sp not in maps:
            continue

        seq2aln = maps[sp]["seq2aln"]

        try:
            p = float(r["p_value"])
        except (TypeError, ValueError):
            continue
        if p <= 0:
            continue

        s = int(r["start"])  # 1-based
        e = int(r["end"])    # 1-based inclusive from FIMO
        # convert to 0-based half-open
        s0 = s - 1
        e0 = e               # since e is inclusive, e0 is already the right bound

        val = -np.log10(p)
        if cap_logp is not None:
            val = min(val, cap_logp)

        # For each seq position in the motif, map to alignment column
        for seq_pos in range(s0, e0):
            if seq_pos < 0 or seq_pos >= len(seq2aln):
                continue
            aln_col = seq2aln[seq_pos]
            if np.isnan(scores[sp][aln_col]):
                scores[sp][aln_col] = val
            else:
                scores[sp][aln_col] = max(scores[sp][aln_col], val)

    return scores

def load_tf_calls_from_merge_all(fimo_dir, species_list):
    rows = []
    for sp in species_list:
        path = os.path.join(fimo_dir, f"{sp}__TF_promoter_calls.tsv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue

        required = {"TF", "start", "stop", "best_p"}
        if not required.issubset(df.columns):
            raise ValueError(f"{path} missing required columns: {required - set(df.columns)}")

        for _, r in df.iterrows():
            try:
                rows.append({
                    "species": sp,
                    "TF": str(r["TF"]),
                    "start": int(r["start"]),   # 1-based
                    "end": int(r["stop"]),      # 1-based inclusive
                    "p_value": float(r["best_p"])
                })
            except Exception:
                continue
    return rows

def get_tree_leaf_order_biophylo(newick_path, species_in_data):
    newick_path = os.path.abspath(os.path.expanduser(newick_path))
    print(f"[tree] Reading Newick tree from: {newick_path}")

    tree = Phylo.read(newick_path, "newick")

    # prune to mammalian clade using available species IDs
    tree = prune_tree_to_mammals(tree, set(species_in_data))

    tree_leaves = [leaf.name for leaf in tree.get_terminals()]
    print(f"[tree] Number of leaves in (pruned) tree: {len(tree_leaves)}")
    print("[tree] First 20 leaf names:", tree_leaves[:20])

    species_in_data = list(species_in_data)
    data_set = set(species_in_data)

    ordered = []
    for leaf in tree.get_terminals():
        if leaf.name in data_set:
            ordered.append(leaf.name)

    print(f"[tree] Number of species found in both tree and data: {len(ordered)}")
    if not ordered:
        print("[tree][WARN] No overlap between tree leaves and species_in_data!")

    return ordered, tree

def dedupe_fimo_sites(rows):
    """
    Collapse FIMO rows that represent the same physical site:
    (species, motif_file, sequence_name, start, end), keeping the lowest p-value.
    Strand is ignored for deduplication.
    """
    best = {}
    for r in rows:
        key = (
            r["species"],
            r["motif_file"],
            r.get("sequence_name", ""),
            r["start"],
            r["end"],
        )
        try:
            p = float(r["p_value"])
        except (TypeError, ValueError):
            p = 1.0
        if key not in best or p < float(best[key]["p_value"]):
            best[key] = r
    return list(best.values())

with open("/home/nleathlm/workspace/ifne/data/trees/ucsc_genomes_official.json") as f:
    data = json.load(f)

genomes = data["ucscGenomes"]   # <-- this is the dict with ailMel1, hg38, etc.

id_to_common = {}
for db, meta in genomes.items():
    if not isinstance(meta, dict):
        continue
    # UCSC uses "organism" for the common name, "genome" is often the same string
    name = meta.get("organism") or meta.get("genome") or db
    id_to_common[db] = name

print(f"Loaded {len(id_to_common)} assembly → common-name mappings")
for test_id in ["hg38", "panTro4", "otoGar3", "ailMel1", "cavPor3"]:
    print(test_id, "→", id_to_common.get(test_id))

def plot_tree_and_heatmap_biophylo(tree_newick, ordered_species, M, out_png, id_to_common=None):
    tree = Phylo.read(tree_newick, "newick")
    tree = prune_tree_to_mammals(tree, set(ordered_species))
    
    # Prune unwanted leaves
    keep = set(ordered_species)
    for leaf in list(tree.get_terminals()):
        if leaf.name not in keep:
            tree.prune(leaf)

    # Rename leaves to common names
    if id_to_common is not None:
        for leaf in tree.get_terminals():
            if leaf.name in id_to_common:
                leaf.name = id_to_common[leaf.name]

    common_labels = [id_to_common.get(sp, sp) for sp in ordered_species] if id_to_common else ordered_species

    fig = plt.figure(figsize=(14, 0.35 * len(ordered_species)))
    gs = GridSpec(1, 2, width_ratios=[1.2, 4])

    # Tree panel
    ax_tree = fig.add_subplot(gs[0, 0])
    Phylo.draw(tree, axes=ax_tree, do_show=False)

    # Remove the box/frame and all axis elements
    ax_tree.set_axis_off()
    for spine in ax_tree.spines.values():
        spine.set_visible(False)
    ax_tree.patch.set_visible(False)

    # Heatmap panel
    ax_heat = fig.add_subplot(gs[0, 1])
    masked = np.ma.masked_invalid(M)
    im = ax_heat.imshow(masked, aspect="auto", interpolation="nearest")
    ax_heat.set_yticks(range(len(ordered_species)))
    ax_heat.set_yticklabels(common_labels)
    ax_heat.set_xlabel("Alignment column")
    fig.colorbar(im, ax=ax_heat, label=f"-log10(p) ({TF_NAME} motif)")
    #fig.colorbar(im, ax=ax_heat, label="-log10(p) (PGR motif)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[tree+heatmap] wrote {out_png}")

def resolve_common_name(sp, mapping):
    if sp in mapping:
        return mapping[sp]
    # try prefix match (bosTau -> bosTau9)
    for k, v in mapping.items():
        if k.startswith(sp) or sp.startswith(k):
            return v
    return sp

def build_presence_matrix(fimo_dir, ordered_species, id_to_common=None):
    """
    Returns:
      P: (n_species x n_TF) matrix of 0/1
      tf_list: list of TFs (columns)
      ylabels: common names
    """
    dfs = []
    for sp in ordered_species:
        path = os.path.join(fimo_dir, f"{sp}__TF_promoter_calls.tsv")
        if not os.path.exists(path):
            continue
        try:
            df = pd.read_csv(path, sep="\t")
        except pd.errors.EmptyDataError:
            continue
        if df.empty:
            continue
        # only keep present calls
        if "present" in df.columns:
            df = df[df["present"] == True]
        dfs.append(df[["species", "TF"]])

    if not dfs:
        return None, None, None

    all_df = pd.concat(dfs, ignore_index=True)

    # presence per species×TF (any hit => present)
    pres = (
        all_df.assign(present=1)
        .drop_duplicates(["species", "TF"])
        .pivot_table(index="species", columns="TF", values="present", fill_value=0, aggfunc="max")
    )

    # Ensure all ordered_species are present (missing => 0)
    pres = pres.reindex(index=ordered_species, fill_value=0)

    tf_list = list(pres.columns)
    ylabels = [id_to_common.get(sp, sp) for sp in ordered_species] if id_to_common else ordered_species
    return pres.values.astype(int), tf_list, ylabels

from collections import defaultdict
import numpy as np

def map_hits_to_alignment_by_tf(fimo_rows, maps, aln_len, cap_logp=20.0):
    """
    Returns:
      scores_by_tf: dict[TF] -> dict[species] -> np.array(aln_len) with max -log10(p) per column
    """
    scores_by_tf = defaultdict(lambda: {
        sp: np.full(aln_len, np.nan, dtype=float) for sp in maps.keys()
    })

    for r in fimo_rows:
        sp = r["species"]
        tf = r["TF"]
        if sp not in maps:
            continue

        try:
            p = float(r["p_value"])
        except (TypeError, ValueError):
            continue
        if p <= 0:
            continue

        s0 = int(r["start"]) - 1
        e0 = int(r["end"])      # inclusive -> half-open bound already correct in your original logic

        val = -np.log10(p)
        if cap_logp is not None:
            val = min(val, cap_logp)

        seq2aln = maps[sp]["seq2aln"]
        for seq_pos in range(s0, e0):
            if seq_pos < 0 or seq_pos >= len(seq2aln):
                continue
            aln_col = seq2aln[seq_pos]
            cur = scores_by_tf[tf][sp][aln_col]
            scores_by_tf[tf][sp][aln_col] = val if np.isnan(cur) else max(cur, val)

    return scores_by_tf

def plot_presence_heatmap(P, tf_list, ylabels, out_png):
    # Discrete purple / yellow palette. Can be changed!
    cmap = ListedColormap(["#440154", "#FDE725"])  # viridis endpoints
    norm = BoundaryNorm([-0.5, 0.5, 1.5], cmap.N)

    plt.figure(figsize=(0.35 * len(tf_list) + 4, 0.32 * len(ylabels)))
    im = plt.imshow(P, aspect="auto", cmap=cmap, norm=norm)
    
    plt.yticks(range(len(ylabels)), ylabels)
    plt.xticks(range(len(tf_list)), tf_list, rotation=90)
    plt.xlabel("TF")
    plt.title("Motif presence across species")

    cbar = plt.colorbar(im, ticks=[0, 1])
    cbar.set_ticklabels(["Absent", "Present"])
    cbar.set_label("Motif presence")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[presence] wrote {out_png}")

# -------------------------------
# 1. LOAD SPECIES + ALIGNMENT
# -------------------------------

maf_path = os.path.expanduser("~/workspace/ifne/data/maf/IFNE_100way.maf")
print("Using MAF:", maf_path)

species_list = get_species_from_maf(maf_path)
print(f"Found {len(species_list)} species in MAF")

aln_seqs = load_alignment(maf_path, species_list)
aln_len = len(next(iter(aln_seqs.values())))
print("Alignment length:", aln_len)

maps = {}
for sp in species_list:
    aln2seq, seq2aln = build_maps(aln_seqs[sp])
    maps[sp] = {"aln2seq": aln2seq, "seq2aln": seq2aln}

# -------------------------------
# 2. MAP FIMO HITS TO ALIGNMENT
# -------------------------------

FIMO_DIR = "/home/nleathlm/workspace/ifne/analysis/fimo_out/"   # wherever merge_tf outputs live
TF_NAME = "PGR"  # or whichever TF you want

fimo_rows = load_tf_calls_from_merge_all(FIMO_DIR, species_list)
print("Example fimo_row:", fimo_rows[0] if fimo_rows else None)
scores = map_hits_to_alignment(fimo_rows, maps, aln_len)

print(f"Loaded {len(fimo_rows)} TF-call rows for {TF_NAME}")
print("Example row:", fimo_rows[0] if fimo_rows else "None")

# -------------------------------
# 3. GET TREE ORDER 
# -------------------------------

tree_path = "/home/nleathlm/workspace/ifne/data/trees/hg38_100way.nh"

ordered_species, tree = get_tree_leaf_order_biophylo(
    tree_path,
    species_list
)

print("Tree-ordered species:", ordered_species)

if len(ordered_species) == 0:
    raise RuntimeError("No overlap between tree leaves and MAF species IDs!")

# -------------------------------
# 4. BUILD TREE-ORDERED MATRIX
# -------------------------------

M = np.vstack([scores[sp] for sp in ordered_species])

# -------------------------------
# 5. SIMPLE TREE-ORDERED HEATMAP
# -------------------------------

masked = np.ma.masked_invalid(M)

common_labels = [resolve_common_name(sp, id_to_common) for sp in ordered_species]

# Hard-coded label overrides for IDs not covered (or wrong) in the JSON mapping

COMMON_NAME_OVERRIDES = {
    "tupChi1": "Chinese tree shrew",
    "jacJac1": "Lesser Egyptian jerboa",
    "micOch1": "Prairie vole",
    "chiLan1": "Long-tailed chinchilla",
    "octDeg1": "Degu",
    "camFer1": "Wild Bactrian camel",
    "orcOrc1": "Killer whale",
    "panHod1": "Tibetan antelope",
    "capHir1": "Goat",
    "odoRosDiv1": "Walrus",
    "lepWed1": "Weddell seal",
    "pteAle1": "Black flying-fox",
    "eptFus1": "Big brown bat",
    "myoDav1": "David's mouse-eared bat",
    "myoLuc2": "Little brown bat",
    "conCri1": "Star-nosed mole",
    "eleEdw1": "Cape elephant shrew",
    "chrAsi1": "Cape golden mole",
    "oryAfe1": "Aardvark"
}

# Apply overrides (override JSON mapping)
id_to_common.update(COMMON_NAME_OVERRIDES)

print("[heatmap] wrote tree-ordered alignment heatmap")

# -------------------------------
# 6. TREE + HEATMAP COMBINED PANEL
# -------------------------------

# Build M matrices for the TFs you want
tfs_to_plot = ["ESR1", "PGR","ELF1","ELF3"]

# Map all TF hits to the alignment
scores_by_tf = map_hits_to_alignment_by_tf(
    fimo_rows,
    maps,
    aln_len,
    cap_logp=20.0
)

M_by_tf = {}
for tf in tfs_to_plot:
    M_by_tf[tf] = np.vstack([scores_by_tf[tf][sp] for sp in ordered_species])

plot_tree_and_multi_heatmaps(
    tree_newick=tree_path,
    ordered_species=ordered_species,
    M_by_tf=M_by_tf,
    out_png="/home/nleathlm/workspace/ifne/analysis/figures/ifne_tree_panels.png",
    id_to_common=id_to_common,
    tf_order=tfs_to_plot,
    max_tfs=len(tfs_to_plot),
    width_tree=1.0,
    width_tf=1.4,               #  TF panel size 
    show_tree_tip_labels=False,
    global_vmax=None            # or set e.g. 20 for consistent scaling
)

P, tf_list, ylabels = build_presence_matrix(
    fimo_dir="/home/nleathlm/workspace/ifne/analysis/fimo_out",
    ordered_species=ordered_species,
    id_to_common=id_to_common
)
if P is not None:
    plot_presence_heatmap(
        P, tf_list, ylabels,
        "/home/nleathlm/workspace/ifne/analysis/figures/ifne_presence_species_by_TF.png"
    )
