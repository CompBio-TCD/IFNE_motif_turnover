import os,csv
import matplotlib.pyplot as plt
import numpy as np
import json

from Bio import AlignIO 
from Bio import Phylo
from matplotlib.gridspec import GridSpec

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

def get_tree_leaf_order_biophylo(newick_path, species_in_data):
    """
    Extract leaf order from UCSC Newick tree using Biopython.
    Only returns species that are present in your data.

    Adds debug prints so you can see what leaf names the tree actually has.
    """
    newick_path = os.path.abspath(os.path.expanduser(newick_path))
    print(f"[tree] Reading Newick tree from: {newick_path}")

    tree = Phylo.read(newick_path, "newick")

    # Leaves as they appear in the tree
    tree_leaves = [leaf.name for leaf in tree.get_terminals()]
    print(f"[tree] Number of leaves in tree: {len(tree_leaves)}")
    print("[tree] First 20 leaf names:", tree_leaves[:20])

    # Species you have from MAF/FIMO
    species_in_data = list(species_in_data)
    print(f"[tree] Number of species in data: {len(species_in_data)}")
    print("[tree] First 20 species in data:", species_in_data[:20])

    data_set = set(species_in_data)

    ordered = []
    for leaf in tree.get_terminals():
        if leaf.name in data_set:
            ordered.append(leaf.name)

    print(f"[tree] Number of species found in both tree and data: {len(ordered)}")
    if not ordered:
        print("[tree][WARN] No overlap between tree leaves and species_in_data!")
        # Helpful suggestion:
        # You may need to map between IDs (e.g. hg38 vs Homo_sapiens, or hg19 vs hg38)

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
    ax_tree.set_xticks([])
    ax_tree.set_yticks([])

    # Heatmap panel
    ax_heat = fig.add_subplot(gs[0, 1])
    masked = np.ma.masked_invalid(M)
    im = ax_heat.imshow(masked, aspect="auto", interpolation="nearest")
    ax_heat.set_yticks(range(len(ordered_species)))
    ax_heat.set_yticklabels(common_labels)
    ax_heat.set_xlabel("Alignment column")
    fig.colorbar(im, ax=ax_heat, label="-log10(p) (PGR motif)")

    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()
    print(f"[tree+heatmap] wrote {out_png}")

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

fimo_rows = read_fimo_combined(
    "/home/nleathlm/workspace/ifne/analysis/tables/fimo_combined.tsv"
)

#fimo_rows = dedupe_fimo_sites(fimo_rows_raw)
scores = map_hits_to_alignment(fimo_rows, maps, aln_len)
# scores[sp] is now a vector of length aln_len


# -------------------------------
# 3. GET TREE ORDER (THIS WAS MISSING BEFORE)
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

common_labels = [id_to_common.get(sp, sp) for sp in ordered_species]

plt.figure(figsize=(12, 0.35 * len(ordered_species)))
masked = np.ma.masked_invalid(M)
plt.imshow(masked, aspect="auto", interpolation="nearest")
plt.yticks(range(len(ordered_species)), common_labels)
plt.xlabel("Alignment column (IFNE promoter alignment)")
plt.colorbar(label="-log10(p) (PGR motif)")
plt.tight_layout()
plt.savefig("/home/nleathlm/workspace/ifne/analysis/figures/ifne_pgr_alignment_heatmap_tree_order.png", dpi=200)
plt.close()

print("[heatmap] wrote tree-ordered alignment heatmap")


# -------------------------------
# 6. TREE + HEATMAP COMBINED PANEL
# -------------------------------

plot_tree_and_heatmap_biophylo(
    "/home/nleathlm/workspace/ifne/data/trees/hg38_100way.nh",
    ordered_species,
    M,
    "/home/nleathlm/workspace/ifne/analysis/figures/ifne_pgr_tree_alignment_common.png",
    id_to_common=id_to_common
)

