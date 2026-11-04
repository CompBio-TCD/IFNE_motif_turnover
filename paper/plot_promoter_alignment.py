from pathlib import Path
import textwrap, py_compile

script = r'''#!/usr/bin/env python3
"""
make_promoter_figure.py

Publication 3-panel comparative promoter figure from the same
MAF / tree / FIMO-style outputs used by the IFNE promoter pipeline.

Panels
------
A  Motif prevalence across IFNE, IFNA1 and IFNB.
B  Phylogeny + species-by-motif presence/absence matrices for IFNE, IFNA1, IFNB.
C  Phylogeny + positional maps for selected IFNE motifs (default: ELF3, PGR, IRF1).

Input expectations
------------------
1) MAF files:
   UCSC/Multiz-style MAFs such as:
       data/maf/IFNE_100way.maf
       data/maf/IFNA1_100way.maf
       data/maf/IFNB_100way.maf

2) Tree:
   Newick tree whose leaf names match the MAF species prefixes, e.g.
       data/trees/hg38_100way.nh

3) UCSC genome metadata JSON (optional):
       data/trees/ucsc_genomes_official.json
   Used to convert UCSC assembly IDs to common names.

4) FIMO combined tables:
   The script accepts the column structure used by the plotting pipeline.
   Required information:
       species
       motif identity (one of: TF, motif_alt_id, motif_file, motif_id)
       hit start/end (for panel C)
       p-value (one of: p_value, p-value, best_p)
   Optional:
       q-value / q_used / q_fimo / best_q
       sequence_name

   If your combined TSV is already filtered, use --no-filter.
   Otherwise use --p-thresh and optionally --q-thresh.

Path handling
-------------
All CLI paths may be absolute or relative. Relative paths are resolved against
the CURRENT WORKING DIRECTORY.

If --config is supplied, it is only used for optional defaults:
    output.figure_dir
    fimo.p_value
    fimo.q_value

Example
-------
python make_promoter_figure.py \
  --config ../config/config.yaml \
  --ifne-maf ../data/maf/IFNE_100way.maf \
  --ifne-fimo ../analysis/tables/IFNE_fimo_combined.tsv \
  --ifna-maf ../data/maf/IFNA1_100way.maf \
  --ifna-fimo ../analysis/tables/IFNA1_fimo_combined.tsv \
  --ifnb-maf ../data/maf/IFNB_100way.maf \
  --ifnb-fimo ../analysis/tables/IFNB_fimo_combined.tsv \
  --tree ../data/trees/hg38_100way.nh \
  --ucsc-json ../data/trees/ucsc_genomes_official.json \
  --out-prefix ../analysis/figures/IFNE_promoter_figure

Outputs
-------
<out-prefix>.pdf
<out-prefix>.svg
<out-prefix>.png
<out-prefix>_panelA_prevalence.tsv
<out-prefix>_panelB_presence_IFNE.tsv
<out-prefix>_panelB_presence_IFNA1.tsv
<out-prefix>_panelB_presence_IFNB.tsv

Notes
-----
- The plotting code maps FIMO hits back to alignment columns using the same
  gap-stripped-sequence -> MAF-alignment coordinate logic as the original script.

"""

from __future__ import annotations

import argparse
import copy
import csv
import json
import os
import re
import sys
import warnings
from collections import defaultdict
from pathlib import Path

import numpy as np
import pandas as pd
import yaml
import matplotlib.pyplot as plt
from matplotlib import gridspec
from matplotlib.colors import LinearSegmentedColormap, ListedColormap
from Bio import AlignIO, Phylo


# ---------------------------------------------------------------------
# Defaults
# ---------------------------------------------------------------------

DEFAULT_MOTIFS = [
    "ELF1", "ELF3", "ELF5", "ETS1",
    "PGR", "ESR1", "ESR2", "AR",
    "IRF1", "IRF7", "IRF9", "NF-κB",
]

DEFAULT_POSITIONAL_MOTIFS = ["ELF3", "PGR", "IRF1"]

FAMILY_BREAKS = [3.5, 7.5]
FAMILY_LABELS = [
    (1.5, "ETS/ELF"),
    (5.5, "Steroid receptor"),
    (9.5, "Antiviral"),
]

# Complementary, colour-blind-conscious-ish palette.
GENE_COLORS = {
    "IFNE": "#5B3F8C",   # violet
    "IFNA1": "#2F80B7",  # blue
    "IFNB": "#2AA876",   # green/teal
}
ACCENT = "#E07A5F"       # warm coral
DARK = "#2B2B2B"
GRID = "#D9D9D9"
ABSENT = "#F2F2F2"


# ---------------------------------------------------------------------
# CLI / paths
# ---------------------------------------------------------------------

def parse_args():
    p = argparse.ArgumentParser(
        description="Create a 3-panel publication-ready comparative IFN promoter figure."
    )

    p.add_argument("--config", default=None, help="Optional YAML config file.")

    p.add_argument("--ifne-maf", required=True)
    p.add_argument("--ifne-fimo", required=True)

    p.add_argument("--ifna-maf", required=True)
    p.add_argument("--ifna-fimo", required=True)
    p.add_argument("--ifna-label", default="IFNA1",
                   help="Display label for the IFNA dataset (default: IFNA1).")

    p.add_argument("--ifnb-maf", required=True)
    p.add_argument("--ifnb-fimo", required=True)

    p.add_argument("--tree", required=True, help="Newick species tree.")
    p.add_argument("--ucsc-json", default=None,
                   help="Optional UCSC genome metadata JSON for common names.")

    p.add_argument(
        "--out-prefix",
        default=None,
        help="Output prefix. If omitted, uses config output.figure_dir or ./IFNE_promoter_figure."
    )

    p.add_argument(
        "--motifs",
        nargs="+",
        default=DEFAULT_MOTIFS,
        help="Motifs to show in panels A/B."
    )
    p.add_argument(
        "--positional-motifs",
        nargs="+",
        default=DEFAULT_POSITIONAL_MOTIFS,
        help="IFNE motifs to show in panel C."
    )

    p.add_argument("--p-thresh", type=float, default=None)
    p.add_argument("--q-thresh", type=float, default=None)
    p.add_argument(
        "--no-filter",
        action="store_true",
        help="Assume combined FIMO tables are already filtered."
    )

    p.add_argument("--dpi", type=int, default=600)
    p.add_argument(
        "--max-logp",
        type=float,
        default=20.0,
        help="Cap for -log10(p) colour scale in panel C."
    )
    p.add_argument(
        "--panel-c-hit-width",
        type=int,
        default=1,
        help="Visual-only half-width in alignment columns for panel C marks. "
             "Use 0 for exact columns, 1-2 if marks are hard to see."
    )
    p.add_argument(
        "--strict",
        action="store_true",
        help="Treat ambiguous species/multiple sequence records as errors rather than warnings."
    )

    return p.parse_args()


def abspath(path):
    return Path(path).expanduser().resolve()


def load_config(path):
    if path is None:
        return {}
    path = abspath(path)
    with open(path) as fh:
        return yaml.safe_load(fh) or {}


def get_nested(d, keys, default=None):
    cur = d
    for key in keys:
        if not isinstance(cur, dict) or key not in cur:
            return default
        cur = cur[key]
    return cur


# ---------------------------------------------------------------------
# MAF handling — mirrors the original plotting script
# ---------------------------------------------------------------------

def get_species_from_maf(maf_path):
    species = set()
    with open(maf_path) as f:
        for line in f:
            if not line.startswith("s "):
                continue
            fields = line.split()
            if len(fields) < 2:
                continue
            src = fields[1]
            species.add(src.split(".")[0])
    return sorted(species)


def load_alignment(maf_path, species_prefixes=None):
    """
    Same logic as the user's original script:
    concatenate aligned MAF chunks per species, preserving gap columns.
    """
    if species_prefixes is None:
        species_prefixes = get_species_from_maf(maf_path)

    aln_seqs = {sp: [] for sp in species_prefixes}

    for block in AlignIO.parse(str(maf_path), "maf"):
        for rec in block:
            sp = rec.id.split(".")[0]
            if sp in aln_seqs:
                aln_seqs[sp].append(str(rec.seq))

    for sp in aln_seqs:
        aln_seqs[sp] = "".join(aln_seqs[sp])

    lengths = {len(v) for v in aln_seqs.values() if v}
    if len(lengths) > 1:
        warnings.warn(
            f"{maf_path}: concatenated aligned sequences do not all have the same length: "
            f"{sorted(lengths)[:10]}"
        )

    return aln_seqs


def build_maps(aligned_seq):
    """
    Map gap-stripped sequence coordinates back to MAF alignment columns.
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


# ---------------------------------------------------------------------
# Tree / UCSC names
# ---------------------------------------------------------------------

def load_ucsc_common_names(json_path):
    if json_path is None:
        return {}

    with open(json_path) as f:
        data = json.load(f)

    genomes = data.get("ucscGenomes", {})
    out = {}
    for db, meta in genomes.items():
        if not isinstance(meta, dict):
            continue
        out[db] = meta.get("organism") or meta.get("genome") or db
    return out


def get_tree_leaf_order(newick_path, species_in_data):
    tree = Phylo.read(str(newick_path), "newick")
    data_set = set(species_in_data)

    ordered = [
        leaf.name for leaf in tree.get_terminals()
        if leaf.name in data_set
    ]

    if not ordered:
        raise RuntimeError(
            f"No overlap between tree leaves and species IDs in data.\n"
            f"Tree: {newick_path}"
        )

    return ordered


def prune_and_relabel_tree(tree_path, ordered_species, id_to_common):
    tr = Phylo.read(str(tree_path), "newick")
    keep = set(ordered_species)

    for leaf in list(tr.get_terminals()):
        if leaf.name not in keep:
            tr.prune(leaf)

    for leaf in tr.get_terminals():
        if leaf.name in id_to_common:
            leaf.name = id_to_common[leaf.name]

    return tr


# ---------------------------------------------------------------------
# FIMO input parsing
# ---------------------------------------------------------------------

P_COLS = ["p_value", "p-value", "best_p", "pvalue"]
Q_COLS = ["q_used", "q-value", "q_value", "q_fimo", "best_q"]
START_COLS = ["start"]
END_COLS = ["end", "stop"]
MOTIF_COLS = ["TF", "motif_alt_id", "motif_file", "motif_id"]
SPECIES_COLS = ["species"]


def first_existing(columns, candidates, required=False, label="column"):
    for c in candidates:
        if c in columns:
            return c
    if required:
        raise ValueError(
            f"Could not find required {label}. Tried: {', '.join(candidates)}. "
            f"Available columns: {', '.join(columns)}"
        )
    return None


def canonical_motif_name(value):
    """
    Convert common motif-file labels to a compact TF label.

    Examples:
      /path/ELF3.meme -> ELF3
      ELF3__jaspar -> ELF3
      MA0640.2_ELF3 -> ELF3 (best effort)
    """
    if pd.isna(value):
        return ""

    s = str(value).strip()
    s = Path(s).stem

    # common filename/source suffixes
    s = re.sub(r"__(jaspar|homococo|hocomoco)$", "", s, flags=re.I)
    s = re.sub(r"_(jaspar|homococo|hocomoco)$", "", s, flags=re.I)

    # Prefer a recognisable TF token if present.
    known = [
        "ELF1","ELF2","ELF3","ELF4","ELF5","ETS1",
        "PGR","ESR1","ESR2","AR",
        "IRF1","IRF7","IRF9",
        "NFKB","NF-KB","NFΚB","RELA",
    ]

    upper = s.upper().replace("κ", "K")
    for tf in known:
        if tf in upper:
            if tf in {"NFKB","NF-KB","NFΚB","RELA"}:
                return "NF-κB"
            return tf

    # Otherwise use final readable token.
    toks = re.split(r"[|:;,/ _-]+", s)
    toks = [t for t in toks if t]
    return toks[-1] if toks else s


def read_fimo_table(path, p_thresh=None, q_thresh=None, no_filter=False):
    df = pd.read_csv(path, sep="\t")

    species_col = first_existing(df.columns, SPECIES_COLS, required=True, label="species column")
    motif_col = first_existing(df.columns, MOTIF_COLS, required=True, label="motif column")
    p_col = first_existing(df.columns, P_COLS, required=True, label="p-value column")
    q_col = first_existing(df.columns, Q_COLS, required=False, label="q-value column")
    start_col = first_existing(df.columns, START_COLS, required=False)
    end_col = first_existing(df.columns, END_COLS, required=False)

    out = pd.DataFrame({
        "species": df[species_col].astype(str),
        "motif_raw": df[motif_col].astype(str),
        "p": pd.to_numeric(df[p_col], errors="coerce"),
    })

    out["motif"] = out["motif_raw"].map(canonical_motif_name)

    if q_col is not None:
        out["q"] = pd.to_numeric(df[q_col], errors="coerce")
    else:
        out["q"] = np.nan

    if start_col is not None:
        out["start"] = pd.to_numeric(df[start_col], errors="coerce")
    else:
        out["start"] = np.nan

    if end_col is not None:
        out["end"] = pd.to_numeric(df[end_col], errors="coerce")
    else:
        out["end"] = np.nan

    if "sequence_name" in df.columns:
        out["sequence_name"] = df["sequence_name"].astype(str)
    else:
        out["sequence_name"] = ""

    out = out.dropna(subset=["p"]).copy()

    if not no_filter:
        if p_thresh is not None:
            out = out[out["p"] <= p_thresh].copy()
        if q_thresh is not None and out["q"].notna().any():
            out = out[out["q"] <= q_thresh].copy()

    return out


def dedupe_hits(df):
    """
    Keep the smallest-p hit for duplicate species/motif/position records.
    """
    key = ["species", "motif", "start", "end"]
    use = [k for k in key if k in df.columns]
    if not use:
        return df

    return (
        df.sort_values("p")
          .drop_duplicates(subset=use, keep="first")
          .reset_index(drop=True)
    )


# ---------------------------------------------------------------------
# Summaries for panels A and B
# ---------------------------------------------------------------------

def species_presence_matrix(hits, species_order, motifs):
    mat = pd.DataFrame(0, index=species_order, columns=motifs, dtype=int)

    present = (
        hits[hits["motif"].isin(motifs)]
        .drop_duplicates(["species", "motif"])
    )

    for row in present.itertuples(index=False):
        if row.species in mat.index and row.motif in mat.columns:
            mat.loc[row.species, row.motif] = 1

    return mat


def prevalence_from_matrix(mat):
    if len(mat) == 0:
        return pd.Series(index=mat.columns, dtype=float)
    return mat.mean(axis=0)


# ---------------------------------------------------------------------
# Panel C mapping
# ---------------------------------------------------------------------

def map_hits_to_alignment(hits, maps, aln_len, motif, cap_logp=20.0):
    scores = {
        sp: np.full(aln_len, np.nan, dtype=float)
        for sp in maps
    }

    rows = hits[hits["motif"] == motif]

    for r in rows.itertuples(index=False):
        sp = r.species
        if sp not in maps:
            continue

        if pd.isna(r.start) or pd.isna(r.end):
            continue

        s = int(r.start)
        e = int(r.end)
        p = float(r.p)

        if p <= 0:
            continue

        seq2aln = maps[sp]["seq2aln"]

        # FIMO uses 1-based inclusive coordinates.
        s0 = s - 1
        e0 = e

        val = min(-np.log10(p), cap_logp)

        for seq_pos in range(s0, e0):
            if 0 <= seq_pos < len(seq2aln):
                col = seq2aln[seq_pos]
                if np.isnan(scores[sp][col]):
                    scores[sp][col] = val
                else:
                    scores[sp][col] = max(scores[sp][col], val)

    return scores


def widen_for_display(M, half_width=1):
    """
    VISUAL ONLY. Does not change analytical calls.
    """
    if half_width <= 0:
        return M

    out = np.full_like(M, np.nan, dtype=float)

    for i in range(M.shape[0]):
        cols = np.where(~np.isnan(M[i]))[0]
        for col in cols:
            lo = max(0, col - half_width)
            hi = min(M.shape[1], col + half_width + 1)
            value = M[i, col]

            current = out[i, lo:hi]
            replacement = np.full(hi - lo, value)

            if np.isnan(current).all():
                out[i, lo:hi] = replacement
            else:
                out[i, lo:hi] = np.nanmax(
                    np.vstack([current, replacement]), axis=0
                )

    return out


# ---------------------------------------------------------------------
# Plot helpers
# ---------------------------------------------------------------------

def display_labels(species_order, id_to_common):
    return [id_to_common.get(sp, sp) for sp in species_order]


def draw_tree(ax, tree_path, species_order, id_to_common):
    tr = prune_and_relabel_tree(tree_path, species_order, id_to_common)

    Phylo.draw(
        tr,
        axes=ax,
        do_show=False,
        show_confidence=False,
    )

    ax.set_xlabel("")
    ax.set_ylabel("")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_frame_on(False)


def italic_gene_label(ax, text):
    ax.text(
        0.5, 1.025, text,
        transform=ax.transAxes,
        ha="center", va="bottom",
        fontsize=10.5,
        fontstyle="italic",
        fontweight="normal",
        color=DARK,
    )


def clean_axis(ax):
    for s in ax.spines.values():
        s.set_visible(False)
    ax.tick_params(length=0)


# ---------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------

def main():
    args = parse_args()
    cfg = load_config(args.config)

    # Resolve paths from the caller's working directory.
    tree_path = abspath(args.tree)
    ucsc_json = abspath(args.ucsc_json) if args.ucsc_json else None

    datasets = {
        "IFNE": {
            "maf": abspath(args.ifne_maf),
            "fimo": abspath(args.ifne_fimo),
        },
        args.ifna_label: {
            "maf": abspath(args.ifna_maf),
            "fimo": abspath(args.ifna_fimo),
        },
        "IFNB": {
            "maf": abspath(args.ifnb_maf),
            "fimo": abspath(args.ifnb_fimo),
        },
    }

    if args.ifna_label != "IFNA1":
        # Reuse IFNA1 plotting colour for any IFNA display label.
        GENE_COLORS[args.ifna_label] = GENE_COLORS["IFNA1"]

    # thresholds
    p_thresh = args.p_thresh
    q_thresh = args.q_thresh

    if p_thresh is None:
        p_thresh = get_nested(cfg, ["fimo", "p_value"], default=1e-4)
        if p_thresh is not None:
            p_thresh = float(p_thresh)

    if q_thresh is None:
        q_thresh = get_nested(cfg, ["fimo", "q_value"], default=None)
        if q_thresh is not None:
            q_thresh = float(q_thresh)

    # output prefix
    if args.out_prefix:
        out_prefix = abspath(args.out_prefix)
    else:
        figdir = get_nested(cfg, ["output", "figure_dir"], default=None)
        if figdir:
            out_prefix = abspath(Path(figdir) / "IFNE_promoter_figure")
        else:
            out_prefix = abspath("IFNE_promoter_figure")

    out_prefix.parent.mkdir(parents=True, exist_ok=True)

    motifs = args.motifs
    positional_motifs = args.positional_motifs

    id_to_common = load_ucsc_common_names(ucsc_json) if ucsc_json else {}

    # ------------------------------------------------------------
    # Load all datasets
    # ------------------------------------------------------------
    gene_data = {}

    for gene, d in datasets.items():
        print(f"[load] {gene}")
        print(f"  MAF : {d['maf']}")
        print(f"  FIMO: {d['fimo']}")

        species_maf = get_species_from_maf(d["maf"])
        hits = read_fimo_table(
            d["fimo"],
            p_thresh=p_thresh,
            q_thresh=q_thresh,
            no_filter=args.no_filter,
        )
        hits = dedupe_hits(hits)

        species_hits = set(hits["species"])
        species_common = [sp for sp in species_maf if sp in species_hits]

        if not species_common:
            raise RuntimeError(
                f"{gene}: no overlap between MAF species and FIMO species."
            )

        order = get_tree_leaf_order(tree_path, species_common)

        missing_tree = sorted(set(species_common) - set(order))
        if missing_tree:
            msg = (
                f"{gene}: {len(missing_tree)} MAF/FIMO species are absent from the tree "
                f"and will be omitted: {missing_tree[:10]}"
            )
            if args.strict:
                raise RuntimeError(msg)
            warnings.warn(msg)

        # Warn on >1 sequence_name per species because original mapping is species-level.
        if "sequence_name" in hits.columns:
            multi = (
                hits[hits["sequence_name"].astype(str) != ""]
                .groupby("species")["sequence_name"]
                .nunique()
            )
            multi = multi[multi > 1]
            if len(multi):
                msg = (
                    f"{gene}: {len(multi)} species have >1 sequence_name in FIMO hits. "
                    "Panel C uses the original species-level concatenated-MAF coordinate mapping. "
                    "For targeted promoter MAFs, one sequence record per species is preferred."
                )
                if args.strict:
                    raise RuntimeError(msg)
                warnings.warn(msg)

        presence = species_presence_matrix(hits, order, motifs)

        gene_data[gene] = {
            "maf": d["maf"],
            "hits": hits,
            "order": order,
            "presence": presence,
        }

    # ------------------------------------------------------------
    # Shared species order for panel B
    # ------------------------------------------------------------
    # Use species present in all three datasets so rows are directly comparable.
    gene_names = list(gene_data.keys())

    shared = set(gene_data[gene_names[0]]["order"])
    for gene in gene_names[1:]:
        shared &= set(gene_data[gene]["order"])

    if not shared:
        raise RuntimeError(
            "No species are shared across IFNE, IFNA and IFNB datasets. "
            "Panel B requires a common species set."
        )

    shared_order = [
        sp for sp in get_tree_leaf_order(tree_path, shared)
        if sp in shared
    ]

    # Reindex B matrices to shared species.
    for gene in gene_names:
        gene_data[gene]["presence_shared"] = (
            gene_data[gene]["presence"]
            .reindex(shared_order)
            .fillna(0)
            .astype(int)
        )

    # ------------------------------------------------------------
    # Panel A prevalence
    # ------------------------------------------------------------
    prevalence = pd.DataFrame({"motif": motifs})

    for gene in gene_names:
        prevalence[gene] = (
            gene_data[gene]["presence_shared"]
            .mean(axis=0)
            .reindex(motifs)
            .values
        )

    prevalence.to_csv(
        f"{out_prefix}_panelA_prevalence.tsv",
        sep="\t", index=False
    )

    for gene in gene_names:
        gene_data[gene]["presence_shared"].to_csv(
            f"{out_prefix}_panelB_presence_{gene}.tsv",
            sep="\t"
        )

    # ------------------------------------------------------------
    # IFNE MAF mapping for C
    # ------------------------------------------------------------
    ifne = gene_data["IFNE"]
    aln_seqs = load_alignment(ifne["maf"], ifne["order"])

    nonempty = [v for v in aln_seqs.values() if v]
    if not nonempty:
        raise RuntimeError("IFNE MAF produced no aligned sequence.")

    aln_len = len(nonempty[0])

    maps = {}
    for sp in ifne["order"]:
        seq = aln_seqs.get(sp, "")
        if not seq:
            continue
        aln2seq, seq2aln = build_maps(seq)
        maps[sp] = {"aln2seq": aln2seq, "seq2aln": seq2aln}

    positional_matrices = {}

    for motif in positional_motifs:
        scores = map_hits_to_alignment(
            ifne["hits"],
            maps,
            aln_len,
            motif=motif,
            cap_logp=args.max_logp,
        )

        M = np.vstack([
            scores.get(sp, np.full(aln_len, np.nan))
            for sp in ifne["order"]
        ])

        positional_matrices[motif] = widen_for_display(
            M, half_width=args.panel_c_hit_width
        )

    # ------------------------------------------------------------
    # Figure
    # ------------------------------------------------------------
    n_pos = len(positional_motifs)

    fig = plt.figure(figsize=(13.2, 11.0))

    outer = gridspec.GridSpec(
        3, 1,
        height_ratios=[1.0, 2.1, 2.35],
        hspace=0.34,
        left=0.055, right=0.975,
        top=0.965, bottom=0.06,
    )

    # ===================== A =====================
    axA = fig.add_subplot(outer[0])

    x = np.arange(len(motifs))
    n_genes = len(gene_names)
    total_width = 0.72
    bar_w = total_width / n_genes
    offsets = (
        np.arange(n_genes) - (n_genes - 1) / 2
    ) * bar_w

    for offset, gene in zip(offsets, gene_names):
        color_key = gene if gene in GENE_COLORS else "IFNA1"
        axA.bar(
            x + offset,
            prevalence[gene].values,
            width=bar_w * 0.92,
            color=GENE_COLORS[color_key],
            label=gene,
        )

    axA.set_ylim(0, 1)
    axA.set_ylabel(
        "Proportion of species\nwith ≥1 significant motif",
        fontsize=9.5,
    )
    axA.set_xticks(x)
    axA.set_xticklabels(motifs, rotation=45, ha="right", fontsize=8.5)
    axA.tick_params(axis="y", labelsize=8.5)
    axA.spines["top"].set_visible(False)
    axA.spines["right"].set_visible(False)
    axA.grid(axis="y", color=GRID, linewidth=0.6)

    leg = axA.legend(
        frameon=False, ncol=n_genes,
        loc="upper right", fontsize=9
    )
    for t in leg.get_texts():
        t.set_fontstyle("italic")
        t.set_fontweight("normal")

    # Family separators/labels only when using the default 12-motif order.
    if motifs == DEFAULT_MOTIFS:
        for xpos in FAMILY_BREAKS:
            axA.axvline(xpos, color=GRID, linewidth=1.0)
        for xpos, label in FAMILY_LABELS:
            axA.text(
                xpos, 0.97, label,
                ha="center", va="top",
                fontsize=8.5, color=DARK
            )

    fig.text(
        0.012, 0.965, "A",
        fontsize=14, fontweight="bold",
        ha="left", va="top"
    )

    # ===================== B =====================
    gsB = gridspec.GridSpecFromSubplotSpec(
        1, 1 + len(gene_names),
        subplot_spec=outer[1],
        width_ratios=[1.6] + [1.8] * len(gene_names),
        wspace=0.07,
    )

    ax_tree_B = fig.add_subplot(gsB[0])
    draw_tree(ax_tree_B, tree_path, shared_order, id_to_common)

    for j, gene in enumerate(gene_names, start=1):
        ax = fig.add_subplot(gsB[j])

        color_key = gene if gene in GENE_COLORS else "IFNA1"
        cmap = ListedColormap([ABSENT, GENE_COLORS[color_key]])

        M = gene_data[gene]["presence_shared"].values

        ax.imshow(
            M,
            aspect="auto",
            interpolation="nearest",
            cmap=cmap,
            vmin=0, vmax=1,
        )

        ax.set_xticks(range(len(motifs)))
        ax.set_xticklabels(
            motifs, rotation=60, ha="right",
            fontsize=7.2
        )

        ax.set_yticks(range(len(shared_order)))
        ax.set_yticklabels([])
        clean_axis(ax)

        if motifs == DEFAULT_MOTIFS:
            for xpos in FAMILY_BREAKS:
                ax.axvline(xpos, color="white", linewidth=1.3)

        italic_gene_label(ax, gene)

    fig.text(
        0.012, 0.705, "B",
        fontsize=14, fontweight="bold",
        ha="left", va="top"
    )

    # ===================== C =====================
    gsC = gridspec.GridSpecFromSubplotSpec(
        1, 1 + n_pos,
        subplot_spec=outer[2],
        width_ratios=[1.6] + [1.8] * n_pos,
        wspace=0.08,
    )

    ax_tree_C = fig.add_subplot(gsC[0])
    draw_tree(
        ax_tree_C,
        tree_path,
        ifne["order"],
        id_to_common
    )

    sig_cmap = LinearSegmentedColormap.from_list(
        "ifne_sig",
        ["#F7F3FA", "#C9B8E4", GENE_COLORS["IFNE"], ACCENT]
    )
    sig_cmap.set_bad("white")

    last_im = None

    for j, motif in enumerate(positional_motifs, start=1):
        ax = fig.add_subplot(gsC[j])

        M = positional_matrices[motif]
        masked = np.ma.masked_invalid(M)

        last_im = ax.imshow(
            masked,
            aspect="auto",
            interpolation="nearest",
            cmap=sig_cmap,
            vmin=0,
            vmax=args.max_logp,
        )

        # Alignment columns, intentionally matching the original plotting script.
        ax.set_xlabel(
            "Alignment column\n(IFNE orthologous interval)",
            fontsize=7.8,
        )
        ax.set_yticks(range(len(ifne["order"])))
        ax.set_yticklabels([])
        ax.tick_params(length=0, labelsize=7)
        clean_axis(ax)

        ax.text(
            0.5, 1.025, motif,
            transform=ax.transAxes,
            ha="center", va="bottom",
            fontsize=9.5,
            fontweight="normal",
            color=DARK,
        )

    # IFNE label centred above C data panels.
    # This uses figure coordinates so it spans the motif subpanels.
    fig.text(
        0.705, 0.385, "IFNE",
        fontsize=10.5,
        fontstyle="italic",
        fontweight="normal",
        ha="center", va="bottom",
        color=DARK,
    )

    fig.text(
        0.012, 0.39, "C",
        fontsize=14, fontweight="bold",
        ha="left", va="top",
    )

    if last_im is not None:
        cax = fig.add_axes([0.922, 0.085, 0.012, 0.16])
        cb = fig.colorbar(last_im, cax=cax)
        cb.set_label("−log10(p)", fontsize=7.8)
        cb.ax.tick_params(labelsize=7)

    # ------------------------------------------------------------
    # Export
    # ------------------------------------------------------------
    pdf = f"{out_prefix}.pdf"
    svg = f"{out_prefix}.svg"
    png = f"{out_prefix}.png"

    fig.savefig(pdf, bbox_inches="tight")
    fig.savefig(svg, bbox_inches="tight")
    fig.savefig(png, dpi=args.dpi, bbox_inches="tight")
    plt.close(fig)

    print("[done]")
    print("  ", pdf)
    print("  ", svg)
    print("  ", png)


if __name__ == "__main__":
    main()
'''

path = Path("/mnt/data/make_promoter_figure.py")
path.write_text(script)
py_compile.compile(str(path), doraise=True)
print(path)
