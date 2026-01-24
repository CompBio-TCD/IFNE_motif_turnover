#!/usr/bin/env python3
import os
import re
import glob
import yaml
import pandas as pd
import numpy as np
from pathlib import Path

# -----------------------------
# Load config
# -----------------------------
CONFIG_PATH = Path("../config/config.yaml").resolve()
with open(CONFIG_PATH) as fh:
    cfg = yaml.safe_load(fh)

# resolve fimo_dir relative to config file directory
FIMO_DIR = (CONFIG_PATH.parent / cfg["output"]["fimo_dir"]).resolve()
#P_THRESH = float(cfg["fimo"]["p_value"]) # if reading in from config
#Q_THRESH = float(cfg["fimo"]["q_value"])

Q_THRESH = 0.1
P_THRESH = 0.0001

# -----------------------------
# Helpers
# -----------------------------

def bh_qvalues(pvals: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg q-values for an array of p-values."""
    p = np.asarray(pvals, dtype=float)
    n = p.size
    order = np.argsort(p)
    ranks = np.empty(n, dtype=int)
    ranks[order] = np.arange(1, n+1)
    q = p * n / ranks
    # enforce monotonicity
    q_sorted = q[order]
    q_sorted = np.minimum.accumulate(q_sorted[::-1])[::-1]
    q_final = np.empty_like(q_sorted)
    q_final[order] = np.clip(q_sorted, 0, 1)
    return q_final

def load_and_filter_fimo(tsv_path, p_thresh, q_thresh):
    try:
        df = pd.read_csv(tsv_path, sep="\t")
    except pd.errors.EmptyDataError:
        print(f"[merge][WARN] Empty TSV (skipping): {tsv_path}")
        return pd.DataFrame()
    df = pd.read_csv(tsv_path, sep="\t")
    df["p-value"] = pd.to_numeric(df["p-value"], errors="coerce")
    df = df.dropna(subset=["p-value"])

    # Always compute BH q-values so you can inspect them
    df["q_bh"] = bh_qvalues(df["p-value"].values)

    # Parse FIMO q-values if present
    if "q-value" in df.columns:
        df["q_fimo"] = pd.to_numeric(df["q-value"], errors="coerce")
    else:
        df["q_fimo"] = np.nan

    # Decide which q to use
    if df["q_fimo"].notna().any():
        df["q_used"] = df["q_fimo"]
        df["q_source"] = "fimo"
        # Fill missing q_fimo values with BH so they aren't auto-dropped
        na_mask = df["q_used"].isna()
        if na_mask.any():
            df.loc[na_mask, "q_used"] = df.loc[na_mask, "q_bh"]
            df.loc[na_mask, "q_source"] = "fimo_or_bh_fill"
    else:
        df["q_used"] = df["q_bh"]
        df["q_source"] = "bh_from_p"

    if df["q_source"].eq("bh_from_p").any():
        print(
            os.path.basename(tsv_path),
            "BH q summary:",
            df["q_bh"].describe()
        )

    df_filt = df[(df["p-value"] <= p_thresh) & (df["q_used"] <= q_thresh)].copy()
    df_filt["evidence"] = "q_pass"
    return df_filt

def parse_filename(path):
    """
    Expect: <species>__<TF>__<source>.tsv
    """
    name = Path(path).name
    m = re.match(
        r"(?P<species>[^_]+)__(?P<tf>.+)__(?P<source>jaspar|homococo)\.tsv$",
        name
    )
    if not m:
        raise ValueError(f"Unexpected filename format: {name}")
    return m.group("species"), m.group("tf"), m.group("source")


# -----------------------------
# Merge logic
# -----------------------------
def make_tf_promoter_calls(species):
    paths = glob.glob(f"{FIMO_DIR}/{species}__*__*.tsv")

    all_hits = []
    for path in paths:
        sp, tf, source = parse_filename(path)
        if sp != species:
            continue

        hits = load_and_filter_fimo(path, p_thresh=P_THRESH, q_thresh=Q_THRESH)
        if hits.empty:
            continue

        hits["species"] = sp
        hits["TF"] = tf
        hits["source"] = source

        keep = [
            "species", "TF", "source", "evidence",
            "sequence_name", "start", "stop", "strand",
            "score", "p-value",
            "q_bh", "q_fimo", "q_used", "q_source",
            "matched_sequence"
        ]
        if "q-value" in hits.columns:
            keep.append("q-value")
        keep.append("matched_sequence")

        all_hits.append(hits[keep])

    if not all_hits:
        return pd.DataFrame()

    hits_df = pd.concat(all_hits, ignore_index=True)

    # Prefer q_pass over p_only, then smallest p-value
    hits_df["evidence_rank"] = hits_df["evidence"].map(
        {"q_pass": 0, "p_only": 1}
    )

    hits_df = hits_df.sort_values(
        ["species", "sequence_name", "TF", "evidence_rank", "p-value"]
    )

    best = hits_df.groupby(
        ["species", "sequence_name", "TF"],
        as_index=False
    ).first()

    counts = hits_df.groupby(
        ["species", "sequence_name", "TF"]
    ).size().reset_index(name="n_hits")

    out = best.merge(
        counts, on=["species", "sequence_name", "TF"], how="left"
    )

    out = out.rename(columns={
        "source": "best_source",
        "evidence": "best_evidence",
        "p-value": "best_p",
        "q-value": "best_q"
    })

    if "best_q" not in out.columns:
        out["best_q"] = pd.NA

    out["present"] = True

    # flag whether ANY source gave q-pass support
    any_q = (
        hits_df["evidence"].eq("q_pass")
        .groupby([hits_df["species"], hits_df["sequence_name"], hits_df["TF"]])
        .any()
        .reset_index(name="any_q_support")
    )

    out = out.merge(
        any_q, on=["species", "sequence_name", "TF"], how="left"
    )

    return out

# -----------------------------
# Run for all species present
# -----------------------------
species_set = {
    parse_filename(p)[0]
    for p in glob.glob(f"{FIMO_DIR}/*__*__*.tsv")
}

for sp in sorted(species_set):
    calls = make_tf_promoter_calls(sp)

    # If no TF–promoter calls, write nothing
    if calls.empty:
        print(f"[merge] {sp}: no TF–promoter calls (skipping)")
        continue

    out_path = f"{FIMO_DIR}/{sp}__TF_promoter_calls.tsv"
    calls.to_csv(out_path, sep="\t", index=False)
    print(f"[merge] wrote {out_path} ({len(calls)} TF–promoter calls)")
