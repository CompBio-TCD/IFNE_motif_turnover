#!/usr/bin/env python
import os, re, glob, math, csv, json, argparse
from collections import defaultdict, namedtuple

def read_fimo_tsv(path):
    rows = []
    with open(path, "r") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

def load_all(fimo_dir):
    all_rows = []
    for tsv in glob.glob(os.path.join(fimo_dir, "*.tsv")):
        rows = read_fimo_tsv(tsv)
        for r in rows:
            r["_file"] = os.path.basename(tsv)
        all_rows.extend(rows)
    return all_rows

def to_table(rows, q_thresh=0.05):
    """Convert FIMO rows to a normalized table with species, motif, score, p, q, start, stop, strand."""
    out = []
    for r in rows:
        sp = r["_file"].split("__")[0]
        motif = os.path.splitext(r["_file"].split("__")[1])[0]
        out.append({
            "species": sp,
            "motif_file": motif,
            "motif": r.get("motif_id") or r.get("motif_alt_id") or motif,
            "sequence_name": r.get("sequence_name","NA"),
            "strand": r.get("strand","."),
            "start": int(r.get("start","0")),
            "end": int(r.get("stop") or r.get("end") or r.get("stop","0")),
            "p_value": float(r.get("p-value","1") or 1),
            "q_value": float(r.get("q-value","1") or 1),
            "matched_sequence": r.get("matched_sequence","").upper(),
            "score": float(r.get("score","0")),
        })
    return out

def write_tsv(rows, path):
    if not rows:
        return
    keys = list(rows[0].keys())
    with open(path, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, delimiter="\t")
        w.writeheader()
        for r in rows:
            w.writerow(r)

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--fimo_dir", default="data/fimo_out")
    ap.add_argument("--out", default="analysis/tables/fimo_combined.tsv")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    rows = load_all(args.fimo_dir)
    tbl = to_table(rows)
    write_tsv(tbl, args.out)
    print(f"[parse_fimo] wrote {args.out} with {len(tbl)} rows")
