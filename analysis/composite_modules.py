#!/usr/bin/env python
import os, csv, re, argparse, math
from collections import defaultdict

def load_table(path):
    rows = []
    with open(path, "r") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            r["start"] = int(r["start"])
            r["end"] = int(r["end"])
            r["q_value"] = float(r["q_value"])
            r["p_value"] = float(r["p_value"])
            rows.append(r)
    return rows

def find_half_pre(seq):
    # Return list of half-PRE matches with start positions (1-based)
    seq = seq.upper()
    hits = []
    for m in re.finditer(r"TGTTCT|AGAACA", seq):
        hits.append(("PGR_HALF", "+", m.start()+1, m.group()))
    # reverse complement
    comp = str.maketrans("ACGT","TGCA")
    rc = seq.translate(comp)[::-1]
    L = len(seq)
    for m in re.finditer(r"TGTTCT|AGAACA", rc):
        # map back: rc pos -> original strand coordinates (approximate)
        start_rc = m.start()+1
        start = L - (start_rc - 1) - 5
        hits.append(("PGR_HALF", "-", start, m.group()))
    return hits

def build_index(rows, q_thr):
    idx = defaultdict(list)
    for r in rows:
        if r["q_value"] <= q_thr:
            idx[(r["species"], r["sequence_name"])].append(r)
    return idx

def within(a, b, min_d, max_d):
    return min_d <= abs(a-b) <= max_d

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="analysis/tables/fimo_combined.tsv")
    ap.add_argument("--fasta_dir", default="data/fasta")
    ap.add_argument("--q", type=float, default=0.05)
    ap.add_argument("--min_spacing", type=int, default=5)
    ap.add_argument("--max_spacing", type=int, default=60)
    ap.add_argument("--out", default="analysis/tables/composites.tsv")
    args = ap.parse_args()

    rows = load_table(args.combined)
    idx = build_index(rows, args.q)

    # Load FASTAs (simple loader, one record per file)
    fastas = {}
    for fn in os.listdir(args.fasta_dir):
        if not fn.endswith((".fa",".fasta")): continue
        with open(os.path.join(args.fasta_dir, fn)) as f:
            lines = [l.strip() for l in f if l.strip()]
        seq = "".join([l for l in lines if not l.startswith(">")])
        name = fn.split(".")[0]
        fastas[(name, name)] = seq  # sequence_name defaults to species name

    out_rows = []
    for key, feats in idx.items():
        sp, seqname = key
        seq = fastas.get((sp, seqname))
        if not seq:
            continue
        halves = find_half_pre(seq)
        # Extract partner motifs (ETS/SP1/AP1) by file name heuristic
        ets_sites = [r for r in feats if "ETS" in r["motif_file"].upper()]
        sp1_sites = [r for r in feats if "SP1" in r["motif_file"].upper()]
        ap1_sites = [r for r in feats if "AP1" in r["motif_file"].upper()]

        for anchor in halves:
            _, strand, pos, kmer = anchor
            for fam, sites in [("ETS", ets_sites), ("SP1", sp1_sites), ("AP1", ap1_sites)]:
                for s in sites:
                    # distance computed on start coordinates
                    d = s["start"] - pos
                    if within(s["start"], pos, args.min_spacing, args.max_spacing):
                        out_rows.append({
                            "species": sp,
                            "sequence_name": seqname,
                            "anchor_pos": pos,
                            "anchor_strand": strand,
                            "anchor_seq": kmer,
                            "partner_family": fam,
                            "partner_motif": s["motif"],
                            "partner_pos": s["start"],
                            "partner_q": s["q_value"],
                            "distance_bp": d
                        })

    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    with open(args.out, "w", newline="") as f:
        w = csv.DictWriter(f, delimiter="\t", fieldnames=list(out_rows[0].keys()) if out_rows else
                           ["species","sequence_name","anchor_pos","anchor_strand","anchor_seq","partner_family","partner_motif","partner_pos","partner_q","distance_bp"])
        w.writeheader()
        for r in out_rows:
            w.writerow(r)

    print(f"[composites] wrote {args.out} with {len(out_rows)} rows")
