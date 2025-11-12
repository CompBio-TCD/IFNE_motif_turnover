#!/usr/bin/env python
import os, csv, argparse
import numpy as np
import matplotlib.pyplot as plt

def read_table(path):
    rows = []
    with open(path, "r") as f:
        rdr = csv.DictReader(f, delimiter="\t")
        for r in rdr:
            rows.append(r)
    return rows

def heatmap_scores(fimo_combined, out_png):
    # Build species x motif max -log10(p)
    rows = read_table(fimo_combined)
    species = sorted(set(r["species"] for r in rows))
    motifs = sorted(set(r["motif_file"] for r in rows))
    M = np.zeros((len(species), len(motifs)))
    for r in rows:
        i = species.index(r["species"])
        j = motifs.index(r["motif_file"])
        p = float(r["p_value"])
        if p > 0:
            val = -np.log10(p)
            M[i,j] = max(M[i,j], val)
    plt.figure()
    plt.imshow(M, aspect="auto")
    plt.xticks(range(len(motifs)), motifs, rotation=90)
    plt.yticks(range(len(species)), species)
    plt.colorbar(label="-log10(p) (max per species × motif)")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

def lollipop_composites(composites_tsv, out_png):
    rows = read_table(composites_tsv)
    # simple scatter of distances by partner family
    fams = sorted(set(r["partner_family"] for r in rows))
    x = []
    y = []
    labels = []
    for r in rows:
        x.append(float(r["distance_bp"]))
        y.append(fams.index(r["partner_family"])+1)
        labels.append(r["partner_family"])
    plt.figure()
    plt.scatter(x, y)
    plt.yticks(range(1, len(fams)+1), fams)
    plt.xlabel("Distance from half-PRE (bp)")
    plt.ylabel("Partner family")
    plt.tight_layout()
    plt.savefig(out_png, dpi=200)
    plt.close()

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--combined", default="analysis/tables/fimo_combined.tsv")
    ap.add_argument("--composites", default="analysis/tables/composites.tsv")
    ap.add_argument("--figdir", default="analysis/figures")
    args = ap.parse_args()
    os.makedirs(args.figdir, exist_ok=True)
    heatmap_scores(args.combined, os.path.join(args.figdir, "heatmap_scores.png"))
    lollipop_composites(args.composites, os.path.join(args.figdir, "composite_distances.png"))
    print("[plots] wrote figures")
