#!/usr/bin/bash

mkdir -p ../analysis/qc

python qc_ifn_maf.py \
  --maf ../data/maf/IFNE_100way_vertebrates_hg38.maf \ # or whatever IFN MAF you're using
  --out-prefix ../analysis/qc/ifne_promoter_qc \ # update prefix based on IFN of interest
  --write-aligned-fasta
