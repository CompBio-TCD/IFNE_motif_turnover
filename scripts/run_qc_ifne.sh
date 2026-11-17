#!/usr/bin/bash

mkdir -p ../analysis/qc

python qc_ifne_maf.py \
  --maf ../data/maf/IFNE_100way_vertebrates_hg38.maf \
  --out-prefix ../analysis/qc/ifne_promoter_qc \
  --write-aligned-fasta
