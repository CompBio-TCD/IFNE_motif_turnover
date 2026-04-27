#!/bin/bash
python select_and_classify_ifne.py \
  --protein-aln aligned_proteins.fas \
  --cds-fasta cds_sequences.fasta \
  --metadata-tsv metadata.tsv \
  --human-id Homo_sapiens \
  --conserved-tsv ifne_conserved_sites.tsv \
  --signalp-tsv signalp_results.tsv \
  --selected-output selected_representatives.tsv \
  --classification-output ifne_classification.tsv
