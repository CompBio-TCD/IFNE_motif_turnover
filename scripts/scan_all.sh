#!/usr/bin/env bash
set -euo pipefail

# Iterate species × motif files
FIMO_DIR="$(python -c 'import yaml; print(yaml.safe_load(open("config/config.yaml"))["output"]["fimo_dir"])')"
BG_DIR="$(python -c 'import yaml; print(yaml.safe_load(open("config/config.yaml"))["fimo"]["background_dir"])')"

for fasta in data/fasta/*.fa data/fasta/*.fasta; do
  [[ -e "$fasta" ]] || continue
  species="$(basename "$fasta")"
  species="${species%%.*}"
  for motif in motifs/*.meme; do
    bg_file="${BG_DIR}/${species}.bg"
    if [[ -f "${bg_file}" ]]; then
      bash scripts/run_fimo.sh "${species}" "${motif}" "${fasta}" "${FIMO_DIR}" "${bg_file}"
    else
      bash scripts/run_fimo.sh "${species}" "${motif}" "${fasta}" "${FIMO_DIR}"
    fi
  done
done
