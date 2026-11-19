#!/usr/bin/env bash
set -euo pipefail
if [[ $# -lt 2 ]]; then
  echo "Usage: $0 <fasta_dir> <fimo_out_dir>" >&2
  exit 2
fi
FASTA_DIR="$1"
FIMO_DIR="$2"
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
BG_DIR="${FIMO_DIR}/background"
mkdir -p "${FIMO_DIR}" "${BG_DIR}"
for fasta in "${FASTA_DIR}"/*.fa "${FASTA_DIR}"/*.fasta; do
  [[ -e "$fasta" ]] || continue
  species="$(basename "$fasta")"; species="${species%%.*}"
  for source in jaspar homococo; do
    for motif in "${SCRIPT_DIR}/../motifs/${source}"/*.meme; do
      [[ -e "$motif" ]] || continue
      bg_file="${BG_DIR}/${species}.bg"
      bash "${SCRIPT_DIR}/run_fimo.sh" "${species}" "${source}" "${motif}" "${fasta}" "${FIMO_DIR}" "${bg_file}"
    done
  done
done
