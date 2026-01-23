#!/usr/bin/env bash
set -euo pipefail

# Usage: bash scripts/run_fimo.sh <species> <motif_file> <fasta> <outdir> [background]
species="${1}"
motifs="${2}"
fasta="${3}"
outdir="${4}"
bg="${5:-}"

mkdir -p "${outdir}"

# FIMO parameters
PVAL="$(python -c 'import yaml,sys; print(yaml.safe_load(open("config/config.yaml"))["fimo"]["p_value"])')"
# If you don't have python-yaml, set PVAL manually above.

cmd=(fimo --verbosity 1 --thresh "${PVAL}")
if [[ -n "${bg}" && -f "${bg}" ]]; then
  cmd+=( --bgfile "${bg}" )
fi
cmd+=( "${motifs}" "${fasta}" )

# Write TSV
"${cmd[@]}" > "${outdir}/${species}__$(basename "${motifs%.*}").tsv"
echo "[run_fimo] Wrote ${outdir}/${species}__$(basename "${motifs%.*}").tsv"
