#!/usr/bin/env bash
set -euo pipefail

if [[ $# -lt 5 ]]; then
  echo "Usage: $0 <species> <source> <motif_file> <fasta> <outdir> [background]" >&2
  exit 2
fi

species="$1"
source="$2"
motifs="$3"
fasta="$4"
outdir="$5"
bg="${6:-}"

mkdir -p "${outdir}"

# ---- thresholds from config ----
PVAL="$(python - <<'PY'
import yaml
print(yaml.safe_load(open("../config/config.yaml"))["fimo"].get("p_value","1e-3"))
PY
 2>/dev/null || true)"
PVAL="${PVAL:-1e-3}"

# ---- background: auto-generate if not supplied ----
# cache backgrounds here
BG_DIR="../analysis/bg0"
mkdir -p "${BG_DIR}"

if [[ -z "${bg}" ]]; then
  bg="${BG_DIR}/${species}.bg"
fi

if [[ ! -f "${bg}" ]]; then
  if ! command -v fasta-get-markov >/dev/null 2>&1; then
    echo "[run_fimo] ERROR: fasta-get-markov not found, cannot generate background ${bg}" >&2
    echo "[run_fimo] Install MEME-suite (or add fasta-get-markov to PATH), or pass an existing bg file as arg6." >&2
    exit 1
  fi
  #echo "[run_fimo] Generating 0th-order background: ${bg} from ${fasta}"
  #fasta-get-markov -m 0 "${fasta}" "${bg}"
  echo "[run_fimo] Generating 1st-order background: ${bg} from ${fasta}"
  fasta-get-markov -m 1 "${fasta}" "${bg}"
fi

motif_base="$(basename "${motifs%.*}")"
run_dir="${outdir}/.fimo_tmp_${species}__${source}__${motif_base}"

rm -rf "${run_dir}"
mkdir -p "${run_dir}"

cmd=(fimo --verbosity 1 --thresh "${PVAL}" --oc "${run_dir}" --bgfile "${bg}")
cmd+=( "${motifs}" "${fasta}" )

"${cmd[@]}"

if [[ ! -f "${run_dir}/fimo.tsv" ]]; then
  echo "[run_fimo] ERROR: fimo.tsv not found in ${run_dir}" >&2
  exit 1
fi

# Strip comment lines starting with '#'
clean_tmp="${run_dir}/fimo.clean.tsv"
grep -v '^#' "${run_dir}/fimo.tsv" > "${clean_tmp}"

# Determine TF label:
# Prefer motif_alt_id from the first hit row (column 2).
# If there are no hits, fall back to motif_base.
tf_label="$(awk 'BEGIN{FS="\t"} NR==2 {print $2}' "${clean_tmp}" || true)"
tf_label="${tf_label:-$motif_base}"

# Sanitize TF label for filenames
tf_safe="$(echo "${tf_label}" | sed 's/[^A-Za-z0-9_.-]/_/g')"

out_tsv="${outdir}/${species}__${tf_safe}__${source}.tsv"
mv "${clean_tmp}" "${out_tsv}"

rm -rf "${run_dir}"

echo "[run_fimo] Wrote ${out_tsv}"
echo "[run_fimo] Used bg: ${bg}"
