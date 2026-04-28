#!/usr/bin/env bash

# ==============================================================================
# Script: 02_annotate_te_de.sh
# Package: Transposable element analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for TE catalog construction, TE differential-expression annotation, TE-gene overlap, IAP summaries, methylation metaprofiles, and heatmaps.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# Annotate SalmonTE differential expression CSV files with TE metadata.

# Usage:
#   bash 02_annotate_te_de.sh --te-family-summary TE_family_summary.tsv --de-dir SalmonTE_DE --pattern '*TE_DE_results.csv'
set -euo pipefail
CAT=""; DE_DIR=""; PATTERN="*TE_DE_results.csv"
usage(){ echo "Usage: $0 --te-family-summary <TE_family_summary.tsv> --de-dir <dir> [--pattern '*TE_DE_results.csv']"; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --te-family-summary) CAT="$2"; shift 2;;
    --de-dir) DE_DIR="$2"; shift 2;;
    --pattern) PATTERN="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument $1" >&2; usage;;
  esac
done
[[ -f "$CAT" ]] || { echo "ERROR: TE family summary not found" >&2; exit 1; }
[[ -d "$DE_DIR" ]] || { echo "ERROR: DE directory not found" >&2; exit 1; }

ALIAS=$(mktemp)
trap 'rm -f "$ALIAS"' EXIT
cat > "$ALIAS" <<'ALIASES'
from	to
IAPA_MM	IAPA-int
IAPEZI	IAPEZ-int
IAPEYI	IAPEY-int
IAPEY3_I	IAPEY3-int
IAPEY4_I	IAPEY4-int
IAPEY5_I	IAPEY5-int
MERVL_LTR	MT2_Mm
IAP1-MM_I	IAP1_MM-int
MMERVK10D3_I	MMERVK10D3-int
RLTR4I_MM	RLTR4_MM
RLTR6I_MM	RLTR6_MM
ERVB4_2B-LTR_MM	ERVB4_2-LTR_MM
RLTR1IAP_MM	IAPLTR1_Mm
ALIASES

shopt -s nullglob
FILES=("$DE_DIR"/$PATTERN)
[[ ${#FILES[@]} -gt 0 ]] || { echo "ERROR: no files matched $DE_DIR/$PATTERN" >&2; exit 1; }

for de in "${FILES[@]}"; do
  out="${de%/*}/TE_DE_results_annotated.tsv"
  echo "Annotating $de -> $out"
  awk -v OFS='\t' -v cat="$CAT" -v alias="$ALIAS" '
    BEGIN{
      FS=",";
      while((getline l<cat)>0){if(l ~ /^family\t/)continue; split(l,f,"\t"); fam=f[1]; rc[fam]=f[2]; rf[fam]=f[3]; nc[fam]=f[4]; mpd[fam]=f[6]; mlen[fam]=f[7]}
      while((getline a<alias)>0){if(a ~ /^from/)continue; split(a,x,"\t"); if(x[1] && x[2]) al[x[1]]=x[2]}
      print "family","baseMean","log2FoldChange","lfcSE","pvalue","padj","repClass","repFamily","n_copies","mean_percDiv","mean_length_bp"
    }
    NR>1{
      gsub(/\r/,""); fam=$1; gsub(/^"|"$/, "", fam); gsub(/^[ \t]+|[ \t]+$/, "", fam); key=(fam in al ? al[fam] : fam)
      print fam,$2,$3,$4,$5,$6,(key in rc?rc[key]:"NA"),(key in rf?rf[key]:"NA"),(key in nc?nc[key]:"NA"),(key in mpd?mpd[key]:"NA"),(key in mlen?mlen[key]:"NA")
    }' "$de" > "$out"
done

echo "QC: rows with missing TE metadata"
for f in "$DE_DIR"/*TE_DE_results_annotated.tsv; do
  [[ -f "$f" ]] || continue
  printf "%s\tNA_rows=" "$(basename "$(dirname "$f")")"
  awk -F'\t' 'NR>1 && ($7=="NA" || $8=="NA"){c++} END{print c+0}' "$f"
done
