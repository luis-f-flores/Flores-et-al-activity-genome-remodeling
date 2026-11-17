#!/bin/bash
set -euo pipefail

[[ $# -eq 5 ]] || { echo "Usage: gene_body_metagene_profile.sh genes.bed output_dir bigwigs_1h bigwigs_6h bigwigs_recovery" >&2; exit 1; }
GENES_BED="$1"
OUT_DIR="$2"
BW_1H="$3"
BW_6H="$4"
BW_RECOVERY="$5"
[[ -s "$GENES_BED" ]] || { echo "Missing gene-body BED: $GENES_BED" >&2; exit 1; }
N_GENES=$(awk 'NF && $1 !~ /^#/ {n++} END {print n}' "$GENES_BED")
N_LABEL=$(python3 -c 'import sys; print(f"{int(sys.argv[1]):,}")' "$N_GENES")
mkdir -p "$OUT_DIR"

for tp in 1hr 6hr 24hr; do
  case "$tp" in
    1hr) BW_DIR="$BW_1H" ;;
    6hr) BW_DIR="$BW_6H" ;;
    24hr) BW_DIR="$BW_RECOVERY" ;;
  esac
  KCL_BW="$BW_DIR/16_genomewide_KCl_${tp}_mean_pct_meth.bw"
  MANN_BW="$BW_DIR/16_genomewide_Mannitol_${tp}_mean_pct_meth.bw"
  MATRIX="$OUT_DIR/gene_body_matrix_${tp}.gz"
  OUT_PDF="$OUT_DIR/gene_body_profile_${tp}.pdf"
  case "$tp" in
    1hr) TITLE="1 hr" ;;
    6hr) TITLE="6 hr" ;;
    24hr) TITLE="Recovery" ;;
  esac

  [[ -f "$KCL_BW" && -f "$MANN_BW" ]] || {
    echo "Missing $tp methylation bigWigs; run condition_methylation_bigwigs.R $tp first." >&2
    exit 1
  }
  [[ ! -e "$MATRIX" && ! -e "$OUT_PDF" ]] || {
    echo "Choose a new output directory: $OUT_DIR" >&2
    exit 1
  }
  computeMatrix scale-regions \
    -R "$GENES_BED" -S "$KCL_BW" "$MANN_BW" \
    --regionBodyLength 10000 --binSize 100 -b 5000 -a 5000 \
    --samplesLabel "KCl" "Mannitol" -p 4 -o "$MATRIX"

  plotProfile -m "$MATRIX" -o "$OUT_PDF" --perGroup \
    --regionsLabel "n=$N_LABEL stable loci >=2 kb" --plotTitle "$TITLE" \
    --yAxisLabel "Mean %methylation" --startLabel "TSS" --endLabel "TES" \
    --colors "#d73027" "#4575b4" --yMin 0 --yMax 100 \
    --plotHeight 7 --plotWidth 8
done
echo "Gene-body methylation profiles written to $OUT_DIR"
