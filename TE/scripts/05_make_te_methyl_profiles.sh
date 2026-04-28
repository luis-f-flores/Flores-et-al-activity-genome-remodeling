#!/usr/bin/env bash

# ==============================================================================
# Script: 05_make_te_methyl_profiles.sh
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

# Generate WGBS methylation metaprofiles over TE copies using deepTools.

# Usage:
#   bash 05_make_te_methyl_profiles.sh --bigwig-dir group_mean_bigwigs --beds-dir TE_family_beds --outdir TE_methyl_profiles --times '1 6 24'
set -euo pipefail
BIGWIG_DIR=""; BEDS_DIR=""; OUT_DIR=""; TIMES="1 6 24"; BIN_SIZE=50; N_BINS=50; FLANK_BINS=4; AVERAGE_TYPE="mean"; LOCK_Y=1; PAD=2
usage(){ echo "Usage: $0 --bigwig-dir <dir> --beds-dir <dir> --outdir <dir> [--times '1 6 24']"; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --bigwig-dir) BIGWIG_DIR="$2"; shift 2;;
    --beds-dir) BEDS_DIR="$2"; shift 2;;
    --outdir) OUT_DIR="$2"; shift 2;;
    --times) TIMES="$2"; shift 2;;
    --bin-size) BIN_SIZE="$2"; shift 2;;
    --n-bins) N_BINS="$2"; shift 2;;
    --flank-bins) FLANK_BINS="$2"; shift 2;;
    --average-type) AVERAGE_TYPE="$2"; shift 2;;
    --lock-y) LOCK_Y="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument $1" >&2; usage;;
  esac
done
[[ -d "$BIGWIG_DIR" ]] || { echo "ERROR: --bigwig-dir not found" >&2; exit 1; }
[[ -d "$BEDS_DIR" ]] || { echo "ERROR: --beds-dir not found" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || usage
mkdir -p "$OUT_DIR/matrices" "$OUT_DIR/plots"
command -v computeMatrix >/dev/null || { echo "ERROR: computeMatrix not found" >&2; exit 1; }
command -v plotProfile >/dev/null || { echo "ERROR: plotProfile not found" >&2; exit 1; }
read -r -a TIME_ARRAY <<< "$TIMES"
for tp in "${TIME_ARRAY[@]}"; do
  [[ -f "$BIGWIG_DIR/${tp}K_mean.bw" ]] || { echo "Missing ${tp}K_mean.bw" >&2; exit 1; }
  [[ -f "$BIGWIG_DIR/${tp}M_mean.bw" ]] || { echo "Missing ${tp}M_mean.bw" >&2; exit 1; }
done
BODY_LEN=$((N_BINS * BIN_SIZE)); FLANK_LEN=$((FLANK_BINS * BIN_SIZE))
shopt -s nullglob
BEDS=("$BEDS_DIR"/*.bed)
[[ ${#BEDS[@]} -gt 0 ]] || { echo "ERROR: no BED files in $BEDS_DIR" >&2; exit 1; }
for BED in "${BEDS[@]}"; do
  FAM=$(basename "$BED" .bed)
  case "$FAM" in L1MdGf_I|L1MdTf_I) echo "Skipping ALT-only family $FAM"; continue;; esac
  echo "[family] $FAM"
  PLOT_DIR="$OUT_DIR/plots/$FAM"; mkdir -p "$PLOT_DIR"; TSVS=()
  for TP in "${TIME_ARRAY[@]}"; do
    MAT="$OUT_DIR/matrices/${FAM}_${TP}h.mat.gz"; TSV="$OUT_DIR/matrices/${FAM}_${TP}h_profile.tsv"
    computeMatrix scale-regions -R "$BED" -S "$BIGWIG_DIR/${TP}K_mean.bw" "$BIGWIG_DIR/${TP}M_mean.bw" \
      --beforeRegionStartLength "$FLANK_LEN" --regionBodyLength "$BODY_LEN" --afterRegionStartLength "$FLANK_LEN" \
      --binSize "$BIN_SIZE" --smartLabels -o "$MAT" --outFileNameMatrix "${MAT%.gz}.tab" --numberOfProcessors max
    plotProfile -m "$MAT" --samplesLabel "${TP}h KCl" "${TP}h Mannitol" --plotTitle "${FAM} - ${TP}h" \
      --yAxisLabel "%mCG" --perGroup --averageType "$AVERAGE_TYPE" --regionsLabel "" --outFileName "$PLOT_DIR/${FAM}_${TP}h_profile.pdf" \
      --outFileNameData "$TSV" --plotFileFormat pdf --numPlotsPerRow 1
    plotProfile -m "$MAT" --samplesLabel "${TP}h KCl" "${TP}h Mannitol" --yAxisLabel "%mCG" --perGroup \
      --averageType "$AVERAGE_TYPE" --regionsLabel "" --outFileName "$PLOT_DIR/${FAM}_${TP}h_profile.png" --plotFileFormat png --numPlotsPerRow 1
    TSVS+=("$TSV")
  done
  if [[ "$LOCK_Y" -eq 1 ]]; then
    read -r YMIN YMAX < <(awk -v pad="$PAD" 'FNR==1{next}{for(i=2;i<=NF;i++) if($i!="nan"){sum[i]+=$i; cnt[i]++}} ENDFILE{for(i in sum){m=sum[i]/cnt[i]; if(min==""||m<min)min=m; if(max==""||m>max)max=m}; delete sum; delete cnt} END{if(min=="")min=0;if(max=="")max=0;ymin=min-pad;if(ymin<0)ymin=0;ymax=max+pad;if(ymax>100)ymax=100;if(ymax<=ymin)ymax=ymin+1;printf "%.6f %.6f\n",ymin,ymax}' "${TSVS[@]}")
    echo "  [locked y-axis] $YMIN $YMAX"
    for TP in "${TIME_ARRAY[@]}"; do
      MAT="$OUT_DIR/matrices/${FAM}_${TP}h.mat.gz"
      plotProfile -m "$MAT" --samplesLabel "${TP}h KCl" "${TP}h Mannitol" --plotTitle "${FAM} - ${TP}h" --yAxisLabel "%mCG" \
        --perGroup --averageType "$AVERAGE_TYPE" --regionsLabel "" --plotFileFormat pdf --numPlotsPerRow 1 --yMin "$YMIN" --yMax "$YMAX" \
        --outFileName "$PLOT_DIR/${FAM}_${TP}h_profile.pdf"
      plotProfile -m "$MAT" --samplesLabel "${TP}h KCl" "${TP}h Mannitol" --yAxisLabel "%mCG" --perGroup \
        --averageType "$AVERAGE_TYPE" --regionsLabel "" --plotFileFormat png --numPlotsPerRow 1 --yMin "$YMIN" --yMax "$YMAX" \
        --outFileName "$PLOT_DIR/${FAM}_${TP}h_profile.png"
    done
  fi
done
echo "Done. Outputs written to $OUT_DIR"
