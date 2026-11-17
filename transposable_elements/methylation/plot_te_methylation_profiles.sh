#!/usr/bin/env bash
set -euo pipefail

OUT_DIR="${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_methyl_profiles"
TIMES=(1 6 24)
AVERAGE_TYPE="mean"
COLOR_KCL="#1f77b4"
COLOR_MANN="#d62728"

declare -A yMin yMax
yMin[RLTR13B3]=25;         yMax[RLTR13B3]=75
yMin[IAPEY4_I]=40;         yMax[IAPEY4_I]=65
yMin[IAPEZI]=40;           yMax[IAPEZI]=60
yMin[IAPA_MM]=20;          yMax[IAPA_MM]=100
yMin[IAPLTR1a_I_MM]=30;    yMax[IAPLTR1a_I_MM]=85

for FAM in "${!yMin[@]}"; do
  echo "[replot] $FAM  y=(${yMin[$FAM]}, ${yMax[$FAM]})"
  PLOT_DIR_FAM="$OUT_DIR/plots/$FAM"
  mkdir -p "$PLOT_DIR_FAM"

  for TP in "${TIMES[@]}"; do
    MAT="$OUT_DIR/matrices/${FAM}_${TP}h.mat.gz"
    [[ -f "$MAT" ]] || { echo "  [skip] $TPh: matrix not found -> $MAT"; continue; }

    PDF="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.pdf"
    PNG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.png"
    JPG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.jpg"
    TSV="$OUT_DIR/matrices/${FAM}_${TP}h_profile.tsv"

    plotProfile -m "$MAT" \
      --samplesLabel "${TP}h KCl" "${TP}h Mann" \
      --plotTitle "${FAM} — ${TP}h (KCl vs Mann)" \
      --yAxisLabel "%mCG" \
      --perGroup --averageType "$AVERAGE_TYPE" --regionsLabel "" \
      --colors "$COLOR_KCL" "$COLOR_MANN" \
      --numPlotsPerRow 1 --plotFileFormat pdf \
      --yMin "${yMin[$FAM]}" --yMax "${yMax[$FAM]}" \
      --outFileName "$PDF" --outFileNameData "$TSV"

    plotProfile -m "$MAT" \
      --samplesLabel "${TP}h KCl" "${TP}h Mann" \
      --yAxisLabel "%mCG" \
      --perGroup --averageType "$AVERAGE_TYPE" --regionsLabel "" \
      --colors "$COLOR_KCL" "$COLOR_MANN" \
      --numPlotsPerRow 1 --plotFileFormat png \
      --yMin "${yMin[$FAM]}" --yMax "${yMax[$FAM]}" \
      --outFileName "$PNG"

    if command -v convert >/dev/null 2>&1; then convert "$PNG" "$JPG" || true; fi
  done
done

echo "[done] Replotted selected families with fixed y-limits."
