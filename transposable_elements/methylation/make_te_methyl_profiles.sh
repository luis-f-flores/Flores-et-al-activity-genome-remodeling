#!/usr/bin/env bash

BIGWIG_DIR="${BIGWIG_DIR:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/wgbs/methylKit_bedgraph_bigwig}"
OUT_DIR="${OUT_DIR:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_methyl_profiles}"
BEDS_DIR="${BEDS_DIR:-$OUT_DIR/beds}"

BIN_SIZE="${BIN_SIZE:-50}"
N_BINS="${N_BINS:-50}"
FLANK_BINS="${FLANK_BINS:-4}"

TIMES=(${TIMES:-1 6 24})

AVERAGE_TYPE="${AVERAGE_TYPE:-mean}"

COLOR_KCL="${COLOR_KCL:-#1f77b4}"
COLOR_MANN="${COLOR_MANN:-#d62728}"

LOCK_Y_PER_FAMILY=${LOCK_Y_PER_FAMILY:-1}   # 1 = lock across timepoints within family
PAD_DOWN_PP=${PAD_DOWN_PP:-2}               # subtract from family-wide min mean
PAD_UP_PP=${PAD_UP_PP:-2}                   # add to family-wide max mean
YMIN_FALLBACK=${YMIN_FALLBACK:-0}           # used if no data
YMAX_FALLBACK=${YMAX_FALLBACK:-100}         # used if no data
MAKE_JPG=${MAKE_JPG:-1}                     # convert PNG -> JPG if ImageMagick is present

set -euo pipefail
mkdir -p "$OUT_DIR/matrices" "$OUT_DIR/plots"

for t in computeMatrix plotProfile; do
  command -v "$t" >/dev/null || { echo "ERROR: deepTools '$t' not found"; exit 1; }
done

for tp in "${TIMES[@]}"; do
  [[ -f "$BIGWIG_DIR/${tp}K_mean.bw" ]] || { echo "Missing $BIGWIG_DIR/${tp}K_mean.bw"; exit 1; }
  [[ -f "$BIGWIG_DIR/${tp}M_mean.bw" ]] || { echo "Missing $BIGWIG_DIR/${tp}M_mean.bw"; exit 1; }
done

BODY_LEN=$(( N_BINS * BIN_SIZE ))
FLANK_LEN=$(( FLANK_BINS * BIN_SIZE ))
echo "[info] body=${BODY_LEN} bp (${N_BINS}×${BIN_SIZE}), flanks=±${FLANK_LEN} bp (${FLANK_BINS}×${BIN_SIZE}), keep true zeros"

shopt -s nullglob
for BED in "$BEDS_DIR"/*.bed; do
  FAM=$(basename "$BED" .bed)

  case "$FAM" in
    L1MdGf_I|L1MdTf_I) echo "[skip ALT_ONLY] $FAM"; continue;;
  esac

  echo "[family] $FAM"
  PLOT_DIR_FAM="$OUT_DIR/plots/$FAM"
  mkdir -p "$PLOT_DIR_FAM"

  declare -a TSVS=()

  for TP in "${TIMES[@]}"; do
    K="$BIGWIG_DIR/${TP}K_mean.bw"
    M="$BIGWIG_DIR/${TP}M_mean.bw"

    MAT="$OUT_DIR/matrices/${FAM}_${TP}h.mat.gz"
    TAB="${MAT%.gz}.tab"
    PDF="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.pdf"
    PNG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.png"
    JPG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.jpg"
    TSV="$OUT_DIR/matrices/${FAM}_${TP}h_profile.tsv"

    computeMatrix scale-regions \
      -R "$BED" \
      -S "$K" "$M" \
      --beforeRegionStartLength "$FLANK_LEN" \
      --regionBodyLength "$BODY_LEN" \
      --afterRegionStartLength "$FLANK_LEN" \
      --binSize "$BIN_SIZE" \
      --smartLabels \
      -o "$MAT" \
      --outFileNameMatrix "$TAB" \
      --numberOfProcessors max

    plotProfile \
      -m "$MAT" \
      --samplesLabel "${TP}h KCl" "${TP}h Mann" \
      --plotTitle "${FAM} — ${TP}h (KCl vs Mann)" \
      --yAxisLabel "%mCG" \
      --perGroup \
      --averageType "$AVERAGE_TYPE" \
      --regionsLabel "" \
      --colors "$COLOR_KCL" "$COLOR_MANN" \
      --outFileName "$PDF" \
      --outFileNameData "$TSV" \
      --plotFileFormat pdf \
      --numPlotsPerRow 1

    plotProfile \
      -m "$MAT" \
      --samplesLabel "${TP}h KCl" "${TP}h Mann" \
      --yAxisLabel "%mCG" \
      --perGroup \
      --averageType "$AVERAGE_TYPE" \
      --regionsLabel "" \
      --colors "$COLOR_KCL" "$COLOR_MANN" \
      --outFileName "$PNG" \
      --plotFileFormat png \
      --numPlotsPerRow 1

    if [[ "${MAKE_JPG:-1}" -eq 1 ]] && command -v convert >/dev/null 2>&1; then
      convert "$PNG" "$JPG" || true
    fi

    if command -v awk >/dev/null; then
      NAN_K=$(awk 'NR>1{c+=($2=="nan")} END{print (NR>1?100*c/(NR-1):0)}' "$TSV")
      NAN_M=$(awk 'NR>1{c+=($3=="nan")} END{print (NR>1?100*c/(NR-1):0)}' "$TSV")
      printf "  QC: %s %sh  NaN  KCl=%.1f%%  Mann=%.1f%%\n" "$FAM" "$TP" "$NAN_K" "$NAN_M"
    fi

    TSVS+=("$TSV")
  done

  if [[ "$LOCK_Y_PER_FAMILY" -eq 1 && ${#TSVS[@]} -gt 0 ]]; then
    read -r YMIN_FAM YMAX_FAM < <(
      awk -v pad_dn="$PAD_DOWN_PP" -v pad_up="$PAD_UP_PP" '
        FNR==1 { next }  # skip header row
        {
          for(i=2;i<=NF;i++){
            if($i!="nan"){ sum[i]+=$i; cnt[i]++ }
          }
        }
        ENDFILE {
          for(i in sum){
            if(cnt[i]>0){
              m = sum[i]/cnt[i]
              if(minm=="" || m<minm) minm=m
              if(maxm=="" || m>maxm) maxm=m
            }
          }
          delete sum; delete cnt
        }
        END {
          if(minm==""){minm=0} if(maxm==""){maxm=0}
          ymin = minm - pad_dn; if(ymin < 0) ymin = 0
          ymax = maxm + pad_up; if(ymax > 100) ymax = 100
          if(ymax <= ymin){ ymax = (ymin < 100 ? ymin+1 : 100); if(ymin >= ymax) ymin = 0 }
          printf "%.6f %.6f\n", ymin, ymax
        }
      ' "${TSVS[@]}"
    )

    echo "  [lock y] ${FAM} -> yMin=${YMIN_FAM}  yMax=${YMAX_FAM}"

    for TP in "${TIMES[@]}"; do
      MAT="$OUT_DIR/matrices/${FAM}_${TP}h.mat.gz"
      PDF="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.pdf"
      PNG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.png"
      JPG="$PLOT_DIR_FAM/${FAM}_${TP}h_profile.jpg"

      plotProfile \
        -m "$MAT" \
        --samplesLabel "${TP}h KCl" "${TP}h Mann" \
        --plotTitle "${FAM} — ${TP}h (KCl vs Mann)" \
        --yAxisLabel "%mCG" \
        --perGroup \
        --averageType "$AVERAGE_TYPE" \
        --regionsLabel "" \
        --colors "$COLOR_KCL" "$COLOR_MANN" \
        --plotFileFormat pdf \
        --numPlotsPerRow 1 \
        --yMin "$YMIN_FAM" --yMax "$YMAX_FAM" \
        --outFileName "$PDF"

      plotProfile \
        -m "$MAT" \
        --samplesLabel "${TP}h KCl" "${TP}h Mann" \
        --yAxisLabel "%mCG" \
        --perGroup \
        --averageType "$AVERAGE_TYPE" \
        --regionsLabel "" \
        --colors "$COLOR_KCL" "$COLOR_MANN" \
        --plotFileFormat png \
        --numPlotsPerRow 1 \
        --yMin "$YMIN_FAM" --yMax "$YMAX_FAM" \
        --outFileName "$PNG"

      if [[ "${MAKE_JPG:-1}" -eq 1 ]] && command -v convert >/dev/null 2>&1; then
        convert "$PNG" "$JPG" || true
      fi
    done
  fi

done

echo "[done] Plots -> $OUT_DIR/plots/<FAM> ; Matrices -> $OUT_DIR/matrices"
