#!/usr/bin/env bash

# ==============================================================================
# Script: 07_pick_and_extract_TE_families_for_plot_profile.sh
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

# 07_pick_and_extract_TE_families_for_plot_profile.sh
#
# Select TE families from TE differential-expression heatmap matrices and extract
# per-family BED files from a TE catalog for downstream methylation metaprofiles.
#

# Inputs:
#   --te-summary-dir   Directory containing:
#                        TE_DE_heatmap_SIGONLY_log2FC_matrix.csv
#                        TE_DE_heatmap_SIGONLY_stars_matrix.csv
#                      or the corresponding FULL matrices.
#   --te-bed           TE catalog BED file, usually TE_mm10.bed.
#                      Expected BED6 format:
#                        chr start end copyID score strand
#                      where copyID contains the TE family label.
#   --outdir           Output directory.
#

# Outputs:
#   top_families.txt
#   beds/<FAMILY>.bed
#   match_report.tsv
#
# Example:
#   bash 07_pick_and_extract_TE_families_for_plot_profile.sh \
#     --te-summary-dir /path/to/TE_DE_summary \
#     --te-bed /path/to/TE_catalog_mm10/TE_mm10.bed \
#     --outdir /path/to/TE_methyl_profiles \
#     --n-up 12 \
#     --n-down 12 \
#     --require-sig 1

set -euo pipefail

usage() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

TE_SUMMARY_DIR=""
TE_BED=""
OUT_DIR=""
FAMILY_COL=4
N_UP=12
N_DOWN=12
REQUIRE_SIG=1
OVERRIDES_CSV=""
VALIDATE_BED=""
SAMPLE_VALIDATE_N=10

while [[ $# -gt 0 ]]; do
  case "$1" in
    --te-summary-dir) TE_SUMMARY_DIR="$2"; shift 2;;
    --te-bed) TE_BED="$2"; shift 2;;
    --outdir) OUT_DIR="$2"; shift 2;;
    --family-col) FAMILY_COL="$2"; shift 2;;
    --n-up) N_UP="$2"; shift 2;;
    --n-down) N_DOWN="$2"; shift 2;;
    --require-sig) REQUIRE_SIG="$2"; shift 2;;
    --overrides-csv) OVERRIDES_CSV="$2"; shift 2;;
    --validate-bed) VALIDATE_BED="$2"; shift 2;;
    --sample-validate-n) SAMPLE_VALIDATE_N="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument: $1" >&2; usage;;
  esac
done

[[ -d "$TE_SUMMARY_DIR" ]] || { echo "ERROR: --te-summary-dir not found" >&2; exit 1; }
[[ -f "$TE_BED" ]] || { echo "ERROR: --te-bed not found" >&2; exit 1; }
[[ -n "$OUT_DIR" ]] || { echo "ERROR: --outdir required" >&2; exit 1; }
[[ "$FAMILY_COL" =~ ^[0-9]+$ ]] || { echo "ERROR: --family-col must be an integer" >&2; exit 1; }

command -v Rscript >/dev/null || { echo "ERROR: Rscript not found" >&2; exit 1; }

mkdir -p "$OUT_DIR/beds" "$OUT_DIR/validate"

LFC_SIGONLY="$TE_SUMMARY_DIR/TE_DE_heatmap_SIGONLY_log2FC_matrix.csv"
STAR_SIGONLY="$TE_SUMMARY_DIR/TE_DE_heatmap_SIGONLY_stars_matrix.csv"
LFC_FULL="$TE_SUMMARY_DIR/TE_DE_heatmap_FULL_log2FC_matrix.csv"
STAR_FULL="$TE_SUMMARY_DIR/TE_DE_heatmap_FULL_stars_matrix.csv"

if [[ -f "$LFC_SIGONLY" && -f "$STAR_SIGONLY" ]]; then
  LFC="$LFC_SIGONLY"
  STAR="$STAR_SIGONLY"
else
  LFC="$LFC_FULL"
  STAR="$STAR_FULL"
fi

[[ -f "$LFC" && -f "$STAR" ]] || {
  echo "ERROR: heatmap matrices not found in $TE_SUMMARY_DIR" >&2
  exit 1
}

TOP_LIST="$OUT_DIR/top_families.txt"

Rscript - <<RS
suppressPackageStartupMessages({
  library(readr)
  library(dplyr)
  library(tidyr)
})

lfc <- read_csv("$LFC", show_col_types = FALSE)
stars <- read_csv("$STAR", show_col_types = FALSE)

lfc_long <- pivot_longer(lfc, -c(family, type), names_to = "contrast", values_to = "lfc")
star_long <- pivot_longer(stars, -c(family, type), names_to = "contrast", values_to = "star")

x <- inner_join(lfc_long, star_long, by = c("family", "type", "contrast"))

if ($REQUIRE_SIG == 1) {
  x <- filter(x, !is.na(star), star != "")
}

x <- filter(x, !is.na(lfc))

ranked <- x %>%
  group_by(family) %>%
  summarize(
    max_up = max(lfc, na.rm = TRUE),
    max_down = min(lfc, na.rm = TRUE),
    .groups = "drop"
  )

top <- unique(c(
  ranked %>% arrange(desc(max_up)) %>% slice_head(n = $N_UP) %>% pull(family),
  ranked %>% arrange(max_down) %>% slice_head(n = $N_DOWN) %>% pull(family)
))

writeLines(top, "$TOP_LIST")
cat("Selected", length(top), "families -> $TOP_LIST\n")
RS

[[ -s "$TOP_LIST" ]] || { echo "ERROR: no TE families selected" >&2; exit 1; }

CAT_IDX="$OUT_DIR/.catalog_index.tsv"

awk -v c="$FAMILY_COL" 'BEGIN{FS=OFS="\t"}
{
  chrom=$1
  start=$2
  end=$3
  name=$c
  score=(NF>=5?$5:0)
  strand=(NF>=6?$6:"+")
  rep=name

  # The catalog may use FAMILY|copyID, repName@coord, or just FAMILY.
  sub(/\|.*/, "", rep)
  sub(/@.*/, "", rep)

  rep_norm=toupper(rep)
  gsub(/[^A-Z0-9]/, "", rep_norm)
  sub(/INT$/, "I", rep_norm)

  print chrom,start,end,rep,rep_norm,score,strand,name
}' "$TE_BED" > "$CAT_IDX"

OVR_MAP="$OUT_DIR/.overrides.tsv"

if [[ -n "$OVERRIDES_CSV" && -f "$OVERRIDES_CSV" ]]; then
  awk -F',' 'NR==1{next} NF>=2{
    from=$1; to=$2
    gsub(/\r/,"",from); gsub(/\r/,"",to)
    fn=toupper(from); gsub(/[^A-Z0-9]/,"",fn); sub(/INT$/,"I",fn)
    tn=toupper(to);   gsub(/[^A-Z0-9]/,"",tn); sub(/INT$/,"I",tn)
    print fn "\t" tn
  }' "$OVERRIDES_CSV" > "$OVR_MAP"
else
  : > "$OVR_MAP"
fi

REPORT="$OUT_DIR/match_report.tsv"
echo -e "family\tstrategy\tmatched_token_example\tcopies" > "$REPORT"

while IFS= read -r FAM; do
  [[ -z "$FAM" ]] && continue

  OUT_BED="$OUT_DIR/beds/${FAM}.bed"
  FAM_NORM="$(printf '%s' "$FAM" | tr '[:lower:]' '[:upper:]' | sed -E 's/[^A-Z0-9]//g' | sed -E 's/INT$/I/')"

  OVR_TO="$(awk -v fn="$FAM_NORM" 'BEGIN{FS=OFS="\t"} $1==fn{print $2; found=1} END{if(!found)print ""}' "$OVR_MAP")"

  if [[ -n "$OVR_TO" ]]; then
    awk -v t="$OVR_TO" 'BEGIN{FS=OFS="\t"} $5==t{print $1,$2,$3,$8,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="override"
  else
    awk -v t="$FAM_NORM" 'BEGIN{FS=OFS="\t"} $5==t{print $1,$2,$3,$8,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="norm_equal"
  fi

  if [[ ! -s "$OUT_BED" ]]; then
    BASE="$(printf '%s' "$FAM_NORM" | sed -E 's/I$//')"
    awk -v b="$BASE" 'BEGIN{FS=OFS="\t"} ($5 ~ "^" b){print $1,$2,$3,$8,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="base_prefix"
  fi

  if [[ -s "$OUT_BED" ]]; then
    N=$(wc -l < "$OUT_BED" | tr -d ' ')
    EX=$(head -n1 "$OUT_BED" | awk '{print $4}')
    echo -e "${FAM}\t${STRAT}\t${EX}\t${N}" >> "$REPORT"
    echo "[match] $FAM -> $N copies [$STRAT]"
  else
    echo -e "${FAM}\tno_match\t.\t0" >> "$REPORT"
    echo "[no match] $FAM; check naming or provide --overrides-csv" >&2
    rm -f "$OUT_BED"
  fi
done < "$TOP_LIST"

if [[ -n "$VALIDATE_BED" && -f "$VALIDATE_BED" ]]; then
  if command -v bedtools >/dev/null 2>&1; then
    SUMMARY="$OUT_DIR/validate/validation_summary.tsv"
    echo -e "family\tn_matched\tn_tested" > "$SUMMARY"

    while IFS= read -r FAM; do
      BED="$OUT_DIR/beds/${FAM}.bed"
      [[ -s "$BED" ]] || continue

      SAMPLE="$OUT_DIR/validate/${FAM}.sample.bed"
      head -n "$SAMPLE_VALIDATE_N" "$BED" > "$SAMPLE" || true
      [[ -s "$SAMPLE" ]] || continue

      HIT=$(bedtools intersect -u -f 0.9 -r -a "$SAMPLE" -b "$VALIDATE_BED" | wc -l | tr -d ' ')
      TOT=$(wc -l < "$SAMPLE" | tr -d ' ')
      echo -e "${FAM}\t${HIT}\t${TOT}" >> "$SUMMARY"
    done < "$TOP_LIST"

    echo "[validate] $SUMMARY"
  else
    echo "[validate] bedtools not found; skipping validation"
  fi
fi

rm -f "$CAT_IDX" "$OVR_MAP"

echo "[done] Families: $TOP_LIST"
echo "[done] BEDs: $OUT_DIR/beds"
echo "[done] Report: $REPORT"
