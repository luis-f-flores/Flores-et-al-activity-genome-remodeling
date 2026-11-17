#!/usr/bin/env bash

set -euo pipefail

TE_SUMMARY_DIR="${TE_SUMMARY_DIR:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_DE_summary}"
TE_CATALOG_DIR="${TE_CATALOG_DIR:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_catalog_mm10}"
TE_BED="${TE_BED:-$TE_CATALOG_DIR/TE_mm10.bed}"        # 6-column BED; col4 is repName@coord
FAMILY_COL="${FAMILY_COL:-4}"                           # 1-based col with repName@coord (we parse repName from it)
OUT_DIR="${OUT_DIR:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_methyl_profiles}"

N_UP="${N_UP:-12}"                  # top N most up
N_DOWN="${N_DOWN:-12}"              # top N most down
REQUIRE_SIG="${REQUIRE_SIG:-1}"     # 1=require starred cells in DE; 0=ignore significance

OVERRIDES_CSV="${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna/TE_name_overrides.csv"  # e.g., /path/TE_name_overrides.csv

VALIDATE_BED="${VALIDATE_BED:-}"    # e.g., /path/RepeatMasker_mm10.bed
SAMPLE_VALIDATE_N="${SAMPLE_VALIDATE_N:-10}"  # per family, how many copies to sample for check

LFC_SIGONLY="$TE_SUMMARY_DIR/TE_DE_heatmap_SIGONLY_log2FC_matrix.csv"
STAR_SIGONLY="$TE_SUMMARY_DIR/TE_DE_heatmap_SIGONLY_stars_matrix.csv"
LFC_FULL="$TE_SUMMARY_DIR/TE_DE_heatmap_FULL_log2FC_matrix.csv"
STAR_FULL="$TE_SUMMARY_DIR/TE_DE_heatmap_FULL_stars_matrix.csv"

command -v Rscript >/dev/null || { echo "ERROR: Rscript not found"; exit 1; }
[[ -f "$TE_BED" ]] || { echo "ERROR: TE_BED not found: $TE_BED"; exit 1; }
[[ "$FAMILY_COL" =~ ^[0-9]+$ ]] || { echo "ERROR: FAMILY_COL must be an integer"; exit 1; }
mkdir -p "$OUT_DIR/beds" "$OUT_DIR/validate"

if [[ -f "$LFC_SIGONLY" && -f "$STAR_SIGONLY" ]]; then
  LFC="$LFC_SIGONLY"; STAR="$STAR_SIGONLY"
else
  LFC="$LFC_FULL"; STAR="$STAR_FULL"
fi
[[ -f "$LFC" && -f "$STAR" ]] || { echo "ERROR: DE matrices not found in $TE_SUMMARY_DIR"; exit 1; }

TOP_LIST="$OUT_DIR/top_families.txt"
Rscript - <<RS
suppressPackageStartupMessages({library(readr);library(dplyr);library(tidyr)})
lfc   <- read_csv("$LFC", show_col_types=FALSE)
stars <- read_csv("$STAR", show_col_types=FALSE)
lfc_long <- pivot_longer(lfc,   -c(family,type), names_to="contrast", values_to="lfc")
st_long  <- pivot_longer(stars, -c(family,type), names_to="contrast", values_to="star")
x <- inner_join(lfc_long, st_long, by=c("family","type","contrast"))
if ($REQUIRE_SIG == 1) x <- filter(x, !is.na(star), star!="")
x <- filter(x, !is.na(lfc))
ranked <- x %>% group_by(family) %>%
  summarize(max_up=max(lfc, na.rm=TRUE), max_down=min(lfc, na.rm=TRUE), .groups="drop")
top <- unique(c(
  arrange(ranked, desc(max_up)) %>% slice_head(n=$N_UP) %>% pull(family),
  arrange(ranked, max_down)     %>% slice_head(n=$N_DOWN) %>% pull(family)
))
writeLines(top, "$TOP_LIST")
cat("Selected", length(top), "families -> $TOP_LIST\n")
RS
[[ -s "$TOP_LIST" ]] || { echo "ERROR: no families selected (empty $TOP_LIST)"; exit 1; }

CAT_IDX="$OUT_DIR/.catalog_index.tsv"
awk -v c="$FAMILY_COL" 'BEGIN{
  FS=OFS="\t"
}
{
  chrom=$1; start=$2; end=$3; name=$c; score=(NF>=5?$5:0); strand=(NF>=6?$6:"+");
  rep=name; sub(/@.*/,"",rep);                      # strip @chr:start-end
  rep_norm=toupper(rep);
  gsub(/[^A-Z0-9]/,"",rep_norm);                    # drop punctuation
  sub(/INT$/,"I",rep_norm);                         # INT == I
  print chrom,start,end,rep,rep_norm,score,strand
}' "$TE_BED" > "$CAT_IDX"

OVR_MAP="$OUT_DIR/.overrides.tsv"
if [[ -n "$OVERRIDES_CSV" && -f "$OVERRIDES_CSV" ]]; then
  awk -F',' 'NR==1{next} NF>=2{from=$1;to=$2; gsub(/\r/,"",from); gsub(/\r/,"",to);
    fn=toupper(from); gsub(/[^A-Z0-9]/,"",fn); sub(/INT$/,"I",fn);
    tn=toupper(to);   gsub(/[^A-Z0-9]/,"",tn); sub(/INT$/,"I",tn);
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
    awk -v t="$OVR_TO" 'BEGIN{FS=OFS="\t"} $5==t{print $1,$2,$3,$4,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="override"
  else
    awk -v t="$FAM_NORM" 'BEGIN{FS=OFS="\t"} $5==t{print $1,$2,$3,$4,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="norm_equal"
  fi

  if [[ ! -s "$OUT_BED" ]]; then
    BASE="$(printf '%s' "$FAM_NORM" | sed -E 's/I$//')"   # drop trailing I if present
    awk -v b="$BASE" 'BEGIN{FS=OFS="\t"} ($5 ~ "^" b){print $1,$2,$3,$4,$6,$7}' "$CAT_IDX" > "$OUT_BED"
    STRAT="base_prefix"
  fi

  if [[ -s "$OUT_BED" ]]; then
    N=$(wc -l < "$OUT_BED")
    EX=$(head -n1 "$OUT_BED" | awk '{print $4}')
    echo -e "${FAM}\t${STRAT}\t${EX}\t${N}" >> "$REPORT"
    echo "  $FAM -> $OUT_BED ($N copies) [$STRAT]"
  else
    echo -e "${FAM}\tno_match\t.\t0" >> "$REPORT"
    echo "  $FAM -> no matches; check naming or add override" >&2
    rm -f "$OUT_BED"
  fi
done < "$TOP_LIST"

echo "[report] $(wc -l < "$REPORT") lines -> $REPORT"

if [[ -n "$VALIDATE_BED" && -f "$VALIDATE_BED" ]]; then
  if command -v bedtools >/dev/null 2>&1; then
    echo "[validate] using bedtools intersect vs $VALIDATE_BED"
    while IFS= read -r FAM; do
      BED="$OUT_DIR/beds/${FAM}.bed"
      [[ -s "$BED" ]] || continue
      SAMPLE="$OUT_DIR/validate/${FAM}.sample.bed"
      head -n "$SAMPLE_VALIDATE_N" "$BED" > "$SAMPLE" || true
      [[ -s "$SAMPLE" ]] || continue
      HIT=$(bedtools intersect -u -f 0.9 -r -a "$SAMPLE" -b "$VALIDATE_BED" | wc -l | tr -d ' ')
      TOT=$(wc -l < "$SAMPLE")
      printf "%s\t%d/%d matched (>=90%% recip overlap)\n" "$FAM" "$HIT" "$TOT" \
        >> "$OUT_DIR/validate/validation_summary.tsv"
    done < "$TOP_LIST"
    echo "[validate] summary -> $OUT_DIR/validate/validation_summary.tsv"
  else
    echo "[validate] bedtools not found; skipping coordinate overlap check"
  fi
fi

echo "[done] Families: $TOP_LIST"
echo "[done] BEDs in:  $OUT_DIR/beds/"
echo "[done] Report:   $REPORT"
