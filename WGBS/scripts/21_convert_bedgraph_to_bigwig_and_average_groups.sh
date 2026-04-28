#!/usr/bin/env bash

# ==============================================================================
# Script: 21_convert_bedgraph_to_bigwig_and_average_groups.sh
# Package: WGBS / methylation analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for WGBS DMC/DMR/DMG analysis, annotation, overlap analysis, methylation-expression integration, signal tracks, and OLS modeling.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# 21_convert_bedgraph_to_bigwig_and_average_groups.sh
#

# Convert WGBS bedGraph files to bigWig and average replicate bigWigs by group.
#

# Requirements:
#   bedGraphToBigWig or wigToBigWig from UCSC utilities
#   wiggletools
#   sort
#

# Usage:
#   bash 21_convert_bedgraph_to_bigwig_and_average_groups.sh \
#     --input-dir methylKit_bedgraph_bigwig \
#     --chrom-sizes mm10.chrom.sizes \
#     --groups groups.tsv \
#     --outdir methylKit_bigwig_group_means
#
# groups.tsv format:
#   group    sample_id
#   0h       A0
#   0h       B0
#   0h       C0
#   1K       A1K
#   ...
#
# Expected bedGraph names:
#   <sample_id>_normalized_cov4.bedGraph

set -euo pipefail

INPUT_DIR=""
CHROM_SIZES=""
GROUPS_TSV=""
OUTDIR="methylKit_bigwig_group_means"
COV="4"

usage() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --input-dir) INPUT_DIR="$2"; shift 2;;
    --chrom-sizes) CHROM_SIZES="$2"; shift 2;;
    --groups) GROUPS_TSV="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --cov) COV="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument: $1" >&2; usage;;
  esac
done

[[ -d "$INPUT_DIR" ]] || { echo "ERROR: --input-dir not found" >&2; exit 1; }
[[ -f "$CHROM_SIZES" ]] || { echo "ERROR: --chrom-sizes not found" >&2; exit 1; }
[[ -f "$GROUPS_TSV" ]] || { echo "ERROR: --groups TSV not found" >&2; exit 1; }

mkdir -p "$OUTDIR/bigwigs" "$OUTDIR/group_means" "$OUTDIR/tmp"

if command -v bedGraphToBigWig >/dev/null 2>&1; then
  BG2BW="bedGraphToBigWig"
else
  echo "ERROR: bedGraphToBigWig not found. Install UCSC utilities." >&2
  exit 1
fi

command -v wiggletools >/dev/null 2>&1 || { echo "ERROR: wiggletools not found" >&2; exit 1; }
command -v wigToBigWig >/dev/null 2>&1 || { echo "ERROR: wigToBigWig not found" >&2; exit 1; }

# Convert individual bedGraphs to bigWigs.
awk 'BEGIN{FS=OFS="\t"} NR==1 && tolower($1)=="group"{next} NF>=2{print $2}' "$GROUPS_TSV" | sort -u |
while read -r SAMPLE; do
  BG="$INPUT_DIR/${SAMPLE}_normalized_cov${COV}.bedGraph"
  SORTED="$OUTDIR/tmp/${SAMPLE}.sorted.bedGraph"
  BW="$OUTDIR/bigwigs/${SAMPLE}.bw"

  [[ -f "$BG" ]] || { echo "ERROR: missing bedGraph for sample $SAMPLE: $BG" >&2; exit 1; }

  sort -k1,1 -k2,2n "$BG" > "$SORTED"
  "$BG2BW" "$SORTED" "$CHROM_SIZES" "$BW"
  echo "[bigWig] $BW"
done

# Average replicate bigWigs by group.
cut -f1 "$GROUPS_TSV" | tail -n +2 | sort -u |
while read -r GROUP; do
  mapfile -t SAMPLES < <(awk -v g="$GROUP" 'BEGIN{FS=OFS="\t"} NR==1{next} $1==g{print $2}' "$GROUPS_TSV")
  BWS=()
  for S in "${SAMPLES[@]}"; do
    BWS+=("$OUTDIR/bigwigs/${S}.bw")
  done

  WIG="$OUTDIR/tmp/${GROUP}_mean.wig"
  OUTBW="$OUTDIR/group_means/${GROUP}_mean.bw"

  wiggletools mean "${BWS[@]}" > "$WIG"
  wigToBigWig "$WIG" "$CHROM_SIZES" "$OUTBW"
  rm -f "$WIG"
  echo "[mean] $OUTBW"
done

rm -rf "$OUTDIR/tmp"
echo "[done] BigWig conversion and group averaging complete."
