#!/usr/bin/env bash

# ==============================================================================
# Script: 05_run_damMer_batch.sh
# Package: DamID / LAG analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for DamID preprocessing, DiffBind/GTF annotation, LAG/DEG overlap analysis, and GO enrichment.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

# 05_run_damMer_batch.sh
#
# Generic wrapper for damMer.py using Dam-LaminB1 experimental FASTQs and Dam-only control FASTQs.
#

# Usage:
#   bash 05_run_damMer_batch.sh \
#     --workdir /path/to/NanoDam_analysis \
#     --dammer /path/to/damMer.py \
#     --pipeline /path/to/damidseq_pipeline_vR.1.pl \
#     --bowtie2-index /path/to/mm10_bowtie2_index/mm10 \
#     --gatc /path/to/mm10.GATC.gff \
#     --bowtie2 /path/to/bowtie2 \
#     --samtools /path/to/samtools \
#     --experimental exp_rep1.fq.gz,exp_rep2.fq.gz \
#     --control dam_rep1.fq.gz,dam_rep2.fq.gz \
#     --email user@example.edu
#
# Notes:
#   - This is a wrapper around the published damMer/NanoDam workflow.
#   - The experimental files correspond to Dam-LaminB1 libraries.
#   - The control files correspond to Dam-only libraries.

set -euo pipefail

WORKDIR=""
DAMMER="damMer.py"
PIPELINE="damidseq_pipeline_vR.1.pl"
BOWTIE2_INDEX=""
GATC=""
BOWTIE2="bowtie2"
SAMTOOLS="samtools"
EXPERIMENTAL=""
CONTROL=""
EMAIL=""

usage() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --workdir) WORKDIR="$2"; shift 2;;
    --dammer) DAMMER="$2"; shift 2;;
    --pipeline) PIPELINE="$2"; shift 2;;
    --bowtie2-index) BOWTIE2_INDEX="$2"; shift 2;;
    --gatc) GATC="$2"; shift 2;;
    --bowtie2) BOWTIE2="$2"; shift 2;;
    --samtools) SAMTOOLS="$2"; shift 2;;
    --experimental) EXPERIMENTAL="$2"; shift 2;;
    --control) CONTROL="$2"; shift 2;;
    --email) EMAIL="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument: $1" >&2; usage;;
  esac
done

[[ -n "$WORKDIR" && -d "$WORKDIR" ]] || { echo "ERROR: --workdir required" >&2; exit 1; }
[[ -f "$DAMMER" || "$(command -v "$DAMMER" || true)" ]] || { echo "ERROR: damMer.py not found" >&2; exit 1; }
[[ -f "$PIPELINE" ]] || { echo "ERROR: damidseq pipeline not found" >&2; exit 1; }
[[ -n "$BOWTIE2_INDEX" ]] || { echo "ERROR: --bowtie2-index required" >&2; exit 1; }
[[ -f "$GATC" ]] || { echo "ERROR: --gatc GFF not found" >&2; exit 1; }
[[ -n "$EXPERIMENTAL" && -n "$CONTROL" ]] || { echo "ERROR: --experimental and --control required" >&2; exit 1; }

IFS=',' read -r -a EXP_FILES <<< "$EXPERIMENTAL"
IFS=',' read -r -a CTRL_FILES <<< "$CONTROL"

cd "$WORKDIR"

cmd=(python3 "$DAMMER" -e "${EXP_FILES[@]}" -c "${CTRL_FILES[@]}" \
  -i "$BOWTIE2_INDEX" \
  -g "$GATC" \
  -b "$BOWTIE2" \
  -s "$SAMTOOLS" \
  -q "$PIPELINE")

if [[ -n "$EMAIL" ]]; then
  cmd+=(-f "$EMAIL")
fi

printf '[run]'; printf ' %q' "${cmd[@]}"; printf '\n'
"${cmd[@]}"
