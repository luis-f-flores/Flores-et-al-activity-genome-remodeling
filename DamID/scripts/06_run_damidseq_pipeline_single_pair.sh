#!/usr/bin/env bash

# ==============================================================================
# Script: 06_run_damidseq_pipeline_single_pair.sh
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

# 06_run_damidseq_pipeline_single_pair.sh
#
# Generic wrapper for damidseq_pipeline_vR.1.pl for a Dam-LaminB1 vs Dam-only pair.
#

# Usage:
#   bash 06_run_damidseq_pipeline_single_pair.sh \
#     --pipeline /path/to/damidseq_pipeline_vR.1.pl \
#     --gatc /path/to/mm10.GATC.gff \
#     --bowtie2-index /path/to/mm10_bowtie2_index/mm10 \
#     --samtools-dir /path/to/samtools/bin \
#     --bowtie2-dir /path/to/bowtie2/bin \
#     --dam-only Dam_only.fq.gz \
#     --experimental DamLaminB1.fq.gz \
#     --outdir output \
#     --bins 300
#
# Output:
#   damidseq_pipeline bedGraph/BAM outputs copied to --outdir.

set -euo pipefail

PIPELINE=""
GATC=""
BOWTIE2_INDEX=""
SAMTOOLS_DIR=""
BOWTIE2_DIR=""
DAM_ONLY=""
EXPERIMENTAL=""
OUTDIR="output"
BINS="300"

usage() {
  grep '^#' "$0" | sed 's/^# \{0,1\}//'
  exit 1
}

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pipeline) PIPELINE="$2"; shift 2;;
    --gatc) GATC="$2"; shift 2;;
    --bowtie2-index) BOWTIE2_INDEX="$2"; shift 2;;
    --samtools-dir) SAMTOOLS_DIR="$2"; shift 2;;
    --bowtie2-dir) BOWTIE2_DIR="$2"; shift 2;;
    --dam-only) DAM_ONLY="$2"; shift 2;;
    --experimental) EXPERIMENTAL="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --bins) BINS="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument: $1" >&2; usage;;
  esac
done

[[ -f "$PIPELINE" ]] || { echo "ERROR: --pipeline not found" >&2; exit 1; }
[[ -f "$GATC" ]] || { echo "ERROR: --gatc not found" >&2; exit 1; }
[[ -n "$BOWTIE2_INDEX" ]] || { echo "ERROR: --bowtie2-index required" >&2; exit 1; }
[[ -d "$SAMTOOLS_DIR" ]] || { echo "ERROR: --samtools-dir not found" >&2; exit 1; }
[[ -d "$BOWTIE2_DIR" ]] || { echo "ERROR: --bowtie2-dir not found" >&2; exit 1; }
[[ -f "$DAM_ONLY" ]] || { echo "ERROR: --dam-only FASTQ not found" >&2; exit 1; }
[[ -f "$EXPERIMENTAL" ]] || { echo "ERROR: --experimental FASTQ not found" >&2; exit 1; }

mkdir -p "$OUTDIR"

perl "$PIPELINE" \
  --bins="$BINS" \
  --gatc_frag_file="$GATC" \
  --bowtie2_genome_dir="$BOWTIE2_INDEX" \
  --samtools_path="$SAMTOOLS_DIR/" \
  --bowtie2_path="$BOWTIE2_DIR/" \
  --dam="$DAM_ONLY" \
  "$EXPERIMENTAL"

# Move common pipeline outputs if generated in current directory.
find . -maxdepth 1 \( -name "*.bedgraph" -o -name "*.bedGraph" -o -name "*.bam" \) -exec mv -f {} "$OUTDIR"/ \;

echo "[done] damidseq pipeline outputs moved to: $OUTDIR"
