#!/usr/bin/env bash
# Align trimmed Read 1 reads using the manuscript Bismark settings.
set -euo pipefail
if [[ $# -ne 3 ]]; then echo "Usage: $0 genome_directory trimmed_read1.fastq.gz output_directory" >&2; exit 2; fi
mkdir -p "$3"
bismark --non_directional --hisat2 --local --genome "$1" --output_dir "$3" "$2"
