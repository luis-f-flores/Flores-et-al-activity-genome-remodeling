#!/usr/bin/env bash
# Run the GPL-licensed R counter with a strand-neutral GATC annotation.
# Site and 400-kb count outputs feed prepare_window_counts.py for Fig. S01G-N.
set -euo pipefail
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
bam_dir="${1:?Usage: count_gatc_batch.sh BAM_DIRECTORY GATC_GFF CHROM_SIZES OUTPUT_DIRECTORY}"
gff_file="${2:?Provide the strand-neutral mm10 GATC GFF}"
genome_size_file="${3:?Provide the mm10 chromosome-size table}"
output_dir="${4:?Provide a new count-output directory}"
if [[ -e "$output_dir" ]]; then
    echo "Output directory already exists: $output_dir" >&2
    exit 1
fi
# Require strand-neutral GATC annotations.
awk '!/^#/ && NF >= 8 {n++; if ($7 != ".") bad=1} END {exit (bad || !n)}' "$gff_file"
for hours in 6 24; do
    for rep in A B; do
        for sample in DBK DK DBM DM; do
            test -f "$bam_dir/${hours}${rep}-${sample}.bam"
        done
    done
done
mkdir -p "$output_dir/logs"
for hours in 6 24; do
    for rep in A B; do
        for sample in DBK DK DBM DM; do
            name="${hours}${rep}-${sample}"
            Rscript "$script_dir/gatc_counting_mm10.R" -b "$bam_dir/$name.bam" -f "$gff_file" -o "$output_dir" -n "$name" -B 400 --genome_size "$genome_size_file" --mapqual 5 --maxgap 2 --log "$output_dir/logs/$name.log"
        done
    done
done
