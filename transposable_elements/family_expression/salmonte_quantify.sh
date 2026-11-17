#!/usr/bin/env bash
# The recorded KCl family-quantification command, parameterized for the three cohorts.
set -euo pipefail
if [[ $# -ne 4 ]]; then echo "Usage: $0 1|6|24 fastq_directory SalmonTE_mouse_index output_directory" >&2; exit 2; fi
case "$1" in 1|6|24) ;; *) echo "Timepoint must be 1, 6 or 24" >&2; exit 2;; esac
for rep in A B C; do
  for condition in M K; do
    sample="${rep}${1}${condition}"
    salmon quant -q -i "$3" -l A -1 "$2/${sample}_1.fq.gz" -2 "$2/${sample}_2.fq.gz" -o "$4/$sample" -p 4
  done
done
