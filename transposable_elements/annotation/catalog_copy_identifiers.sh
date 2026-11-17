#!/usr/bin/env bash
set -euo pipefail
# TE catalogue identifiers use family@chr:start-end.
TE="${1:?Supply the generated TE catalog directory}"
test -f "$TE/TE_mm10.bed"
test -f "$TE/TE_meta.tsv"
test ! -e "$TE/TE_mm10.bed.bak"
test ! -e "$TE/TE_meta.tsv.bak"
cp "$TE/TE_mm10.bed" "$TE/TE_mm10.bed.bak"
cp "$TE/TE_meta.tsv" "$TE/TE_meta.tsv.bak"
awk -v OFS='	' '{split($4,a,"|"); $4=a[1] "@" $1 ":" $2 "-" $3; print}' "$TE/TE_mm10.bed.bak" > "$TE/TE_mm10.bed"
awk -v OFS='	' 'NR==1{print; next}{$1=$2 "@" $8 ":" $9 "-" $10; print}' "$TE/TE_meta.tsv.bak" > "$TE/TE_meta.tsv"
awk -v OFS='	' '{mid=int(($2+$3)/2); print $1,mid,mid+1,$4,0,$6}' "$TE/TE_mm10.bed" > "$TE/TE_midpoints.bed"
for name in TE_mm10.bed TE_midpoints.bed; do
  LC_ALL=C sort -k1,1V -k2,2n -k3,3n "$TE/$name" -o "$TE/$name"
done
