#!/usr/bin/env bash

# ==============================================================================
# Script: 01_build_te_catalog_mm10.sh
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

# Build a per-copy mm10 TE catalog with family labels matched to SalmonTE.

# Usage:
#   bash 01_build_te_catalog_mm10.sh --rmsk mm10_rmsk.txt.gz --salmonte-ref mm.fa --outdir TE_catalog_mm10
set -euo pipefail

usage(){ echo "Usage: $0 --rmsk <rmsk.txt.gz> --salmonte-ref <mm.fa> --outdir <TE_catalog_mm10> [--alias-map aliases.tsv]"; exit 1; }
RMSK=""; SALMONTE_REF=""; OUTDIR=""; ALIAS_MAP=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --rmsk) RMSK="$2"; shift 2;;
    --salmonte-ref) SALMONTE_REF="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --alias-map) ALIAS_MAP="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument $1" >&2; usage;;
  esac
done
[[ -f "$RMSK" ]] || { echo "ERROR: rmsk file not found" >&2; exit 1; }
[[ -f "$SALMONTE_REF" ]] || { echo "ERROR: SalmonTE reference not found" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || usage
mkdir -p "$OUTDIR"

grep -E '^>' "$SALMONTE_REF" | sed 's/^>//; s/[[:space:]].*$//' | awk '{print toupper($1)}' | sort -u > "$OUTDIR/SalmonTE_mm_families.txt"

if [[ -z "$ALIAS_MAP" ]]; then
  ALIAS_MAP="$OUTDIR/te_name_aliases.tsv"
  cat > "$ALIAS_MAP" <<'ALIASES'
from	to
IAPA_MM	IAPA-int
IAPEZI	IAPEZ-int
IAPEYI	IAPEY-int
IAPEY3_I	IAPEY3-int
IAPEY4_I	IAPEY4-int
IAPEY5_I	IAPEY5-int
MERVL_LTR	MT2_Mm
IAP1-MM_I	IAP1_MM-int
MMERVK10D3_I	MMERVK10D3-int
RLTR4I_MM	RLTR4_MM
RLTR6I_MM	RLTR6_MM
ERVB4_2B-LTR_MM	ERVB4_2-LTR_MM
RLTR1IAP_MM	IAPLTR1_Mm
ALIASES
fi

zcat "$RMSK" | awk -v OFS='\t' -v fams="$OUTDIR/SalmonTE_mm_families.txt" -v aliases="$ALIAS_MAP" '
BEGIN{
  FS="\t";
  while((getline f<fams)>0){keep[toupper(f)]=1}
  while((getline a<aliases)>0){if(a ~ /^from/) continue; split(a,x,"\t"); if(x[1] && x[2]) alias[toupper(x[1])]=x[2]}
}
function map_family(rep, key,tmp,base){
  key=toupper(rep)
  if(key in alias) return alias[key]
  if(key in keep) return key
  tmp=key; sub(/-INT$/,"_I",tmp); if(tmp in keep) return tmp
  tmp=key; sub(/-INT$/,"I",tmp);  if(tmp in keep) return tmp
  tmp=key; sub(/-INT$/,"",tmp);   if(tmp in keep) return tmp
  if((tmp "_MM") in keep) return tmp "_MM"
  if(key ~ /^B1/ && ("B1" in keep)) return "B1"
  if(key ~ /^B2/ && ("B2" in keep)) return "B2"
  base=key; sub(/_MM$/,"",base); if(base in keep) return base; if((base "_MM") in keep) return base "_MM"
  return ""
}
{
  fam=map_family($11)
  if(fam!="") print $6,$7,$8,fam"|"$17,0,$10
}' | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$OUTDIR/TE_mm10.bed"

awk -v OFS='\t' '{mid=int(($2+$3)/2); print $1,mid,mid+1,$4,0,$6}' "$OUTDIR/TE_mm10.bed" | LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$OUTDIR/TE_midpoints.bed"

{
  echo -e "copyID\tfamily\trepClass\trepFamily\tpercDiv\tlength\tstrand\tchr\tstart\tend"
  zcat "$RMSK" | awk -v OFS='\t' -v fams="$OUTDIR/SalmonTE_mm_families.txt" -v aliases="$ALIAS_MAP" '
  BEGIN{FS="\t"; while((getline f<fams)>0){keep[toupper(f)]=1}; while((getline a<aliases)>0){if(a ~ /^from/) continue; split(a,x,"\t"); if(x[1] && x[2]) alias[toupper(x[1])]=x[2]}}
  function map_family(rep, key,tmp,base){key=toupper(rep); if(key in alias)return alias[key]; if(key in keep)return key; tmp=key; sub(/-INT$/,"_I",tmp); if(tmp in keep)return tmp; tmp=key; sub(/-INT$/,"I",tmp); if(tmp in keep)return tmp; tmp=key; sub(/-INT$/,"",tmp); if(tmp in keep)return tmp; if((tmp "_MM") in keep)return tmp "_MM"; if(key ~ /^B1/ && ("B1" in keep))return "B1"; if(key ~ /^B2/ && ("B2" in keep))return "B2"; base=key; sub(/_MM$/,"",base); if(base in keep)return base; if((base "_MM") in keep)return base "_MM"; return ""}
  {fam=map_family($11); if(fam!=""){len=$8-$7; pd=$3/10.0; print fam"|"$17,fam,$12,$13,pd,len,$10,$6,$7,$8}}
  '
} > "$OUTDIR/TE_meta.tsv"

awk -F'\t' 'NR>1{fam=$2; n[fam]++; bp[fam]+=$6; div[fam]+=$5; cls[fam]=$3; grp[fam]=$4} END{OFS="\t"; print "family","repClass","repFamily","n_copies","total_bp","mean_percDiv","mean_length_bp"; for(f in n) print f,cls[f],grp[f],n[f],bp[f],div[f]/n[f],bp[f]/n[f]}' "$OUTDIR/TE_meta.tsv" | LC_ALL=C sort -k4,4nr > "$OUTDIR/TE_family_summary.tsv"

echo "Done. Outputs written to $OUTDIR"
