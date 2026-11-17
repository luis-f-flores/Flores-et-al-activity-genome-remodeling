#!/usr/bin/env bash
set -euo pipefail

[[ $# -eq 3 ]] || { echo "Usage: build_te_catalog_mm10.sh rmsk.txt.gz salmonte_mm.fa output_directory"; exit 1; }
RMSK="$1"
SALMONTE_REF="$2"
CATALOG="$3"
for name in TE_mm10.bed TE_meta.tsv; do
  [[ ! -e "$CATALOG/$name" ]] || { echo "Output already exists: $CATALOG/$name"; exit 1; }
done

[[ -f "$RMSK" ]] || { echo "ERROR: RepeatMasker table not found: $RMSK"; exit 1; }
[[ -f "$SALMONTE_REF" ]] || { echo "ERROR: SalmonTE mm.fa not found: $SALMONTE_REF"; exit 1; }
mkdir -p "$CATALOG"

echo "==> RMSK:        $RMSK"
echo "==> SalmonTE:    $SALMONTE_REF"
echo "==> Output dir:  $CATALOG"

echo "[1/4] Extracting SalmonTE family list ..."
grep -E '^>' "$SALMONTE_REF" | sed 's/^>//; s/[[:space:]].*$//' | gawk '{print toupper($1)}' | sort -u > "$CATALOG/SalmonTE_mm_families.txt"
echo "  families: $(wc -l < "$CATALOG/SalmonTE_mm_families.txt")"

alias_awk='
  function fam_key(rep,   repU,tmp,k1,k2,k2b,k0,kPlusI,kMinusI,baseAlpha,baseNum) {
    repU = toupper(rep)

if (keep["IAPLTR1A_I_MM"] && repU ~ /^IAPLTR1A(_MM)?$/) return "IAPLTR1A_I_MM";
if (keep["IAPLTR1A_I_MM"] && repU == "IAPLTR1A_I")      return "IAPLTR1A_I_MM";

if (keep["LTRIAPEY"] && repU=="IAPEY_LTR") return "LTRIAPEY";

if (keep["RLTR1IAP_MM"] && repU ~ /^RLTR1[BCDEF]_MM$/) return "RLTR1IAP_MM";

if (keep["MERVL_LTR"] && repU=="MERVL") return "MERVL_LTR";

    k1   = repU                          # exact
    k2   = repU; sub(/-INT$/, "_I", k2)  # ...-INT -> _I (IAPEZ-INT -> IAPEZ_I)
    k2b  = repU; sub(/-INT$/, "I", k2b)  # ...-INT -> I  (IAPEZ-INT -> IAPEZI)
    k0   = repU; sub(/-INT$/, "", k0)    # ...-INT -> "" (MERVL-INT -> MERVL)

    if (keep[k1])            return k1
    if (keep[k2])            return k2
    if (keep[k2b])           return k2b
    if (keep[k0])            return k0
    if (keep[k0 "_MM"])      return k0 "_MM"

    kPlusI  = repU
    if (kPlusI ~ /^IAPLTR/ && kPlusI !~ /_I(_MM)?$/) sub(/_MM$/, "_I_MM", kPlusI)
    if (keep[kPlusI])        return kPlusI

    kMinusI = repU
    if (kMinusI ~ /^IAPLTR/ && kMinusI ~ /_I(_MM)?$/) sub(/_I(_MM)?$/, "_MM", kMinusI)
    if (keep[kMinusI])       return kMinusI

    if (repU ~ /^B1F/ && keep["B1F"]) return "B1F"
    if (repU ~ /^B1/  && keep["B1"])  return "B1"
    if (repU ~ /^B2/  && keep["B2"])  return "B2"

    if (MAP_RLTR6BC_TO_6I && repU ~ /^RLTR6[BC](_MM)?$/ && keep["RLTR6I_MM"]) return "RLTR6I_MM"

    tmp = repU
    sub(/_MM$/, "", tmp)
    baseAlpha = tmp
    if (match(tmp, /^(RLTR[0-9]+[A-Z]+)/, a)) baseAlpha = a[1]   # RLTR13B4 -> RLTR13B ; RLTR12A -> RLTR12A
    baseNum = tmp
    if (match(tmp, /^(RLTR[0-9]+)/, b)) baseNum = b[1]           # RLTR20C1 -> RLTR20
    if (ALLOW_UMBRELLA) {
      if (keep[baseAlpha])         return baseAlpha
      if (keep[baseAlpha "_MM"])   return baseAlpha "_MM"
      if (keep[baseNum])           return baseNum
      if (keep[baseNum "_MM"])     return baseNum "_MM"
    }

    return ""
  }
'

echo "[2/4] Building TE_mm10.bed ..."
zcat "$RMSK" | gawk -v OFS="\t" \
    -v ALLOW_UMBRELLA="${ALLOW_UMBRELLA:-1}" \
    -v MAP_RLTR6BC_TO_6I="${MAP_RLTR6BC_TO_6I:-0}" '
  BEGIN{FS="\t"} NR==FNR {keep[$1]=1; next}
  '"$alias_awk"'
  {
    fam = fam_key($11)
    if (fam != "") print $6,$7,$8,fam"|" $17,0,$10
  }
' "$CATALOG/SalmonTE_mm_families.txt" - \
| LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$CATALOG/TE_mm10.bed"

echo "[3/4] Building TE_midpoints.bed ..."
gawk -v OFS="\t" '{mid=int(($2+$3)/2); print $1,mid,mid+1,$4,0,$6}' "$CATALOG/TE_mm10.bed" \
| LC_ALL=C sort -k1,1V -k2,2n -k3,3n > "$CATALOG/TE_midpoints.bed"

echo "[4/4] Writing TE_meta.tsv and TE_family_summary.tsv ..."
{
  echo -e "copyID\tfamily\trepClass\trepFamily\tpercDiv\tlength\tstrand\tchr\tstart\tend\trepStart\trepEnd"
  zcat "$RMSK" | gawk -v OFS="\t" \
      -v ALLOW_UMBRELLA="${ALLOW_UMBRELLA:-1}" \
      -v MAP_RLTR6BC_TO_6I="${MAP_RLTR6BC_TO_6I:-0}" '
    BEGIN{FS="\t"} NR==FNR {keep[$1]=1; next}
    '"$alias_awk"'
    {
      fam = fam_key($11)
      if (fam != "") {
        len=$8-$7; pd=$3/10.0
        printf "%s|%s\t%s\t%s\t%s\t%.3f\t%d\t%s\t%s\t%d\t%d\t%d\t%d\n",
               fam,$17,
               fam,
               $12,
               $13,
               pd,
               len,
               $10,
               $6,
               $7,
               $8,
               $14,
               $15
      }
    }
  ' "$CATALOG/SalmonTE_mm_families.txt" - 
} > "$CATALOG/TE_meta.tsv"

gawk -F"\t" 'NR>1{
  fam=$2; cls=$3; grp=$4; pd=$5+0; len=$6+0;
  n[fam]++; bp[fam]+=len; spd[fam]+=pd; slen[fam]+=len; C[fam]=cls; G[fam]=grp
}
END{
  OFS="\t"; print "family","repClass","repFamily","n_copies","total_bp","mean_percDiv","mean_length_bp";
  for(f in n) printf "%s\t%s\t%s\t%d\t%d\t%.3f\t%.1f\n", f,C[f],G[f],n[f],bp[f],(spd[f]/n[f]),(slen[f]/n[f]);
}' "$CATALOG/TE_meta.tsv" \
| LC_ALL=C sort -k4,4nr > "$CATALOG/TE_family_summary.tsv"

cut -f1 "$CATALOG/SalmonTE_mm_families.txt" | sort -u > "$CATALOG/.salmonte_fams.tmp"
tail -n +2 "$CATALOG/TE_meta.tsv" | cut -f2 | sort -u > "$CATALOG/.matched_fams.tmp"
comm -23 "$CATALOG/.salmonte_fams.tmp" "$CATALOG/.matched_fams.tmp" > "$CATALOG/unmatched_families.txt"
rm -f "$CATALOG/.salmonte_fams.tmp" "$CATALOG/.matched_fams.tmp"
echo "Unmatched SalmonTE families: $(wc -l < "$CATALOG/unmatched_families.txt") (see $CATALOG/unmatched_families.txt)"

echo "Catalog written to $CATALOG"
bash "$(dirname "$0")/catalog_copy_identifiers.sh" "$CATALOG"
