BASE="${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna"
CAT="$BASE/TE_catalog_mm10/TE_family_summary.tsv"
ALIAS="$BASE/TE_catalog_mm10/te_name_aliases.tsv"

rm -f "$BASE"/SalmonTE_DE_*/*TE_DE_results_annotated.tsv

cat > "$ALIAS" <<'EOF'
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
EOF

for de in "$BASE"/SalmonTE_DE_*/*TE_DE_results.csv; do
  out="${de%/*}/TE_DE_results_annotated.tsv"
  echo "Annotating: $de"
  awk -v OFS="\t" -v cat="$CAT" -v alias="$ALIAS" '
    BEGIN{
      FS=",";                             # DE files are CSV
      while ((getline l < cat) > 0) {
        if(l ~ /^family\t/) continue
        split(l,f,"\t")
        fam=f[1]
        rc[fam]=f[2]; rf[fam]=f[3]; nc[fam]=f[4]; mpd[fam]=f[6]; mlen[fam]=f[7]
      }
      close(cat)
      while ((getline a < alias) > 0) {
        if(a ~ /^from[ \t]/) continue
        split(a,aa,"\t")
        if (aa[1] != "" && aa[2] != "") al[aa[1]]=aa[2]
      }
      close(alias)

      print "family","baseMean","log2FoldChange","lfcSE","pvalue","padj",
            "repClass","repFamily","n_copies","mean_percDiv","mean_length_bp"
    }
    NR>1{
      gsub(/\r/,"")
      fam=$1; gsub(/^"|"$/,"",fam); gsub(/^[ \t]+|[ \t]+$/,"",fam)
      key=(fam in al ? al[fam] : fam)
      print fam,$2,$3,$4,$5,$6,
            (key in rc?rc[key]:"NA"),
            (key in rf?rf[key]:"NA"),
            (key in nc?nc[key]:"NA"),
            (key in mpd?mpd[key]:"NA"),
            (key in mlen?mlen[key]:"NA")
    }
  ' "$de" > "$out"
  echo "  -> $out"
done

for f in "$BASE"/SalmonTE_DE_*/*TE_DE_results_annotated.tsv; do
  printf "%s\tNA_rows=" "$(basename "$(dirname "$f")")"
  awk -F'\t' 'NR>1 && ($7=="NA" || $8=="NA"){c++} END{print (c+0)}' "$f"
done
