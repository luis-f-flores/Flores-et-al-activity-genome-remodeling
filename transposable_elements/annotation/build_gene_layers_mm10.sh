#!/usr/bin/env bash

set -euo pipefail

BASE="${1:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/rna}"
GTF="${2:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/references/mm10/annot/mm10.refGene.gtf.gz}"
CHROMS="${3:-${GB_WORK_ROOT:?Set GB_WORK_ROOT}/references/mm10/chrom/mm10.chrom.sizes}"
OUT="$BASE/gene_layers_mm10"
TMP="$OUT/.tmp"

mkdir -p "$OUT" "$TMP"

command -v gawk >/dev/null || { echo "ERROR: GNU awk (gawk) is required." >&2; exit 1; }
for b in genes.bed exons.bed introns.bed tss.bed promoters_-1kb_+100bp.bed; do
  [[ ! -e "$OUT/$b" ]] || { echo "ERROR: Output already exists: $OUT/$b. Choose an empty output directory." >&2; exit 1; }
done

[[ -s "$GTF" ]]   || { echo "ERROR: GTF not found: $GTF" >&2; exit 1; }
[[ -s "$CHROMS" ]]|| { echo "ERROR: chrom.sizes not found: $CHROMS" >&2; exit 1; }

echo ">> Building gene layers from:"
echo "   GTF    : $GTF"
echo "   CHROMS : $CHROMS"
echo "   OUT    : $OUT"

echo ">> Extracting and sorting exon records…"
gzip -dc "$GTF" \
| gawk -F'\t' -v OFS='\t' '
  BEGIN{ }
  $0 ~ /^#/ { next }
  $3 != "exon" { next }
  {
    chr=$1; start=$4-1; end=$5; strand=$7; attr=$9
    gname=""; n=split(attr,a,";")
    for(i=1;i<=n;i++){
      gsub(/^ +| +$/,"",a[i]); if(a[i]=="") continue
      split(a[i],kv," "); k=kv[1]; v=kv[2]; gsub(/"/,"",v)
      if(k=="gene_name" && gname=="") gname=v
      if(k=="gene"      && gname=="") gname=v
      if(k=="gene_id"   && gname=="") gname=v
      if(k=="transcript_id" && gname=="") gname=v
    }
    if(gname=="") next
    print chr, start, end, gname, 0, strand
  }
' \
| sort -k4,4 -k1,1 -k2,2n > "$TMP/exons.sorted.bed"

echo ">> Merging exons per gene, computing genes/TSS/promoters/introns…"
gawk -v OFS='\t' -v CHR="$CHROMS" '
  BEGIN{
    while ((getline < CHR) > 0){
      if($1=="") continue; chrlen[$1]=$2
    }
    close(CHR)
    prevG=""; prevChr=""; prevStrand="";
    ms=-1; me=-1; k=0; haveGene=0
  }
  function flush_gene(   i,chr,len,tss,pstart,pend,start0,end0){
    if(!haveGene) return

    for(i=1;i<=k;i++){
      print prevChr, ExS[i], ExE[i], prevG, 0, prevStrand >> exons
    }

    print prevChr, gStart, gEnd, prevG, 0, prevStrand >> genes

    if(prevStrand=="+"){ tss = minTss } else { tss = maxTss }
    print prevChr, tss, tss+1, prevG, 0, prevStrand >> tssbed

    len = chrlen[prevChr]
    if(prevStrand=="+"){ pstart=tss-1000; pend=tss+100 }
    else               { pstart=tss-100;  pend=tss+1000 }
    if(pstart < 0) pstart=0
    if(len>0 && pend > len) pend=len
    if(pend > pstart) print prevChr, pstart, pend, prevG, 0, prevStrand >> prom

    for(i=1;i<k;i++){
      start0 = ExE[i]
      end0   = ExS[i+1]
      if(end0 > start0){
        print prevChr, start0, end0, prevG, 0, prevStrand >> intr
      }
    }

    k=0; ms=-1; me=-1; haveGene=0
  }
  function start_gene(chr,g,strand, s,e){
    prevChr=chr; prevG=g; prevStrand=strand
    gStart=s; gEnd=e
    minTss=(strand=="+")? s : e
    maxTss=minTss
    k=1; ExS[1]=s; ExE[1]=e
    ms=s; me=e
    haveGene=1
  }
  function append_exon(s,e){
    if(s<gStart) gStart=s
    if(e>gEnd)   gEnd=e
    if(prevStrand=="+"){ if(s<minTss) minTss=s }
    else               { if(e>maxTss) maxTss=e }
    if(s <= me){ if(e>me) me=e; ExE[k]=me }
    else{ k++; ExS[k]=s; ExE[k]=e; ms=s; me=e }
  }

  BEGINFILE{
    exons="'"$OUT"'/exons.bed"
    genes="'"$OUT"'/genes.bed"
    tssbed="'"$OUT"'/tss.bed"
    prom="'"$OUT"'/promoters_-1kb_+100bp.bed"
    intr="'"$OUT"'/introns.bed"
  }

  {
    chr=$1; s=$2; e=$3; g=$4; strand=$6
    if(g!=prevG || chr!=prevChr || strand!=prevStrand){
      flush_gene()
      start_gene(chr,g,strand,s,e)
    } else {
      append_exon(s,e)
    }
  }
  END{
    flush_gene()
  }
' "$TMP/exons.sorted.bed"

echo ">> Sorting outputs…"
for b in genes.bed exons.bed introns.bed tss.bed promoters_-1kb_+100bp.bed; do
  sort -k1,1 -k2,2n "$OUT/$b" > "$TMP/$b.sorted"
  mv "$TMP/$b.sorted" "$OUT/$b"
done

echo ">> Done."
for b in genes.bed exons.bed introns.bed tss.bed promoters_-1kb_+100bp.bed; do
  n=$(wc -l < "$OUT/$b")
  echo "   $(printf '%-28s' "$b") rows: $n"
done

rm -rf "$TMP"
