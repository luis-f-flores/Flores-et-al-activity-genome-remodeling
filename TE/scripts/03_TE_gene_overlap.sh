#!/usr/bin/env bash

# ==============================================================================
# Script: 03_TE_gene_overlap.sh
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

# Intersect TE copies with DEG gene-feature layers and summarize overlaps.

# Usage:
#   bash 03_TE_gene_overlap.sh --gene-layers gene_layers_mm10 --te-catalog TE_catalog_mm10 --deg-dir deg_tables --outdir overlaps_mm10 --pattern 'Hour_*_KClvsHour_*_Mann_deg.xls'
set -euo pipefail
LAY=""; TE=""; DEG_DIR=""; OUTDIR=""; PATTERN="Hour_*_KClvsHour_*_Mann_deg.xls"
usage(){ echo "Usage: $0 --gene-layers <dir> --te-catalog <dir> --deg-dir <dir> --outdir <dir> [--pattern pattern]"; exit 1; }
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene-layers) LAY="$2"; shift 2;;
    --te-catalog) TE="$2"; shift 2;;
    --deg-dir) DEG_DIR="$2"; shift 2;;
    --outdir) OUTDIR="$2"; shift 2;;
    --pattern) PATTERN="$2"; shift 2;;
    -h|--help) usage;;
    *) echo "ERROR: unknown argument $1" >&2; usage;;
  esac
done
for f in "$LAY/genes.bed" "$LAY/exons.bed" "$LAY/introns.bed" "$LAY/tss.bed" "$LAY/promoters_-1kb_+100bp.bed" "$TE/TE_mm10.bed" "$TE/TE_midpoints.bed" "$TE/TE_meta.tsv"; do
  [[ -s "$f" ]] || { echo "ERROR: missing required file: $f" >&2; exit 1; }
done
[[ -d "$DEG_DIR" ]] || { echo "ERROR: DEG directory not found" >&2; exit 1; }
[[ -n "$OUTDIR" ]] || usage
mkdir -p "$OUTDIR"
command -v bedtools >/dev/null || { echo "ERROR: bedtools not found" >&2; exit 1; }
python3 - <<'PYCHK' >/dev/null || { echo "ERROR: Python packages pandas/openpyxl missing" >&2; exit 1; }
import pandas, openpyxl
PYCHK

filter_genes_exact(){ local list="$1" inbed="$2" outbed="$3"; awk 'NR==FNR{a[$1]=1;next} ($4 in a)' "$list" "$inbed" > "$outbed"; }
add_stats(){ local stats="$1" inbed="$2" outbed="$3"; awk -v OFS='\t' -v stats="$stats" 'BEGIN{while((getline s<stats)>0){split(s,a,"\t"); if(a[1]!=""){L2[a[1]]=a[2]; P[a[1]]=a[3]}}} {g=$4; print $1,$2,$3,$4,$5,$6,(g in L2?L2[g]:"NA"),(g in P?P[g]:"NA")}' "$inbed" > "$outbed"; }

shopt -s nullglob
DEG_FILES=("$DEG_DIR"/$PATTERN)
[[ ${#DEG_FILES[@]} -gt 0 ]] || { echo "ERROR: no DEG files matched" >&2; exit 1; }

for DEG_TABLE in "${DEG_FILES[@]}"; do
  LABEL=$(basename "$DEG_TABLE" .xls)
  OUT="$OUTDIR/$LABEL"
  mkdir -p "$OUT"
  echo ">>> $LABEL"

  DEG_STATS="$OUT/deg_stats.tsv"
  awk -v OFS='\t' 'NR==1{for(i=1;i<=NF;i++){gsub(/\r/,"",$i); H[$i]=i}; if(!("gene_name" in H && "log2FoldChange" in H && "padj" in H)){print "ERROR: DEG table missing gene_name/log2FoldChange/padj" > "/dev/stderr"; exit 2}; next} {g=$H["gene_name"]; if(g!="") print g,$H["log2FoldChange"],$H["padj"]}' FS='\t' "$DEG_TABLE" | sort -u > "$DEG_STATS"
  cut -f1 "$DEG_STATS" > "$OUT/deg_genes.txt"

  filter_genes_exact "$OUT/deg_genes.txt" "$LAY/genes.bed" "$OUT/deg_genes.bed"
  filter_genes_exact "$OUT/deg_genes.txt" "$LAY/exons.bed" "$OUT/deg_exons.bed"
  filter_genes_exact "$OUT/deg_genes.txt" "$LAY/introns.bed" "$OUT/deg_introns.bed"
  filter_genes_exact "$OUT/deg_genes.txt" "$LAY/tss.bed" "$OUT/deg_tss.bed"
  filter_genes_exact "$OUT/deg_genes.txt" "$LAY/promoters_-1kb_+100bp.bed" "$OUT/deg_promoters.bed"

  add_stats "$DEG_STATS" "$OUT/deg_genes.bed" "$OUT/deg_genes_with_stats.bed"
  add_stats "$DEG_STATS" "$OUT/deg_exons.bed" "$OUT/deg_exons_with_stats.bed"
  add_stats "$DEG_STATS" "$OUT/deg_introns.bed" "$OUT/deg_introns_with_stats.bed"
  add_stats "$DEG_STATS" "$OUT/deg_tss.bed" "$OUT/deg_tss_with_stats.bed"
  add_stats "$DEG_STATS" "$OUT/deg_promoters.bed" "$OUT/deg_promoters_with_stats.bed"
  rm -f "$OUT"/deg_{genes,exons,introns,tss,promoters}.bed
  for bed in "$OUT"/*_with_stats.bed; do LC_ALL=C sort -k1,1V -k2,2n -k3,3n "$bed" -o "$bed"; done

  bedtools intersect -wa -wb -a "$OUT/deg_genes_with_stats.bed" -b "$TE/TE_mm10.bed" > "$OUT/degGenes__overlap_TEcopies.tsv" || true
  bedtools intersect -wa -wb -a "$TE/TE_mm10.bed" -b "$OUT/deg_promoters_with_stats.bed" > "$OUT/TE_in_promoters_deg.tsv" || true
  bedtools intersect -wa -wb -a "$TE/TE_mm10.bed" -b "$OUT/deg_exons_with_stats.bed" > "$OUT/TE_in_exons_deg.tsv" || true
  bedtools intersect -wa -wb -a "$TE/TE_mm10.bed" -b "$OUT/deg_introns_with_stats.bed" > "$OUT/TE_in_introns_deg.tsv" || true
  bedtools closest -s -D ref -sorted -a "$TE/TE_midpoints.bed" -b "$OUT/deg_tss_with_stats.bed" > "$OUT/TEmidpoint__to_nearestTSS_deg.tsv" || true

  python3 - "$OUT" "$TE" <<'PY'
import os, sys, pandas as pd
out, te = sys.argv[1], sys.argv[2]
meta = pd.read_csv(os.path.join(te,'TE_meta.tsv'), sep='\t', usecols=['copyID','family','repClass','repFamily']).drop_duplicates()
def read(path, cols):
    return pd.DataFrame(columns=cols) if (not os.path.exists(path) or os.path.getsize(path)==0) else pd.read_csv(path, sep='\t', header=None, names=cols)
gb_cols=['g_chr','g_start','g_end','gene','g_score','g_strand','log2FC','padj','te_chr','te_start','te_end','copyID','te_score','te_strand']
feat_cols=['te_chr','te_start','te_end','copyID','te_score','te_strand','feat_chr','feat_start','feat_end','gene','feat_score','feat_strand','log2FC','padj']
gb=read(os.path.join(out,'degGenes__overlap_TEcopies.tsv'),gb_cols).merge(meta,on='copyID',how='left')
prom=read(os.path.join(out,'TE_in_promoters_deg.tsv'),feat_cols).merge(meta,on='copyID',how='left')
exon=read(os.path.join(out,'TE_in_exons_deg.tsv'),feat_cols).merge(meta,on='copyID',how='left')
intron=read(os.path.join(out,'TE_in_introns_deg.tsv'),feat_cols).merge(meta,on='copyID',how='left')
def famsum(df):
    if df.empty: return pd.DataFrame(columns=['family','repClass','repFamily','n_overlaps','n_unique_copies','n_genes','pct_overlaps'])
    total=len(df)
    s=df.groupby(['family','repClass','repFamily'],dropna=False).agg(n_overlaps=('copyID','size'), n_unique_copies=('copyID','nunique'), n_genes=('gene','nunique')).reset_index()
    s['pct_overlaps']=(s['n_overlaps']/total*100).round(2)
    return s.sort_values(['n_overlaps','n_unique_copies'],ascending=[False,False])
with pd.ExcelWriter(os.path.join(out,'TE_gene_overlap_master.xlsx')) as xw:
    for name, df in [('geneBodies',gb),('promoters',prom),('exons',exon),('introns',intron)]:
        df.to_excel(xw,index=False,sheet_name=f'overlaps_{name}')
        famsum(df).to_excel(xw,index=False,sheet_name=f'summary_family_{name}')
PY
  echo ">>> wrote $OUT/TE_gene_overlap_master.xlsx"
done
