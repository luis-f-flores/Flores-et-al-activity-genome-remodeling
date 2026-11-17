#!/usr/bin/env python3

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)

from pathlib import Path
import gzip, hashlib, json, re
from datetime import datetime
from zoneinfo import ZoneInfo
import numpy as np
import pandas as pd
import argparse
ap=argparse.ArgumentParser(description='Match methylated gene-body loci to RNA records by name, chromosome, strand and overlapping coordinates.')
ap.add_argument('--rna-gtf',type=Path,required=True)
ap.add_argument('--gene-bed',type=Path,required=True)
ap.add_argument('--gene-body',type=Path,nargs=3,required=True)
ap.add_argument('--rna',type=Path,nargs=3,required=True)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();TABLE=args.output;TABLE.mkdir(parents=True,exist_ok=True)
W=GB_WORK_ROOT/'wgbs_package';TP=['1hr','6hr','24hr']
inputs=set();outputs=[];rows=[];excluded=[];ledger=[];reference_checks=[]
def read(p,**kw):
 inputs.add(p);return pd.read_csv(p,**kw)
def save_table(df,name):
 p=TABLE/name;df.to_csv(p,index=False);outputs.append(p)
def chrom(ch):
 ch=str(ch).removeprefix('chr')
 if ch=='X_GL456233_random':return 'GL456233.1'
 return ch
def coordinate_ok(t):
 return (t.chr.map(chrom)==t.gene_chr.map(chrom))&(t.strand==t.gene_strand)&(t.start<=t.gene_end)&(t.end>=t.gene_start)
rna_gtf=args.rna_gtf
inputs.add(rna_gtf)
attr_re=re.compile(r'(\w+) "([^"]+)"')
gtf_genes={}
with gzip.open(rna_gtf,'rt') as stream:
 for line in stream:
  if line.startswith('#'):continue
  p=line.rstrip('\n').split('\t')
  if p[2]!='gene':continue
  at=dict(attr_re.findall(p[8]));gid=at['gene_id']
  assert gid not in gtf_genes,gid
  gtf_genes[gid]=dict(gene_id_rna=gid,gene_name=at.get('gene_name'),gene_chr=p[0],gene_start=int(p[3]),gene_end=int(p[4]),gene_strand=p[6])
bed_path=args.gene_bed;inputs.add(bed_path)
bed=pd.read_csv(bed_path,sep='\t',header=None,names=['chr','start0','end','gene_symbol','score','strand'])
bed_keys=set(zip(bed.chr,bed.start0+1,bed.end,bed.strand,bed.gene_symbol))
for i,tp in enumerate(TP):
 n=tp.replace('hr','')
 diff=read(args.gene_body[i])
 de=read(args.rna[i],sep='\t')
 assert de.gene_id.is_unique
 sig=diff[(diff['meth.diff'].abs()>=5)&(diff.qvalue<=.05)].copy()
 sig['body_locus_id']=sig.chr+':'+sig.start.astype(str)+'-'+sig.end.astype(str)+':'+sig.strand
 assert sig.body_locus_id.is_unique
 assert all(k in bed_keys for k in zip(sig.chr,sig.start,sig.end,sig.strand,sig.gene_symbol))
 cols=['gene_id','gene_name','gene_chr','gene_start','gene_end','gene_strand','log2FoldChange','padj']
 legacy=sig.merge(de[cols],left_on='gene_symbol',right_on='gene_name',suffixes=('_meth','_rna'))
 for row in legacy.drop_duplicates('gene_id_rna').itertuples():
  assert row.gene_id_rna in gtf_genes,row.gene_id_rna
  g=gtf_genes[row.gene_id_rna]
  assert (row.gene_name,chrom(row.gene_chr),row.gene_start,row.gene_end,row.gene_strand)==(g['gene_name'],chrom(g['gene_chr']),g['gene_start'],g['gene_end'],g['gene_strand']),row.gene_id_rna
  reference_checks.append(dict(timepoint=tp,**g,delivered_matches_original_GTF=True))
 legacy['same_strand_coordinate_overlap']=coordinate_ok(legacy)
 bad=legacy[~legacy.same_strand_coordinate_overlap].copy()
 bad['reason']=np.where(bad.chr.map(chrom)!=bad.gene_chr.map(chrom),'same_name_different_chromosome',np.where(bad.strand!=bad.gene_strand,'same_name_opposite_strand','same_name_nonoverlapping_coordinates'))
 bad['valid_same_name_alternative_retained']=bad.body_locus_id.isin(legacy.loc[legacy.same_strand_coordinate_overlap,'body_locus_id'])
 excluded.append(bad)
 
 mapping=legacy[legacy.same_strand_coordinate_overlap][['body_locus_id','gene_id_rna']]
 assert mapping.body_locus_id.is_unique and mapping.gene_id_rna.is_unique
 fixed=sig.merge(mapping,on='body_locus_id',validate='one_to_one').merge(de[cols],left_on='gene_id_rna',right_on='gene_id',validate='many_to_one',suffixes=('_meth','_rna_source'))
 fixed['same_strand_coordinate_overlap']=coordinate_ok(fixed)
 assert fixed.same_strand_coordinate_overlap.all()
 fixed['join_scope']='same name, equivalent chromosome, strand and overlapping span; RNA joined by unique stable ID'
 fixed['timepoint']=tp
 rows.append(fixed)
 account=sig[['timepoint','body_locus_id','gene_symbol','chr','start','end','strand']].copy()
 account['RNA_match_status']=np.where(account.body_locus_id.isin(mapping.body_locus_id),'retained',np.where(account.body_locus_id.isin(legacy.body_locus_id),'excluded_no_supported_same_name_locus','no_delivered_RNA_record_with_same_name'))
 ledger.append(account)
matched=pd.concat(rows,ignore_index=True);crosswalk=matched[['body_locus_id','gene_id_rna']].drop_duplicates()
assert crosswalk.body_locus_id.is_unique and crosswalk.gene_id_rna.is_unique
rejected=pd.concat(excluded,ignore_index=True)
save_table(matched,'64_displayed_gene_locus_RNA_matches.csv')
save_table(crosswalk,'64_unique_locus_to_RNA_ID_crosswalk.csv')
save_table(rejected,'64_rejected_wrong_locus_matches.csv')
save_table(pd.concat(ledger,ignore_index=True),'64_all_DMG_match_accounting.csv')
save_table(pd.DataFrame(reference_checks),'64_RNA_original_GTF_verification.csv')
print('Matched gene bodies:',matched.groupby('timepoint').size().to_dict())
