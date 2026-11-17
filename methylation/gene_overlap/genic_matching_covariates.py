"""Calculate genic matching covariates from reference sequence and provider FPKM."""
from pathlib import Path
import pandas as pd,numpy as np,hashlib,json
import argparse
ap=argparse.ArgumentParser(description='Measure gene-span CpG density and baseline RNA for genic matching.')
ap.add_argument('--inputs',type=Path,required=True)
ap.add_argument('--fasta',type=Path,required=True)
ap.add_argument('--master-fpkm',type=Path,required=True)
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();TAB=args.inputs;OUT=args.output;OUT.mkdir(parents=True,exist_ok=True)
cv=pd.read_csv(TAB/'67_existing_covariate_universe_coordinate_review.tsv',sep='\t')
cv=cv[cv.eligible].copy()
body=pd.read_csv(TAB/'67_gene_body_union_intervals.tsv',sep='\t')
body=body[body.gene_symbol.isin(cv.gene_symbol)].copy()
assert cv.gene_symbol.is_unique
seqrows=[]
groups={ch:q for ch,q in body.groupby('chr')}
fasta=args.fasta;sha=hashlib.sha256()
ch=None;chunks=[]
def finish(ch,chunks):
 if ch not in groups:return
 seq=b''.join(chunks).upper()
 for r in groups[ch].itertuples():
  assert r.start>=1 and r.end<=len(seq)
  seqrows.append(dict(gene_symbol=r.gene_symbol,chr=ch,start=r.start,end=r.end,cg_count=seq[r.start-1:r.end].count(b'CG')))
with fasta.open('rb') as f:
 for line in f:
  sha.update(line)
  if line.startswith(b'>'):
   finish(ch,chunks);ch=line[1:].split()[0].decode();chunks=[]
  elif ch in groups:chunks.append(line.strip())
 finish(ch,chunks)
counts=pd.DataFrame(seqrows)
assert len(counts)==len(body)
counts.to_csv(OUT/'67_body_interval_CpG_counts.tsv',sep='\t',index=False)
lengths=body.assign(length=body.end-body.start+1).groupby('gene_symbol').length.sum()
cpg=counts.groupby('gene_symbol').cg_count.sum()/lengths*1000
master=pd.read_csv(args.master_fpkm,sep='\t')
assert master.gene_id.is_unique
master=master.set_index('gene_id');expr=[]
for r in cv.itertuples():
 ids=r.RNA_IDs.split(';');mm=master.reindex(ids)
 assert mm.gene_name.eq(r.gene_symbol).all(),r.gene_symbol
 
 expr.append(mm[['A0','B0','C0']].mean(axis=1).max())
out=cv[['gene_symbol','RNA_IDs','chr','strand']].copy()
out['gene_length']=out.gene_symbol.map(lengths)
out['true_cg_density_per_kb']=out.gene_symbol.map(cpg)
out['baseline_fpkm']=expr
assert out.notna().all().all()
assert (out.gene_length>0).all() and (out.baseline_fpkm>=0).all()
out.to_csv(OUT/'gene_covariates.tsv',sep='\t',index=False)
print('Calculated gene length, CpG density and baseline RNA for',len(out),'genes.')
