"""Build the genic DMR/DEG overlap background from genome annotation and RNA tables."""
from pathlib import Path
import re,math,json,hashlib,gzip
import numpy as np,pandas as pd
import argparse
ap=argparse.ArgumentParser(description='Build coordinate-matched genic DMR and DEG memberships.')
ap.add_argument('--rna',type=Path,nargs=3,required=True,help='Provider DE tables for 1 h, 6 h and recovery.')
ap.add_argument('--master-fpkm',type=Path,required=True)
ap.add_argument('--rna-gtf',type=Path,required=True)
ap.add_argument('--refgene-gtf',type=Path,required=True)
ap.add_argument('--gene-bed',type=Path,required=True)
ap.add_argument('--fai',type=Path,required=True)
ap.add_argument('--dmrs',type=Path,nargs=6,required=True,help='Hyper/hypo merged DMR tables for 1 h, 6 h and recovery, in that order.')
ap.add_argument('--output',type=Path,required=True)
args=ap.parse_args();TAB=args.output;TAB.mkdir(parents=True,exist_ok=True)
refs=[]
for path in args.rna:
 r=pd.read_csv(path,sep='\t')
 refs.append(r)
cols=['gene_id','gene_name','gene_chr','gene_start','gene_end','gene_strand']
gene_rows=[]
rna_gtf=args.rna_gtf
for line in gzip.open(rna_gtf,'rt'):
 if line.startswith('#'):continue
 a=line.rstrip('\n').split('\t')
 if len(a)!=9 or a[2]!='gene':continue
 at=dict(re.findall(r'(\w+) "([^"]*)"',a[8]))
 gene_rows.append([at['gene_id'],at.get('gene_name',''),a[0],int(a[3]),int(a[4]),a[6]])
rna=pd.DataFrame(gene_rows,columns=cols)
assert rna.gene_id.is_unique
for ref in refs:
 joined=ref[cols].merge(rna,on='gene_id',how='left',suffixes=('_delivered','_GTF'),validate='1:1')
 for c in cols[1:]:assert joined[c+'_delivered'].fillna('').astype(str).equals(joined[c+'_GTF'].fillna('').astype(str)),c
print('RNA GTF gene records verified',len(rna),flush=True)
rna['chr']=rna.gene_chr.astype(str).map(lambda c:{'GL456233.1':'chrX_GL456233_random','MT':'chrM'}.get(c,c if c.startswith('chr') else 'chr'+c))
tx=[]
for line in (gzip.open(args.refgene_gtf,'rt') if args.refgene_gtf.suffix=='.gz' else args.refgene_gtf.open()):
 a=line.rstrip('\n').split('\t')
 if len(a)!=9 or a[2]!='transcript':continue
 at=dict(re.findall(r'(\w+) "([^"]*)"',a[8]));tx.append([a[0],int(a[3]),int(a[4]),at['gene_name'],a[6],at['transcript_id']])
tx=pd.DataFrame(tx,columns=['chr','start','end','gene_symbol','strand','transcript'])
standard={f'chr{i}' for i in range(1,20)}|{'chrX','chrY','chrM'}
bed=pd.read_csv(args.gene_bed,sep='\t',header=None,names=['chr','start0','end','gene_symbol','score','strand'])
# Exclude ambiguous symbol and coordinate keys before defining the background.
bed=bed[~bed.gene_symbol.duplicated(keep=False)].copy()
bed=bed[~bed.duplicated(['chr','start0','end'],keep=False)].copy()
chromosomes=set(pd.read_csv(args.fai,sep='\t',header=None,usecols=[0])[0])
master=pd.read_csv(args.master_fpkm,sep='\t')
base=master.assign(baseline=master[['A0','B0','C0']].mean(axis=1)).groupby('gene_name').baseline.max()
bed['baseline_fpkm']=bed.gene_symbol.map(base)
cov=bed[bed.chr.isin(chromosomes)&np.isfinite(bed.baseline_fpkm)].copy()
cov['gene_length']=cov.end-cov.start0
cov=cov[['gene_symbol','gene_length','baseline_fpkm']]

rna_by={g:q for g,q in rna.groupby('gene_name')};tx_by={g:q for g,q in tx.groupby('gene_symbol')}
def union(iv):
 ans=[]
 for a,b in sorted(iv):
  if ans and a<=ans[-1][1]+1:ans[-1][1]=max(ans[-1][1],b)
  else:ans.append([int(a),int(b)])
 return ans
def overlaps(feat,dm):
 mask=np.zeros(len(feat),bool)
 for ch,q in feat.groupby('chr'):
  d=dm[dm.chr==ch].sort_values('start')
  if len(d)==0:continue
  starts=d.start.to_numpy();ends=np.maximum.accumulate(d.end.to_numpy());p=np.searchsorted(starts,q.end.to_numpy(),side='right')-1
  mask[q.index]=(p>=0)&(ends[np.maximum(p,0)]>=q.start.to_numpy())
 return mask

maprows=[];features=[];bodies=[];idfeatures=[]
allnames=sorted(set(tx.gene_symbol)&set(rna.gene_name.dropna()))
for g in allnames:
 ts=tx_by[g];rs=rna_by[g]
 candidates=[]
 for rr in rs.itertuples():
  same=ts[(ts.chr==rr.chr)&(ts.strand==rr.gene_strand)&(ts.start<=rr.gene_end)&(ts.end>=rr.gene_start)&ts.chr.isin(standard)]
  if len(same):candidates.append((rr,same))
 
 
 loci={(c[0].chr,c[0].gene_strand) for c in candidates}
 connected=len(loci)==1 and len(union((c[0].gene_start,c[0].gene_end) for c in candidates))==1
 if not candidates or not connected:
  maprows.append(dict(gene_symbol=g,status='no_standard_same_locus_RNA_match' if len(candidates)==0 else 'multiple_disconnected_RNA_loci',n_matching_RNA_IDs=len(candidates),RNA_IDs=';'.join(c[0].gene_id for c in candidates)))
  continue
 rr=candidates[0][0];same=pd.concat([c[1] for c in candidates]).drop_duplicates('transcript')
 iv=union(zip(same.start,same.end))
 row=dict(gene_symbol=g,status='single_locus_coordinate_match',n_matching_RNA_IDs=len(candidates),RNA_IDs=';'.join(c[0].gene_id for c in candidates),chr=rr.chr,strand=rr.gene_strand,RNA_start=min(c[0].gene_start for c in candidates),RNA_end=max(c[0].gene_end for c in candidates),transcript_count=len(same),body_interval_count=len(iv),body_union_length=sum(b-a+1 for a,b in iv),body_min_start=iv[0][0],body_max_end=iv[-1][1])
 maprows.append(row)
 for a,b in iv:bodies.append(dict(gene_symbol=g,chr=rr.chr,strand=rr.gene_strand,start=a,end=b))
 allivs=list(zip(same.start,same.end))
 for a,b,s in same[['start','end','strand']].itertuples(index=False,name=None):allivs.append((max(1,a-1000),a-1) if s=='+' else (b+1,b+1000))
 for a,b in union(allivs):features.append(dict(gene_symbol=g,chr=rr.chr,strand=rr.gene_strand,start=a,end=b))
 for ri,si in candidates:
  ints=list(zip(si.start,si.end))
  for a,b,s in si[['start','end','strand']].itertuples(index=False,name=None):ints.append((max(1,a-1000),a-1) if s=='+' else (b+1,b+1000))
  for a,b in union(ints):idfeatures.append(dict(gene_symbol=g,RNA_ID=ri.gene_id,chr=ri.chr,strand=ri.gene_strand,start=a,end=b))
mapping=pd.DataFrame(maprows);feat=pd.DataFrame(features);body=pd.DataFrame(bodies)
idfeat=pd.DataFrame(idfeatures)
mapping.to_csv(TAB/'67_gene_RNA_coordinate_crosswalk.tsv',sep='\t',index=False)
feat.to_csv(TAB/'67_promoter_body_union_intervals.tsv',sep='\t',index=False)
body.to_csv(TAB/'67_gene_body_union_intervals.tsv',sep='\t',index=False)
idfeat.to_csv(TAB/'67_RNA_ID_promoter_body_intervals.tsv',sep='\t',index=False)
valid=mapping[mapping.status=='single_locus_coordinate_match']
cv=cov.merge(valid,on='gene_symbol',how='left',validate='1:1')
cv['eligible']=cv.status.eq('single_locus_coordinate_match')
cv['length_changed']=cv.eligible&cv.gene_length.ne(cv.body_union_length)
cv.to_csv(TAB/'67_existing_covariate_universe_coordinate_review.tsv',sep='\t',index=False)
summary=[];members=[]
for ti,(t,r) in enumerate(zip([1,6,24],refs)):
 tp=f'{t}hr';sig=r[r.gene_name.notna()&(r.padj<=.01)&(r.log2FoldChange.abs()>=math.log2(1.3))]
 full_degs=set(sig.gene_name)
 allowed_pairs={(rr.gene_symbol,i) for rr in valid.itertuples() for i in rr.RNA_IDs.split(';')}
 sig=sig[[ (rr.gene_name,rr.gene_id) in allowed_pairs for rr in sig.itertuples() ]]
 degs=set(sig.gene_name)
 dm=pd.concat([pd.read_csv(path,sep='\t') for path in args.dmrs[2*ti:2*ti+2]],ignore_index=True)
 dmrgenes=set(feat.loc[overlaps(feat,dm),'gene_symbol'])
 hitids=set(idfeat.loc[overlaps(idfeat,dm),'RNA_ID'])
 new=set(sig.loc[sig.gene_id.isin(hitids),'gene_name'])
 eligible=set(cv.loc[cv.eligible,'gene_symbol'])
 
 assert new==degs&dmrgenes,(tp,sorted(new^(degs&dmrgenes)))
 row=dict(timepoint=tp,corrected_main=len(new),eligible_observed=len(new&eligible),full_main_outside_eligible=len(new-eligible),eligible_DEGs=len(degs&eligible),all_original_DEGs=len(full_degs))
 summary.append(row)
 for g in sorted(set(cov.gene_symbol)|new):
  members.append(dict(timepoint=tp,gene_symbol=g,is_DEG=g in degs,is_promoter_body_DMR=g in dmrgenes,corrected_main=g in new,permutation_eligible=g in eligible,eligible_observed=g in new&eligible))
pd.DataFrame(summary).to_csv(TAB/'67_preparation_count_summary.tsv',sep='\t',index=False)
pd.DataFrame(members).to_csv(TAB/'67_permutation_gene_membership.tsv',sep='\t',index=False)
print('Built genic memberships for 1 h, 6 h and recovery.')
