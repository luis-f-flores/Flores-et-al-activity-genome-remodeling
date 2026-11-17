"""Measure gene-span B1 changes within pairwise shared-LAD gene sets."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy import sparse
from scipy.stats import fisher_exact

CONS=['Bic_vs_TTX','Bic_Nimo_vs_Bic','Bic_Nimo_vs_TTX']
PAIRS=[(0,1),(1,2),(0,2)]
GROUPS=[[0,1],[2,3,4],[5,6,7]]
SCHEMES=['condition_depth']

def read(path):
    return pd.read_csv(path,sep='\t')

def tab(path, data):
    pd.DataFrame(data).to_csv(path,sep='\t',index=False)

def merge(rows):
 out=[]
 for a,b in sorted(rows):
  if out and a<=out[-1][1]:out[-1]=(out[-1][0],max(b,out[-1][1]))
  else:out.append((int(a),int(b)))
 return out

def build(signals, decoded, membership, autosomal_reference, sex_reference, output):
 R=output; (R/'tables').mkdir(parents=True,exist_ok=True); (R/'GO_RNA').mkdir(exist_ok=True)
 a=dict(np.load(signals));st=np.load(decoded)['states'];eligible=read(membership)
 masks=dict(np.load(autosomal_reference));masks.update(dict(np.load(sex_reference)))
 bad={};idx={}
 for ch in np.unique(a['chromosome']):
  ints=merge([(x,y) for k in ['mask','nonACGT'] for x,y in zip(masks[ch+'_'+k+'_start'],masks[ch+'_'+k+'_end'])]);bad[ch]=(np.array([x for x,y in ints]),np.array([y for x,y in ints]));ii=np.flatnonzero(a['chromosome']==ch);idx[ch]=(ii,a['start'][ii])
 def clean(ch,l,r):
  ss,ee=bad[ch];i=np.searchsorted(ee,l,side='right');j=np.searchsorted(ss,r,side='left');return int(r-l-np.maximum(0,np.minimum(ee[i:j],r)-np.maximum(ss[i:j],l)).sum())
 good=eligible[eligible.TSS_bin_index.ge(0)].drop_duplicates('gene_id').sort_values('gene_id');weights=[];coverage=[];ids=[]
 for row in good.itertuples():
  ii,ss=idx[row.chromosome];bins=ii[np.searchsorted(ss,row.start0-30000,side='right'):np.searchsorted(ss,row.end0,side='left')];ww=[];bb=[]
  for b in bins:
   w=clean(row.chromosome,max(row.start0,int(a['start'][b])),min(row.end0,int(a['start'][b]+30000)))
   if w:bb.append(b);ww.append(w)
  if not ww:continue
  k=len(ids);ids.append(row.gene_id);weights.extend((k,int(b),int(w)) for b,w in zip(bb,ww));coverage.append(dict(gene_id=row.gene_id,assessed_gene_bp=sum(ww),overlapping_bins=len(bb)))
 wm=sparse.csr_matrix(([w for k,b,w in weights],([k for k,b,w in weights],[b for k,b,w in weights])),shape=(len(ids),len(st)))
 den=np.asarray(wm.sum(1)).ravel();rep=[];scores={}
 for s in SCHEMES:
  sc=np.asarray(wm@a[s])/den[:,None];scores[s]=pd.DataFrame(sc,index=ids,columns=a['IP_ids']);z=scores[s].reset_index(names='gene_id');z['scheme']=s;rep.append(z)
 tab(R/'tables/76_gene_replicates.tsv',pd.concat(rep));tab(R/'tables/gene_coverage.tsv',coverage)
 np.savez_compressed(R/'tables/gene_weights.npz',data=wm.data,indices=wm.indices,indptr=wm.indptr,shape=wm.shape,gene_ids=ids)
 allout=[];gom=[];venn=[];clusters=[]
 for con,(base,target) in zip(CONS,PAIRS):
  d=eligible[eligible.contrast==con].copy();bi=d.TSS_bin_index.to_numpy(int);ok=bi>=0
  for name,grp in [('baseline_LAD',base),('target_LAD',target)]:d[name]=ok & (st[np.maximum(bi,0)][:,GROUPS[grp]]>=1).any(1)
  gom.append(d)
  valid=d[d.gene_id.isin(ids)&np.isfinite(d.pvalue)&np.isfinite(d.padj)&np.isfinite(d.log2FoldChange)].copy();valid['LAD']=valid.baseline_LAD|valid.target_LAD
  x=valid.LAD.to_numpy();y=valid.DEG.to_numpy();cells=[int((x&y).sum()),int((x&~y).sum()),int((~x&y).sum()),int((~x&~y).sum())];or_,p=fisher_exact(np.array(cells).reshape(2,2));venn.append(dict(contrast=con,eligible=len(valid),LAD=int(x.sum()),DEG=int(y.sum()),overlap=cells[0],a=cells[0],b=cells[1],c=cells[2],d=cells[3],sample_OR=or_,P=p))
  tab(R/'tables'/f'{con}_eligible.tsv',valid)
  d=valid[valid.LAD].sort_values('gene_id');ci={};cid=0;union=(st[:,GROUPS[base]+GROUPS[target]]>=1).any(1)
  for ch in np.unique(a['chromosome']):
   rows=[(int(a['start'][b]),int(a['start'][b]+30000)) for b in np.flatnonzero(union&(a['chromosome']==ch))];sub=d[d.chromosome==ch];rows+=list(zip(sub.start0,sub.end0))
   for l,r in merge(rows):
    clusters.append(dict(contrast=con,cluster_id=cid,chromosome=ch,start0=l,end0=r))
    for gid in sub.loc[(sub.start0>=l)&(sub.end0<=r),'gene_id']:ci[gid]=cid
    cid+=1
  assert set(ci)==set(d.gene_id)
  for s in SCHEMES:
   v=scores[s].loc[d.gene_id].to_numpy();mu=np.column_stack([v[:,ix].mean(1) for ix in GROUPS]);delta=mu[:,target]-mu[:,base]
   out=d.copy();out['scheme']=s;out['RNA_classification']=out.classification;out['RNA_log2FC']=out.log2FoldChange;out['RNA_padj']=out.padj;out['cluster_id']=out.gene_id.map(ci)
   for j,k in enumerate(['B1_TTX','B1_Bic','B1_Nimo']):out[k]=mu[:,j]
   out['B1_delta_log2']=delta;out['B1_direction']=np.where(delta<0,'decreased',np.where(delta>0,'increased','zero'));allout.append(out)
 tab(R/'tables/76_gene_measurements.tsv',pd.concat(allout));tab(R/'tables/dependence_clusters.tsv',clusters);tab(R/'tables/DEG_LAD_overlap.tsv',venn)
 tab(R/'GO_RNA/92_GO_membership_input.tsv',pd.concat(gom));tab(R/'tables/all_gene_membership.tsv',pd.concat(gom))
 print('Wrote gene-span B1 measurements.',flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['signals','decoded','membership','autosomal-reference','sex-reference','output']:
        p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args()
    build(a.signals,a.decoded,a.membership,a.autosomal_reference,a.sex_reference,a.output)

if __name__=='__main__':
    main()
