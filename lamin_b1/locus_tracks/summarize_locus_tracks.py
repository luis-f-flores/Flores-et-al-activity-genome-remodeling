'Calculate locus changes from per-IP 10-kb ChIP/Input values.'
import numpy as np,pandas as pd
from pathlib import Path
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'locus_tracks')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--fonts',type=Path,default=None)
args=ap.parse_args(); IN=args.inputs; OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)

raw=np.load(IN/'raw_10kb.npz'); lookup={(str(c),int(s)):i for i,(c,s) in enumerate(zip(raw['chromosome'],raw['start']))}
group=[[0,1],[2,3,4],[5,6,7]]; rows=[];computed={}
for source,genes in [('main.npz',['Ccbe1','Grin2a']),('supplement.npz',['Bdnf','Mpp7'])]:
 a=np.load(IN/source)
 for gene in genes:
  meta=pd.read_csv(IN/('main_genes.tsv' if source=='main.npz' else 'supplement_genes.tsv'),sep='	').set_index('gene_name').loc[gene]
  x=a[gene+'_starts'] if source=='main.npz' else a[gene+'__B1_starts'];cur=np.full((len(x),3),np.nan)
  for j,s in enumerate(x):
   ix=[lookup.get((meta.chromosome,int(s)+10000*k)) for k in [-1,0,1]]
   if any(k is None for k in ix):continue
   ip=np.log2(np.mean(2**raw['condition_depth'][ix],axis=0));cur[j]=[np.mean(ip[g]) for g in group]
  delta=cur[:,1:]-cur[:,[0]];computed[gene+'_starts']=x;computed[gene+'_delta']=delta
  rows.append(dict(gene=gene,bins=len(x),supported=int(np.isfinite(delta).all(axis=1).sum())))
np.savez_compressed(OUT/'subtraction_tracks.npz',**computed)
pd.DataFrame(rows).to_csv(OUT/'support.tsv',sep='	',index=False)
