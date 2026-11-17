"""Summarize activity-gene promoter methylation and matched RNA changes across time."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument('--counts',type=Path,required=True);p.add_argument('--activity-membership',type=Path,required=True)
 p.add_argument('--rna',type=Path,nargs=3,required=True);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args();d=pd.read_csv(a.counts,sep='\t');m=pd.read_csv(a.activity_membership,sep='\t')
 if d.duplicated(['gene','sample']).any() or not m.gene_name.is_unique:raise ValueError('Gene/sample keys must be unique.')
 np.testing.assert_allclose(d.mCG_pct,100*d.numCs/d.coverage,rtol=1e-12,atol=1e-12)
 g=m[['gene_name','gene_set','gene_id']].rename(columns={'gene_name':'gene','gene_set':'gene_class','gene_id':'RNA_KCl_gene_id'}).set_index('gene')
 g['gene_class']=g.gene_class.replace({'rapidPRGs':'rPRG','delayedPRGs':'dPRG','SRGs':'SRG'})
 if not g.gene_class.isin(['rPRG','dPRG','SRG']).all():raise ValueError('Unknown activity-response class.')
 wide=d.pivot(index='gene',columns='sample',values='mCG_pct').reindex(g.index)
 for sid in ['A0','B0','C0']+[f'{r}{t}{c}' for t in [1,6,24] for c in ['M','K'] for r in 'ABC']:g['mCG_'+sid]=wide[sid]
 g['mCG_0hr_mean']=g[['mCG_'+r+'0' for r in 'ABC']].mean(axis=1)
 for t,source in zip([1,6,24],a.rna):
  for condition in ['M','K']:
   cols=[f'mCG_{r}{t}{condition}' for r in 'ABC'];key=f'mCG_{t}{condition}'
   g[key+'_n']=g[cols].notna().sum(axis=1);g[key+'_mean']=g[cols].mean(axis=1,skipna=False);g[key+'_SD']=g[cols].std(axis=1,ddof=1,skipna=False)
  g[f'mCG_{t}K_minus_{t}M']=g[f'mCG_{t}K_mean']-g[f'mCG_{t}M_mean']
  rna=pd.read_csv(source,sep='\t').set_index('gene_id')
  if not rna.index.is_unique:raise ValueError('RNA gene IDs must be unique.')
  rows=rna.reindex(g.RNA_KCl_gene_id);rows.index=g.index
  if not np.isfinite(rows.log2FoldChange).all():raise ValueError('Missing activity-gene expression effects.')
  g[f'RNA_{t}h_log2FC']=rows.log2FoldChange;g[f'RNA_{t}h_padj']=rows.padj
 g=g.reset_index();g=pd.concat([g[g.gene_class==cl].sort_values(['mCG_0hr_mean','gene'],kind='stable') for cl in ['rPRG','dPRG','SRG']],ignore_index=True)
 a.output.parent.mkdir(parents=True,exist_ok=True);g.to_csv(a.output,sep='\t',index=False)

if __name__=='__main__':main()
