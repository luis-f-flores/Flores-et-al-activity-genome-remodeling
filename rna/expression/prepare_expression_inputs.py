"""Prepare provider gene counts and gene lengths for TPM calculations."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd

KCL=['A0','B0','C0']+[f'{r}{t}{c}' for t in [1,6,24] for c in ['M','K'] for r in 'ABC']
BIC=[f'{c}_{r}' for c in ['Basal','Bic','Bic_nimo'] for r in 'ABC']

def prepare(count_file,output,name,samples,gene_order=None):
    d=pd.read_csv(count_file,sep='\t',float_precision='round_trip')
    if not d.gene_id.is_unique:raise ValueError('Provider gene IDs must be unique.')
    if gene_order is not None:
        order=pd.read_csv(gene_order,sep='\t',usecols=['gene_id']).gene_id
        if not order.is_unique or set(order)!=set(d.gene_id):raise ValueError('Gene universes differ between provider tables.')
        d=d.set_index('gene_id').loc[order].reset_index()
    values=d[samples].to_numpy(float)
    if not (d.gene_length>0).all() or not np.isfinite(values).all() or (values<0).any():
        raise ValueError('Gene lengths and raw counts must be finite and nonnegative.')
    if not np.array_equal(values,np.rint(values)):raise ValueError('Expected integer raw gene counts.')
    output=Path(output);output.mkdir(parents=True,exist_ok=True)
    d.to_csv(output/name,sep='\t',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--kcl-counts',type=Path,required=True)
    p.add_argument('--kcl-fpkm',type=Path,required=True,help='Provider gene table defining the common gene order.')
    p.add_argument('--bic-counts',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    a=p.parse_args()
    prepare(a.kcl_counts,a.output,'counts_KCl.tsv',KCL,a.kcl_fpkm)
    prepare(a.bic_counts,a.output,'counts_Bic.tsv',BIC)
