'Calculate Figure 5 TPM values using the gene-wise display transformation.'
from pathlib import Path
import numpy as np
import pandas as pd
COLS=['A0','B0','C0']+[f'{rep}{tp}{tr}' for tp in (1,6,24) for tr in ('M','K') for rep in 'ABC']
CLASSES=['rapidPRGs','delayedPRGs','SRGs']
def row_z(a):
 x=np.log2(np.asarray(a,float)+1);sd=x.std(axis=1,ddof=1)
 assert (sd>0).all()
 return (x-x.mean(axis=1,keepdims=True))/sd[:,None]
def calculate(inputs,output):
 I=Path(inputs);Q=Path(output);Q.mkdir(parents=True,exist_ok=True)
 counts=pd.read_csv(I/'TPM_gene_counts_lengths.tsv.gz',sep='\t',float_precision='round_trip').set_index('gene_id')
 order=pd.read_csv(I/'TPM_gene_order.csv');members=pd.read_csv(I/'TPM_activity_membership.tsv',sep='\t')
 assert len(counts)==54532 and counts.index.is_unique and (counts.gene_length>0).all()
 assert len(order)==7257 and order.gene_id.is_unique and order.groupby('Cluster').size().tolist()==[1601,1781,1006,2263,606]
 assert len(members)==173 and members.gene_id.is_unique and [int(members.gene_set.eq(c).sum()) for c in CLASSES]==[19,116,38]
 c=counts[COLS].to_numpy(float);assert np.isfinite(c).all() and (c>=0).all() and np.array_equal(c,np.rint(c))
 rpk=counts[COLS].div(counts.gene_length/1000,axis=0);tpm=rpk.div(rpk.sum(axis=0),axis=1)*1e6
 assert np.allclose(tpm.sum(),1e6,atol=1e-7,rtol=0)
 z=row_z(tpm.loc[order.gene_id,COLS]);az=row_z(tpm.loc[members.gene_id,COLS]);am=np.stack([az[members.gene_set.eq(cl)].mean(axis=0) for cl in CLASSES])
 d=order.copy();d[COLS]=z;d.to_csv(Q/'Fig5B_TPM_Z.tsv',sep='\t',index=False)
 d=members.copy();d[COLS]=az;d.to_csv(Q/'Fig5C_gene_TPM_Z.tsv',sep='\t',index=False)
 pd.DataFrame(am,index=CLASSES,columns=COLS).to_csv(Q/'Fig5C_class_mean_TPM_Z.tsv',sep='\t')
 return z,am
if __name__=='__main__':
 import argparse
 a=argparse.ArgumentParser();a.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'kcl_rna');a.add_argument('--output',type=Path,required=True);p=a.parse_args();calculate(p.inputs,p.output)
