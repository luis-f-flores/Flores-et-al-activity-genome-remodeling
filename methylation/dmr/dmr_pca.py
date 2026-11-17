'Calculate the locus-centered DMR-union PCA from the supplied sample matrix.'
from pathlib import Path
import argparse
import numpy as np,pandas as pd
ap=argparse.ArgumentParser(description=__doc__);ap.add_argument('--inputs',type=Path,required=True);ap.add_argument('--output',type=Path,required=True)
a=ap.parse_args();IN=a.inputs;OUT=a.output;OUT.mkdir(parents=True,exist_ok=True)
z=np.load(IN/'DMR_sample_matrix.npz');ids=list(z['samples']);m=z['mCG_percent']
keep=np.isfinite(m).all(axis=1)&(np.nanstd(m,axis=1)>0);x=m[keep];center=x-x.mean(axis=1,keepdims=True)
gram=center.T@center
vals,vec=np.linalg.eigh(gram);vals=vals[::-1].clip(0);vec=vec[:,::-1]
scores=vec*np.sqrt(vals)[None,:];load=center@vec[:,:2]/np.sqrt(vals[:2])[None,:]
for k in range(2):
 if load[np.argmax(abs(load[:,k])),k]<0:load[:,k]*=-1;scores[:,k]*=-1
explained=100*vals/vals.sum()
pd.DataFrame({'sample':ids,'PC1':scores[:,0],'PC2':scores[:,1]}).to_csv(OUT/'DMR_PCA_scores.tsv',sep='	',index=False)
pd.DataFrame({'PC':np.arange(1,22),'variance_percent':explained}).to_csv(OUT/'DMR_PCA_variance.tsv',sep='	',index=False)
print('Calculated DMR-union PCA:',int(keep.sum()),'regions,',len(ids),'samples')
