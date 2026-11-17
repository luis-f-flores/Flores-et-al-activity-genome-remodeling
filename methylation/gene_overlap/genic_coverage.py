"""Measure CpG-report read support across the genic analysis intervals."""
from pathlib import Path
from concurrent.futures import ProcessPoolExecutor
import argparse
import numpy as np
import pandas as pd
def coverage_job(arg):
    path, intervals, symbols, out = arg
    ix={g:i for i,g in enumerate(symbols)}
    n=len(symbols); sums=np.zeros((n,4),dtype=np.float64)
    bychr={}
    for ch, d in intervals.groupby('chr',sort=False):
        bychr[ch]=(d.start.to_numpy(),d.end.to_numpy(),d.gene_symbol.map(ix).to_numpy())
    # All source rows, including zero-depth cytosines, contribute to denominators.
    for z in pd.read_csv(path,sep='\t',header=None,usecols=[0,1,3,4],
                         names=['chr','pos','m','u'],chunksize=1000000):
        for ch,d in z.groupby('chr',sort=False):
            if ch not in bychr: continue
            s,e,gi=bychr[ch]
            p=d.pos.to_numpy(); depth=(d.m+d.u).to_numpy()
            order=np.argsort(p,kind='stable'); p=p[order]; depth=depth[order]
            left=np.searchsorted(p,s,'left'); right=np.searchsorted(p,e,'right')
            for col,v in enumerate([np.ones(len(p)),depth,(depth>0),(depth>=4)]):
                c=np.r_[0,np.cumsum(v,dtype=np.float64)]
                np.add.at(sums[:,col],gi,c[right]-c[left])
    d=pd.DataFrame(sums,columns=['cytosines','depth_sum','covered','coverage4'])
    d.insert(0,'gene_symbol',symbols); d.to_csv(out,sep='\t',index=False)
    return str(out)


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--covariates',type=Path,required=True)
    p.add_argument('--intervals',type=Path,required=True)
    p.add_argument('--samples',type=Path,required=True,help='TSV with sample and report columns for the nine mannitol libraries.')
    p.add_argument('--output',type=Path,required=True)
    p.add_argument('--workers',type=int,default=2)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    cv=pd.read_csv(a.covariates,sep='\t');symbols=cv.gene_symbol.tolist()
    feat=pd.read_csv(a.intervals,sep='\t');feat=feat[feat.gene_symbol.isin(symbols)].copy()
    samples=pd.read_csv(a.samples,sep='\t').set_index('sample')
    jobs=[]
    for tp in [1,6,24]:
        for rep in 'ABC':
            sample=f'{rep}{tp}M';path=Path(samples.loc[sample,'report'])
            if not path.is_absolute():path=a.samples.parent/path
            jobs.append((path,feat,symbols,a.output/f'{sample}_coverage.tsv'))
    with ProcessPoolExecutor(max_workers=a.workers) as pool:
        for result in pool.map(coverage_job,jobs):print(result,flush=True)

if __name__=='__main__':main()
