"""Aggregate the three supplied baseline CpG reports on the published B1 grid."""
from pathlib import Path
import argparse,gzip,collections,json
import numpy as np,pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True)
parser.add_argument('--reports',type=Path,required=True)
parser.add_argument('--reference-mask',type=Path)
parser.add_argument('--output',type=Path,required=True)
args=parser.parse_args();inputs=args.inputs;output=args.output;output.mkdir(parents=True,exist_ok=True)
bins=np.load(inputs/'B1_bins.npz');mask=np.load(args.reference_mask or inputs/'reference_mask.npz')
bc=bins['chromosome'];bs=bins['start'];n=len(bs)
def dump(path,value):path.write_text(json.dumps(value,indent=2)+'\n')
off={ch:np.full(int(mask['lengths'][i])//30000+1,-1,dtype=np.int32) for i,ch in enumerate(mask['chromosomes'])}
for i,(ch,s) in enumerate(zip(bc,bs)):off[ch][int(s)//30000]=i
meth=np.zeros((n,3));depth=np.zeros((n,3));sites=np.zeros((n,3),dtype=np.int64);equal=np.zeros((n,3))
raws=[];coverage_stats=[]
for j,sample in enumerate(['A0','B0','C0']):
 p=next(args.reports.glob(sample+'_*.txt'));raws.append(p)
 cache=output/f'{sample}_bins.npz'
 if cache.exists():
  z=np.load(cache);meth[:,j]=z['meth'];depth[:,j]=z['depth'];sites[:,j]=z['sites'];equal[:,j]=z['equal']
  coverage_stats.append(json.loads((output/f'{sample}_coverage.json').read_text()));continue
 print('BASELINE coverage histogram',sample,flush=True)
 histogram={}
 for q in pd.read_csv(p,sep='\t',header=None,usecols=[3,4],names=['m','u'],dtype=np.int32,chunksize=1000000):
  cov=q.m.to_numpy(dtype=np.int64)+q.u.to_numpy(dtype=np.int64);vals,cnt=np.unique(cov[cov>=4],return_counts=True)
  for a,b in zip(vals,cnt):histogram[int(a)]=histogram.get(int(a),0)+int(b)
 vv=np.array(sorted(histogram));nn=np.array([histogram[int(v)] for v in vv]);cumul=np.cumsum(nn);N=int(nn.sum())
 pos=(N-1)*.999;lower=int(np.floor(pos));upper=int(np.ceil(pos))
 cutoff=float(vv[np.searchsorted(cumul,lower+1)]*(1-(pos-lower))+vv[np.searchsorted(cumul,upper+1)]*(pos-lower))
 print('BASELINE aggregate',sample,'upper coverage',cutoff,flush=True)
 retained=0
 for q in pd.read_csv(p,sep='\t',header=None,usecols=[0,1,3,4],names=['chrom','pos','m','u'],dtype={'chrom':'category','pos':'int64','m':'int64','u':'int64'},chunksize=1000000):
  cov=q.m.to_numpy()+q.u.to_numpy();q=q[(cov>=4)&(cov<cutoff)]
  for ch,g in q.groupby('chrom',observed=True):
   if ch not in off:continue
   pos0=g.pos.to_numpy()-1;cc=g.m.to_numpy()+g.u.to_numpy();mm=g.m.to_numpy()
   ind=off[ch][pos0//30000];keep=ind>=0
   for typ in ['mask','nonACGT']:
    a=mask[ch+'_'+typ+'_start'];b=mask[ch+'_'+typ+'_end'];k=np.searchsorted(a,pos0,side='right')-1
    keep &= (k<0)|(pos0>=b[np.maximum(k,0)])
   ind=ind[keep];cc=cc[keep];mm=mm[keep];retained+=len(ind)
   meth[:,j]+=np.bincount(ind,weights=mm,minlength=n);depth[:,j]+=np.bincount(ind,weights=cc,minlength=n)
   sites[:,j]+=np.bincount(ind,minlength=n)
   equal[:,j]+=np.bincount(ind,weights=mm/cc,minlength=n)
 covmeta={'sample':sample,'minimum_coverage':4,'upper_coverage_exclusive':cutoff,'original_covered_records':N,'retained_selected_bin_records':retained,'input_bytes':p.stat().st_size,'input_mtime_ns':p.stat().st_mtime_ns}
 np.savez_compressed(cache,meth=meth[:,j],depth=depth[:,j],sites=sites[:,j],equal=equal[:,j]);dump(output/f'{sample}_coverage.json',covmeta);coverage_stats.append(covmeta)

np.savez_compressed(output/'baseline_counts.npz',meth=meth,depth=depth,sites=sites,B1=bins['means'][:,0])
