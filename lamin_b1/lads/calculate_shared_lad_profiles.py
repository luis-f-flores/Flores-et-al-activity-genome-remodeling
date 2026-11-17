"""Calculate B1 profiles across the all-condition LAD union and its boundaries."""
from pathlib import Path
import argparse
import csv
import numpy as np
from scipy import sparse
from lad_profile_windows import profile_queries

GROUPS = [[0,1],[2,3,4],[5,6,7]]
NAMES = ['TTX','Bic','Bic+Nimo']
SCHEMES = ['condition_depth','condition_fly']
W = 30000


def main():
 p=argparse.ArgumentParser(description=__doc__)
 for key in ['signals','decoded','fai','output']:p.add_argument('--'+key,type=Path,required=True)
 args=p.parse_args();RUN=args.output;RUN.mkdir(parents=True,exist_ok=True)
 def table(name,rows):
  with (RUN/name).open('w',newline='') as f:
   writer=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');writer.writeheader();writer.writerows(rows)
 v=dict(np.load(args.signals));states=np.load(args.decoded)['states']
 ch,st,bp=v['chromosome'],v['start'],v['clean_bp'];ips=v['IP_ids'].tolist()
 masks=np.column_stack([(states[:,g]>0).any(1) for g in GROUPS]);union=masks.any(1)
 sizes={r[0]:int(r[1]) for r in csv.reader(args.fai.open(),delimiter='\t')}
 idx={c:np.flatnonzero(ch==c) for c in dict.fromkeys(ch)}
 regions=[]
 for c,ii in idx.items():
  starts=st[ii[union[ii]]]
  for group in np.split(starts,np.flatnonzero(np.diff(starts)!=W)+1):
   if len(group):regions.append((str(c),int(group[0]),min(int(group[-1])+W,sizes[c])))
 table('regions.tsv',[dict(chromosome=c,start0=a,end0=b) for c,a,b in regions])
 z=np.column_stack([v[k] for k in SCHEMES])
 def overlap(c,a,b,masked=False):
  a=max(float(a),0);b=min(float(b),sizes[c]);ii=idx[c]
  if b<=a:return np.array([],int),np.array([],float)
  lo=np.searchsorted(st[ii],a-W,side='right');hi=np.searchsorted(st[ii],b,side='left');ids=ii[lo:hi]
  weights=np.maximum(0,np.minimum(st[ids]+W,b)-np.maximum(st[ids],a))*bp[ids]/W
  keep=(weights>0)&(~union[ids] if masked else True)
  return ids[keep],weights[keep]
 arrays={};rows=[];supports=[];error=0.;checked=0
 for kind in ['body','border']:
  queries,x=profile_queries(regions,kind);pos=np.array([q[3] for q in queries]);outside=np.array([(p<17 or p>=67) if kind=='body' else x[p]<0 for _,_,_,p,_ in queries])
  rr=[];cc=[];ww=[];covered=np.zeros(len(queries));excluded=np.zeros(len(queries))
  for j,(c,a,b,position,rid) in enumerate(queries):
   ids,w=overlap(c,a,b);total=w.sum()
   if outside[j]:
    keep=~union[ids];excluded[j]=w[~keep].sum();ids,w=ids[keep],w[keep]
    assert not union[ids].any()
   covered[j]=w.sum()
   if not len(w) or not w.sum():continue
   rr.extend([j]*len(ids));cc.extend(ids.tolist());ww.extend((w/w.sum()).tolist())
  mat=sparse.csr_matrix((ww,(rr,cc)),shape=(len(queries),len(st)))
  measured=np.asarray(mat@z);measured[covered==0]=np.nan
  n=np.bincount(pos[covered>0],minlength=len(x));assert (n>0).all()
  means=np.array([np.nanmean(measured[pos==i],axis=0) for i in range(len(x))])
  arrays[kind+'_x']=x;arrays[kind+'_means']=means;arrays[kind+'_n']=n
  for p in range(len(x)):
   take=pos==p
   supports.append(dict(profile=kind,position=p,x=float(x[p]),contributing_regions=int(n[p]),query_weighted_bp=float(covered[take].sum()),excluded_LAD_weighted_bp=float(excluded[take].sum())))
   for si,scheme in enumerate(SCHEMES):
    for i,ip in enumerate(ips):rows.append(dict(profile=kind,position=p,x=float(x[p]),scheme=scheme,sample=ip,mean_log2_signal=float(means[p,si*8+i]),contributing_regions=int(n[p])))
 table('profiles.tsv',rows);table('profile_support.tsv',supports)
 
 contrasts=[]
 for ref,mask in zip(NAMES+['Union'],[*masks.T,union]):
  for scheme in SCHEMES:
   inside=np.average(v[scheme][mask],axis=0,weights=bp[mask]);outside=np.average(v[scheme][~mask],axis=0,weights=bp[~mask])
   for i,ip in enumerate(ips):contrasts.append(dict(reference=ref,scheme=scheme,sample=ip,condition=next(n for n,g in zip(NAMES,GROUPS) if i in g),inside_log2=float(inside[i]),outside_log2=float(outside[i]),contrast_log2=float(inside[i]-outside[i]),contrast_ratio=float(2**(inside[i]-outside[i]))))
 table('inside_outside.tsv',contrasts);np.savez_compressed(RUN/'profile_arrays.npz',**arrays)

if __name__=='__main__':main()
