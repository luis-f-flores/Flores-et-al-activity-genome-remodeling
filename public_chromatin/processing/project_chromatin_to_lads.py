"""Project published mm9 ChIP counts onto fixed mm10 TTX LADs."""
from pathlib import Path
import argparse,gzip,collections
import numpy as np,pandas as pd
parser=argparse.ArgumentParser(description=__doc__)
parser.add_argument('--inputs',type=Path,required=True)
parser.add_argument('--output',type=Path,required=True)
parser.add_argument('--lads',type=Path)
parser.add_argument('--reference-mask',type=Path)
parser.add_argument('--b1-grid',type=Path)
parser.add_argument('--chain',type=Path)
args=parser.parse_args();inputs=args.inputs;output=args.output;output.mkdir(parents=True,exist_ok=True)
lads=pd.read_csv(args.lads or inputs/'TTX_LADs.tsv',sep='\t');mask=np.load(args.reference_mask or inputs/'reference_mask.npz');universe=np.load(args.b1_grid or inputs/'B1_bins.npz')
sizes=dict(zip(mask['chromosomes'],mask['lengths']));chroms=['chr'+str(i) for i in range(1,20)]
grid=pd.read_csv(inputs/'mm9_5kb_windows.tsv',sep='\t')
off={c:int(np.flatnonzero(grid.chromosome.eq(c))[0]) for c in chroms};nwin={c:int(grid.chromosome.eq(c).sum()) for c in chroms}
countfiles=[inputs/(n+'_mm9_5kb.npz') for n in ['H3K9me2','H3_1','H3_2']]+[inputs/f'{n}_{r}_mm9_5kb.npz' for r in [1,2,3] for n in ['DNMT3A_ChIP','Input']]
cpms=[np.load(p)['counts']/int(np.load(p)['total_unique_reads'])*1e6 for p in countfiles]
with np.errstate(divide='ignore',invalid='ignore'):
 signal=np.column_stack([np.log2(cpms[0]/((cpms[1]+cpms[2])/2)),np.log2(cpms[0]/cpms[1]),np.log2(cpms[0]/cpms[2])]+[np.log2(cpms[i]/cpms[i+1]) for i in [3,5,7]])
valid=np.isfinite(signal).all(1);D=6
names=['H3K9me2_H3mean','H3K9me2_H3rep1','H3K9me2_H3rep2','DNMT3A_rep1','DNMT3A_rep2','DNMT3A_rep3']
blocks=collections.defaultdict(list);chain=args.chain or inputs/'mm9ToMm10.over.chain.gz'
for line in gzip.open(chain,'rt'):
 a=line.split()
 if not a:continue
 if a[0]=='chain':header=a;tp=int(a[5]);qp=int(a[10]);continue
 width=int(a[0])
 if header[2] in chroms and header[7] in chroms:
  assert header[4]=='+';rev=header[9]=='-';lo=int(header[8])-qp-width if rev else qp
  blocks[header[7]].append((lo,lo+width,chroms.index(header[2]),tp,int(rev)))
 tp+=width;qp+=width
 if len(a)==3:tp+=int(a[1]);qp+=int(a[2])
class Track:
 def __init__(self,st,en,v):
  self.st=st;self.en=en;self.v=v;self.length=en-st
  assert np.all(st[1:]>=en[:-1]) and np.all(en>st)
  self.den=np.r_[0,np.cumsum(self.length,dtype=float)]
  self.num=np.vstack([np.zeros(D),np.cumsum(self.length[:,None]*v,axis=0)])
 def integral(self,x):
  x=np.asarray(x,dtype=float)
  if not len(self.st):return np.zeros(x.shape+(D,)),np.zeros_like(x)
  i=np.searchsorted(self.st,x,side='right')-1;k=np.maximum(i,0);z=np.clip(x-self.st[k],0,self.length[k]);ok=i>=0
  return np.where(ok[...,None],self.num[k]+z[...,None]*self.v[k],0),np.where(ok,self.den[k]+z,0)
 def sums(self,a,b):
  nb,db=self.integral(b);na,da=self.integral(a);return nb-na,db-da
tracks={};gin=np.zeros(D);gout=np.zeros(D);bin_=bout=0.;direct_checks=[]
for c in chroms:
 b=np.asarray(blocks[c],dtype=np.int64);p,ii=np.unique(np.r_[b[:,0],b[:,1]],return_inverse=True);k=len(b)
 cnt=np.zeros(len(p),dtype=np.int64);ids=cnt.copy()
 np.add.at(cnt,ii[:k],1);np.add.at(cnt,ii[k:],-1);np.add.at(ids,ii[:k],np.arange(k)+1);np.add.at(ids,ii[k:],-np.arange(k)-1)
 cnt=np.cumsum(cnt);ids=np.cumsum(ids);keep=cnt[:-1]==1;ql=p[:-1][keep];qr=p[1:][keep];bb=b[ids[:-1][keep]-1]
 sl=bb[:,3]+np.where(bb[:,4]>0,bb[:,1]-qr,ql-bb[:,0]);sr=bb[:,3]+np.where(bb[:,4]>0,bb[:,1]-ql,qr-bb[:,0])
 num=(sr-1)//5000-sl//5000+1;seg=np.repeat(np.arange(len(sl)),num);within=np.arange(num.sum())-np.repeat(np.cumsum(num)-num,num)
 bins=sl[seg]//5000+within;left=np.maximum(sl[seg],bins*5000);right=np.minimum(sr[seg],(bins+1)*5000)
 st=ql[seg]+np.where(bb[seg,4]>0,sr[seg]-right,left-sl[seg]);en=ql[seg]+np.where(bb[seg,4]>0,sr[seg]-left,right-sl[seg]);code=bb[seg,2]
 offset=np.asarray([off[x] for x in chroms]);nn=np.asarray([nwin[x] for x in chroms]);keep=bins<nn[code]
 st,en,idx=st[keep],en[keep],offset[code[keep]]+bins[keep];keep=valid[idx];st,en,idx=st[keep],en[keep],idx[keep]
 order=np.argsort(st,kind='stable');st,en,idx=st[order],en[order],idx[order];assert np.all(st[1:]>=en[:-1])
 merged=[]
 for a,z in sorted((int(a),int(z)) for typ in ['mask','nonACGT'] for a,z in zip(mask[c+'_'+typ+'_start'],mask[c+'_'+typ+'_end'])):
  if merged and a<=merged[-1][1]:merged[-1][1]=max(z,merged[-1][1])
  else:merged.append([a,z])
 m=np.asarray(merged,dtype=np.int64);dom=lads[lads.chromosome==c].sort_values('start0');la=dom.start0.to_numpy();lb=dom.end0.to_numpy()
 ua=universe['start'][universe['chromosome']==c];ub=np.minimum(ua+30000,int(sizes[c]))
 cuts=np.unique(np.r_[st,en,m.ravel(),la,lb,ua,ub]);aa=cuts[:-1];zz=cuts[1:]
 j=np.searchsorted(st,aa,side='right')-1;jj=np.maximum(j,0);mi=np.searchsorted(m[:,0],aa,side='right')-1;ui=np.searchsorted(ua,aa,side='right')-1
 good=(j>=0)&(zz<=en[jj])&((mi<0)|(aa>=m[np.maximum(mi,0),1]))&(zz>aa)&(ui>=0)&(zz<=ub[np.maximum(ui,0)])
 aa,zz,jj=aa[good],zz[good],jj[good];v=signal[idx[jj]];li=np.searchsorted(la,aa,side='right')-1;inside=(li>=0)&(aa<lb[np.maximum(li,0)])
 full=Track(aa,zz,v);inn=Track(aa[inside],zz[inside],v[inside]);out=Track(aa[~inside],zz[~inside],v[~inside]);tracks[c]=(full,out)
 gin+=inn.num[-1];gout+=out.num[-1];bin_+=inn.den[-1];bout+=out.den[-1]
 for row in dom.iloc[:3].itertuples():
  w=np.maximum(0,np.minimum(zz,row.end0)-np.maximum(aa,row.start0));bn,bp=full.sums(np.array([row.start0]),np.array([row.end0]))
  np.testing.assert_allclose(bn[0],np.sum(w[:,None]*v,axis=0),atol=1e-6,rtol=1e-9);assert abs(bp[0]-w.sum())<.01
  direct_checks.append(row.LAD_id)
 print('PROJECTED',c,len(aa),flush=True)

def means(tr,a,b):
 n,d=tr.sums(np.asarray(a),np.asarray(b))
 with np.errstate(divide='ignore',invalid='ignore'):v=n/d[...,None]
 return v,d
N=len(lads);body=np.full((N,84,D),np.nan);bbp=np.zeros((N,84))
whole=np.full((N,D),np.nan);flank=whole.copy()
x=np.r_[np.linspace(-1,0,18)[:-1]+1/34,np.linspace(0,1,51)[:-1]+.01,np.linspace(1,2,18)[:-1]+1/34]
for rid,row in enumerate(lads.itertuples()):
 c,a,b=row.chromosome,int(row.start0),int(row.end0);full,out=tracks[c]
 w,wd=means(full,[a],[b]);whole[rid]=w[0]
 fn,fd=out.sums(np.array([max(0,a-500000),b]),np.array([a,min(int(sizes[c]),b+500000)]))
 if fd.sum()>0:flank[rid]=fn.sum(axis=0)/fd.sum()
 edges=np.r_[np.linspace(a-500000,a,18),np.linspace(a,b,51)[1:],np.linspace(b,b+500000,18)[1:]]
 left=np.clip(edges[:-1],0,sizes[c]);right=np.clip(edges[1:],0,sizes[c])
 for ix,tr in [(np.arange(17,67),full),(np.r_[np.arange(17),np.arange(67,84)],out)]:body[rid,ix],bbp[rid,ix]=means(tr,left[ix],right[ix])

table=lads[['LAD_id','chromosome','start0','end0','B1_TTX_mean','length_bp','mCG_baseline_mean']].copy()
for j,name in enumerate(names):table[name]=whole[:,j];table[name+'_outside']=flank[:,j]
table['DNMT3A_mean']=whole[:,3:6].mean(axis=1)
table.to_csv(output/'LAD_values.tsv',sep='\t',index=False)
np.savez_compressed(output/'body_profile.npz',LAD_ids=lads.LAD_id.to_numpy(dtype=str),body=body[:,:,[0,3,4,5]],body_bp=bbp)
print('Calculated LAD enrichment and scaled-body profiles.')
