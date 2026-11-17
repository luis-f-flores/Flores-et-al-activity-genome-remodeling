#!/usr/bin/env python3
"""Interval, methylation-opportunity and coverage utilities for IAP genomic matching."""
import argparse,json,hashlib,os,sys,platform,subprocess,gzip,csv,time
from pathlib import Path
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
import pyBigWig
import scipy

def sha(p):
 h=hashlib.sha256()
 with open(p,'rb') as f:
  for b in iter(lambda:f.read(8*1024**2),b''):h.update(b)
 return h.hexdigest()
def record(p,full=True):
 p=Path(p);s=p.stat();d=dict(path=str(p),resolved=str(p.resolve()),bytes=s.st_size,mtime_ns=s.st_mtime_ns)
 if full:d['sha256']=sha(p)
 return d
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if isinstance(d,pd.DataFrame):d.to_csv(p,sep='\t',index=False)
 else:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def truth(s):return s.astype(str).str.lower().eq('true')
def union(intervals):
 out=[]
 for a,b in sorted(intervals):
  if out and a<=out[-1][1]:out[-1][1]=max(out[-1][1],int(b))
  else:out.append([int(a),int(b)])
 return out
def groups(d):
 out=[];ch=None;end=-1
 for i,r in d.sort_values(['chrom','start0','end0']).iterrows():
  if r.chrom!=ch or r.start0>end+100:out.append([i]);ch=r.chrom;end=r.end0
  else:out[-1].append(i);end=max(end,r.end0)
 return out
class Index:
 def __init__(self,iv):
  a=np.asarray(sorted(iv),dtype=np.int64).reshape(-1,2)
  self.s=a[:,0];self.e=a[:,1];self.mx=np.maximum.accumulate(self.e) if len(a) else np.array([],dtype=np.int64)
 def hit(self,s,e):
  j=np.searchsorted(self.s,e,'left')-1
  return (j>=0)&(self.mx[np.maximum(j,0)]>s) if len(self.s) else np.zeros(np.shape(s),bool)
 def contains(self,s,e):
  j=np.searchsorted(self.s,s,'right')-1
  return (j>=0)&(self.mx[np.maximum(j,0)]>=e) if len(self.s) else np.zeros(np.shape(s),bool)
class Mask:
 def __init__(self,iv):
  a=np.asarray(union(iv),dtype=np.int64).reshape(-1,2);self.s=a[:,0];self.e=a[:,1];self.pre=np.r_[0,np.cumsum(self.e-self.s)]
 def cdf(self,x):
  j=np.searchsorted(self.s,x,'right')-1
  return np.where(j>=0,self.pre[np.maximum(j,0)]+np.maximum(0,np.minimum(x-self.s[np.maximum(j,0)],self.e[np.maximum(j,0)]-self.s[np.maximum(j,0)])),0) if len(self.s) else np.zeros(np.shape(x))
 def fraction(self,s,e):return (self.cdf(e)-self.cdf(s))/(e-s)
def bed_index(path):
 iv=defaultdict(list)
 with open(path) as f:
  for l in f:
   z=l.rstrip().split('\t');iv[z[0]].append((int(z[1]),int(z[2])))
 return {k:Index(v) for k,v in iv.items()}
def bw_positions(bw,ch,lo,hi):
 v=bw.intervals(ch,lo,hi)
 if not v:return np.array([],np.int64),np.array([],float)
 a=np.asarray(v);s=np.maximum(a[:,0].astype(np.int64),lo);e=np.minimum(a[:,1].astype(np.int64),hi);n=e-s
 p=np.repeat(s,n)+np.arange(n.sum())-np.repeat(np.r_[0,np.cumsum(n)[:-1]],n)
 return p,np.repeat(a[:,2],n)
class Features:
 def __init__(self,cfg,masks):
  self.k=pyBigWig.open(cfg['kcl_bw']);self.m=pyBigWig.open(cfg['mann_bw']);self.map=pyBigWig.open(cfg['bismap'])
  self.sizes=self.k.chroms();self.masks=masks;self.g=bed_index(cfg['genes']);self.x=bed_index(cfg['exons']);self.pr=bed_index(cfg['promoters']);self.key=None
 def load(self,ch,lo,hi):
  hi=min(self.sizes[ch],hi)
  if self.key==(ch,lo,hi):return
  p,k=bw_positions(self.k,ch,lo,hi);q,m=bw_positions(self.m,ch,lo,hi)
  self.pos,ki,mi=np.intersect1d(p,q,assume_unique=True,return_indices=True)
  self.mpre=np.r_[0,np.cumsum(m[mi])]
  z=np.asarray(self.map.values(ch,lo,hi,numpy=True),dtype=np.float64)
  good=np.isfinite(z);self.map_pre=np.r_[0,np.cumsum(np.where(good,z,0))];self.map_cov=np.r_[0,np.cumsum(good,dtype=np.int64)]
  self.lo=lo;self.hi=hi;self.ch=ch;self.key=(ch,lo,hi)
 def calc(self,s,e):
  s=np.atleast_1d(s).astype(np.int64);e=np.atleast_1d(e).astype(np.int64)
  assert np.all(s>=self.lo) and np.all(e<=self.hi) and np.all(e>s)
  a=np.searchsorted(self.pos,s);b=np.searchsorted(self.pos,e);n=b-a
  gene=self.g.get(self.ch,Index([]));exon=self.x.get(self.ch,Index([]));prom=self.pr.get(self.ch,Index([]))
  gh=gene.hit(s,e);xh=exon.hit(s,e);contained=gene.contains(s,e)
  context=np.where(xh,0,np.where(~gh,3,np.where(contained,1,2)))
  L=e-s;aa=s-self.lo;bb=e-self.lo
  return dict(shared=n,baseline=np.divide(self.mpre[b]-self.mpre[a],n,out=np.full(len(n),np.nan),where=n>0),context=context,genic=gh.astype(int),promoter=prom.hit(s,e).astype(int),mappability=(self.map_pre[bb]-self.map_pre[aa])/L,map_recorded=(self.map_cov[bb]-self.map_cov[aa])/L,tested_fraction=self.masks.get(self.ch,Mask([])).fraction(s,e))
 def point(self,ch,s,e,block):
  lo=int(s//block)*block;hi=max(lo+block,int(e));self.load(ch,lo,hi);return {k:v[0].item() for k,v in self.calc([s],[e]).items()}
def load_masks(path):
 iv=defaultdict(list);n=0
 for d in pd.read_csv(path,sep='\t',usecols=['chr','start','end','pvalue'],chunksize=500000):
  d=d[np.isfinite(pd.to_numeric(d.pvalue,errors='coerce'))];n+=len(d)
  for ch,x in d.groupby('chr',sort=False):iv[ch].extend(zip(x.start.to_numpy(np.int64)-1,x.end.to_numpy(np.int64)))
 masks={ch:Mask(v) for ch,v in iv.items()};print('Assessed window mask',n,'rows',flush=True);return masks,n
def add_coverage(regions,paths,cfg,out,source_dir):
 bed=out/'coverage_regions.bed';regions[['chrom','start0','end0','region_id']].to_csv(bed,sep='\t',header=False,index=False)
 binary=out/'coverage_sum';subprocess.run(['/usr/bin/g++','-O3','-std=c++17',str(source_dir/'cpg_interval_coverage.cpp'),'-o',str(binary)],check=True)
 # Boundary test includes zero-depth cytosines, excluded end coordinate and nested intervals.
 toy=out/'toy_report.tsv';toy.write_text('chr1\t2\t+\t2\t3\tCG\tCGA\nchr1\t3\t-\t0\t0\tCG\tCGT\nchr1\t8\t+\t1\t0\tCG\tCGA\n')
 tb=out/'toy_regions.bed';tb.write_text('chr1\t1\t2\t0\nchr1\t1\t7\t1\nchr1\t0\t8\t2\n')
 tt=out/'toy_coverage.tsv';subprocess.run([str(binary),str(toy),str(tb),str(tt)],check=True)
 t=pd.read_csv(tt,sep='\t');assert t.cytosines.tolist()==[1,2,3] and t.depth_sum.tolist()==[5,5,6]
 def job(item):
  sample,path=item;dest=out/f'{sample}_interval_coverage.tsv'
  a=subprocess.run([str(binary),str(path),str(bed),str(dest)],capture_output=True,text=True,check=True)
  (out/f'{sample}_coverage.log').write_text(a.stderr)
  print('Coverage complete',sample,flush=True);return sample,pd.read_csv(dest,sep='\t').sort_values('region_id')
 with ThreadPoolExecutor(max_workers=2) as pool:tables=dict(pool.map(job,paths.items()))
 n=tables['A6M'].cytosines.to_numpy();assert np.all(n>=0), 'Invalid reference CpG count'
 for sample,t in tables.items():assert np.array_equal(n,t.cytosines),sample
 for condition in ['M','K']:
  regions['depth_'+condition]=np.mean([np.divide(tables[r+'6'+condition].depth_sum.to_numpy(),n,out=np.zeros(len(n)),where=n>0) for r in 'ABC'],axis=0)
  regions['coverage4_'+condition]=np.mean([np.divide(tables[r+'6'+condition].ge4.to_numpy(),n,out=np.zeros(len(n)),where=n>0) for r in 'ABC'],axis=0)
 regions['cytosines']=n
 return regions
