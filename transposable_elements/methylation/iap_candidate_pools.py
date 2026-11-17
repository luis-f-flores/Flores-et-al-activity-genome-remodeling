"""Construct IAP genomic-matching candidate pools without using DMR outcomes."""
import numpy as np,pandas as pd
import iap_utils as core
from iap_matching import start_domain
def measured(feat,ch,iv):
 L=sum(e-s for s,e in iv);N=0;ms=0.;mp=0.;mc=0.;tf=0.;genic=False;exon=False;prom=False
 for s,e in iv:
  k=np.asarray(feat.k.values(ch,s,e,numpy=True),float);m=np.asarray(feat.m.values(ch,s,e,numpy=True),float);g=np.isfinite(k)&np.isfinite(m)
  N+=int(g.sum());ms+=float(m[g].sum());z=np.asarray(feat.map.values(ch,s,e,numpy=True),float);v=np.isfinite(z);mp+=float(np.nan_to_num(z,nan=0.).sum());mc+=int(v.sum())
  tf+=float(feat.masks.get(ch,core.Mask([])).fraction(np.array([s]),np.array([e]))[0])*(e-s)
  genic|=bool(feat.g.get(ch,core.Index([])).hit(np.array([s]),np.array([e]))[0]);exon|=bool(feat.x.get(ch,core.Index([])).hit(np.array([s]),np.array([e]))[0]);prom|=bool(feat.pr.get(ch,core.Index([])).hit(np.array([s]),np.array([e]))[0])
 return dict(shared=N,baseline=ms/N if N else np.nan,mappability=mp/L,map_recorded=mc/L,tested_fraction=tf/L,genic=int(genic),exon=int(exon),promoter=int(prom),length_bp=L)

def aggregate_features(feat,starts,offsets):
 L=sum(e-s for s,e in offsets);n=len(starts)
 out={k:np.zeros(n) for k in ['shared','mappability','map_recorded','tested_fraction','genic','exon','promoter','baseline']}
 ms=np.zeros(n)
 for s,e in offsets:
  q=feat.calc(starts+s,starts+e);out['shared']+=q['shared'];ms+=np.nan_to_num(q['baseline'],nan=0.)*q['shared']
  for k in ['mappability','map_recorded','tested_fraction']:out[k]+=(e-s)*q[k]
  out['genic']=np.maximum(out['genic'],q['genic']);out['exon']=np.maximum(out['exon'],q['context']==0);out['promoter']=np.maximum(out['promoter'],q['promoter'])
 for k in ['mappability','map_recorded','tested_fraction']:out[k]/=L
 out['baseline']=np.divide(ms,out['shared'],out=np.full(n,np.nan),where=out['shared']>0)
 return out

def matched_pool(unit,feat,cfg,rng):
 ch=unit['chrom'];origin=unit['start0'];span=unit['end0']-origin;offsets=[(s-origin,e-origin) for s,e in unit['segments']]
 block=cfg['block_bp'];bl=origin//block*block;bh=min(feat.sizes[ch],bl+block)
 for scope,lo,hi,cap in [('local20Mb',bl,bh,cfg['max_proposals']),('chromosome',0,feat.sizes[ch],5000000)]:
  p,k=core.bw_positions(feat.k,ch,lo,hi);q,m=core.bw_positions(feat.m,ch,lo,hi);pos=np.intersect1d(p[np.isfinite(k)],q[np.isfinite(m)],assume_unique=True)
  starts,cumulative=start_domain(pos,offsets,lo,hi,span,int(unit['shared']));del p,q,k,m,pos
  total=int(cumulative[-1]) if len(cumulative) else 0
  chosen={};proposed=0
  while len(chosen)<cfg['pre_pool'] and proposed<cap and total:
   count=(cap-proposed) if scope=='chromosome' else min(50000,cap-proposed);u=rng.integers(0,total,count);j=np.searchsorted(cumulative,u,'right');v=starts[j]+u-np.r_[0,cumulative][j];proposed+=count
   v=np.unique(v);v=v[~np.isin(v,list(chosen))];batch=[]
   for b in np.unique(v//block):
    vv=v[v//block==b];left=int(b)*block;right=min(feat.sizes[ch],left+block+max(100000,span))
    if not(feat.key and feat.ch==ch and feat.lo==left and feat.hi>=right):feat.load(ch,left,right)
    f=aggregate_features(feat,vv,offsets);ok=np.ones(len(vv),bool)
    for key in ['genic','exon','promoter']:ok&=f[key]==unit[key]
    for key in ['mappability','map_recorded']:ok&=np.abs(f[key]-unit[key])<=.1
    ok&=(f['tested_fraction']==0) if unit['tested_fraction']==0 else np.abs(f['tested_fraction']-unit['tested_fraction'])<=.1
    ids=np.flatnonzero(ok)
    for i in ids:batch.append((int(vv[i]),{k:float(v[i]) for k,v in f.items()}))
   rng.shuffle(batch)
   for start,f in batch[:cfg['pre_pool']-len(chosen)]:chosen[start]=f
  if len(chosen)>=cfg['pre_min'] or scope=='chromosome':break
 rows=[]
 for i,(start,f) in enumerate(chosen.items()):rows.append(dict(unit_id=unit['unit_id'],unit_type=unit['unit_type'],placement=i,chrom=ch,start0=start,end0=start+span,segments=[(start+a,start+b) for a,b in offsets],**f))
 assert len(rows)>0,('No matches even after exact conditional-CpG sampling',unit['unit_id'],ch,origin)
 return rows,dict(unit_id=unit['unit_id'],scope=scope,proposals=proposed,pre_placements=len(rows),eligible_CpG_start_bases=total)

