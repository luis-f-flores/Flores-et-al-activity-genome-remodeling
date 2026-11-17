#!/usr/bin/env python3
import os,sys,json,argparse,hashlib,subprocess,platform
from pathlib import Path
from collections import defaultdict
import numpy as np,pandas as pd,pyBigWig,scipy
import iap_utils as core

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
def start_domain(pos,offsets,lo,hi,span,target):
 # Half-open integer start domains. A CpG at p enters [x+a,x+b) at p-b+1
 # and exits at p-a+1. Disjoint union segments prevent double counting.
 ev=[np.array([lo,hi-span+1],dtype=np.int64)];delta=[np.zeros(2,dtype=np.int32)]
 for a,b in offsets:
  ev.extend([pos-b+1,pos-a+1]);delta.extend([np.ones(len(pos),np.int32),-np.ones(len(pos),np.int32)])
 ev=np.concatenate(ev);delta=np.concatenate(delta);order=np.argsort(ev,kind='stable');ev=ev[order];delta=delta[order]
 unique,ix=np.unique(ev,return_index=True);count=np.cumsum(np.add.reduceat(delta,ix))
 a=np.maximum(unique[:-1],lo);b=np.minimum(unique[1:],hi-span+1);tol=0 if target<=10 else .2*target
 good=(b>a)&(np.abs(count[:-1]-target)<=tol)
 a=a[good];b=b[good];return a,np.cumsum(b-a)
_cache=None
def matched_pool(unit,feat,cfg,rng):
 global _cache
 if _cache is None:
  old=pd.read_csv(cfg['reuse_geometry'],sep='\t')
  diags=pd.read_csv(cfg['reuse_diagnostics'],sep='\t').set_index('unit_id')
  sizes=pd.read_csv(cfg['reuse_final_sizes'],sep='\t').set_index('unit_id').placements
  _cache=(old,diags,sizes)
 old,diags,sizes=_cache
 prev=diags.loc[unit['unit_id']].to_dict()
 ch=unit['chrom'];origin=unit['start0'];span=unit['end0']-origin
 offsets=[(s-origin,e-origin)for s,e in unit['segments']]
 if sizes.loc[unit['unit_id']]>=cfg['final_min']:
  rows=old[old.unit_id==unit['unit_id']].drop(columns='parent_id').to_dict('records')
  for r in rows:r['segments']=[(int(r['start0'])+a,int(r['start0'])+b)for a,b in offsets]
  return rows,dict(unit_id=unit['unit_id'],**prev,reused_before_outcomes=True,exhaustive=False)
 block=cfg['block_bp'];scope=prev['scope']
 lo=origin//block*block if scope=='local20Mb' else 0
 hi=min(feat.sizes[ch],lo+block) if scope=='local20Mb' else feat.sizes[ch]
 p,k=core.bw_positions(feat.k,ch,lo,hi);q,m=core.bw_positions(feat.m,ch,lo,hi)
 pos=np.intersect1d(p[np.isfinite(k)],q[np.isfinite(m)],assume_unique=True)
 starts,cumulative=start_domain(pos,offsets,lo,hi,span,int(unit['shared']))
 del p,q,k,m,pos
 lens=np.diff(np.r_[0,cumulative]);ends=starts+lens
 total=int(cumulative[-1])if len(cumulative)else 0
 rows=[];tested=0;context_n=0;accepted=0
 # Scan every integer CpG-eligible start exactly once, chromosome blocks in order.
 # Context queries precede expensive dense mappability loading; no DMR labels used.
 for left in range(lo,hi,block):
  right=min(hi,left+block);aa=np.maximum(starts,left);bb=np.minimum(ends,right)
  keep=bb>aa;aa=aa[keep];bb=bb[keep]
  if not len(aa):continue
  for first in range(left,right,250000):
     a=np.maximum(aa,first);b=np.minimum(bb,first+250000);good=b>a;a=a[good];b=b[good]
     n=b-a
     if not len(n):continue
     vv=np.repeat(a,n)+np.arange(n.sum())-np.repeat(np.r_[0,np.cumsum(n)[:-1]],n);tested+=len(vv)
     ok=np.ones(len(vv),bool)
     for label,index in [('genic',feat.g),('exon',feat.x),('promoter',feat.pr)]:
      hits=np.zeros(len(vv),bool);idx=index.get(ch,core.Index([]))
      for x,y in offsets:hits|=idx.hit(vv+x,vv+y)
      ok&=hits==unit[label]
     vv=vv[ok];context_n+=len(vv)
     if not len(vv):continue
     loadhi=min(feat.sizes[ch],left+block+max(100000,span))
     if not(feat.key and feat.ch==ch and feat.lo==left and feat.hi>=loadhi):feat.load(ch,left,loadhi)
     f=aggregate_features(feat,vv,offsets);ok=np.ones(len(vv),bool)
     for label in ['genic','exon','promoter']:ok&=f[label]==unit[label]
     for label in ['mappability','map_recorded']:ok&=np.abs(f[label]-unit[label])<=.1
     ok&=(f['tested_fraction']==0)if unit['tested_fraction']==0 else np.abs(f['tested_fraction']-unit['tested_fraction'])<=.1
     tol=0 if unit['shared']<=10 else .2*unit['shared'];ok&=np.abs(f['shared']-unit['shared'])<=tol
     for j in np.flatnonzero(ok):
      start=int(vv[j]);accepted+=1
      row=dict(unit_id=unit['unit_id'],unit_type=unit['unit_type'],placement=-1,chrom=ch,start0=start,end0=start+span,segments=[(start+x,start+y)for x,y in offsets],**{key:float(value[j])for key,value in f.items()})
      # Uniform reservoir over the complete admissible geometric population.
      if len(rows)<cfg['exact_pool_cap']:rows.append(row)
      else:
       pick=int(rng.integers(0,accepted))
       if pick<len(rows):rows[pick]=row
 assert tested==total,('Incomplete exhaustive start scan',unit['unit_id'],tested,total)
 for i,r in enumerate(rows):r['placement']=i
 assert rows,('No matched geometric placements in complete allowed scope',unit['unit_id'])
 return rows,dict(unit_id=unit['unit_id'],scope=scope,proposals=total,pre_placements=len(rows),eligible_CpG_start_bases=total,context_eligible_start_bases=context_n,all_geometry_eligible_starts=accepted,reused_before_outcomes=False,exhaustive=True)

def chance_test(cand,units,cfg,out,spec):
 result=[];rng=np.random.default_rng(cfg['seed']+(30 if spec=='primary' else 40))
 for typ in ['fragment','group100']:
  target=[u for u in units if u['unit_type']==typ];counts=np.zeros(cfg['permutations'],int);info=[];draws=[]
  for u in target:
   p=cand[cand.unit_id==u['unit_id']];assert len(p)>0 and (len(p)>=cfg['final_min'] or u['unit_id'] in cfg['complete_pool_units']),(spec,u['unit_id'],len(p))
   y=p.hypo.to_numpy(int);idx=rng.integers(0,len(y),cfg['permutations']);counts+=y[idx];draws.append(idx)
   info.append(dict(unit_id=u['unit_id'],block=u['chrom']+':'+str(u['start0']//cfg['block_bp']),observed=int(u['hypo']),expected=float(y.mean()),placements=len(y)))
  info=pd.DataFrame(info);core.write(out/f'genomic_{spec}_{typ}_pools.tsv',info)
  core.write(out/f'genomic_{spec}_{typ}_null.tsv.gz',pd.DataFrame(dict(permutation=np.arange(1,cfg['permutations']+1),overlaps=counts,fraction=counts/len(target))))
  np.savez_compressed(out/f'genomic_{spec}_{typ}_draw_indices.npz',indices=np.asarray(draws,np.int32))
  blocks=info.groupby('block')[['observed','expected']].sum().to_numpy();boot=[]
  for _ in range(cfg['bootstrap']):
   a=blocks[rng.integers(0,len(blocks),len(blocks))].sum(axis=0);boot.append(a[0]/a[1] if a[1]>0 else np.nan)
  core.write(out/f'genomic_{spec}_{typ}_block_bootstrap.tsv.gz',pd.DataFrame(dict(fold=boot)))
  n=len(target);obs=int(info.observed.sum());exp=float(info.expected.sum());ci=np.nanquantile(boot,[.025,.975]);q=np.quantile(counts,[.025,.975])
  result.append(dict(analysis='genomic_chance',specification=spec,unit=typ,n=n,observed=obs,observed_pct=100*obs/n,expected=exp,expected_pct=100*exp/n,fold=obs/exp if exp>0 else np.nan,absolute_difference_pct=100*(obs-exp)/n,p_enrichment=(1+int((counts>=obs).sum()))/(cfg['permutations']+1),null_count_low=q[0],null_count_high=q[1],null_pct_low=100*q[0]/n,null_pct_high=100*q[1]/n,fold_ci_low=ci[0],fold_ci_high=ci[1],bootstrap_blocks=len(blocks),bootstrap_undefined=int(np.isnan(boot).sum()),minimum_pool=int(info.placements.min()),null_mc_mean=float(counts.mean()),null_mc_SE=float(counts.std()/np.sqrt(len(counts)))))
 return result
def main():
 ap=argparse.ArgumentParser();ap.add_argument('--settings');ap.add_argument('--output');ap.add_argument('--prepare-pools',action='store_true',help='Generate candidate geometry and depth-support tables before exhaustive completion.');a=ap.parse_args();cfg=json.load(open(a.settings));cfg={k:str((Path(a.settings).resolve().parent/Path(v)).resolve()) if isinstance(v,str) and k in ['context','rna','dmrs','tested_windows','genes','exons','promoters','kcl_bw','mann_bw','bismap','raw_reports','reuse_geometry','reuse_diagnostics','reuse_final_sizes'] else v for k,v in cfg.items()};out=Path(a.output).resolve();out.mkdir(parents=True,exist_ok=True);src=Path(__file__).parent
 for offsets in [[(0,3)],[(0,2),(4,6)]]:
  pos=np.array([1,4,7],dtype=np.int64);span=max(e for s,e in offsets)
  for n in [1,2]:
   left,cumulative=start_domain(pos,offsets,0,12,span,n);lens=np.diff(np.r_[0,cumulative]);got=set(x for aa,nn in zip(left,lens)for x in range(int(aa),int(aa+nn)))
   expected={x for x in range(13-span)if sum(sum(x+s<=p<x+e for p in pos)for s,e in offsets)==n}
   assert got==expected,('Conditional proposal boundary test',offsets,n,got,expected)
 raw=pd.read_csv(cfg['context'],sep='\t');d=raw[core.truth(raw.iap_related)&raw.shared_covered_positions.gt(0)].copy();d['RNA_valid']=np.isfinite(pd.to_numeric(d.rna_q_BH,errors='coerce'))&np.isfinite(pd.to_numeric(d.rna_log2fc_shrunk,errors='coerce'));core.write(out/'RNA_eligibility_all172.tsv',d[['locus_id','RNA_valid','rna_q_BH','shared_covered_positions']]);d=d[d.RNA_valid].reset_index(drop=True)
 d['up']=core.truth(d.rna_higher_q_lt_0_05);d['hypo']=core.truth(d.direct_6hr_hypo_dmr).astype(int);assert len(d)==171 and d.up.sum()==61 and d.loc[d.up,'hypo'].sum()==9
 masks,nw=core.load_masks(cfg['tested_windows']);feat=core.Features(cfg,masks);parents=[];units=[];raw_regions=[]
 def add_parent(ch,iv,label):
  pid=len(parents);parents.append(dict(parent_id=pid,label=label))
  for s,e in iv:raw_regions.append(dict(region_id=len(raw_regions),parent_id=pid,chrom=ch,start0=s,end0=e))
  return pid
 observed=[]
 for i,r in d.iterrows():
  q=measured(feat,r.chrom,[(int(r.start0),int(r.end0))]);assert q['shared']==r.shared_covered_positions
  observed.append(dict(**q,parent_id=add_parent(r.chrom,[(int(r.start0),int(r.end0))],r.locus_id)))
 f=pd.DataFrame(observed)
 for col in f:d[col]=f[col]
 core.write(out/'observed_fragment_opportunity_before_outcomes.tsv',d)
 selected=core.groups(d[d.up]);assert len(selected)==38 and sum(d.loc[ix,'hypo'].any() for ix in selected)==6
 for i,r in d[d.up].iterrows():units.append(dict(unit_id='F'+str(i),unit_type='fragment',chrom=r.chrom,start0=int(r.start0),end0=int(r.end0),segments=[(int(r.start0),int(r.end0))],hypo=int(r.hypo),obs_parent=int(r.parent_id),**{k:float(r[k]) for k in ['shared','baseline','mappability','map_recorded','tested_fraction','genic','exon','promoter','length_bp']}))
 for g,ix in enumerate(selected):
  dd=d.loc[ix];ch=dd.chrom.iloc[0];iv=core.union(zip(dd.start0,dd.end0));q=measured(feat,ch,iv)
  units.append(dict(unit_id='G'+str(g),unit_type='group100',chrom=ch,start0=int(dd.start0.min()),end0=int(dd.end0.max()),segments=iv,hypo=int(dd.hypo.any()),obs_parent=add_parent(ch,iv,'selected_group'+str(g)),**q))
 globalgroups=core.groups(d);grows=[];membership=[]
 for g,ix in enumerate(globalgroups):
  dd=d.loc[ix];ch=dd.chrom.iloc[0];iv=core.union(zip(dd.start0,dd.end0));q=measured(feat,ch,iv)
  grows.append(dict(unit_id='population_group'+str(g),chrom=ch,start0=int(dd.start0.min()),end0=int(dd.end0.max()),up=int(dd.up.any()),hypo=int(dd.hypo.any()),n_fragments=len(ix),mixed_up=int(dd.up.nunique()>1),parent_id=add_parent(ch,iv,'population_group'+str(g)),**q))
  for i in ix:membership.append(dict(unit_id='population_group'+str(g),locus_id=d.loc[i,'locus_id'],up=int(d.loc[i,'up']),hypo=int(d.loc[i,'hypo'])))
 candidates=[];diags=[];rng=np.random.default_rng(cfg['seed'])
 pool_builder=matched_pool
 if a.prepare_pools:
  from iap_candidate_pools import matched_pool as pool_builder
 for number,u in enumerate(sorted(units,key=lambda u:(u['chrom'],u['start0'],u['unit_type']))):
  rows,diag=pool_builder(u,feat,cfg,rng)
  for r in rows:r['parent_id']=add_parent(r['chrom'],r['segments'],r['unit_id']+':'+str(r['placement']))
  candidates.extend(rows);diags.append(diag);core.write(out/'pool_search_diagnostics.tsv',pd.DataFrame(diags))
  print('Unit pools',number+1,'of',len(units),diag,flush=True)
  # Candidate geometry saved without hypo-DMR labels; matching remains outcome-blind.
 cand=pd.DataFrame(candidates);core.write(out/'candidate_geometry_before_depth.tsv.gz',cand.drop(columns='segments'))
 paths={}
 for c in ['M','K']:
  for rep in 'ABC':
   files=list(Path(cfg['raw_reports']).glob(rep+'6'+c+'_*.CpG_report.txt'));assert len(files)==1;paths[rep+'6'+c]=files[0]
 rawinput={s:core.record(p,False) for s,p in paths.items()}
 regions=core.add_coverage(pd.DataFrame(raw_regions),paths,cfg,out,src)
 parent_ids=regions.parent_id.to_numpy(int);nr=len(parents);N=np.bincount(parent_ids,weights=regions.cytosines,minlength=nr);assert np.all(N>0)
 pcv=pd.DataFrame(dict(parent_id=np.arange(nr),cytosines=N))
 for c in ['M','K']:
  for key in ['depth_','coverage4_']:pcv[key+c]=np.bincount(parent_ids,weights=regions[key+c]*regions.cytosines,minlength=nr)/N
 d=d.merge(pcv,on='parent_id',validate='one_to_one');gd=pd.DataFrame(grows).merge(pcv,on='parent_id',validate='one_to_one');cand=cand.merge(pcv,on='parent_id',validate='one_to_one')
 target=pd.DataFrame(units).set_index('unit_id').join(pcv.set_index('parent_id'),on='obs_parent',rsuffix='_coverage')
 good=np.ones(len(cand),bool)
 for c in ['M','K']:
  td=cand.unit_id.map(target['depth_'+c]);good&=np.abs(np.log2(cand['depth_'+c]+1)-np.log2(td+1))<=1
  tf=cand.unit_id.map(target['coverage4_'+c]);good&=np.abs(cand['coverage4_'+c]-tf)<=.2
 cand=cand[good].copy();counts=cand.groupby('unit_id').size().reindex(target.index,fill_value=0);core.write(out/'final_pool_sizes.tsv',counts.rename('placements').reset_index())
 if a.prepare_pools:
  print('Prepared candidate geometry and measured depth support.',flush=True)
  return
 cfg['complete_pool_units']=[x['unit_id'] for x in diags if x.get('exhaustive',False) and x['pre_placements']==x['all_geometry_eligible_starts']]
 supported=(counts>=cfg['final_min'])|counts.index.isin(cfg['complete_pool_units'])
 assert counts.min()>0 and supported.all(),('Insufficient sampled opportunity support',counts[~supported].to_dict())
 dm=pd.read_csv(cfg['dmrs'],sep='\t');dmidx={ch:core.Index(zip(x.start.to_numpy(np.int64)-1,x.end.to_numpy(np.int64))) for ch,x in dm.groupby('chr')}
 def event(ch,iv):return any(bool(dmidx.get(ch,core.Index([])).hit(np.array([s]),np.array([e]))[0]) for s,e in iv)
 cand['hypo']=[int(event(r.chrom,r.segments)) for r in cand.itertuples()]
 for u in units:assert event(u['chrom'],u['segments'])==bool(u['hypo'])
 core.write(out/'genomic_matched_candidates.tsv.gz',cand.assign(segments=cand.segments.map(json.dumps)))
 core.write(out/'eligible_IAP_fragment_covariates.tsv',d);core.write(out/'eligible_IAP_group100_covariates.tsv',gd);core.write(out/'eligible_IAP_group100_membership.tsv',pd.DataFrame(membership));core.write(out/'selected_genomic_units.tsv',pd.DataFrame(units).assign(segments=lambda x:x.segments.map(json.dumps)))
 bal=[]
 for u in units:
  z=cand[cand.unit_id==u['unit_id']]
  for c in ['shared','mappability','tested_fraction','depth_M','depth_K','coverage4_M','coverage4_K']:
   base=target.loc[u['unit_id'],c];bal.append(dict(unit_id=u['unit_id'],unit_type=u['unit_type'],covariate=c,observed=float(base),expected=float(z[c].mean()),minimum=float(z[c].min()),maximum=float(z[c].max())))
 core.write(out/'genomic_covariate_balance.tsv',pd.DataFrame(bal))
 results=chance_test(cand,units,cfg,out,'primary')
 core.write(out/'genomic_enrichment_results.tsv',pd.DataFrame(results))
if __name__=='__main__':main()
