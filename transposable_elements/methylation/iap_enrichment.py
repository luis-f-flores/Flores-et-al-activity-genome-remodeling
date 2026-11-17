"""Fragment and 100-bp group hypo-DMR enrichment against matched genomic positions."""
import json
import numpy as np,pandas as pd
from pathlib import Path
import argparse
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'iap_matching')
ap.add_argument('--output',type=Path,required=True)
ap.add_argument('--fonts',type=Path,default=None)
args=ap.parse_args(); IN=args.inputs; OUT=args.output; OUT.mkdir(parents=True,exist_ok=True)
def write(p,d):
 p=Path(p);p.parent.mkdir(parents=True,exist_ok=True)
 if isinstance(d,pd.DataFrame):d.to_csv(p,sep='\t',index=False)
 else:p.write_text(json.dumps(d,indent=2,allow_nan=False)+'\n')
def chance_test(cand,units,cfg,out,spec):
 result=[];rng=np.random.default_rng(cfg['seed']+(30 if spec=='primary' else 40))
 for typ in ['fragment','group100']:
  target=[u for u in units if u['unit_type']==typ];counts=np.zeros(cfg['permutations'],int);info=[];draws=[]
  for u in target:
   p=cand[cand.unit_id==u['unit_id']];assert len(p)>0 and (len(p)>=cfg['final_min'] or u['unit_id'] in cfg['complete_pool_units']),(spec,u['unit_id'],len(p))
   y=p.hypo.to_numpy(int);idx=rng.integers(0,len(y),cfg['permutations']);counts+=y[idx];draws.append(idx)
   info.append(dict(unit_id=u['unit_id'],block=u['chrom']+':'+str(u['start0']//cfg['block_bp']),observed=int(u['hypo']),expected=float(y.mean()),placements=len(y)))
  info=pd.DataFrame(info);write(out/f'genomic_{spec}_{typ}_pools.tsv',info)
  write(out/f'genomic_{spec}_{typ}_null.tsv.gz',pd.DataFrame(dict(permutation=np.arange(1,cfg['permutations']+1),overlaps=counts,fraction=counts/len(target))))
  np.savez_compressed(out/f'genomic_{spec}_{typ}_draw_indices.npz',indices=np.asarray(draws,np.int32))
  blocks=info.groupby('block')[['observed','expected']].sum().to_numpy();boot=[]
  for _ in range(cfg['bootstrap']):
   a=blocks[rng.integers(0,len(blocks),len(blocks))].sum(axis=0);boot.append(a[0]/a[1] if a[1]>0 else np.nan)
  write(out/f'genomic_{spec}_{typ}_block_bootstrap.tsv.gz',pd.DataFrame(dict(fold=boot)))
  n=len(target);obs=int(info.observed.sum());exp=float(info.expected.sum());ci=np.nanquantile(boot,[.025,.975]);q=np.quantile(counts,[.025,.975])
  result.append(dict(analysis='genomic_chance',specification=spec,unit=typ,n=n,observed=obs,observed_pct=100*obs/n,expected=exp,expected_pct=100*exp/n,fold=obs/exp if exp>0 else np.nan,absolute_difference_pct=100*(obs-exp)/n,p_enrichment=(1+int((counts>=obs).sum()))/(cfg['permutations']+1),null_count_low=q[0],null_count_high=q[1],null_pct_low=100*q[0]/n,null_pct_high=100*q[1]/n,fold_ci_low=ci[0],fold_ci_high=ci[1],bootstrap_blocks=len(blocks),bootstrap_undefined=int(np.isnan(boot).sum()),minimum_pool=int(info.placements.min()),null_mc_mean=float(counts.mean()),null_mc_SE=float(counts.std()/np.sqrt(len(counts)))))
 return result

cfg=json.loads((IN/'settings.json').read_text())
cand=pd.read_csv(IN/'genomic_matched_candidates.tsv.gz',sep='	')
units=pd.read_csv(IN/'selected_genomic_units.tsv',sep='	').to_dict('records')
results=chance_test(cand,units,cfg,OUT,'primary')
write(OUT/'genomic_enrichment_results.tsv',pd.DataFrame(results))
print(pd.DataFrame(results)[['unit','n','observed','expected','fold','p_enrichment']].to_string(index=False))
