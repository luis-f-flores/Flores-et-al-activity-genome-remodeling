"""Continuous RNA endpoints from complete counts and fixed published cohorts."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from scipy import stats
ap=argparse.ArgumentParser(description=__doc__)
ap.add_argument('--inputs',type=Path,default=Path(__import__('os').environ.get('GB_DATA_ROOT', '../inputs'))/'continuous_rna')
ap.add_argument('--output',type=Path,default=Path(__file__).resolve().parent/'output')
args=ap.parse_args();I=args.inputs;O=args.output;O.mkdir(parents=True,exist_ok=True)
def read(name):return pd.read_csv(I/name,sep='\t',float_precision='round_trip')
def tpm(file,groups):
    d=read(file).set_index('gene_id');samples=sum(groups.values(),[])
    assert d.index.is_unique and len(d)==54532 and (d.gene_length>0).all()
    assert np.isfinite(d[samples]).all().all() and (d[samples]>=0).all().all()
    rpk=d[samples].divide(d.gene_length/1000,axis=0)
    v=rpk.divide(rpk.sum(axis=0),axis=1)*1e6
    assert np.allclose(v.sum(),1e6,rtol=0,atol=1e-7)
    return pd.DataFrame({k:v[s].mean(axis=1) for k,s in groups.items()})
groups={f'{t}{c}':[f'{r}{t}{c}' for r in 'ABC'] for t in [1,6,24] for c in ['K','M']}
km=tpm('counts_KCl.tsv',groups)
ke=pd.DataFrame({t:np.log2((km[f'{t}K']+1)/(km[f'{t}M']+1)) for t in [1,6,24]})
bg={'TTX':['Basal_A','Basal_B','Basal_C'],'Bic':['Bic_A','Bic_B','Bic_C'],'Bic_Nimo':['Bic_nimo_A','Bic_nimo_B','Bic_nimo_C']}
bm=tpm('counts_Bic.tsv',bg)
cons=['Bic_vs_TTX','Bic_Nimo_vs_TTX','Bic_Nimo_vs_Bic'];pairs=[('TTX','Bic'),('TTX','Bic_Nimo'),('Bic','Bic_Nimo')]
be=pd.DataFrame({c:np.log2((bm[t]+1)/(bm[b]+1)) for c,(b,t) in zip(cons,pairs)})
assert np.max(abs(be.iloc[:,1]-be.iloc[:,0]-be.iloc[:,2]))<1e-14
results=[]
def correlation(panel,comparison,x,y,cohort='all'):
    x=np.asarray(x);y=np.asarray(y);assert np.isfinite(x).all() and np.isfinite(y).all()
    p=stats.pearsonr(x,y);s=stats.spearmanr(x,y);z=(x-min(x))/(max(x)-min(x))-.5;fit=stats.linregress(z,y)
    row=dict(panel=panel,comparison=str(comparison),cohort=cohort,n=len(x),pearson_r=p.statistic,pearson_P=p.pvalue,R2=p.statistic**2,spearman_rho=s.statistic,scaled_OLS_slope=fit.slope,scaled_OLS_SE=fit.stderr)
    results.append(row);return row
prom=read('promoter_cohort.tsv');body=read('body_cohorts.tsv')
for t in [1,6,24]:
    q=prom[prom.timepoint==t];correlation('Fig_6',t,q.delta_mCG_percent,ke.loc[q.gene_id,t])
    for g in ['differentially_methylated','DEG_subset']:
        q=body[(body.timepoint==t)&(body.cohort==g)];correlation('Fig_S14',t,q.delta_mCG_percent,ke.loc[q.gene_id,t],g)
links=read('S17_links.tsv');reg=read('S17_regions.tsv');regions=[]
for c in cons:
    q=links[links.contrast==c].copy();q['RNA_change']=be.loc[q.gene_id,c].to_numpy();med=q.groupby('region_index').RNA_change.median()
    z=reg[reg.contrast==c].copy();z['RNA_change']=med.loc[z.region_index].to_numpy();assert len(z)==1920
    correlation('Fig_S09_LADs',c,z.B1_delta,z.RNA_change);regions.append(z)
pd.DataFrame(results).to_csv(O/'correlations.tsv',sep='\t',index=False)
prom['RNA_change']=[ke.loc[r.gene_id,r.timepoint] for r in prom.itertuples()]
body['RNA_change']=[ke.loc[r.gene_id,r.timepoint] for r in body.itertuples()]
prom.to_csv(O/'promoter_points.tsv',sep='\t',index=False)
body.to_csv(O/'gene_body_points.tsv',sep='\t',index=False)
pd.concat(regions).to_csv(O/'lad_points.tsv',sep='\t',index=False)
bm.to_csv(O/'Bic_mean_TPM.tsv',sep='\t')
base=read('baseline_genes_states.tsv');base['mean_TTX_TPM']=bm.loc[base.gene_id,'TTX'].to_numpy();base['log2_TTX_TPM_plus1']=np.log2(base.mean_TTX_TPM+1)
base.to_csv(O/'baseline_TPM.tsv',sep='\t',index=False)
base.groupby('state').agg(n=('gene_id','size'),median_TPM=('mean_TTX_TPM','median'),median_log2_TPM_plus1=('log2_TTX_TPM_plus1','median')).to_csv(O/'baseline_summary.tsv',sep='\t')
# One observation per region, in the same order for all contrasts.
ids=regions[0].region_index.to_numpy();chrom=regions[0].chromosome.to_numpy();pos=(regions[0].start0.to_numpy()+regions[0].end0.to_numpy())/2
assert all(np.array_equal(z.region_index,ids) for z in regions)
X=np.column_stack([stats.rankdata(z.B1_delta) for z in regions]);Y=np.column_stack([stats.rankdata(z.RNA_change) for z in regions]);X-=X.mean(0);Y-=Y.mean(0)
den=np.sqrt((X*X).sum(0)*(Y*Y).sum(0));observed=(X*Y).sum(0)/den
spatial=[]
for width in [5000000]:
    blocks=[]
    for ch in dict.fromkeys(chrom):
        ii=np.flatnonzero(chrom==ch);label=(pos[ii]//width).astype(int);blocks.append([ii[label==k] for k in dict.fromkeys(label)])
    rng=np.random.default_rng(20261012+width);null=np.empty((10000,3))
    for k in range(10000):
        order=np.concatenate([np.concatenate([b[j] for j in rng.permutation(len(b))]) for b in blocks]);null[k]=(X*Y[order]).sum(0)/den
    for j,c in enumerate(cons):
        lower=int((null[:,j]<=observed[j]+1e-14).sum());upper=int((null[:,j]>=observed[j]-1e-14).sum());p=min(1,2*(min(lower,upper)+1)/10001)
        spatial.append(dict(contrast=c,block_bp=width,rho=observed[j],P=p,draws=10000,method='within_chromosome_block_permutation'))
pd.DataFrame(spatial).to_csv(O/'spatial_tests.tsv',sep='\t',index=False)
print('Computed TPM correlations and regressions.')
