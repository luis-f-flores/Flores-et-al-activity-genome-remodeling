"""Construct chromosome-specific propensity strata for the genic DMR/DEG comparison."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
import pyBigWig
from sklearn.preprocessing import SplineTransformer,StandardScaler
from sklearn.linear_model import LogisticRegression
def weights(groups,deg):
    w=np.zeros(len(deg)); pools=[]
    for g in np.unique(groups):
        ids=np.where(groups==g)[0]; k=int(deg[ids].sum())
        w[ids]=k/len(ids); pools.append((g,ids,k))
    return w,pools


def prepare(covariates,intervals,membership,bismap,coverage_dir,out):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    cv=pd.read_csv(covariates,sep='\t');members=pd.read_csv(membership,sep='\t')
    feat=pd.read_csv(intervals,sep='\t');feat=feat[feat.gene_symbol.isin(cv.gene_symbol)].copy()
    for g,d in feat.groupby('gene_symbol'):
        assert d.chr.nunique()==1
        a=d.sort_values('start');assert np.all(a.start.to_numpy()[1:]>a.end.to_numpy()[:-1])
    bw=pyBigWig.open(str(bismap));chromsizes=bw.chroms();mapp={};loci={};covered={}
    for g,d in feat.groupby('gene_symbol',sort=False):
        total=0.;length=0;cover=0.
        for r in d.itertuples():
            start=int(r.start)-1;end=int(r.end);ln=end-start
            assert r.chr in chromsizes and 0<=start<end<=chromsizes[r.chr]
            mean=bw.stats(r.chr,start,end,type='mean',exact=True)[0]
            cov=bw.stats(r.chr,start,end,type='coverage',exact=True)[0]
            total+=(mean or 0.)*(cov or 0.)*ln;cover+=(cov or 0.)*ln;length+=ln
        mapp[g]=total/length;covered[g]=cover/length
        loci[g]=((d.start.min()+d.end.max())/2)/chromsizes[d.chr.iloc[0]]
    bw.close()
    cv['bismap100_mean']=cv.gene_symbol.map(mapp);cv['bismap_recorded_fraction']=cv.gene_symbol.map(covered)
    cv['chromosome_relative_midpoint']=cv.gene_symbol.map(loci)
    assert np.isfinite(cv[['bismap100_mean','chromosome_relative_midpoint']]).all().all()
    cv.to_csv(out/'gene_reference_covariates.tsv',sep='\t',index=False)
    first=pd.read_csv(Path(coverage_dir)/'A1M_coverage.tsv',sep='\t').set_index('gene_symbol')
    missing=first.index[first.cytosines==0].tolist()
    pd.DataFrame({'gene_symbol':missing,'reason':'No CpG-report positions; coverage undefined'}).to_csv(out/'coverage_exclusions.tsv',sep='\t',index=False)
    excluded=members[members.gene_symbol.isin(missing)]
    assert not excluded.is_DEG.any() and not excluded.is_promoter_body_DMR.any()
    cv=cv[~cv.gene_symbol.isin(missing)].reset_index(drop=True);symbols=cv.gene_symbol.tolist()
    for tp in [1,6,24]:
        tables=[pd.read_csv(Path(coverage_dir)/f'{r}{tp}M_coverage.tsv',sep='\t').set_index('gene_symbol').loc[symbols] for r in 'ABC']
        for a in tables:assert np.array_equal(a.cytosines,tables[0].cytosines)
        denominator=tables[0].cytosines.to_numpy();assert np.all(denominator>0)
        depth=np.mean([a.depth_sum.to_numpy()/denominator for a in tables],axis=0)
        frac=np.mean([a.coverage4.to_numpy()/denominator for a in tables],axis=0)
        present=np.mean([a.covered.to_numpy()/denominator for a in tables],axis=0)
        X=pd.DataFrame({
            'log_gene_length':np.log1p(cv.gene_length),
            'log_CpG_density':np.log1p(cv.true_cg_density_per_kb),
            'log_baseline_RNA':np.log1p(cv.baseline_fpkm),
            'log_control_depth':np.log1p(depth),
            'control_fraction_coverage4':frac,
            'control_fraction_covered':present,
            'bismap100_mean':cv.bismap100_mean,
            'chromosome_relative_midpoint':cv.chromosome_relative_midpoint})
        assert np.isfinite(X).all().all()
        X.assign(gene_symbol=symbols).to_csv(out/f'{tp}h_covariates.tsv',sep='\t',index=False)
        m=members[members.timepoint==f'{tp}hr'].set_index('gene_symbol').loc[symbols]
        assert m.permutation_eligible.all()
        deg=m.is_DEG.to_numpy(bool);y=m.is_promoter_body_DMR.to_numpy(bool);ch=cv.chr.to_numpy()
        # Outcome labels do not enter propensity fitting or stratum construction.
        blocks=[]
        for col in X:
            values=X[[col]].to_numpy();knots=np.unique(np.quantile(values,[0,1/3,2/3,1]))
            if len(knots)>1:
                blocks.append(SplineTransformer(degree=3,knots=knots[:,None],include_bias=False).fit_transform(values))
        Z=StandardScaler().fit_transform(np.column_stack(blocks))
        Z=np.c_[Z,pd.get_dummies(cv.chr,dtype=float).to_numpy()]
        fit=LogisticRegression(C=1.0,max_iter=3000,tol=1e-8,solver='lbfgs').fit(Z,deg)
        assert np.max(fit.n_iter_)<3000
        ps=fit.predict_proba(Z)[:,1];groups=np.empty(len(cv),dtype=object)
        for chrom in np.unique(ch):
            ids=np.where(ch==chrom)[0];order=ids[np.argsort(ps[ids],kind='stable')]
            for j,block in enumerate(np.array_split(order,20)):groups[block]=f'{chrom}_{j:02d}'
        ew,epools=weights(groups.astype(str),deg)
        account=[dict(stratum=group,N=len(ids),n=k,DMR_genes=int(y[ids].sum())) for group,ids,k in epools]
        pd.DataFrame(account).to_csv(out/f'{tp}h_strata.tsv',sep='\t',index=False)
        pd.DataFrame({'gene_symbol':symbols,'DEG':deg,'DMR':y,'propensity':ps,'stratum':groups,'expected_inclusion':ew}).to_csv(out/f'{tp}h_membership.tsv',sep='\t',index=False)

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['covariates','intervals','membership','bismap','coverage-dir','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();prepare(a.covariates,a.intervals,a.membership,a.bismap,a.coverage_dir,a.output)
