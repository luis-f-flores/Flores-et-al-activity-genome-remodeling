"""Calculate conditional promoter-DMG/DEG matching from tested regions and baseline measurements."""

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)

import argparse, hashlib, json, math, platform
from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import scipy
from scipy.spatial.distance import cdist
from scipy.stats import hypergeom, ks_2samp

def sha(p):
    h=hashlib.sha256()
    with Path(p).open("rb") as f:
        for b in iter(lambda:f.read(8*1024*1024),b""): h.update(b)
    return h.hexdigest()

def read(p, **kw): return pd.read_csv(p,sep="\t",**kw)

def merge(v):
    out=[]
    for s,e in sorted(set(v)):
        if out and s<=out[-1][1]+1: out[-1][1]=max(out[-1][1],e)
        else: out.append([s,e])
    return out

def main(config_path, out_path, source_config=None, union_tables=None, reports=None):
    c=json.loads(Path(config_path).read_text())
    out=Path(out_path)
    out.mkdir(parents=True,exist_ok=True)
    cfg_path=Path(source_config)
    cfg=json.loads(cfg_path.read_text())
    root=(cfg_path.parent/str(cfg.get('package_root','.'))).resolve()
    resolve=lambda p: (root/p).resolve()
    old=Path(union_tables)
    inputs=[cfg_path,Path(__file__).with_name('promoter_gene_id_overrides.tsv'),
            old/'21_original_promoter_intervals.tsv',old/'promoter_union_union_gene_membership.tsv',
            resolve(cfg['master_fpkm'])]
    aliases=read(inputs[1]); alias=dict(zip(aliases.dmg_symbol,aliases.rna_gene_name))
    canon=lambda x:alias.get(x,x)
    membership=read(old/"promoter_union_union_gene_membership.tsv")
    genes=membership.loc[membership.in_union_test].copy().set_index("gene_symbol")
    assert genes.index.is_unique and len(genes)>0
    prom=read(old/"21_original_promoter_intervals.tsv")
    prom["gene_symbol"]=prom.gene_symbol.map(canon)
    prom=prom.loc[prom.gene_symbol.isin(genes.index)].drop_duplicates()
    genes["n_promoters"]=prom.groupby("gene_symbol").size().reindex(genes.index,fill_value=0)
    rna=read(resolve(cfg['master_fpkm']))
    rna["gene_symbol"]=rna.gene_name.map(canon)

    base=rna.groupby("gene_symbol")[["A0","B0","C0"]].mean().mean(axis=1)
    genes["baseline_RNA_FPKM"]=base.reindex(genes.index)
    genes["log_RNA"]=np.log1p(genes.baseline_RNA_FPKM)
    genes["testability"]=""
    windows=defaultdict(lambda:defaultdict(set))
    for tp in ("1","6","24"):
        d=cfg["timepoints"][tp]
        wp=resolve(d["tested_windows"]); rp=resolve(d["rna"])
        inputs.extend([wp,rp])
        w=read(wp,usecols=["chr","start","end","pvalue","qvalue"])
        w=w.loc[np.isfinite(w.pvalue)&np.isfinite(w.qvalue)].drop_duplicates(["chr","start","end"])
        nwin=defaultdict(int)
        for chrom,ps in prom.groupby("chr"):
            wc=w.loc[w.chr.eq(chrom)].sort_values("start")
            starts=wc.start.to_numpy(); ends=wc.end.to_numpy()
            assert np.all(np.diff(ends)>=0)
            for row in ps.itertuples():
                lo=np.searchsorted(ends,row.start,side="left")
                hi=np.searchsorted(starts,row.end,side="right")
                pairs=list(zip(starts[lo:hi].tolist(),ends[lo:hi].tolist()))
                windows[row.gene_symbol][chrom].update(pairs)
            
        wt=defaultdict(set)
        for chrom,ps in prom.groupby("chr"):
            wc=w.loc[w.chr.eq(chrom)].sort_values("start")
            starts=wc.start.to_numpy(); ends=wc.end.to_numpy()
            for row in ps.itertuples():
                lo=np.searchsorted(ends,row.start); hi=np.searchsorted(starts,row.end,side="right")
                wt[row.gene_symbol].update((chrom,int(s),int(e)) for s,e in zip(starts[lo:hi],ends[lo:hi]))
        genes["windows_"+tp]=[len(wt[g]) for g in genes.index]
        r=read(rp)
        valid=np.isfinite(r.padj)&np.isfinite(r.log2FoldChange)&r.gene_name.notna()&~r.gene_name.isin(["","-"])
        rset=set(r.loc[valid,"gene_name"].map(canon))
        genes["testability"] += [str(int(g in rset))+str(int(len(wt[g])>0)) for g in genes.index]
        print("Prepared tested-window opportunity at",tp,"h",flush=True)
    genes["n_tested_windows"]=[sum(len(v) for v in windows[g].values()) for g in genes.index]
    genes["log_windows"]=np.log1p(genes.n_tested_windows)
    intervals=defaultdict(list); all_ranges=defaultdict(list)
    for g in genes.index:
        for chrom,win in windows[g].items():
            for s,e in merge(win):
                intervals[chrom].append((g,s,e)); all_ranges[chrom].append((s,e))
    all_ranges={ch:np.array(merge(v),dtype=np.int64) for ch,v in all_ranges.items()}
    interval_rows=[(g,ch,s,e) for ch,v in intervals.items() for g,s,e in v]
    pd.DataFrame(interval_rows,columns=["gene_symbol","chr","start_1based","end_inclusive"]).to_csv(
        out/"baseline_candidate_regions.tsv",sep="\t",index=False)
    for sample,report_path in zip(('A0','B0','C0'),reports):
        report=Path(report_path); inputs.append(report)
        retained=defaultdict(list)
        for chunk in pd.read_csv(report,sep="\t",header=None,usecols=[0,1,3,4],
                names=["chr","pos","strand","m","u","context","trinucleotide"],chunksize=1000000,
                dtype={"chr":str,"pos":np.int64,"m":np.int64,"u":np.int64}):
            chunk=chunk.loc[(chunk.m+chunk.u)>0]
            for chrom,d in chunk.groupby("chr",sort=False):
                if chrom not in all_ranges: continue
                a=all_ranges[chrom]; pos=d.pos.to_numpy()
                j=np.searchsorted(a[:,0],pos,side="right")-1
                keep=(j>=0)&(pos<=a[np.maximum(j,0),1])
                if keep.any(): retained[chrom].append(d.loc[keep,["pos","m","u"]].to_numpy())
        agg=defaultdict(lambda:np.zeros(3,dtype=np.int64))
        for chrom,iv in intervals.items():
            if not retained[chrom]: continue
            a=np.concatenate(retained[chrom]); a=a[np.argsort(a[:,0])]
            assert len(np.unique(a[:,0]))==len(a),"Duplicate cytosine coordinate"
            pref=np.vstack([np.zeros(2,dtype=np.int64),np.cumsum(a[:,1:],axis=0)])
            for g,s,e in iv:
                lo=np.searchsorted(a[:,0],s); hi=np.searchsorted(a[:,0],e,side="right")
                m,u=pref[hi]-pref[lo]; agg[g] += [m,u,hi-lo]
        vals=np.array([agg[g] for g in genes.index]); total=vals[:,0]+vals[:,1]
        genes[sample+"_mCG"]=np.divide(vals[:,0],total,out=np.full(len(genes),np.nan),where=total>0)
        genes[sample+"_cytosines"]=vals[:,2]
        genes[sample+"_depth"]=np.divide(total,vals[:,2],out=np.full(len(genes),np.nan),where=vals[:,2]>0)
        print("Baseline aggregation finished:",sample,flush=True)
    genes["baseline_mCG"]=genes[["A0_mCG","B0_mCG","C0_mCG"]].mean(axis=1,skipna=False)
    genes["log_cytosines"]=np.log1p(genes[["A0_cytosines","B0_cytosines","C0_cytosines"]].mean(axis=1))
    genes["log_depth"]=np.log1p(genes[["A0_depth","B0_depth","C0_depth"]].mean(axis=1,skipna=False))
    metrics=["log_RNA","log_windows","baseline_mCG","log_cytosines","log_depth"]
    valid=np.isfinite(genes[metrics]).all(axis=1)&genes.n_tested_windows.gt(0)
    genes["eligible_for_matching"]=valid
    genes.to_csv(out/"gene_covariates_and_membership.tsv",sep="\t")
    df=genes.loc[valid].copy(); fg=np.flatnonzero(df.is_DMG.to_numpy())
    y=df.is_DEG.to_numpy(dtype=int); n=len(df); obs=int(y[fg].sum())
    expected=len(fg)*y.mean()
    result={"original":{"universe":18986,"DMGs":202,"DEGs":6817,"overlap":84},
        "eligible_unmatched":{"universe":n,"DMGs":len(fg),"DEGs":int(y.sum()),"overlap":obs,
            "expected":expected,"fold":obs/expected,"p":float(hypergeom.sf(obs-1,n,int(y.sum()),len(fg)))},
        "excluded_DMGs":genes.index[~valid & genes.is_DMG].tolist(),
        "excluded_genes":int((~valid).sum()),"tests":[]}
    allz=(df[metrics]-df[metrics].mean())/df[metrics].std(ddof=0)
    allz.to_csv(out/"standardized_covariates.tsv",sep="\t")
    allz.corr(method="spearman").to_csv(out/"covariate_spearman.tsv",sep="\t")
    
    
    for name,cols in [("primary",metrics[:3])]:
        z=allz[cols].to_numpy()
        pools=[]; support=[]
        for i in fg:
            eligible=np.flatnonzero(df.testability.to_numpy()==df.testability.iloc[i])
            distance=np.sqrt(((z[eligible]-z[i])**2).sum(axis=1))
            ranked=eligible[np.argsort(distance,kind="stable")]
            pool=ranked[:min(20,len(ranked))]
            pools.append(pool)
            support.append({"gene_symbol":df.index[i],"pool_n":len(pool),
                "max_standardized_distance":float(np.linalg.norm(z[pool]-z[i],axis=1).max()),
                "max_single_covariate_difference":float(np.abs(z[pool]-z[i]).max()),
                "includes_self":bool(i in pool)})
        pd.DataFrame(support).to_csv(out/(name+"_support.tsv"),sep="\t",index=False)
        preview=np.mean([allz.iloc[p].mean().to_numpy() for p in pools],axis=0)
        pd.DataFrame({"covariate":metrics,"foreground_mean":allz.iloc[fg].mean(),
            "candidate_pool_mean":preview,
            "approximate_SMD_before_draws":allz.iloc[fg].mean().to_numpy()-preview}).to_csv(
                out/(name+"_pre_draw_balance.tsv"),sep="\t",index=False)
        rng=np.random.default_rng(c["seed"]); draws=c["draws"]
        counts=np.empty(draws,dtype=int); nullmean=np.zeros(len(metrics)); freq=np.zeros(n,dtype=int)
        pool_order=np.argsort([len(p) for p in pools],kind="stable")
        for b in range(draws):
            selected=set()
            
            order=pool_order.copy()
            rng.shuffle(order)
            for j in order:
                available=[int(x) for x in pools[j] if int(x) not in selected]
                if not available: raise RuntimeError("Matching pools cannot provide a unique full-size draw")
                selected.add(int(rng.choice(available)))
            ix=np.fromiter(selected,dtype=int)
            assert len(ix)==len(fg)
            counts[b]=int(y[ix].sum()); freq[ix]+=1
            nullmean+=allz.iloc[ix].mean().to_numpy()
        nullmean/=draws
        smd=allz.iloc[fg].mean().to_numpy()-nullmean
        balance=pd.DataFrame({"covariate":metrics,"foreground_standardized_mean":allz.iloc[fg].mean(),
            "null_standardized_mean":nullmean,"SMD":smd})
        balance.to_csv(out/(name+"_balance.tsv"),sep="\t",index=False)
        pd.DataFrame({"gene_symbol":df.index,"selection_fraction":freq/draws}).to_csv(
            out/(name+"_null_gene_frequencies.tsv"),sep="\t",index=False)
        np.savetxt(out/(name+"_null_overlap_counts.tsv"),counts,fmt="%d",header="overlap",comments="")
        r={"name":name,"matched_covariates":cols,"draws":draws,"seed":c["seed"],
           "observed":obs,"null_mean":float(counts.mean()),"fold":float(obs/counts.mean()),
           "p_greater":float((1+(counts>=obs).sum())/(draws+1)),
           "null_95_range":np.quantile(counts,[.025,.975]).tolist(),
           "max_abs_SMD_matched":float(balance.loc[balance.covariate.isin(cols),"SMD"].abs().max()),
           "max_abs_SMD_all":float(abs(smd).max()),
           "foreground_pool_min":min(x["pool_n"] for x in support),
           "foreground_pool_distance_max":max(x["max_standardized_distance"] for x in support),
           "balance_pass":bool(balance.loc[balance.covariate.isin(cols),"SMD"].abs().max()<=.1)}
        result["tests"].append(r); print(json.dumps(r),flush=True)
    result["definitions"]={
        "endpoint":'Overlap between promoter-DMG and DEG gene unions across the sampled timepoints.',
        "baseline_regions":"Union of full finite-p/q tested 500bp windows touching upstream-1kb transcript promoters; overlapping windows merged per gene, including nonsignificant windows",
        "baseline_mCG":"Mean of three baseline replicate read-weighted m/(m+u), over covered cytosines in those regions; no whole-promoter effect interpretation",
        "RNA":"Mean A0/B0/C0 FPKM; duplicate canonical gene symbols averaged before log1p",
        "matching":"20 nearest eligible genes within exact six-bit RNA/WGBS timepoint testability pattern; Euclidean standardized covariates; foreground eligible in null; each draw unique genes"}
    (out/"RESULTS.json").write_text(json.dumps(result,indent=2)+"\n")
    manifest={"status":"completed","inputs":{str(p):sha(p) for p in dict.fromkeys(inputs)},
        "code_sha256":sha(__file__),"config":c,"python":platform.python_version(),
        "versions":{"numpy":np.__version__,"pandas":pd.__version__,"scipy":scipy.__version__},
        "results":result,"outputs":{p.name:sha(p) for p in out.iterdir() if p.is_file()}}
    (out/"RUN_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")

if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--config',type=Path,default=Path(__file__).with_name('promoter_matching_config.json'))
    ap.add_argument('--source-config',type=Path,required=True,help='Configuration used for promoter_union.py, with source input paths.')
    ap.add_argument('--union-tables',type=Path,required=True)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--reports',type=Path,nargs=3,help='Baseline A0, B0 and C0 CpG reports.')
    group.add_argument('--samples',type=Path,help='Sample/report manifest containing A0, B0 and C0.')
    ap.add_argument('--output',type=Path,required=True)
    a=ap.parse_args()
    if a.samples:
        samples=read(a.samples).set_index('sample')
        if not samples.index.is_unique:raise ValueError('Sample identifiers must be unique.')
        a.reports=[Path(samples.loc[s,'report']) for s in ['A0','B0','C0']]
        a.reports=[p if p.is_absolute() else a.samples.parent/p for p in a.reports]
    main(a.config,a.output,a.source_config,a.union_tables,a.reports)
