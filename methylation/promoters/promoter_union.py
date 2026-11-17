#!/usr/bin/env python3
"""Promoter-DMG/DEG union, assessed background and timing calculations."""

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)

from pathlib import Path
import argparse, hashlib, json, math, subprocess, platform
import pandas as pd
import numpy as np
import scipy
from scipy.stats import hypergeom, fisher_exact, chisquare
def sha(p):
 h=hashlib.sha256()
 with p.open("rb") as f:
  for block in iter(lambda:f.read(2**20),b""):h.update(block)
 return h.hexdigest()
def read(p):return pd.read_csv(p,sep="\t")
def main():
 ap=argparse.ArgumentParser(description=__doc__)
 ap.add_argument("--config",type=Path,default=Path(__file__).with_name("promoter_cohort_config.json"))
 ap.add_argument("--output",type=Path,required=True)
 a=ap.parse_args();c=json.loads(a.config.read_text())
 root=(a.config.parent/c["package_root"]).resolve();resolve=lambda p:(root/p).resolve()
 out=a.output.resolve();td=out/"tables";logs=out/"logs"
 if (out/"RUN_MANIFEST.json").exists():raise FileExistsError("Completed run exists; choose a new output.")
 td.mkdir(parents=True,exist_ok=True);logs.mkdir(exist_ok=True)
 tps=["1","6","24"];helper=Path(__file__).with_name("promoter_testability.R")
 paths={tp:{k:resolve(v) for k,v in c["timepoints"][tp].items()} for tp in tps}
 aliases=Path(__file__).with_name("promoter_gene_id_overrides.tsv");features=resolve(c["promoter_features"])
 inputs=[Path(__file__).resolve(),a.config.resolve(),helper.resolve(),aliases,features]
 for d in paths.values():inputs.extend(d.values())
 before={str(p):sha(p) for p in inputs}
 cmd=[c.get("rscript","Rscript"),str(helper),str(features),str(td)]+[str(paths[tp]["tested_windows"]) for tp in tps]
 r=subprocess.run(cmd,text=True,capture_output=True)
 (logs/"assessed_promoters.log").write_text(r.stdout+r.stderr)
 if r.returncode:raise RuntimeError(r.stderr)
 cross=read(aliases);assert cross.dmg_symbol.is_unique
 alias=dict(zip(cross.dmg_symbol,cross.rna_gene_name));canon=lambda x:alias.get(x,x)
 tested=read(td/"21_promoter_testability.tsv")
 tested["canonical_symbol"]=tested.gene_symbol.map(canon)
 valid_w=tested.wgbs_promoter_tested.astype(str).str.lower().isin(["true","1"])
 rs={};ws={};ds={};ms={};raw_ms={};counts=[]
 for tp in tps:
  rna=read(paths[tp]["rna"]);assert rna.gene_id.is_unique
  good=np.isfinite(rna.padj)&np.isfinite(rna.log2FoldChange)&rna.gene_name.notna()&~rna.gene_name.isin(["","-"])
  sig=good&(rna.padj<=c["deg_padj"])&(rna.log2FoldChange.abs()>=math.log2(c["deg_fold_change"]))
  rs[tp]=set(rna.loc[good,"gene_name"].map(canon));ds[tp]=set(rna.loc[sig,"gene_name"].map(canon))
  m=pd.concat([read(paths[tp]["hyper"]),read(paths[tp]["hypo"])])
  assert m.gene_name.is_unique and (m.qvalue<=c["dmg_qvalue"]).all() and (m.meth_diff.abs()>=c["dmg_delta"]).all()
  raw_ms[tp]=set(m.gene_name);ms[tp]=set(m.gene_name.map(canon));assert len(ms[tp])==len(m)
  ws[tp]=set(tested.loc[valid_w&tested.timepoint.eq(tp+"hr"),"canonical_symbol"])
  assert ms[tp]<=ws[tp],"DMG absent from assessed promoter set"
  counts.append(dict(timepoint=tp+"hr",RNA_rows=len(rna),RNA_symbols=len(rs[tp]),DEGs=len(ds[tp]),DMGs=len(ms[tp])))
 R=set.union(*rs.values());W=set.union(*ws.values());D=set.union(*ds.values());M=set.union(*ms.values())
 U=R&W;Uall=set.intersection(*rs.values())&set.intersection(*ws.values())
 def stats(stage,background):
  d=D&background;m=M&background;k=len(d&m);n=len(background)
  cells=[k,len(m)-k,len(d)-k,n-len(m)-len(d)+k];assert min(cells)>=0
  p=float(hypergeom.sf(k-1,n,len(d),len(m)))
  fp=float(fisher_exact([cells[:2],cells[2:]],alternative="greater")[1]);assert abs(p-fp)<1e-12
  expected=len(d)*len(m)/n
  return dict(stage=stage,background=n,DEGs=len(d),DMGs=len(m),shared=k,expected=expected,fold_enrichment=k/expected,p_enrichment=p,fisher_p=fp)
 summary=[stats("assessed_at_any_time",U)]
 rawM=set.union(*raw_ms.values());gene_rows=[]
 for g in sorted(D&M&U):
  dt=next(tp for tp in tps if g in ds[tp]);mt=next(tp for tp in tps if g in ms[tp])
  category="expression-first" if tps.index(dt)<tps.index(mt) else "methylation-first" if tps.index(dt)>tps.index(mt) else "coincident"
  gene_rows.append(dict(gene_symbol=g,original_DMG_names=";".join(sorted(x for x in rawM if canon(x)==g)),first_DEG=dt+"hr",first_DMG=mt+"hr",timing_class=category,assessed_in_both_assays_at_all_times=g in Uall))
 genes=pd.DataFrame(gene_rows);assert genes.assessed_in_both_assays_at_all_times.all(),"First-significance interpretation requires assessment at each sampled time"
 timing=genes.groupby("timing_class").size().rename("n").reset_index()
 timing["percentage"]=100*timing.n/timing.n.sum()
 chi=float(chisquare(timing.n)[1])
 primary=summary[0]
 display=dict(all_DEGs=len(D),all_DMGs=len(M),shared=primary["shared"],DEG_only=len(D)-primary["shared"],DMG_only=len(M)-primary["shared"],background=primary["background"],eligible_DEGs=primary["DEGs"],eligible_DMGs=primary["DMGs"],expected=primary["expected"],fold_enrichment=primary["fold_enrichment"],p_enrichment=primary["p_enrichment"],timing_chisquare_P=chi)
 pd.DataFrame(summary).to_csv(td/"promoter_union_background_check.tsv",sep="\t",index=False)
 genes.to_csv(td/"promoter_union_gene_timing_checked.tsv",sep="\t",index=False)
 timing.to_csv(td/"promoter_union_timing_counts_checked.tsv",sep="\t",index=False)
 pd.DataFrame(counts).to_csv(td/"promoter_union_source_counts.tsv",sep="\t",index=False)
 membership=[dict(gene_symbol=g,RNA_any=g in R,WGBS_any=g in W,in_union_test=g in U,is_DEG=g in D,is_DMG=g in M) for g in sorted(R|M)]
 pd.DataFrame(membership).to_csv(td/"promoter_union_union_gene_membership.tsv",sep="\t",index=False)
 pd.DataFrame([dict(original_DMG_symbol=g,canonical_symbol=canon(g),reason="no_eligible_RNA_comparison_at_any_time" if canon(g) not in R else "no_tested_WGBS_promoter") for g in sorted(rawM) if canon(g) not in U]).to_csv(td/"promoter_union_DMG_exclusions.tsv",sep="\t",index=False)
 (td/"plot_values.json").write_text(json.dumps(dict(overlap=display,timing=timing.to_dict("records")),indent=2)+"\n")
 assert all(sha(Path(p))==h for p,h in before.items()),"Input changed during run"
 manifest=dict(scope="Promoter-DMG union and sampled significance timing",definitions=dict(primary="Unique gene symbols with eligible RNA statistics at any sampled time and a tested WGBS promoter at any sampled time; foreground unions restricted to this same background.",timing_test="Pearson chi-square goodness-of-fit against equal frequencies in the three timing categories"),inputs=before,result=summary,display=display,timing=timing.to_dict("records"),checks=dict(hypergeom_equals_fisher=True,all_shared_genes_assessed_at_all_three_times=True,unchanged_significance_thresholds=True,input_hashes_unchanged=True),versions=dict(python=platform.python_version(),pandas=pd.__version__,numpy=np.__version__,scipy=scipy.__version__),outputs={str(p.relative_to(out)):sha(p) for p in td.iterdir() if p.is_file()})
 (out/"RUN_MANIFEST.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print(json.dumps(display,indent=2))
if __name__=="__main__":main()

