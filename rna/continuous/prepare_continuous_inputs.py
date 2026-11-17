"""Assemble the generated gene counts and chromatin cohorts for continuous RNA analyses."""
from pathlib import Path
import argparse
import shutil

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['counts','methylation','b1','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    groups=[(a.counts,['counts_KCl.tsv','counts_Bic.tsv']),
            (a.methylation,['promoter_cohort.tsv','body_cohorts.tsv']),
            (a.b1,['baseline_genes_states.tsv','S17_links.tsv','S17_regions.tsv'])]
    for root,names in groups:
        for name in names:
            src=root/name
            if not src.is_file():raise FileNotFoundError(src)
            shutil.copyfile(src,a.output/name)
