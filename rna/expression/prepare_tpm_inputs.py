"""Assemble generated counts, cluster order and activity membership for the TPM display."""
from pathlib import Path
import argparse
import pandas as pd

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    for name in ['counts','gene-order','activity-membership','output']:p.add_argument('--'+name,type=Path,required=True)
    a=p.parse_args();a.output.mkdir(parents=True,exist_ok=True)
    pd.read_csv(a.counts,sep='\t',float_precision='round_trip').to_csv(a.output/'TPM_gene_counts_lengths.tsv.gz',sep='\t',index=False)
    pd.read_csv(a.gene_order).to_csv(a.output/'TPM_gene_order.csv',index=False)
    pd.read_csv(a.activity_membership,sep='\t').to_csv(a.output/'TPM_activity_membership.tsv',sep='\t',index=False)
