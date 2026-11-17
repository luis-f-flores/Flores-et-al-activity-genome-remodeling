"""Summarize B1 and RNA responses to nimodipine within fixed bicuculline DEG sets."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    x = pd.read_csv(args.input, sep='\t', low_memory=False)
    if not x.gene_id.is_unique:
        raise ValueError('Response table must contain one row per gene.')
    cols = ['RNA_initial_DEG', 'RNA_direct_DEG', 'RNA_toward', 'B1_toward', 'LAD_Bic_vs_TTX']
    if not all(x[c].dtype == bool for c in cols):
        raise ValueError('DEG, LAD and response-direction fields must be Boolean.')
    d = x[x.RNA_initial_DEG].copy()
    records=[]
    for scope,y in [('All_jointly_assessable_Bic_DEGs',d),('Pairwise_LAD_Bic_DEGs',d[d.LAD_Bic_vs_TTX])]:
     for direction,z in [('All',y),('Up',y[y.fc_Bic_vs_TTX>0]),('Down',y[y.fc_Bic_vs_TTX<0])]:
      records.append(dict(scope=scope,direction=direction,n=len(z),RNA_toward_TTX=int(z.RNA_toward.sum()),
       B1_toward_TTX=int(z.B1_toward.sum()),both_toward_TTX=int((z.RNA_toward&z.B1_toward).sum()),
       RNA_only_toward=int((z.RNA_toward&~z.B1_toward).sum()),B1_only_toward=int((~z.RNA_toward&z.B1_toward).sum()),
       neither_toward=int((~z.RNA_toward&~z.B1_toward).sum()),
       direct_RNA_DEG_and_toward=int((z.RNA_toward&z.RNA_direct_DEG).sum()),
       both_toward_and_direct_RNA_DEG=int((z.RNA_toward&z.B1_toward&z.RNA_direct_DEG).sum())))
    summary=pd.DataFrame(records);summary.to_csv(args.output/'fixed_Bic_DEG_return_summary_125.tsv',sep='\t',index=False)
    d[['gene_id','gene_name','class_Bic_vs_TTX']+cols+['b','bn','n','fc_Bic_vs_TTX','fc_Bic_Nimo_vs_TTX','fc_Bic_Nimo_vs_Bic','q_Bic_Nimo_vs_Bic']].to_csv(args.output/'fixed_Bic_DEG_return_details_125.tsv',sep='\t',index=False)
    print(summary.to_string(index=False))


if __name__ == '__main__':
    main()
