"""Convert coverage summaries and RNA-cluster GO results to plotting tables."""
from pathlib import Path
import argparse
import pandas as pd

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--coverage',type=Path);p.add_argument('--go-dir',type=Path);p.add_argument('--output',type=Path,required=True)
 a=p.parse_args()
 if a.coverage is None and a.go_dir is None:p.error('Provide a coverage summary or GO results directory.')
 a.output.mkdir(parents=True,exist_ok=True)
 if a.coverage is not None:
  d=pd.read_csv(a.coverage)
  names=dict(Sample='sample',Min_Coverage='min',Max_Coverage='max',Mean_Coverage='mean',Median_Coverage='median',Percent_5X='pct5',Percent_10X='pct10',Percent_20X='pct20')
  d.rename(columns=names)[list(names.values())].to_csv(a.output/'coverage.tsv',sep='\t',index=False)
 if a.go_dir is not None:
  for cl in range(1,6):
   d=pd.read_excel(a.go_dir/f'Cluster_{cl}/GO_Biological_Process_Top15.xlsx')
   d.to_csv(a.output/f'Cluster_{cl}_GO_Biological_Process_Top15.xlsx.csv',index=False)

if __name__=='__main__':main()
