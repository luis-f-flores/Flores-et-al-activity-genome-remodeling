"""Prepare public-LAD gene intervals from the supplied UCSC RefSeq coordinate table."""
from pathlib import Path
import argparse
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--coordinates', type=Path, required=True)
    p.add_argument('--output', type=Path, required=True)
    a = p.parse_args()
    d = pd.read_csv(a.coordinates, sep='\t', header=None,
                    names=['chrom', 'strand', 'geneStart', 'geneEnd', 'name'])
    d['geneStart'] = d.geneStart.astype(int)
    d['geneEnd'] = d.geneEnd.astype(int)
    d = d[~d.chrom.str.contains('_') & d.chrom.str.startswith('chr')]
    a.output.parent.mkdir(parents=True, exist_ok=True)
    d[['chrom', 'geneStart', 'geneEnd', 'name']].to_csv(a.output, sep='\t', header=False, index=False)


if __name__ == '__main__':
    main()
