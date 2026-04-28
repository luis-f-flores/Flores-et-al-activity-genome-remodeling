#!/usr/bin/env python3
"""
00_process_gtf_coordinates.py

Convert a simplified mm10 gene-coordinate table into BED format.

Expected input columns, no header:
    chrom, strand, geneStart, geneEnd, name

Output:
    cleaned_gene_coordinates.bed
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input", required=True, help="Gene coordinate table.")
    parser.add_argument("--output", default="cleaned_gene_coordinates.bed")

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_path = Path(args.input)
    output_path = Path(args.output)

    genes = pd.read_csv(
        input_path,
        sep="\t",
        header=None,
        names=["chrom", "strand", "geneStart", "geneEnd", "name"],
    )

    genes["geneStart"] = genes["geneStart"].astype(int)
    genes["geneEnd"] = genes["geneEnd"].astype(int)

    genes = genes[~genes["chrom"].str.contains("_", regex=False)]
    genes = genes[genes["chrom"].str.startswith("chr")]

    bed = genes[["chrom", "geneStart", "geneEnd", "name"]].copy()
    bed.columns = ["chrom", "start", "end", "name"]

    bed.to_csv(output_path, sep="\t", header=False, index=False)

    print(f"Wrote: {output_path}")


if __name__ == "__main__":
    main()
