#!/usr/bin/env python3
"""
01_summarize_lad.py

Summarize total LAD bases and fraction of genome covered by LAD intervals.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--lad-map", required=True, help="CSV/TSV with columns: tissue,file")
    parser.add_argument("--genome-size", type=int, default=2725521370)
    parser.add_argument("--out", default="LAD_summary.txt")

    return parser.parse_args()


def read_mapping(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    return pd.read_csv(path, sep="\t")


def total_bed_bases(path: Path) -> int:
    total = 0

    with path.open() as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")

            if len(parts) < 3:
                continue

            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue

            total += max(0, end - start)

    return total


def main() -> None:
    args = parse_args()

    mapping = read_mapping(Path(args.lad_map))
    required = {"tissue", "file"}
    missing = required - set(mapping.columns)

    if missing:
        raise ValueError(f"LAD map missing columns: {sorted(missing)}")

    rows = []

    for _, row in mapping.iterrows():
        tissue = str(row["tissue"])
        bed_file = Path(str(row["file"]))

        total_lad_bases = total_bed_bases(bed_file)
        genome_fraction = total_lad_bases / args.genome_size

        rows.append(
            {
                "Tissue": tissue,
                "Total_LAD_bases": total_lad_bases,
                "Total_Genome_bases": args.genome_size,
                "Genome_fraction_in_LADs": genome_fraction,
            }
        )

    out = pd.DataFrame(rows)
    out.to_csv(args.out, sep="\t", index=False)

    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
