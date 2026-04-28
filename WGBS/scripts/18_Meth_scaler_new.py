#!/usr/bin/env python3

# ==============================================================================
# Script: 18_Meth_scaler_new.py
# Package: WGBS / methylation analysis scripts
#
# Purpose:
#   Cleaned author-generated scripts for WGBS DMC/DMR/DMG analysis, annotation,
#   overlap analysis, methylation-expression integration, signal tracks, and
#   OLS modeling.
#
# Notes:
#   - Paths are supplied by command-line arguments where possible.
#   - Standard external tools/packages are described in the Methods.
#   - This script documents the manuscript-facing analysis logic.
# ==============================================================================

"""
18_Meth_scaler_new.py

Scale methylation difference values to [-0.5, 0.5] range for DMG tables.

This enables direct comparison of methylation changes across different
magnitude ranges in methylation-expression correlation analyses.

Usage:
    python3 18_Meth_scaler_new.py \
        --input-dir DMG_tables \
        --pattern 'DMG_meth25_q01_with_FPKM_*hr.tsv' \
        --outdir DMG_scaled_tables
"""

from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd
from sklearn.preprocessing import minmax_scale


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Scale methylation difference values in DMG tables"
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing DMG tables with FPKM values",
    )

    parser.add_argument(
        "--pattern",
        default="DMG_meth25_q01_with_FPKM_*hr.tsv",
        help="File pattern to match (default: DMG_meth25_q01_with_FPKM_*hr.tsv)",
    )

    parser.add_argument(
        "--outdir",
        default=None,
        help="Output directory (default: same as input-dir)",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_dir = Path(args.outdir) if args.outdir else input_dir
    output_dir.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(args.pattern))

    if not files:
        raise FileNotFoundError(
            f"No files matched pattern '{args.pattern}' in {input_dir}"
        )

    for file_path in files:
        df = pd.read_csv(file_path, sep="\t")

        # Find methylation difference and expression columns (case-insensitive)
        meth_col = next(
            (c for c in df.columns if c.lower() == "meth_diff"),
            None,
        )

        expr_col = next(
            (c for c in df.columns if c.lower() == "log2fc_fpkm"),
            None,
        )

        if meth_col is None or expr_col is None:
            print(f"[skip] {file_path.name}: missing meth_diff or log2FC_fpkm column")
            continue

        # Filter to rows with both methylation and expression data
        df_filtered = df[df[expr_col].notna() & df[meth_col].notna()].copy()

        if len(df_filtered) < 2:
            print(f"[skip] {file_path.name}: fewer than 2 valid rows")
            continue

        # Scale methylation difference to [-0.5, 0.5]
        df_filtered["meth_diff_scaled"] = minmax_scale(
            df_filtered[meth_col].astype(float),
            feature_range=(-0.5, 0.5),
        )

        # Write output
        output_path = output_dir / f"{file_path.stem}_scaled.tsv"
        df_filtered.to_csv(output_path, sep="\t", index=False)

        print(f"[done] {output_path}")


if __name__ == "__main__":
    main()
