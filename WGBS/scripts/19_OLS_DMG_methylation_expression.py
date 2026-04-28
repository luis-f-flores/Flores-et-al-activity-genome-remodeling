#!/usr/bin/env python3

# ==============================================================================
# Script: 19_OLS_DMG_methylation_expression.py
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
19_OLS_DMG_methylation_expression.py

Fit ordinary least squares (OLS) regression models to test the relationship
between promoter methylation changes and gene expression changes for DMGs.

Model: log2FC_fpkm ~ meth_diff_scaled

Outputs:
    - Model summary statistics
    - Residuals
    - Scatter plot with regression line
    - Summary table across all timepoints

Usage:
    python3 19_OLS_DMG_methylation_expression.py \
        --input-dir DMG_scaled_tables \
        --pattern 'DMG_meth25_q01_with_FPKM_*hr_scaled.tsv' \
        --outdir OLS_results
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
import statsmodels.api as sm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Fit OLS models for methylation-expression correlation"
    )

    parser.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing scaled DMG tables",
    )

    parser.add_argument(
        "--pattern",
        default="DMG_meth25_q01_with_FPKM_*hr_scaled.tsv",
        help="File pattern to match (default: DMG_meth25_q01_with_FPKM_*hr_scaled.tsv)",
    )

    parser.add_argument(
        "--outdir",
        required=True,
        help="Output directory for OLS results",
    )

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    output_root = Path(args.outdir)
    output_root.mkdir(parents=True, exist_ok=True)

    files = sorted(input_dir.glob(args.pattern))

    if not files:
        raise FileNotFoundError(
            f"No files matched pattern '{args.pattern}' in {input_dir}"
        )

    summary_rows = []

    for file_path in files:
        # Extract timepoint from filename
        timepoint = (
            file_path.stem.replace("DMG_meth25_q01_with_FPKM_", "")
            .replace("_scaled", "")
        )

        # Create timepoint-specific output directory
        output_dir = output_root / timepoint
        output_dir.mkdir(parents=True, exist_ok=True)

        # Load data
        df = pd.read_csv(file_path, sep="\t")

        # Extract columns needed for regression
        subset = df[["gene_name", "log2FC_fpkm", "meth_diff_scaled"]].dropna()

        if len(subset) < 3:
            print(f"[skip] {timepoint}: fewer than 3 valid rows for regression")
            continue

        # Prepare regression variables
        y = subset["log2FC_fpkm"].astype(float)
        X = sm.add_constant(
            subset[["meth_diff_scaled"]].astype(float),
            has_constant="add",
        )

        # Fit OLS model
        model = sm.OLS(y, X).fit()

        # Save model summary
        summary_path = output_dir / "M0_meth_only_summary.txt"
        summary_path.write_text(model.summary().as_text())

        # Calculate and save residuals
        predictions = model.predict(X)
        residuals_df = pd.DataFrame(
            {
                "gene_name": subset["gene_name"],
                "log2FC_fpkm": y,
                "pred": predictions,
                "residual": y - predictions,
            }
        )

        residuals_path = output_dir / "M0_meth_only_residuals.tsv"
        residuals_df.to_csv(residuals_path, sep="\t", index=False)

        # Create scatter plot with regression line
        fig, ax = plt.subplots(figsize=(8, 6))

        ax.scatter(
            subset["meth_diff_scaled"],
            y,
            alpha=0.6,
            s=20,
            label="DMGs",
        )

        ax.plot(
            subset["meth_diff_scaled"],
            predictions,
            color="red",
            linewidth=2,
            label=f"OLS fit (R² = {model.rsquared:.3f})",
        )

        ax.set_xlabel("Methylation difference (scaled)")
        ax.set_ylabel("Gene expression change (log2FC FPKM)")
        ax.set_title(f"{timepoint}: Methylation-Expression Correlation")
        ax.legend()
        ax.grid(True, alpha=0.3)

        plt.tight_layout()

        plt.savefig(output_dir / "M0_meth_only_feature_fit.pdf")
        plt.savefig(output_dir / "M0_meth_only_feature_fit.png", dpi=300)
        plt.close()

        # Collect summary statistics
        summary_rows.append(
            {
                "timepoint": timepoint,
                "n": int(model.nobs),
                "r_squared": model.rsquared,
                "slope": model.params.get("meth_diff_scaled", float("nan")),
                "p_value": model.pvalues.get("meth_diff_scaled", float("nan")),
            }
        )

        print(f"[done] {timepoint}: n={model.nobs}, R²={model.rsquared:.3f}")

    # Save summary table across all timepoints
    summary_df = pd.DataFrame(summary_rows)
    summary_path = output_root / "meth_only_models_summary_meth25_q01.tsv"
    summary_df.to_csv(summary_path, sep="\t", index=False)

    print(f"\n[summary] {summary_path}")


if __name__ == "__main__":
    main()
