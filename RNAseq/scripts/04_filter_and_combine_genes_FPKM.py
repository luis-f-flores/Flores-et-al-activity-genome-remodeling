#!/usr/bin/env python3
"""
04_filter_and_combine_genes_FPKM.py

Build a unified RNA-seq trajectory matrix for the union of DEGs across timepoints.

DEGs are selected using:
    padj <= --padj-threshold
    abs(log2FoldChange) >= log2(--fc-threshold)

For every DEG in the global union, FPKM-derived trajectory values are computed as:
    log2((FPKM_KCl + pseudocount) / (FPKM_Mannitol + pseudocount))
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd


def read_table(path: Path) -> pd.DataFrame:
    suffix = path.suffix.lower()

    if suffix in {".xlsx", ".xls"}:
        return pd.read_excel(path)

    if suffix == ".csv":
        return pd.read_csv(path)

    return pd.read_csv(path, sep="\t")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--mapping", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--fc-threshold", type=float, default=1.3)
    parser.add_argument("--padj-threshold", type=float, default=0.01)
    parser.add_argument("--pseudocount", type=float, default=1.0)

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    mapping = pd.read_csv(args.mapping)

    required_mapping_cols = {"label", "file", "kcl_col", "mann_col"}
    missing = required_mapping_cols - set(mapping.columns)

    if missing:
        raise ValueError(f"Mapping file missing required columns: {sorted(missing)}")

    log2fc_threshold = np.log2(args.fc_threshold)

    all_degs: set[str] = set()
    raw_tables: dict[str, pd.DataFrame] = {}

    for _, row in mapping.iterrows():
        label = str(row["label"])
        file_path = input_dir / str(row["file"])

        if not file_path.exists():
            raise FileNotFoundError(file_path)

        df = read_table(file_path)
        raw_tables[label] = df

        required = {
            "gene_id",
            "gene_name",
            "log2FoldChange",
            "padj",
            str(row["kcl_col"]),
            str(row["mann_col"]),
        }

        missing_cols = required - set(df.columns)

        if missing_cols:
            raise ValueError(f"{file_path} missing columns: {sorted(missing_cols)}")

        passing = df[
            (pd.to_numeric(df["padj"], errors="coerce") <= args.padj_threshold)
            & (pd.to_numeric(df["log2FoldChange"], errors="coerce").abs() >= log2fc_threshold)
        ].copy()

        all_degs.update(passing["gene_id"].astype(str).tolist())

    final_frames: list[pd.DataFrame] = []

    for _, row in mapping.iterrows():
        label = str(row["label"])
        kcl_col = str(row["kcl_col"])
        mann_col = str(row["mann_col"])

        df = raw_tables[label].copy()
        df["gene_id"] = df["gene_id"].astype(str)

        subset = df[df["gene_id"].isin(all_degs)].copy()

        subset[kcl_col] = pd.to_numeric(subset[kcl_col], errors="coerce").fillna(0)
        subset[mann_col] = pd.to_numeric(subset[mann_col], errors="coerce").fillna(0)

        subset[f"log2FC_{label}"] = np.log2(
            (subset[kcl_col] + args.pseudocount)
            / (subset[mann_col] + args.pseudocount)
        )

        trimmed = (
            subset[["gene_id", "gene_name", f"log2FC_{label}"]]
            .drop_duplicates(subset=["gene_id"])
            .set_index("gene_id")
        )

        final_frames.append(trimmed)

    combined = pd.concat(final_frames, axis=1, join="outer")

    gene_name_cols = [c for c in combined.columns if c == "gene_name"]

    if gene_name_cols:
        combined["gene_name"] = combined[gene_name_cols].bfill(axis=1).iloc[:, 0]
        combined = combined.loc[
            :,
            [c for c in combined.columns if c != "gene_name"] + ["gene_name"],
        ]

    log2fc_cols = [c for c in combined.columns if c.startswith("log2FC_")]

    combined[log2fc_cols] = combined[log2fc_cols].fillna(0)
    combined = combined[~combined[log2fc_cols].eq(0).all(axis=1)]

    output_file = outdir / "combined_log2FC_FPKM_GlobalDEG.csv"
    combined.to_csv(output_file)

    params = outdir / "analysis_parameters.txt"
    params.write_text(
        "\n".join(
            [
                f"fc_threshold={args.fc_threshold}",
                f"log2fc_threshold=log2({args.fc_threshold})={log2fc_threshold}",
                f"padj_threshold={args.padj_threshold}",
                f"pseudocount={args.pseudocount}",
                f"n_global_degs={len(all_degs)}",
                f"n_output_genes={combined.shape[0]}",
            ]
        )
        + "\n"
    )

    print(f"Wrote: {output_file}")
    print(f"Wrote: {params}")


if __name__ == "__main__":
    main()
