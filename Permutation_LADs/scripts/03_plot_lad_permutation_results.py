#!/usr/bin/env python3
"""
03_plot_lad_permutation_results.py

Plot LAD overlap permutation results for one or more gene sets.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("pdf")

import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd


mpl.rcParams["pdf.fonttype"] = 42
mpl.rcParams["ps.fonttype"] = 42
mpl.rcParams["font.family"] = "DejaVu Sans"
mpl.rcParams["text.usetex"] = False


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--input-dir", required=True)
    parser.add_argument("--gene-set", required=True, help="Example: dPRG, SRG, rPRG")
    parser.add_argument("--threshold", type=int, default=75)
    parser.add_argument("--out", required=True)
    parser.add_argument("--title", default=None)

    return parser.parse_args()


def main() -> None:
    args = parse_args()

    input_dir = Path(args.input_dir)
    pattern = f"{args.gene_set}_*_permutation_test_results_{args.threshold}_percent.csv"
    files = sorted(input_dir.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No files matched: {input_dir / pattern}")

    rows = []

    for file in files:
        df = pd.read_csv(file)
        row = df.iloc[0].to_dict()
        lad_name = str(row.get("lad_name", file.name))
        row["label"] = lad_name
        rows.append(row)

    plot_df = pd.DataFrame(rows)

    labels = plot_df["label"].tolist()
    observed = plot_df["observed_overlap"].astype(float).tolist()
    expected = plot_df["mean_permutations"].astype(float).tolist()
    upper = plot_df["two_std_above"].astype(float).tolist()
    lower = plot_df["two_std_below"].astype(float).tolist()
    p_values = plot_df["p_value (two-tailed empirical permutation test)"].astype(float).tolist()

    lower_error = [m - b for m, b in zip(expected, lower)]
    upper_error = [a - m for a, m in zip(upper, expected)]

    fig, ax = plt.subplots(figsize=(14, 6))

    ax.bar(
        labels,
        observed,
        width=0.55,
        edgecolor="black",
        linewidth=1.0,
        label="Observed overlap",
    )

    ax.errorbar(
        labels,
        expected,
        yerr=[lower_error, upper_error],
        fmt="o",
        color="black",
        ecolor="black",
        elinewidth=1.5,
        capsize=4,
        label="Expected overlap ± 2 SD",
    )

    for i, p_value in enumerate(p_values):
        y_pos = max(observed[i], upper[i]) + 1
        ax.text(i, y_pos, f"p={p_value:.4f}", ha="center", va="bottom", fontsize=9)

    ax.set_ylabel("% overlap")
    ax.set_title(
        args.title
        if args.title is not None
        else f"{args.gene_set} overlap with LADs ({args.threshold}% gene-overlap threshold)"
    )
    ax.legend()
    ax.grid(False)

    plt.xticks(rotation=25, ha="right")
    plt.tight_layout()
    plt.savefig(args.out, format="pdf", bbox_inches="tight", dpi=300)

    print(f"Wrote: {args.out}")


if __name__ == "__main__":
    main()
