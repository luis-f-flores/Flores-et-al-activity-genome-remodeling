#!/usr/bin/env python3
"""
02_lad_gene_set_permutation.py

Run empirical permutation tests for overlap between gene sets and LAD intervals.

Overlap definition:
    A gene is considered LAD-overlapping if the sum of LAD-overlapping bases across
    that gene is at least --threshold of the gene length.

Permutation test:
    Random gene sets of the same size are sampled from the gene universe.
    A two-tailed empirical p value is computed relative to the permutation mean.
"""

from __future__ import annotations

import argparse
from pathlib import Path
import random
import warnings

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from pybedtools import BedTool
from tqdm import tqdm


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()

    parser.add_argument("--genes-bed", required=True)
    parser.add_argument("--gene-lists", required=True, help="CSV/TSV with columns: gene_set,file")
    parser.add_argument("--lad-map", required=True, help="CSV/TSV with columns: lad_name,file")
    parser.add_argument("--outdir", default="permutation_results")
    parser.add_argument("--thresholds", default="0.75,1.0")
    parser.add_argument("--iterations", type=int, default=10000)
    parser.add_argument("--n-jobs", type=int, default=6)
    parser.add_argument("--seed", type=int, default=1)

    return parser.parse_args()


def read_map(path: Path) -> pd.DataFrame:
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)

    return pd.read_csv(path, sep="\t")


def read_gene_list(path: Path) -> list[str]:
    return pd.read_csv(path, header=None)[0].dropna().astype(str).tolist()


def calculate_overlap(
    gene_list: list[str],
    lads: BedTool,
    genes_df: pd.DataFrame,
    genes_bed: BedTool,
    threshold: float,
) -> tuple[float, list[str]]:
    gene_set = set(gene_list)
    gene_bed = genes_bed.filter(lambda x: x.name in gene_set)

    overlap = gene_bed.intersect(lads, wao=True)

    overlap_df = overlap.to_dataframe(
        names=[
            "chrom",
            "start",
            "end",
            "name",
            "lad_chrom",
            "lad_start",
            "lad_end",
            "overlap",
        ]
    )

    if overlap_df.empty:
        return 0.0, []

    overlap_df["start"] = pd.to_numeric(overlap_df["start"])
    overlap_df["end"] = pd.to_numeric(overlap_df["end"])
    overlap_df["overlap"] = pd.to_numeric(overlap_df["overlap"])

    overlap_df["gene_length"] = overlap_df["end"] - overlap_df["start"]

    gene_overlap = (
        overlap_df.groupby("name")
        .agg(overlap_bases=("overlap", "sum"), gene_length=("gene_length", "first"))
        .reset_index()
    )

    gene_overlap["overlap_fraction"] = gene_overlap["overlap_bases"] / gene_overlap["gene_length"]
    gene_overlap["significant_overlap"] = gene_overlap["overlap_fraction"] >= threshold

    significant_genes = gene_overlap.loc[gene_overlap["significant_overlap"], "name"].tolist()
    observed_percent = 100 * len(significant_genes) / len(gene_list)

    return observed_percent, significant_genes


def permutation_single(
    gene_count: int,
    all_gene_names: list[str],
    lads: BedTool,
    genes_df: pd.DataFrame,
    genes_bed: BedTool,
    threshold: float,
) -> float:
    random_genes = random.sample(all_gene_names, gene_count)
    overlap_percent, _ = calculate_overlap(random_genes, lads, genes_df, genes_bed, threshold)

    return overlap_percent


def run_permutation(
    gene_list: list[str],
    gene_set_name: str,
    lad_name: str,
    lads: BedTool,
    genes_df: pd.DataFrame,
    genes_bed: BedTool,
    threshold: float,
    iterations: int,
    n_jobs: int,
    outdir: Path,
) -> None:
    observed_overlap, significant_genes = calculate_overlap(
        gene_list,
        lads,
        genes_df,
        genes_bed,
        threshold,
    )

    all_gene_names = genes_df["name"].astype(str).tolist()

    permutations = Parallel(n_jobs=n_jobs)(
        delayed(permutation_single)(
            len(gene_list),
            all_gene_names,
            lads,
            genes_df,
            genes_bed,
            threshold,
        )
        for _ in tqdm(range(iterations), desc=f"{gene_set_name} vs {lad_name} ({threshold})")
    )

    permutations = np.asarray(permutations)
    perm_mean = float(np.mean(permutations))
    perm_std = float(np.std(permutations))

    two_std_above = perm_mean + 2 * perm_std
    two_std_below = perm_mean - 2 * perm_std

    observed_distance = abs(observed_overlap - perm_mean)

    p_value = (
        np.sum(np.abs(permutations - perm_mean) >= observed_distance) + 1
    ) / (len(permutations) + 1)

    result = pd.DataFrame(
        [
            {
                "gene_set": gene_set_name,
                "lad_name": lad_name,
                "threshold": threshold,
                "observed_overlap": observed_overlap,
                "p_value (two-tailed empirical permutation test)": p_value,
                "mean_permutations": perm_mean,
                "two_std_above": two_std_above,
                "two_std_below": two_std_below,
                "significant_genes": ",".join(significant_genes),
            }
        ]
    )

    outfile = outdir / f"{gene_set_name}_{lad_name}_permutation_test_results_{int(threshold * 100)}_percent.csv"
    result.to_csv(outfile, index=False)

    print(f"Wrote: {outfile}")


def main() -> None:
    args = parse_args()

    warnings.filterwarnings("ignore", category=DeprecationWarning)

    random.seed(args.seed)
    np.random.seed(args.seed)

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    genes_df = pd.read_csv(
        args.genes_bed,
        sep="\t",
        header=None,
        names=["chrom", "start", "end", "name"],
    )

    genes_bed = BedTool.from_dataframe(genes_df).sort()

    gene_lists = read_map(Path(args.gene_lists))
    lad_map = read_map(Path(args.lad_map))

    thresholds = [float(x) for x in args.thresholds.split(",")]

    for _, lad_row in lad_map.iterrows():
        lad_name = str(lad_row["lad_name"])
        lad_file = Path(str(lad_row["file"]))
        lads = BedTool(str(lad_file)).sort()

        for _, gene_row in gene_lists.iterrows():
            gene_set_name = str(gene_row["gene_set"])
            gene_list_file = Path(str(gene_row["file"]))
            gene_list = read_gene_list(gene_list_file)

            for threshold in thresholds:
                run_permutation(
                    gene_list,
                    gene_set_name,
                    lad_name,
                    lads,
                    genes_df,
                    genes_bed,
                    threshold,
                    args.iterations,
                    args.n_jobs,
                    outdir,
                )


if __name__ == "__main__":
    main()
