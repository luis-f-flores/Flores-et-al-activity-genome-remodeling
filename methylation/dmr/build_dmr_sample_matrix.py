"""Count methylated and unmethylated CpG reads across the merged DMR union."""
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import argparse
import subprocess

import numpy as np
import pandas as pd


def merge_calls(paths):
    calls = []
    for path in paths:
        table = pd.read_csv(path, sep="\t")
        table["source"] = Path(path).name
        table["source_row"] = np.arange(len(table)) + 1
        calls.append(table)
    table = pd.concat(calls, ignore_index=True).sort_values(
        ["chr", "start", "end"], kind="stable")
    union, membership = [], []
    for row in table.itertuples(index=False):
        if not union or row.chr != union[-1][0] or row.start > union[-1][2]:
            union.append([row.chr, int(row.start), int(row.end)])
        else:
            union[-1][2] = max(union[-1][2], int(row.end))
        membership.append((row.source, int(row.source_row), f"DMRU{len(union):06d}"))
    regions = pd.DataFrame(union, columns=["chr", "start", "end"])
    regions["region_id"] = [f"DMRU{i + 1:06d}" for i in range(len(regions))]
    return regions, pd.DataFrame(membership, columns=["source_file", "source_row", "region_id"])


def build(dmrs, samples_file, output, workers=2, compiler="g++"):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    samples = pd.read_csv(samples_file, sep="\t").set_index("sample")
    ids = [f"{rep}{time}{treatment}"
           for treatment, times in [("", [0]), ("M", [1, 6, 24]), ("K", [1, 6, 24])]
           for time in times for rep in "ABC"]
    if not samples.index.is_unique or set(samples.index) != set(ids):
        raise ValueError("The sample table must contain the 21 baseline, mannitol and KCl libraries.")
    regions, membership = merge_calls(dmrs)
    regions.to_csv(output / "DMR_union_1based.tsv", sep="\t", index=False)
    membership.to_csv(output / "DMR_source_membership.tsv.gz", sep="\t", index=False)
    region_file = output / "count_regions_1based.tsv"
    regions.to_csv(region_file, sep="\t", index=False, header=False)
    executable = output / "count_dmr_regions"
    subprocess.run([compiler, "-O3", "-std=c++17",
                    str(Path(__file__).with_name("count_dmr_regions.cpp")),
                    "-o", str(executable)], check=True)

    def count(sample):
        row = samples.loc[sample]
        report = Path(row.report)
        if not report.is_absolute():
            report = Path(samples_file).resolve().parent / report
        prefix = output / sample
        subprocess.run([str(executable), str(report), str(region_file),
                        str(row["upper_99.9"]), str(prefix),
                        str(output / f"{sample}_track_CpGs.tsv")], check=True)
        counts = pd.read_csv(prefix.with_suffix(".tsv"), sep="\t")
        if counts.region_id.tolist() != regions.region_id.tolist():
            raise ValueError(f"Unexpected region order for {sample}")
        return counts

    with ThreadPoolExecutor(max_workers=workers) as pool:
        counts = list(pool.map(count, ids))
    methylated = np.column_stack([d.numCs.to_numpy() for d in counts])
    unmethylated = np.column_stack([d.numTs.to_numpy() for d in counts])
    positions = np.column_stack([d.covered_CpG_records.to_numpy() for d in counts])
    coverage = methylated + unmethylated
    mcg = np.divide(100 * methylated, coverage,
                    out=np.full(methylated.shape, np.nan), where=coverage > 0)
    np.savez_compressed(output / "DMR_sample_matrix.npz", samples=np.array(ids),
                        region_ids=regions.region_id.to_numpy(dtype=str),
                        numCs=methylated, numTs=unmethylated, coverage=coverage,
                        covered_CpG_records=positions, mCG_percent=mcg)
    pd.concat([regions, pd.DataFrame(mcg, columns=ids)], axis=1).to_csv(
        output / "DMR_sample_mCG_percent.tsv.gz", sep="\t", index=False,
        na_rep="NA", float_format="%.6f")
    print(f"Counted {len(regions):,} DMR-union regions in {len(ids)} samples.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dmrs", type=Path, nargs=6, required=True,
                        help="The six merged hypo/hyper DMR tables, with 1-based inclusive chr/start/end columns.")
    parser.add_argument("--samples", type=Path, required=True,
                        help="Coverage-cutoff table with sample, report and upper_99.9 columns.")
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--workers", type=int, default=2)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    build(args.dmrs, args.samples, args.output, args.workers, args.compiler)
