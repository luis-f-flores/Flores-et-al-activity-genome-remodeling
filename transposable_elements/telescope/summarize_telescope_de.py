#!/usr/bin/env python3

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)


import argparse
import csv
import gzip
import json
import math
from pathlib import Path
from statistics import median


def read_table(path):
    with gzip.open(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def number(value):
    if value in {"", "NA", None}:
        return None
    parsed = float(value)
    return parsed if math.isfinite(parsed) else None


def summarize(root):
    ROOT = Path(root).resolve()
    RESULTS = {"primary": ROOT/"results/DE_KCl_vs_Mannitol_final_count.tsv.gz",
               "final_conf": ROOT/"results/DE_final_conf.tsv.gz"}
    (ROOT/"qc").mkdir(parents=True,exist_ok=True)
    tables = {name: read_table(path) for name, path in RESULTS.items()}
    primary = tables["primary"]
    expected_rows = len(primary)
    if expected_rows == 0:
        raise SystemExit("Primary result table is empty")
    if len({row["locus_id"] for row in primary}) != expected_rows:
        raise SystemExit("Primary result table contains duplicate locus IDs")
    for name, rows in tables.items():
        if len(rows) != expected_rows:
            raise SystemExit(f"{name} has {len(rows)} rows; expected {expected_rows}")
        if [row["locus_id"] for row in rows] != [row["locus_id"] for row in primary]:
            raise SystemExit(f"{name} locus order differs from the primary table")

    pvalues = [number(row["pvalue"]) for row in primary]
    padj_all = [number(row["padj_BH_all_retained"]) for row in primary]
    padj_if = [number(row["padj_DESeq2_independent_filter"]) for row in primary]
    shrunk = [number(row["log2FoldChange_shrunk"]) for row in primary]
    base_means = [number(row["baseMean"]) for row in primary]

    significant = [
        (row, effect)
        for row, adjusted, effect in zip(primary, padj_all, shrunk)
        if adjusted is not None and adjusted < 0.05 and effect is not None
    ]
    summary = {
        "status": "pass",
        "scope": "Result-table integrity and DE summary",
        "loci_tested": expected_rows,
        "finite_pvalues": sum(value is not None for value in pvalues),
        "cooks_flagged_loci": sum(row["cooks_flag"].upper() == "TRUE" for row in primary),
        "padj_BH_all_retained_lt_0_05": len(significant),
        "padj_BH_all_retained_lt_0_05_positive_shrunken_lfc": sum(effect > 0 for _, effect in significant),
        "padj_BH_all_retained_lt_0_05_negative_shrunken_lfc": sum(effect < 0 for _, effect in significant),
        "DESeq2_independent_filter_padj_lt_0_05": sum(value is not None and value < 0.05 for value in padj_if),
        "median_baseMean": median(value for value in base_means if value is not None),
        "median_absolute_shrunken_log2FC": median(abs(value) for value in shrunk if value is not None),
        "hard_log2FC_cutoff_used": False,
    }

    (ROOT / "qc" / "result_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")

    print(json.dumps(summary, sort_keys=True))

if __name__ == '__main__':
    parser=argparse.ArgumentParser(description='Summarize Telescope locus differential-expression tables.')
    parser.add_argument('--results-root',type=Path,required=True)
    summarize(parser.parse_args().results_root)
