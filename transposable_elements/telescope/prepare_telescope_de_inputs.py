#!/usr/bin/env python3

from __future__ import annotations

from pathlib import Path as _GB_Path
import os as _GB_os
GB_WORK_ROOT = _GB_Path(_GB_os.environ.get('GB_WORK_ROOT', _GB_Path('../work').resolve()))
def gb_path(name):
    if name.startswith('@tool/'):
        return name.split('/')[-1]
    if name.startswith('@fonts'):
        return str(_GB_Path(_GB_os.environ.get('GB_FONT_DIR', '/usr/share/fonts/truetype/msttcorefonts')) / name.removeprefix('@fonts').lstrip('/'))
    return str(GB_WORK_ROOT / name)


import shutil
import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


REQUIRED_REPORT_COLUMNS = {
    "transcript",
    "transcript_length",
    "final_count",
    "final_conf",
    "unique_count",
    "init_best",
}


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def require_integer(series: pd.Series, label: str) -> pd.Series:
    values = pd.to_numeric(series, errors="raise")
    if not np.isfinite(values).all() or (values < 0).any():
        raise ValueError(f"{label} contains non-finite or negative values")
    rounded = np.rint(values)
    if not np.allclose(values, rounded, rtol=0, atol=1e-8):
        raise ValueError(f"{label} contains non-integer values")
    return pd.Series(rounded.astype(np.int64), index=series.index)


def write_tsv_gz(frame: pd.DataFrame, path: Path, index_label: str) -> None:
    frame.to_csv(path, sep="\t", index=True, index_label=index_label, compression="gzip")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--source-package",
        type=Path,
        required=True,
        help="Directory containing samples.tsv, tables/sample_qc.tsv and source_reports/.",
    )
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    package_root = args.output.resolve()

    source = args.source_package.resolve()
    required_source_files = [source / "samples.tsv", source / "tables" / "sample_qc.tsv"]
    for path in required_source_files:
        if not path.is_file():
            raise FileNotFoundError(path)

    samples = pd.read_csv(source / "samples.tsv", sep="\t", dtype=str)
    if samples["sample"].duplicated().any():
        raise ValueError("Duplicate sample names")
    sample_order = samples["sample"].tolist()
    expected = {"A6M", "A6K", "B6M", "B6K", "C6M", "C6K"}
    if set(sample_order) != expected:
        raise ValueError(f"Unexpected sample set: {sample_order}")

    sample_qc = pd.read_csv(source / "tables" / "sample_qc.tsv", sep="\t")
    sample_qc = sample_qc.set_index("sample").loc[sample_order]
    if sample_qc["version"].nunique() != 1 or sample_qc["annotated_features"].nunique() != 1:
        raise ValueError("Samples do not share one Telescope version and annotation size")

    matrices: dict[str, list[pd.Series]] = {
        "final_count": [],
        "final_conf": [],
        "unique_count": [],
        "init_best": [],
    }
    lengths: list[pd.Series] = []
    source_paths = required_source_files.copy()

    for sample in sample_order:
        path = source / "source_reports" / f"primary_{sample}.telescope_report.tsv.gz"
        if not path.is_file():
            raise FileNotFoundError(path)
        source_paths.append(path)
        report = pd.read_csv(path, sep="\t", comment="#")
        missing = REQUIRED_REPORT_COLUMNS.difference(report.columns)
        if missing:
            raise ValueError(f"{path.name} is missing columns: {sorted(missing)}")
        if report["transcript"].duplicated().any():
            raise ValueError(f"Duplicate transcript IDs in {path.name}")
        
        report = report.loc[report["transcript"] != "__no_feature"].copy()
        report = report.set_index("transcript")
        for field in matrices:
            values = require_integer(report[field], f"{sample}:{field}")
            values.name = sample
            matrices[field].append(values)
        length = pd.to_numeric(report["transcript_length"], errors="raise")
        length.name = sample
        lengths.append(length)

    built: dict[str, pd.DataFrame] = {}
    for field, series_list in matrices.items():
        built[field] = pd.concat(series_list, axis=1).fillna(0).loc[:, sample_order].astype(np.int64)

    if (built["final_conf"] > built["final_count"]).any().any():
        raise ValueError("final_conf exceeds final_count for at least one locus/sample")

    length_matrix = pd.concat(lengths, axis=1)
    conflicts = length_matrix.nunique(axis=1, dropna=True) > 1
    if conflicts.any():
        raise ValueError(f"Conflicting transcript lengths for {int(conflicts.sum())} loci")
    locus_length = length_matrix.bfill(axis=1).iloc[:, 0].astype(np.int64)

    if not all(built[field].index.equals(built["final_count"].index) for field in built):
        raise ValueError("Count/support matrices do not share the same locus order")

    package_root.mkdir(parents=True, exist_ok=True)
    write_tsv_gz(built["final_count"], package_root / "counts_final_count.tsv.gz", "locus_id")
    write_tsv_gz(built["final_conf"], package_root / "counts_final_conf.tsv.gz", "locus_id")

    count_sum = built["final_count"].sum(axis=1)
    conf_sum = built["final_conf"].sum(axis=1)
    support = pd.DataFrame(index=built["final_count"].index)
    support["transcript_length"] = locus_length
    support["sum_final_count"] = count_sum
    support["sum_final_conf"] = conf_sum
    support["pooled_confident_fraction"] = np.where(count_sum > 0, conf_sum / count_sum, np.nan)
    support["samples_positive_final_count"] = (built["final_count"] > 0).sum(axis=1)
    support["samples_final_count_ge_10"] = (built["final_count"] >= 10).sum(axis=1)
    support["samples_final_conf_ge_10"] = (built["final_conf"] >= 10).sum(axis=1)
    support["sum_unique_count"] = built["unique_count"].sum(axis=1)
    support["sum_init_best"] = built["init_best"].sum(axis=1)
    support.to_csv(
        package_root / "locus_support.tsv.gz",
        sep="\t",
        index=True,
        index_label="locus_id",
        compression="gzip",
        na_rep="NA",
    )

    samples.to_csv(package_root / "samples.tsv", sep="\t", index=False)
    numeric_qc = sample_qc.copy()
    mapped = numeric_qc[["pair_mapped", "pair_mixed", "single_mapped"]].sum(axis=1).astype(np.int64)
    geometric_mean = math.exp(np.log(mapped.to_numpy(dtype=float)).mean())
    normalization = pd.DataFrame(
        {
            "sample": sample_order,
            "total_fragments": numeric_qc["total_fragments"].astype(np.int64).to_numpy(),
            "mapped_fragments": mapped.to_numpy(),
            "mapped_fragment_size_factor": mapped.to_numpy(dtype=float) / geometric_mean,
            "te_final_count_sum": built["final_count"].sum(axis=0).loc[sample_order].to_numpy(),
            "te_final_conf_sum": built["final_conf"].sum(axis=0).loc[sample_order].to_numpy(),
        }
    )
    normalization.to_csv(package_root / "normalization_inputs.tsv", sep="\t", index=False)

    source_manifest = []
    for path in source_paths:
        source_manifest.append(
            {
                "source_path": path.relative_to(source).as_posix(),
                "bytes": path.stat().st_size,
                "sha256": sha256(path),
            }
        )
    pd.DataFrame(source_manifest).to_csv(package_root / "SOURCE_MANIFEST.tsv", sep="\t", index=False)

    filter_counts = {
        "reported_union": int(len(built["final_count"])),
        "positive_in_at_least_1": int((built["final_count"] > 0).any(axis=1).sum()),
        "final_count_ge_10_in_at_least_1": int(((built["final_count"] >= 10).sum(axis=1) >= 1).sum()),
        "final_count_ge_10_in_at_least_3": int(((built["final_count"] >= 10).sum(axis=1) >= 3).sum()),
    }
    preparation = {
        "status": "prepared",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "source_package": str(source),
        "samples": sample_order,
        "telescope_version": str(sample_qc["version"].iloc[0]),
        "annotated_features": int(sample_qc["annotated_features"].iloc[0]),
        "matrix_rows": int(len(built["final_count"])),
        "filter_counts": filter_counts,
        "fill_rule": "union of biological locus rows after excluding __no_feature; absent report row filled with zero after shared annotation checks",
        "normalization_input": "mapped fragments = pair_mapped + pair_mixed + single_mapped",
        "script": Path(__file__).name,
    }
    (package_root / "INPUT_PREPARATION.json").write_text(
        json.dumps(preparation, indent=2, sort_keys=True) + "\n"
    )
    config_path = package_root / 'analysis_config.tsv'
    if not config_path.exists():
        shutil.copyfile(Path(__file__).with_name('analysis_config.tsv'), config_path)
    print(json.dumps(preparation, sort_keys=True))


if __name__ == "__main__":
    main()
