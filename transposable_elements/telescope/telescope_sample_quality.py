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
import hashlib
import json
import platform
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd


SAMPLE_ORDER = ["A6M", "A6K", "B6M", "B6K", "C6M", "C6K"]


def sha256(path):
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for chunk in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()




def pca_scores(matrix):
    
    
    x = matrix.T.astype(float)
    x = x - x.mean(axis=0, keepdims=True)
    u, s, _ = np.linalg.svd(x, full_matrices=False)
    scores = u * s
    variance = (s ** 2) / np.sum(s ** 2)
    return scores, variance


def log2_cpm(counts, library_totals):
    return np.log2(counts.divide(library_totals, axis=1) * 1_000_000 + 1)




def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--package-root", type=Path, required=True)
    args = parser.parse_args()
    root = args.package_root.resolve()
    qc = root / "qc"
    qc.mkdir(exist_ok=True)

    samples = pd.read_csv(root / "samples.tsv", sep="\t", dtype=str).set_index("sample").loc[SAMPLE_ORDER]
    sample_qc = pd.read_csv(root / "tables/sample_qc.tsv", sep="\t").set_index("sample").loc[SAMPLE_ORDER]

    matrices = {"final_count": [], "final_conf": []}
    input_paths = [root / "samples.tsv", root / "tables/sample_qc.tsv"]
    for sample in SAMPLE_ORDER:
        path = root / "source_reports" / f"primary_{sample}.telescope_report.tsv.gz"
        input_paths.append(path)
        frame = pd.read_csv(path, sep="\t", comment="#", usecols=["transcript", "final_count", "final_conf"])
        if frame["transcript"].duplicated().any():
            raise ValueError(f"Duplicate transcript IDs in {path}")
        frame = frame.loc[frame["transcript"] != "__no_feature"].set_index("transcript")
        for metric in matrices:
            matrices[metric].append(frame[metric].rename(sample))

    counts = pd.concat(matrices["final_count"], axis=1).fillna(0).loc[:, SAMPLE_ORDER].astype(np.int64)
    confident = pd.concat(matrices["final_conf"], axis=1).fillna(0).loc[:, SAMPLE_ORDER].astype(np.int64)
    if (counts < 0).any().any() or (confident < 0).any().any() or (confident > counts).any().any():
        raise ValueError("Invalid final_count/final_conf relationship")

    keep = (counts >= 10).sum(axis=1) >= 3
    if keep.sum() < 2:
        raise ValueError("Detection filter retained too few loci")
    count_lib = sample_qc["te_final_count_sum"].astype(float)
    conf_lib = sample_qc["te_final_conf_sum"].astype(float)
    count_log = log2_cpm(counts.loc[keep], count_lib)
    conf_log = log2_cpm(confident.loc[keep], conf_lib)
    count_scores, count_var = pca_scores(count_log.to_numpy())
    conf_scores, conf_var = pca_scores(conf_log.to_numpy())
    correlation = count_log.corr(method="pearson")

    pca_rows = []
    for i, sample in enumerate(SAMPLE_ORDER):
        pca_rows.append({
            "sample": sample,
            "replicate_label": samples.loc[sample, "replicate_label"],
            "condition": samples.loc[sample, "condition"],
            "final_count_PC1": count_scores[i, 0],
            "final_count_PC2": count_scores[i, 1],
            "final_conf_PC1": conf_scores[i, 0],
            "final_conf_PC2": conf_scores[i, 1],
        })
    pd.DataFrame(pca_rows).to_csv(qc / "pca_sample_coordinates.tsv", sep="\t", index=False)
    correlation.to_csv(qc / "sample_correlation_final_count.tsv", sep="\t", index=True)

    filter_rows = [
        ["union_of_reported_loci", int(len(counts))],
        ["positive_final_count_in_at_least_1_sample", int(((counts > 0).sum(axis=1) >= 1).sum())],
        ["final_count_ge_10_in_at_least_1_sample", int(((counts >= 10).sum(axis=1) >= 1).sum())],
        ["final_count_ge_10_in_at_least_3_samples_primary_PCA", int(keep.sum())],
        ["final_conf_ge_10_in_at_least_3_samples", int(((confident >= 10).sum(axis=1) >= 3).sum())],
    ]
    pd.DataFrame(filter_rows, columns=["filter_stage", "loci"]).to_csv(qc / "pca_filter_summary.tsv", sep="\t", index=False)

    code_path = Path(__file__).resolve()
    output_paths = [
        qc / "pca_sample_coordinates.tsv", qc / "sample_correlation_final_count.tsv", qc / "pca_filter_summary.tsv",
    ]
    panel_manifest = {
        "status": "pass",
        "created_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Primary and assignment-confidence sample PCA, correlation and filter summaries",
        "samples": SAMPLE_ORDER,
        "primary_pca_filter": "final_count >= 10 in at least 3 of 6 samples",
        "pca_transform": "log2(count / sample TE-assigned total * 1e6 + 1)",
        "pca_center_features": True,
        "pca_scale_features_to_unit_variance": False,
        "confidence_pca_uses_same_loci_as_final_count_pca": True,
        "retained_loci": int(keep.sum()),
        "final_count_explained_variance_pc1_pc2": [float(count_var[0]), float(count_var[1])],
        "final_conf_explained_variance_pc1_pc2": [float(conf_var[0]), float(conf_var[1])],
        "software": {
            "python": platform.python_version(),
            "numpy": np.__version__,
            "pandas": pd.__version__,
        },
        "script": {"path": str(code_path), "sha256": sha256(code_path)},
        "inputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in input_paths],
        "outputs": [{"path": p.relative_to(root).as_posix(), "bytes": p.stat().st_size, "sha256": sha256(p)} for p in output_paths],
    }
    (qc / "QC_PANEL_MANIFEST.json").write_text(json.dumps(panel_manifest, indent=2, sort_keys=True) + "\n")

    print(json.dumps({
        "status": "complete",
        "retained_loci": int(keep.sum()),
        "count_pc1_pc2_percent": [round(float(count_var[0] * 100), 3), round(float(count_var[1] * 100), 3)],
        "conf_pc1_pc2_percent": [round(float(conf_var[0] * 100), 3), round(float(conf_var[1] * 100), 3)],
        "outputs": [p.relative_to(root).as_posix() for p in output_paths],
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
