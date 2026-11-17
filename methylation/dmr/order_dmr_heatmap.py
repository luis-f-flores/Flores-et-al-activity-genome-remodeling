"""Group complete DMR-union methylation profiles and order them for display."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd
from scipy.cluster.vq import kmeans2
from scipy.cluster.hierarchy import linkage, leaves_list


def calculate(matrix_file, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    matrix = np.load(matrix_file)
    values = matrix["mCG_percent"]
    keep = np.isfinite(values).all(axis=1) & (np.nanstd(values, axis=1) > 0)
    centered = values[keep] - values[keep].mean(axis=1, keepdims=True)
    eigenvalues, eigenvectors = np.linalg.eigh(centered.T @ centered)
    eigenvalues = eigenvalues[::-1].clip(0)
    eigenvectors = eigenvectors[:, ::-1]
    loading = centered @ eigenvectors[:, :1] / np.sqrt(eigenvalues[:1])[None, :]
    if loading[np.argmax(abs(loading[:, 0])), 0] < 0:
        loading[:, 0] *= -1
    centers, labels = kmeans2(centered, 6, iter=60, minit="++", missing="raise", seed=20261013)
    group_order = leaves_list(linkage(centers, method="average", metric="euclidean", optimal_ordering=True))
    display = np.empty(6, dtype=int)
    display[group_order] = np.arange(6)
    order = np.lexsort((loading[:, 0], display[labels]))
    pd.DataFrame({"region_id": matrix["region_ids"][keep],
                  "display_group": display[labels] + 1,
                  "PC1_loading": loading[:, 0],
                  "display_order": np.argsort(order)}).to_csv(
        output / "heatmap_cluster_order.tsv.gz", sep="\t", index=False)
    pd.DataFrame(centers[group_order], columns=matrix["samples"], index=np.arange(1, 7)).rename_axis(
        "display_group").to_csv(output / "heatmap_cluster_centers.tsv", sep="\t")
    print(f"Ordered {int(keep.sum()):,} complete, nonconstant DMR profiles.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--matrix", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    calculate(args.matrix, args.output)
