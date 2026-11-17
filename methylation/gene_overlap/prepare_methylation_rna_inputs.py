"""Build promoter and gene-body cohorts for continuous methylation–RNA analyses."""
from pathlib import Path
import argparse
import numpy as np
import pandas as pd


def prepare(promoter_files, gene_body_file, output):
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    promoter_rows = []
    for timepoint, path in zip([1, 6, 24], promoter_files):
        table = pd.read_csv(path, sep="\t", float_precision="round_trip")
        if not table.gene_id.is_unique:
            raise ValueError(f"Nonunique promoter RNA gene IDs at {timepoint} h")
        promoter_rows.append(pd.DataFrame({"timepoint": timepoint,
                                          "gene_id": table.gene_id,
                                          "delta_mCG_percent": table.meth_diff}))
    pd.concat(promoter_rows, ignore_index=True).to_csv(
        output / "promoter_cohort.tsv", sep="\t", index=False)
    body = pd.read_csv(gene_body_file, float_precision="round_trip")
    body_rows = []
    for timepoint in [1, 6, 24]:
        table = body[body.timepoint.eq(f"{timepoint}hr")].copy()
        if not table.gene_id_rna.is_unique or not table.body_locus_id.is_unique:
            raise ValueError(f"Nonunique gene-body/RNA matches at {timepoint} h")
        deg = table[(table.padj <= .01) & (table.log2FoldChange.abs() >= np.log2(1.3))]
        for name, selected in [("differentially_methylated", table), ("DEG_subset", deg)]:
            body_rows.append(pd.DataFrame({"timepoint": timepoint, "cohort": name,
                                          "gene_id": selected.gene_id_rna,
                                          "delta_mCG_percent": selected["meth.diff"]}))
    pd.concat(body_rows, ignore_index=True).to_csv(
        output / "body_cohorts.tsv", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--promoters", type=Path, nargs=3, required=True,
                        help="The promoter_rna_cohorts model-input tables at 1 h, 6 h and recovery.")
    parser.add_argument("--gene-body", type=Path, required=True,
                        help="The CSV of coordinate-supported matches from gene_body_matching.py.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.promoters, args.gene_body, args.output)
