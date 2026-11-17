"""Identify genes tested in each KCl comparison and the bicuculline–TTX comparison."""
from pathlib import Path
import argparse
import pandas as pd


def prepare(kcl_files, bic_file, output):
    bic = pd.read_csv(bic_file, sep="\t", float_precision="round_trip").set_index("gene_id")
    if not bic.index.is_unique:
        raise ValueError("Bicuculline RNA gene IDs are not unique.")
    rows = []
    for path, label in zip(kcl_files, ["1 h", "6 h", "Recovery"]):
        kcl = pd.read_csv(path, sep="\t", float_precision="round_trip").set_index("gene_id")
        if not kcl.index.is_unique:
            raise ValueError(f"KCl RNA gene IDs are not unique: {label}")
        genes = kcl.index.intersection(bic.index).sort_values()
        for field in ["gene_chr", "gene_start", "gene_end", "gene_strand", "gene_length"]:
            if not kcl.loc[genes, field].astype(str).equals(bic.loc[genes, field].astype(str)):
                raise ValueError(f"Gene annotations differ between experiments: {field}")
        rows.append(pd.DataFrame({"gene_id": genes, "timepoint": label,
                                  "KCl_pvalue": kcl.loc[genes, "pvalue"].to_numpy(),
                                  "KCl_log2FoldChange": kcl.loc[genes, "log2FoldChange"].to_numpy(),
                                  "Bic_pvalue": bic.loc[genes, "pvalue"].to_numpy(),
                                  "Bic_log2FoldChange": bic.loc[genes, "log2FoldChange"].to_numpy()}))
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    pd.concat(rows, ignore_index=True).to_csv(output / "testability.tsv.gz", sep="\t", index=False)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--kcl", type=Path, nargs=3, required=True,
                        help="Provider DE tables for 1 h, 6 h and recovery, in that order.")
    parser.add_argument("--bic", type=Path, required=True,
                        help="Provider Bic-versus-TTX DE table.")
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    prepare(args.kcl, args.bic, args.output)
