"""Prepare IAP example tracks, component RNA changes and host-gene expression."""
from __future__ import annotations

from pathlib import Path

import argparse
import pyBigWig
import numpy as np
import pandas as pd


LOCI = [
    {
        "role": "genic",
        "locus_id": "IAPEZI@chr16:32148444-32149702",
        "chrom": "chr16",
        "panel_start0": 32148076,
        "panel_end0": 32153339,
    },
    {
        "role": "intergenic",
        "locus_id": "IAPEZI@chr15:39173390-39174662",
        "chrom": "chr15",
        "panel_start0": 39173040,
        "panel_end0": 39178257,
    },
]


def bin_mcg(frame: pd.DataFrame, start: int, end: int, width: int = 250) -> pd.DataFrame:
    nbin = int(np.ceil((end - start) / width))
    x = frame.copy()
    x["bin"] = ((x.position0 - start) // width).astype(int)
    out = x.groupby("bin").agg(
        mcg_delta_pp=("delta_mcg_pp", "mean"),
        n_shared=("position0", "size"),
    )
    out = out.reindex(range(nbin))
    out["n_shared"] = out.n_shared.fillna(0).astype(int)
    out.loc[out.n_shared < 2, "mcg_delta_pp"] = np.nan
    out["start0"] = start + np.arange(nbin) * width
    out["end0"] = np.minimum(out.start0 + width, end)
    return out.reset_index(drop=True)


def empty_record(locus: dict) -> dict:
    return {
        "record_type": "",
        "role": locus["role"],
        "selected_locus": locus["locus_id"],
        "chrom": locus["chrom"],
        "panel_start0": locus["panel_start0"],
        "panel_end0": locus["panel_end0"],
        "start0": np.nan,
        "end0": np.nan,
        "component": "",
        "rna_log2fc": np.nan,
        "mcg_delta_pp": np.nan,
        "n_shared": np.nan,
        "is_selected_component": False,
        "dmr_direction": "",
        "dmr_meth_diff": np.nan,
        "dmr_qvalue": np.nan,
        "host_gene": "",
        "host_mean_log2_ratio": np.nan,
        "host_sem_log2_ratio": np.nan,
        "host_mean_fpkm": np.nan,
        "host_padj": np.nan,
        "host_n_pairs": np.nan,
    }


def main() -> None:
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ['context','dmr-overlaps','repeat-annotation','kcl-mcg','mannitol-mcg','fpkm','gene-de','output']:
        parser.add_argument('--'+name,type=Path,required=True)
    args=parser.parse_args()
    context=pd.read_csv(args.context,sep='\t')
    dmr=pd.read_csv(args.dmr_overlaps,sep='\t')
    master=pd.read_csv(args.repeat_annotation,sep='\t',low_memory=False)
    fpkm=pd.read_csv(args.fpkm,sep='\t')
    deg=pd.read_csv(args.gene_de,sep='\t')
    shared=[]
    with pyBigWig.open(str(args.kcl_mcg)) as kcl, pyBigWig.open(str(args.mannitol_mcg)) as mann:
        for locus in LOCI:
            lo,hi=locus['panel_start0'],locus['panel_end0']
            k=np.asarray(kcl.values(locus['chrom'],lo,hi))
            m=np.asarray(mann.values(locus['chrom'],lo,hi))
            valid=np.isfinite(k)&np.isfinite(m)
            shared.extend(dict(selected_locus=locus['locus_id'],position0=lo+int(i),
                               delta_mcg_pp=float(k[i]-m[i])) for i in np.flatnonzero(valid))
    mcg=pd.DataFrame(shared)

    rna_lookup = context.set_index("locus_id")["rna_log2fc_shrunk"].to_dict()
    rows: list[dict] = []

    for locus in LOCI:
        values = mcg[mcg.selected_locus.eq(locus["locus_id"])]
        bins = bin_mcg(values, locus["panel_start0"], locus["panel_end0"])
        for b in bins.itertuples(index=False):
            row = empty_record(locus)
            row.update(
                record_type="mcg_bin",
                start0=int(b.start0),
                end0=int(b.end0),
                mcg_delta_pp=float(b.mcg_delta_pp) if np.isfinite(b.mcg_delta_pp) else np.nan,
                n_shared=int(b.n_shared),
            )
            rows.append(row)

        comps = master[
            master.chrom.eq(locus["chrom"])
            & master.start_0based.lt(locus["panel_end0"])
            & master.end_0based_exclusive.gt(locus["panel_start0"])
            & master.is_iap.astype(str).str.upper().eq("TRUE")
        ].sort_values(["start_0based", "end_0based_exclusive"])
        for c in comps.itertuples(index=False):
            row = empty_record(locus)
            row.update(
                record_type="component",
                start0=int(c.start_0based),
                end0=int(c.end_0based_exclusive),
                component="LTR" if c.iap_component == "IAP_LTR" else "IAP internal",
                rna_log2fc=rna_lookup.get(c.locus_id, np.nan),
                is_selected_component=(c.locus_id == locus["locus_id"]),
            )
            rows.append(row)

        hits = dmr[dmr.locus_id.eq(locus["locus_id"])]
        for d in hits.itertuples(index=False):
            row = empty_record(locus)
            row.update(
                record_type="dmr",
                start0=int(d.dmr_start0),
                end0=int(d.dmr_end0),
                dmr_direction=str(d.direction),
                dmr_meth_diff=float(d.meth_diff),
                dmr_qvalue=float(d.qvalue),
            )
            rows.append(row)

    nrros_deg = deg[deg.gene_name.eq("Nrros")]
    nrros_fpkm = fpkm[fpkm.gene_name.eq("Nrros")]
    if len(nrros_deg) != 1 or len(nrros_fpkm) != 1:
        raise RuntimeError("Expected one Nrros record in both delivered tables")
    d = nrros_deg.iloc[0]
    f = nrros_fpkm.iloc[0]
    ratios = np.array([np.log2(float(d[f"{rep}6K"]) / float(d[f"{rep}6M"])) for rep in "ABC"])
    condition_means = [
        np.mean([float(f[f"{rep}6M"]) for rep in "ABC"]),
        np.mean([float(f[f"{rep}6K"]) for rep in "ABC"]),
    ]
    row = empty_record(LOCI[0])
    row.update(
        record_type="host_summary",
        host_gene="Nrros",
        host_mean_log2_ratio=float(ratios.mean()),
        host_sem_log2_ratio=float(ratios.std(ddof=1) / np.sqrt(len(ratios))),
        host_mean_fpkm=float(max(condition_means)),
        host_padj=float(d.padj),
        host_n_pairs=len(ratios),
    )
    rows.append(row)

    out = pd.DataFrame(rows)
    args.output.parent.mkdir(parents=True,exist_ok=True)
    out.to_csv(args.output, sep="\t", index=False, na_rep="NA")
    print(f"Wrote {len(out)} locus-example records to {args.output}")


if __name__ == "__main__":
    main()
