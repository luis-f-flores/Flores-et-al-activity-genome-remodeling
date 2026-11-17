#!/usr/bin/env python3
"""Combine six-hour TE expression, CpG methylation, DMR overlap and gene context.

TE identifiers and gene-layer BED files use 0-based half-open coordinates.
Merged DMR tables use 1-based inclusive coordinates. Condition-mean mCG
bigWigs contain methylation percentages at assayed CpG positions.
"""
from __future__ import annotations
import argparse
import bisect
import csv
import gzip
from collections import Counter, defaultdict
from pathlib import Path
import numpy as np
import pyBigWig

def read_tsv(path: Path) -> list[dict[str, str]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def parse_locus(locus_id: str):
    subfamily, coordinate = locus_id.split("@", 1)
    chrom, bounds = coordinate.split(":", 1)
    start0, end0 = map(int, bounds.split("-", 1))
    if not 0 <= start0 < end0:
        raise ValueError(f"Invalid locus coordinates: {locus_id}")
    return subfamily, chrom, start0, end0


class IntervalIndex:
    def __init__(self, path: Path, bin_size: int = 1_000_000):
        self.bin_size = bin_size
        self.bins = defaultdict(list)
        with path.open() as handle:
            for line in handle:
                if not line.strip() or line.startswith("#"):
                    continue
                f = line.rstrip("\n").split("\t")
                start0, end0 = int(f[1]), int(f[2])
                item = (start0, end0, f[3], f[6] if len(f) > 6 else f[3], f[5] if len(f) > 5 else ".")
                for bin_id in range(start0 // bin_size, (end0 - 1) // bin_size + 1):
                    self.bins[(f[0], bin_id)].append(item)

    def query(self, chrom: str, start0: int, end0: int):
        found = {}
        for bin_id in range(start0 // self.bin_size, (end0 - 1) // self.bin_size + 1):
            for item in self.bins.get((chrom, bin_id), []):
                if item[0] < end0 and item[1] > start0:
                    found[item] = item
        return sorted(found.values())


class TssIndex:
    def __init__(self, path: Path):
        self.rows = defaultdict(list)
        with path.open() as handle:
            for line in handle:
                f = line.rstrip("\n").split("\t")
                self.rows[f[0]].append((int(f[1]), f[3], f[6] if len(f) > 6 else f[3]))
        self.positions = {}
        for chrom in self.rows:
            self.rows[chrom].sort()
            self.positions[chrom] = [x[0] for x in self.rows[chrom]]

    def nearest(self, chrom: str, start0: int, end0: int):
        values = self.rows.get(chrom, [])
        if not values:
            return None, None, None
        pos = self.positions[chrom]
        i = bisect.bisect_left(pos, start0)
        candidates = values[max(0, i - 2): min(len(values), i + 3)]
        best = None
        for tss, locus_id, name in candidates:
            if tss < start0:
                distance = start0 - tss
            elif tss >= end0:
                distance = tss - (end0 - 1)
            else:
                distance = 0
            key = (distance, name, locus_id)
            if best is None or key < best[0]:
                best = (key, locus_id, name, distance)
        return best[1], best[2], best[3]


def load_dmrs(path: Path, direction: str):
    out = defaultdict(list)
    for row in read_tsv(path):
        item = (int(row["start"]) - 1, int(row["end"]), direction,
                float(row["meth_diff"]), float(row["qvalue"]))
        out[row["chr"]].append(item)
    for chrom in out:
        out[chrom].sort()
    return out


def interval_hits(values, chrom: str, start0: int, end0: int):
    rows = values.get(chrom, [])
    stop = bisect.bisect_left([x[0] for x in rows], end0)
    return [x for x in rows[:stop] if x[1] > start0]


def coverage_tier(n: int) -> str:
    if n == 0:
        return "0"
    if n <= 2:
        return "1-2"
    if n <= 4:
        return "3-4"
    if n <= 9:
        return "5-9"
    return "10+"


def semicolon(items) -> str:
    return ";".join(sorted(set(items)))


def write_tsv(path: Path, rows: list[dict], fields=None, gzip_output=False):
    if not rows and fields is None:
        raise ValueError(f"No rows for {path}")
    fields = fields or list(rows[0])
    opener = gzip.open if gzip_output else open
    with opener(path, "xt", newline="") if gzip_output else opener(path, "x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    for name in ["rna-de", "rna-confidence-de", "repeat-annotation", "kcl-mcg", "mannitol-mcg",
                 "hypo-dmrs", "hyper-dmrs", "gene-de", "genes", "exons", "introns", "promoters", "tss", "output"]:
        parser.add_argument("--"+name, type=Path, required=True)
    a=parser.parse_args()
    OUT_DIR=a.output
    OUT_DIR.mkdir(parents=True,exist_ok=True)
    PRIMARY=a.rna_de
    CONFIDENT=a.rna_confidence_de
    CROSSWALK=a.repeat_annotation
    KCL_BW=a.kcl_mcg
    MANN_BW=a.mannitol_mcg
    DMR_HYPO=a.hypo_dmrs
    DMR_HYPER=a.hyper_dmrs
    GENE_DE=a.gene_de
    layer_paths={name:getattr(a,name) for name in ["genes","exons","introns","promoters","tss"]}
    genes = IntervalIndex(layer_paths["genes"])
    exons = IntervalIndex(layer_paths["exons"])
    introns = IntervalIndex(layer_paths["introns"])
    promoters = IntervalIndex(layer_paths["promoters"])
    tss = TssIndex(layer_paths["tss"])

    primary = read_tsv(PRIMARY)
    confident = {x["locus_id"]: x for x in read_tsv(CONFIDENT)}
    crosswalk = {x["locus_id"]: x for x in read_tsv(CROSSWALK)}
    if set(confident) != {x["locus_id"] for x in primary}:
        raise AssertionError("Telescope count and confidence tables must contain the same loci")
    dmrs_hypo = load_dmrs(DMR_HYPO, "hypo")
    dmrs_hyper = load_dmrs(DMR_HYPER, "hyper")

    rows, dmr_rows = [], []
    with pyBigWig.open(str(KCL_BW)) as kcl, pyBigWig.open(str(MANN_BW)) as mann:
        for source in primary:
            locus_id = source["locus_id"]
            subfamily, chrom, start0, end0 = parse_locus(locus_id)
            gene_hits = genes.query(chrom, start0, end0)
            exon_hits = exons.query(chrom, start0, end0)
            intron_hits = introns.query(chrom, start0, end0)
            promoter_hits = promoters.query(chrom, start0, end0)
            nearest_id, nearest_name, nearest_distance = tss.nearest(chrom, start0, end0)
            if exon_hits:
                detail = "exon_overlapping"
            elif not gene_hits:
                detail = "outside_gene_bodies"
            elif any(left <= start0 and right >= end0 for left, right, *_ in gene_hits):
                detail = "fully_intronic"
            else:
                detail = "gene_boundary_overlap"

            k = np.asarray(kcl.values(chrom, start0, end0), dtype=float)
            m = np.asarray(mann.values(chrom, start0, end0), dtype=float)
            shared = np.isfinite(k) & np.isfinite(m)
            shared_n = int(shared.sum())
            if shared_n:
                k_mean = float(np.mean(k[shared])); m_mean = float(np.mean(m[shared])); delta = k_mean - m_mean
            else:
                k_mean = m_mean = delta = None

            hypo_hits = interval_hits(dmrs_hypo, chrom, start0, end0)
            hyper_hits = interval_hits(dmrs_hyper, chrom, start0, end0)
            for hit in hypo_hits + hyper_hits:
                dmr_rows.append({
                    "locus_id": locus_id, "direction": hit[2], "dmr_start0": hit[0], "dmr_end0": hit[1],
                    "overlap_bp": min(end0, hit[1]) - max(start0, hit[0]), "meth_diff": hit[3], "qvalue": hit[4],
                })

            xwalk, conf = crosswalk[locus_id], confident[locus_id]
            conf_lfc = conf["log2FoldChange_shrunk"]
            conf_lfc = None if conf_lfc in {"", "NA"} else float(conf_lfc)
            rna_lfc = float(source["log2FoldChange_shrunk"])
            rna_q = source["padj_BH_all_retained"]
            rna_q = None if rna_q in {"", "NA"} else float(rna_q)
            rows.append({
                "locus_id": locus_id, "chrom": chrom, "start0": start0, "end0": end0,
                "length_bp": end0 - start0, "repeat_subfamily": xwalk["repeat_subfamily"],
                "repeat_class": xwalk["repeat_class"], "repeat_family": xwalk["repeat_family"],
                "iap_related": subfamily.startswith("IAP"),
                "gene_context": "genic" if gene_hits else "intergenic",
                "gene_context_detail": detail,
                "overlapping_gene_ids": semicolon(x[2] for x in gene_hits),
                "overlapping_genes": semicolon(x[3] for x in gene_hits),
                "overlapping_exon_gene_ids": semicolon(x[2] for x in exon_hits),
                "overlapping_exon_genes": semicolon(x[3] for x in exon_hits),
                "overlapping_intron_gene_ids": semicolon(x[2] for x in intron_hits),
                "overlapping_intron_genes": semicolon(x[3] for x in intron_hits),
                "overlapping_promoter_gene_ids": semicolon(x[2] for x in promoter_hits),
                "overlapping_promoter_genes": semicolon(x[3] for x in promoter_hits),
                "nearest_tss_gene_id": nearest_id, "nearest_tss_gene": nearest_name,
                "nearest_tss_distance_bp": nearest_distance,
                "rna_baseMean": float(source["baseMean"]), "rna_log2fc_shrunk": rna_lfc,
                "rna_log2fc_shrunk_final_conf": conf_lfc, "rna_q_BH": rna_q,
                "rna_q_lt_0_05": rna_q is not None and rna_q < 0.05,
                "rna_higher_q_lt_0_05": rna_q is not None and rna_q < 0.05 and rna_lfc > 0,
                "kcl_covered_positions": int(np.isfinite(k).sum()),
                "mannitol_covered_positions": int(np.isfinite(m).sum()),
                "shared_covered_positions": shared_n, "coverage_tier": coverage_tier(shared_n),
                "kcl_mean_pct_mcg_shared": k_mean, "mannitol_mean_pct_mcg_shared": m_mean,
                "mean_delta_pp_kcl_minus_mannitol": delta,
                "direct_6hr_hypo_dmr": bool(hypo_hits), "direct_6hr_hyper_dmr": bool(hyper_hits),
                "direct_6hr_any_dmr": bool(hypo_hits or hyper_hits),
            })

    table_path = OUT_DIR / "te_methylation_context.tsv.gz"
    dmr_path = OUT_DIR / "te_dmr_overlaps.tsv.gz"
    write_tsv(table_path, rows, gzip_output=True)
    write_tsv(dmr_path, dmr_rows,
              fields=["locus_id", "direction", "dmr_start0", "dmr_end0", "overlap_bp", "meth_diff", "qvalue"],
              gzip_output=True)

    gene_de = {x["gene_id"]: x for x in read_tsv(GENE_DE)}
    host_rows = []
    for row in rows:
        if not (row["iap_related"] and row["rna_higher_q_lt_0_05"] and row["direct_6hr_hypo_dmr"]):
            continue
        for gene_id in filter(None, row["overlapping_gene_ids"].split(";")):
            de = gene_de.get(gene_id)
            host_rows.append({
                "locus_id": row["locus_id"], "gene_id": gene_id,
                "gene_name": de["gene_name"] if de else "",
                "gene_log2FoldChange": de["log2FoldChange"] if de else "",
                "gene_padj": de["padj"] if de else "",
                "gene_significant_up_padj_lt_0_05": bool(de and float(de["padj"]) < 0.05 and float(de["log2FoldChange"]) > 0),
                "gene_significant_down_padj_lt_0_05": bool(de and float(de["padj"]) < 0.05 and float(de["log2FoldChange"]) < 0),
            })
    write_tsv(OUT_DIR / "iap_host_gene_expression.tsv", host_rows,
              fields=["locus_id", "gene_id", "gene_name", "gene_log2FoldChange", "gene_padj",
                      "gene_significant_up_padj_lt_0_05", "gene_significant_down_padj_lt_0_05"])

    covered = [x for x in rows if x["shared_covered_positions"] > 0]
    iap = [x for x in covered if x["iap_related"]]
    iap_up = [x for x in iap if x["rna_higher_q_lt_0_05"]]
    context_counts = []
    for label, subset in [("all_tested", rows), ("wgbs_covered", covered),
                          ("iap_wgbs_covered", iap), ("iap_rna_up_wgbs_covered", iap_up),
                          ("iap_rna_up_hypodmr", [x for x in iap_up if x["direct_6hr_hypo_dmr"]])]:
        counts = Counter(x["gene_context_detail"] for x in subset)
        context_counts.append({"subset": label, "n": len(subset), **{k: counts[k] for k in
            ("exon_overlapping", "fully_intronic", "gene_boundary_overlap", "outside_gene_bodies")}})
    write_tsv(OUT_DIR / "gene_context_counts.tsv", context_counts)

    print(f"Wrote methylation and gene context for {len(rows):,} TE loci to {OUT_DIR}")


if __name__ == "__main__":
    main()
