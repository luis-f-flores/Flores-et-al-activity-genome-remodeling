"""Build stable gene loci, exons and TSS coordinates from a GTF annotation."""
from __future__ import annotations
from pathlib import Path
import argparse
import collections
import gzip
import re

ATTR_RE = re.compile(r'(\S+)\s+"([^"]*)"')

def chrom_alias(value: str) -> str:
    if value.startswith("chr"):
        return value
    if value in {"MT", "M"}:
        return "chrM"
    return "chr" + value

def load_fai(path: Path) -> dict[str, int]:
    lengths = {}
    with path.open() as handle:
        for line in handle:
            fields = line.rstrip("\n").split("\t")
            if len(fields) >= 2:
                lengths[chrom_alias(fields[0])] = int(fields[1])
    if not lengths:
        raise ValueError(f"No chromosome lengths in {path}")
    return lengths

def merge_intervals(values: list[tuple[int, int]]) -> list[tuple[int, int]]:
    merged: list[list[int]] = []
    for start, end in sorted(set(values)):
        if not merged or start > merged[-1][1]:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    return [(start, end) for start, end in merged]

def parse_gtf(gtf: Path, chrom_lengths: dict[str, int]):
    opener = gzip.open if gtf.suffix == ".gz" else open
    transcripts: dict[tuple[str, str, str, str], dict] = {}
    gene_rows: dict[tuple[str, str, str], list[tuple[int, int, str]]] = collections.defaultdict(list)
    counters = collections.Counter()

    with opener(gtf, "rt") as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"Malformed GTF row {line_number}")
            feature = fields[2]
            if feature not in {"gene", "transcript", "exon"}:
                continue
            chrom = chrom_alias(fields[0])
            if chrom not in chrom_lengths:
                counters["rows_unknown_chrom"] += 1
                continue
            start0, end0 = int(fields[3]) - 1, int(fields[4])
            strand = fields[6]
            if strand not in {"+", "-"} or not (0 <= start0 < end0 <= chrom_lengths[chrom]):
                raise ValueError(f"Invalid GTF coordinates at row {line_number}")
            attrs = dict(ATTR_RE.findall(fields[8]))
            raw_gid = attrs.get("gene_id")
            if not raw_gid:
                counters["rows_missing_gene_id"] += 1
                continue
            gene_name = attrs.get("gene_name") or attrs.get("gene") or raw_gid
            group_key = (chrom, strand, raw_gid)
            if feature == "gene":
                gene_rows[group_key].append((start0, end0, gene_name))
                counters["gene_rows"] += 1
                continue
            transcript_id = attrs.get("transcript_id")
            if not transcript_id:
                counters["rows_missing_transcript_id"] += 1
                continue
            key = (chrom, strand, raw_gid, transcript_id)
            item = transcripts.setdefault(
                key,
                {"chrom": chrom, "strand": strand, "raw_gene_id": raw_gid,
                 "transcript_id": transcript_id, "gene_names": [], "span": None, "exons": []},
            )
            item["gene_names"].append(gene_name)
            if item["span"] is None:
                item["span"] = [start0, end0]
            else:
                item["span"][0] = min(item["span"][0], start0)
                item["span"][1] = max(item["span"][1], end0)
            if feature == "exon":
                item["exons"].append((start0, end0))
                counters["exon_rows"] += 1
            else:
                counters["transcript_rows"] += 1

    grouped: dict[tuple[str, str, str], list[dict]] = collections.defaultdict(list)
    for (chrom, strand, raw_gid, _tid), tx in transcripts.items():
        if tx["span"] is not None:
            grouped[(chrom, strand, raw_gid)].append(tx)

    all_components = []
    component_counts = collections.Counter()
    for group_key, txs in grouped.items():
        chrom, strand, raw_gid = group_key
        txs.sort(key=lambda x: (x["span"][0], x["span"][1], x["transcript_id"]))
        components: list[list[dict]] = []
        current: list[dict] = []
        current_end = -1
        for tx in txs:
            start0, end0 = tx["span"]
            if current and start0 > current_end:
                components.append(current)
                current = []
                current_end = -1
            current.append(tx)
            current_end = max(current_end, end0)
        if current:
            components.append(current)
        component_counts[raw_gid] += len(components)
        for component in components:
            start0 = min(tx["span"][0] for tx in component)
            end0 = max(tx["span"][1] for tx in component)
            # A source gene row is safe only when the raw ID has one connected
            # component on this chromosome/strand.
            if len(components) == 1 and gene_rows.get(group_key):
                start0 = min(start0, min(row[0] for row in gene_rows[group_key]))
                end0 = max(end0, max(row[1] for row in gene_rows[group_key]))
            names = [name for tx in component for name in tx["gene_names"]]
            if gene_rows.get(group_key):
                names.extend(row[2] for row in gene_rows[group_key])
            name_counts = collections.Counter(names)
            gene_name = sorted(name_counts, key=lambda name: (-name_counts[name], name))[0]
            exons = merge_intervals([interval for tx in component for interval in tx["exons"]])
            all_components.append({
                "chrom": chrom, "strand": strand, "raw_gene_id": raw_gid,
                "gene_name": gene_name, "start0": start0, "end0": end0,
                "transcript_ids": sorted({tx["transcript_id"] for tx in component}),
                "exons": exons,
            })

    raw_id_total_components = collections.Counter(x["raw_gene_id"] for x in all_components)
    for item in all_components:
        if raw_id_total_components[item["raw_gene_id"]] == 1:
            item["locus_id"] = item["raw_gene_id"]
        else:
            item["locus_id"] = (
                f'{item["raw_gene_id"]}|{item["chrom"]}:'
                f'{item["start0"]}-{item["end0"]}:{item["strand"]}'
            )
    return all_components, counters

def write_layers(out_dir: Path, components: list[dict], chrom_lengths: dict[str, int]):
    paths = {
        "genes": out_dir / "genes.bed",
        "exons": out_dir / "exons.bed",
        "introns": out_dir / "introns.bed",
        "tss": out_dir / "tss.bed",
        "promoters": out_dir / "promoters_-1kb_+100bp.bed",
        "metadata": out_dir / "gene_metadata.tsv",
    }
    rows = {name: [] for name in ("genes", "exons", "introns", "tss", "promoters")}
    metadata = []
    for item in components:
        chrom, start0, end0 = item["chrom"], item["start0"], item["end0"]
        locus_id, strand, gene_name = item["locus_id"], item["strand"], item["gene_name"]
        base = (chrom, start0, end0, locus_id, 0, strand, gene_name)
        rows["genes"].append(base)
        for left, right in item["exons"]:
            rows["exons"].append((chrom, left, right, locus_id, 0, strand, gene_name))
        for (_, left_end), (right_start, _) in zip(item["exons"], item["exons"][1:]):
            if right_start > left_end:
                rows["introns"].append((chrom, left_end, right_start, locus_id, 0, strand, gene_name))
        if strand == "+":
            tss_start = start0
            promoter_start, promoter_end = start0 - 1000, start0 + 100
        else:
            tss_start = end0 - 1
            # Strand-aware promoter coordinates use the transcript boundary.
            promoter_start, promoter_end = end0 - 100, end0 + 1000
        rows["tss"].append((chrom, tss_start, tss_start + 1, locus_id, 0, strand, gene_name))
        promoter_start = max(0, promoter_start)
        promoter_end = min(chrom_lengths[chrom], promoter_end)
        if promoter_end > promoter_start:
            rows["promoters"].append((chrom, promoter_start, promoter_end, locus_id, 0, strand, gene_name))
        metadata.append((
            locus_id, item["raw_gene_id"], gene_name, chrom, start0, end0, strand,
            len(item["transcript_ids"]), len(item["exons"]), ",".join(item["transcript_ids"]),
        ))

    for name, values in rows.items():
        values.sort(key=lambda x: (x[0], x[1], x[2], x[3]))
        with paths[name].open("x") as handle:
            for value in values:
                handle.write("\t".join(map(str, value)) + "\n")
    metadata.sort(key=lambda x: (x[3], x[4], x[5], x[0]))
    with paths["metadata"].open("x") as handle:
        handle.write("locus_id\traw_gene_id\tgene_name\tchrom\tstart0\tend0\tstrand\tn_transcripts\tn_merged_exons\ttranscript_ids\n")
        for value in metadata:
            handle.write("\t".join(map(str, value)) + "\n")
    return paths, {name: len(values) for name, values in rows.items()}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--gtf',type=Path,required=True)
    parser.add_argument('--fai',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    args.output.mkdir(parents=True,exist_ok=True)
    lengths=load_fai(args.fai)
    components,_=parse_gtf(args.gtf,lengths)
    write_layers(args.output,components,lengths)
    print(f'Wrote {len(components):,} gene loci.')

if __name__=='__main__':
    main()
