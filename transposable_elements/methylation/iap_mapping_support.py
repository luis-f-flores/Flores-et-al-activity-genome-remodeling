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
import collections
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import subprocess
import sys
import time
import xml.etree.ElementTree as ET
import zipfile

NS = {"s": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}

def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()

def metadata(path, full_hash=False):
    path = Path(path)
    stat = path.stat()
    result = {"path": str(path.resolve()), "bytes": stat.st_size,
              "mtime_ns": stat.st_mtime_ns, "sha256": sha256(path) if full_hash else None}
    if not full_hash:
        result["hash_scope"] = "File metadata only; full BAM was not hashed."
    return result

def read_families(path):
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = ["".join(node.itertext()) for node in
                       ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("s:si", NS)]
        relationships = ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        targets = {node.attrib["Id"]: node.attrib["Target"] for node in relationships}
        workbook = ET.fromstring(archive.read("xl/workbook.xml"))
        target = None
        for sheet in workbook.findall("s:sheets/s:sheet", NS):
            if sheet.attrib["name"] == "02_IAP_families":
                rel_id = sheet.attrib["{http://schemas.openxmlformats.org/officeDocument/2006/relationships}id"]
                target = targets[rel_id]
        if target is None:
            raise ValueError("Existing workbook lacks sheet 02_IAP_families")
        target = target.lstrip("/") if target.startswith("/") else "xl/" + target
        rows = []
        for row in ET.fromstring(archive.read(target)).findall(".//s:row", NS):
            values = {}
            for cell in row.findall("s:c", NS):
                column = re.sub(r"\d", "", cell.attrib["r"])
                value = cell.find("s:v", NS)
                if cell.get("t") == "s" and value is not None:
                    text = strings[int(value.text)]
                elif cell.get("t") == "inlineStr":
                    text = "".join(cell.find("s:is", NS).itertext())
                else:
                    text = value.text if value is not None else ""
                values[column] = text
            rows.append(values)
    column = next(key for key, value in rows[0].items() if value == "family")
    families = [row[column] for row in rows[1:] if row.get(column)]
    if not families or len(families) != len(set(families)):
        raise ValueError("Empty or duplicate family identifiers in existing workbook")
    return families

def merge_intervals(items):
    result = []
    for chrom, start, end in sorted(items):
        if result and result[-1][0] == chrom and start <= result[-1][2]:
            result[-1] = (chrom, result[-1][1], max(end, result[-1][2]))
        else:
            result.append((chrom, start, end))
    return result

def select_regions(bed, families, count):
    groups = {family: [] for family in families}
    with bed.open() as stream:
        for line in stream:
            fields = line.rstrip().split("\t")
            family = fields[3].split("@", 1)[0]
            if family in groups and re.fullmatch(r"chr(?:[1-9]|1[0-9]|X|Y)", fields[0]):
                groups[family].append((fields[0][3:], int(fields[1]), int(fields[2])))
    selected, eligible = [], {}
    for family, items in groups.items():
        merged = merge_intervals(items)
        eligible[family] = len(merged)
        if not merged:
            raise ValueError("No primary chromosome intervals for " + family)
        n = min(count, len(merged))
        indices = sorted({round(i * (len(merged) - 1) / (n - 1)) for i in range(n)}) if n > 1 else [0]
        selected.extend((family, *merged[index]) for index in indices)
    union = merge_intervals([(chrom, start, end) for _, chrom, start, end in selected])
    intervals = collections.defaultdict(list)
    for chrom, start, end in union:
        intervals[chrom].append((start, end))
    return selected, eligible, intervals

def inspect_bam(path, intervals, overlap_minimum, timeout):
    regions = [chrom + ":" + str(start + 1) + "-" + str(end)
               for chrom in sorted(intervals) for start, end in intervals[chrom]]
    command = ["samtools", "view", "-F", "2308", str(path), *regions]
    run = subprocess.run(command, capture_output=True, text=True, check=True, timeout=timeout)
    fragments, seen, spliced = {}, set(), set()
    missing = 0
    for line in run.stdout.splitlines():
        fields = line.split("\t")
        key = (fields[0], int(fields[1]) & 192, fields[2], fields[3], fields[5])
        if key in seen:
            continue
        seen.add(key)
        position = int(fields[3]) - 1
        overlap = 0
        for number, operation in re.findall(r"(\d+)([MIDNSHP=X])", fields[5]):
            length = int(number)
            if operation in "M=X":
                overlap += sum(max(0, min(position + length, end) - max(position, start))
                               for start, end in intervals[fields[2]])
            if operation in "MDN=X":
                position += length
        if overlap < overlap_minimum:
            continue
        tags = {tag.split(":", 2)[0]: tag.split(":", 2)[2]
                for tag in fields[11:] if tag.count(":") >= 2}
        if "NH" not in tags:
            missing += 1
            continue
        previous = fragments.get(fields[0], {"NH": 0, "mapq": 0})
        fragments[fields[0]] = {"NH": max(int(tags["NH"]), previous["NH"]),
                                "mapq": max(int(fields[4]), previous["mapq"])}
        if "N" in fields[5]:
            spliced.add(fields[0])
    histogram = collections.Counter(item["NH"] for item in fragments.values())
    n = len(fragments)
    at_least_five = sum(value for nh, value in histogram.items() if nh >= 5)
    return {"sample": path.stem, "represented_primary_fragments": n,
            "NH_histogram": dict(sorted(histogram.items())),
            "NH_ge5": at_least_five,
            "NH_ge5_percent": round(100 * at_least_five / n, 2) if n else None,
            "spliced_fragments": len(spliced), "mapped_records_missing_NH": missing}

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rna-root", type=Path,
                        default=Path(gb_path('rna')))
    parser.add_argument("--bam-dir", type=Path,
                        default=Path(gb_path('rna/BAM')))
    parser.add_argument("--run-dir", type=Path, required=True)
    parser.add_argument("--regions-per-family", type=int, default=40)
    parser.add_argument("--min-aligned-overlap", type=int, default=30)
    parser.add_argument("--timeout-seconds", type=int, default=40)
    args = parser.parse_args()
    if args.run_dir.exists():
        parser.error("--run-dir must be a new directory; completed evidence is preserved")
    if args.regions_per_family < 1 or args.min_aligned_overlap < 1:
        parser.error("Sampling and overlap settings must be positive")
    started = time.monotonic()
    workbook = args.rna_root / "SalmonTE_DE_6h/IAP_up_6h_location_3part.xlsx"
    bed = args.rna_root / "TE_catalog_mm10/TE_mm10.bed"
    families = read_families(workbook)
    selected, eligible, intervals = select_regions(bed, families, args.regions_per_family)
    samples = ["A0", "B0", "C0"] + [rep + timepoint + condition
        for timepoint in ["1", "6", "24"] for condition in ["K", "M"] for rep in ["A", "B", "C"]]
    rows, inputs = [], [metadata(workbook, True), metadata(bed, True)]
    for sample in samples:
        bam = args.bam_dir / (sample + ".bam")
        index = Path(str(bam) + ".bai")
        if not index.exists():
            index = bam.with_suffix(".bai")
        header = subprocess.run(["samtools", "view", "-H", str(bam)],
                                capture_output=True, check=True).stdout
        bam_input = metadata(bam)
        bam_input["header_sha256"] = hashlib.sha256(header).hexdigest()
        inputs.extend([bam_input, metadata(index, True)])
        rows.append(inspect_bam(bam, intervals, args.min_aligned_overlap, args.timeout_seconds))
    histogram = collections.Counter()
    for row in rows:
        histogram.update(row["NH_histogram"])
    total = sum(row["represented_primary_fragments"] for row in rows)
    output = {
        "purpose": 'Summarize read mapping over sampled IAP intervals in the supplied BAMs.',
        "family_source": str(workbook) + " sheet 02_IAP_families",
        "families": families, "eligible_merged_primary_intervals": eligible,
        "sampled_intervals": len(selected),
        "selected_intervals": [{"family": family, "bam_contig": chrom, "start_0based": start,
                               "end_exclusive": end} for family, chrom, start, end in selected],
        "sampling": "Evenly ranked merged intervals per existing family, independent of RNA signal; all 21 libraries",
        "read_rule": "Primary mapped read overlaps at least the specified number of aligned M/= /X bases; reference gaps skipped; fragments deduplicated by name within sample",
        "limits": 'Summaries describe primary alignments overlapping the sampled intervals; NH values reflect the alignments stored in the BAM.',
        "samples": rows,
        "pooled_summary": {"represented_fragments": total, "NH_histogram": dict(sorted(histogram.items())),
                          "reported_multiple": total - histogram[1],
                          "reported_NH5": histogram[5],
                          "spliced_fragments": sum(row["spliced_fragments"] for row in rows)}
    }
    args.run_dir.mkdir(parents=True)
    result_path = args.run_dir / "bam_mapping_check.json"
    result_path.write_text(json.dumps(output, indent=2) + "\n")
    version = subprocess.run(["samtools", "--version"], capture_output=True,
                             text=True, errors="replace", check=True).stdout.splitlines()[:2]
    manifest = {
        "created_utc": datetime.now(timezone.utc).isoformat(), "status": "read_only_preflight_complete",
        "script": metadata(Path(__file__), True), "command": [sys.executable, *sys.argv],
        "parameters": {key: str(value) if isinstance(value, Path) else value
                       for key, value in vars(args).items()},
        "versions": {"python": sys.version, "samtools": version},
        "inputs": inputs, "outputs": [metadata(result_path, True)],
        "elapsed_seconds": round(time.monotonic() - started, 3),
        "scientific_quantification_run": False, "realignment_run": False,
        "scope": "IAP mapping summaries for sampled intervals; BAM file metadata recorded without full-file checksums",
        "publication_export": False
    }
    (args.run_dir / "RUN_MANIFEST.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(json.dumps(output["pooled_summary"]))
    print("Evidence saved in " + str(args.run_dir))

if __name__ == "__main__":
    main()
