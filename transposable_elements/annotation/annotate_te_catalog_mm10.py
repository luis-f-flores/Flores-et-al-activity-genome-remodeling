"""Join the TE locus catalog to UCSC RepeatMasker names, classes and families."""
from __future__ import annotations
import argparse,csv,gzip,io,json,math,os
from collections import Counter,defaultdict
from pathlib import Path
def normalized_strand(value: str) -> str:
    return "-" if value in {"-", "C"} else "+"

def coord_key(chrom: str, start: str | int, end: str | int, strand: str) -> tuple[str, int, int, str]:
    return chrom, int(start), int(end), normalized_strand(strand)

def target_mapping_rule(raw_name: str, target: str) -> str:
    raw = raw_name.upper()
    tgt = target.upper()
    if raw == tgt:
        return "exact_name"
    if raw.startswith("B1") and tgt in {"B1", "B1F"}:
        return "SINE_B1_harmonization"
    if raw.startswith("B2") and tgt == "B2":
        return "SINE_B2_harmonization"
    if "IAP" in raw or "IAP" in tgt:
        return "IAP_alias_or_component_harmonization"
    if raw.startswith("RLTR") and tgt.startswith("RLTR"):
        return "RLTR_umbrella_harmonization"
    if raw == "MERVL" and tgt == "MERVL_LTR":
        return "MERVL_LTR_alias"
    return "other_SalmonTE_harmonization"

def iap_component(raw_name: str, target: str) -> str:
    raw = raw_name.upper()
    tgt = target.upper()
    if "IAP" not in raw and "IAP" not in tgt:
        return "not_IAP"
    if "LTR" in raw or "LTR" in tgt:
        return "IAP_LTR"
    if raw.endswith("-INT") or tgt.endswith("_I") or tgt.endswith("I"):
        return "IAP_internal"
    return "IAP_other"

def open_deterministic_gzip_text(path: Path):
    raw = path.open("wb")
    gz = gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0)
    return io.TextIOWrapper(gz, encoding="utf-8", newline="")

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--catalog',type=Path,required=True,help='Directory containing TE_meta.tsv and SalmonTE_mm_families.txt.')
    p.add_argument('--rmsk',type=Path,required=True)
    p.add_argument('--output',type=Path,required=True)
    args=p.parse_args();catalog=args.catalog
    te_meta=catalog/'TE_meta.tsv';salmonte_families=catalog/'SalmonTE_mm_families.txt'
    output=args.output/'TE_annotation_master_mm10.tsv.gz'
    if args.output.exists():raise FileExistsError(f'Output directory already exists: {args.output}')
    for path in [te_meta,salmonte_families,args.rmsk]:
        if not path.is_file():raise FileNotFoundError(path)
    args.output.mkdir(parents=True)
    temporary_output=output.with_suffix(output.suffix+f'.tmp.{os.getpid()}')
    salmonte_targets={line.strip().upper() for line in salmonte_families.read_text().splitlines() if line.strip()}
    wanted_keys: set[tuple[str, int, int, str]] = set()
    meta_rows = 0
    meta_ids_seen: set[str] = set()
    duplicate_meta_ids = 0
    duplicate_coordinate_keys = 0
    coord_counts: Counter = Counter()
    with te_meta.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required_columns = {"copyID", "family", "repClass", "repFamily", "percDiv", "length", "strand", "chr", "start", "end", "repStart", "repEnd"}
        if not required_columns.issubset(reader.fieldnames or []):
            raise RuntimeError(f"Unexpected TE_meta columns: {reader.fieldnames}")
        for row in reader:
            meta_rows += 1
            locus_id = row["copyID"]
            if locus_id in meta_ids_seen:
                duplicate_meta_ids += 1
            meta_ids_seen.add(locus_id)
            key = coord_key(row["chr"], row["start"], row["end"], row["strand"])
            coord_counts[key] += 1
            wanted_keys.add(key)
    duplicate_coordinate_keys = sum(1 for count in coord_counts.values() if count > 1)

    raw_by_key: dict[tuple[str, int, int, str], list[dict]] = defaultdict(list)
    raw_rows_scanned = 0
    with gzip.open(args.rmsk, "rt", encoding="utf-8", newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for fields in reader:
            if len(fields) < 17:
                continue
            raw_rows_scanned += 1
            key = coord_key(fields[5], fields[6], fields[7], fields[9])
            if key not in wanted_keys:
                continue
            raw_by_key[key].append({
                "rmsk_id": fields[16],
                "raw_repeat_name": fields[10],
                "raw_repeat_class": fields[11],
                "raw_repeat_family": fields[12],
                "raw_repeat_start": fields[13],
                "raw_repeat_end": fields[14],
                "raw_repeat_left": fields[15],
                "sw_score": fields[1],
                "milli_div": fields[2],
                "milli_del": fields[3],
                "milli_ins": fields[4],
            })

    fields_out = [
        "locus_id", "chrom", "start_0based", "end_0based_exclusive", "strand",
        "rmsk_id", "raw_repeat_name", "raw_repeat_class", "raw_repeat_family",
        "salmonte_target", "target_mapping_rule", "target_name_changed",
        "target_in_salmonte_reference", "is_iap", "iap_component",
        "percent_divergence", "length_bp", "raw_repeat_start", "raw_repeat_end",
        "mapping_status", "coordinate_match_count",
    ]
    status_counts: Counter = Counter()
    mapping_rule_counts: Counter = Counter()
    raw_to_target_counts: Counter = Counter()
    target_meta: dict[str, dict[str, set | int]] = {}
    master_ids: set[str] = set()
    unresolved_rows: list[dict] = []

    with te_meta.open("r", encoding="utf-8", newline="") as source, open_deterministic_gzip_text(temporary_output) as target_handle:
        reader = csv.DictReader(source, delimiter="\t")
        writer = csv.DictWriter(target_handle, fieldnames=fields_out, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for row in reader:
            locus_id = row["copyID"]
            key = coord_key(row["chr"], row["start"], row["end"], row["strand"])
            candidates = raw_by_key.get(key, [])

            def score(candidate: dict) -> int:
                value = 0
                value += 10 if candidate["raw_repeat_class"] == row["repClass"] else 0
                value += 10 if candidate["raw_repeat_family"] == row["repFamily"] else 0
                value += 5 if int(candidate["raw_repeat_start"]) == int(row["repStart"]) else 0
                value += 5 if int(candidate["raw_repeat_end"]) == int(row["repEnd"]) else 0
                value += 5 if math.isclose(int(candidate["milli_div"]) / 10.0, float(row["percDiv"]), abs_tol=1e-6) else 0
                return value

            chosen = None
            mapping_status = "unmatched"
            if candidates:
                scored = [(score(candidate), candidate) for candidate in candidates]
                best_score = max(value for value, _ in scored)
                best = [candidate for value, candidate in scored if value == best_score]
                if len(best) == 1:
                    chosen = best[0]
                    mapping_status = "exact_unique" if best_score == 35 else "best_unique"
                else:
                    signatures = {
                        (x["raw_repeat_name"], x["raw_repeat_class"], x["raw_repeat_family"], x["raw_repeat_start"], x["raw_repeat_end"])
                        for x in best
                    }
                    if len(signatures) == 1:
                        chosen = sorted(best, key=lambda x: int(x["rmsk_id"]))[0]
                        mapping_status = "equivalent_duplicate"
                    else:
                        mapping_status = "ambiguous"

            if chosen is None:
                chosen = {
                    "rmsk_id": "", "raw_repeat_name": "", "raw_repeat_class": "",
                    "raw_repeat_family": "", "raw_repeat_start": "", "raw_repeat_end": "",
                }
                unresolved_rows.append({
                    "locus_id": locus_id,
                    "chrom": row["chr"], "start": row["start"], "end": row["end"],
                    "strand": row["strand"], "salmonte_target": row["family"],
                    "status": mapping_status, "coordinate_match_count": len(candidates),
                })

            raw_name = chosen["raw_repeat_name"]
            rule = target_mapping_rule(raw_name, row["family"]) if raw_name else "unresolved"
            is_iap_value = "TRUE" if ("IAP" in raw_name.upper() or "IAP" in row["family"].upper()) else "FALSE"
            out = {
                "locus_id": locus_id,
                "chrom": row["chr"],
                "start_0based": row["start"],
                "end_0based_exclusive": row["end"],
                "strand": normalized_strand(row["strand"]),
                "rmsk_id": chosen["rmsk_id"],
                "raw_repeat_name": raw_name,
                "raw_repeat_class": chosen["raw_repeat_class"],
                "raw_repeat_family": chosen["raw_repeat_family"],
                "salmonte_target": row["family"],
                "target_mapping_rule": rule,
                "target_name_changed": "TRUE" if raw_name.upper() != row["family"].upper() else "FALSE",
                "target_in_salmonte_reference": "TRUE" if row["family"].upper() in salmonte_targets else "FALSE",
                "is_iap": is_iap_value,
                "iap_component": iap_component(raw_name, row["family"]) if raw_name else "unresolved",
                "percent_divergence": row["percDiv"],
                "length_bp": row["length"],
                "raw_repeat_start": chosen["raw_repeat_start"],
                "raw_repeat_end": chosen["raw_repeat_end"],
                "mapping_status": mapping_status,
                "coordinate_match_count": len(candidates),
            }
            writer.writerow(out)
            master_ids.add(locus_id)
            status_counts[mapping_status] += 1
            mapping_rule_counts[rule] += 1
            raw_to_target_counts[(raw_name, chosen["raw_repeat_class"], chosen["raw_repeat_family"], row["family"], rule)] += 1

            target_key = row["family"]
            if target_key not in target_meta:
                target_meta[target_key] = {"raw_names": set(), "raw_classes": set(), "raw_families": set(), "loci": 0, "iap_loci": 0}
            target_meta[target_key]["raw_names"].add(raw_name)
            target_meta[target_key]["raw_classes"].add(chosen["raw_repeat_class"])
            target_meta[target_key]["raw_families"].add(chosen["raw_repeat_family"])
            target_meta[target_key]["loci"] += 1
            target_meta[target_key]["iap_loci"] += int(is_iap_value == "TRUE")


    temporary_output.replace(output)

    tables = args.output / "tables"
    tables.mkdir()
    mapping_path = tables / "RAW_REPEATMASKER_TO_SALMONTE_TARGET.tsv"
    with mapping_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["raw_repeat_name", "raw_repeat_class", "raw_repeat_family", "salmonte_target", "mapping_rule", "n_loci"])
        for key, count in sorted(raw_to_target_counts.items(), key=lambda item: (item[0][3], item[0][0])):
            writer.writerow([*key, count])

    target_summary_path = tables / "SALMONTE_TARGET_ANNOTATION_SUMMARY.tsv"
    with target_summary_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["salmonte_target", "n_loci", "n_iap_loci", "n_raw_repeat_names", "raw_repeat_names", "raw_classes", "raw_families", "annotation_ambiguous"])
        for target_name in sorted(target_meta):
            item = target_meta[target_name]
            names = sorted(x for x in item["raw_names"] if x)
            classes = sorted(x for x in item["raw_classes"] if x)
            families = sorted(x for x in item["raw_families"] if x)
            ambiguous = len(classes) != 1 or len(families) != 1
            writer.writerow([target_name, item["loci"], item["iap_loci"], len(names), ";".join(names), ";".join(classes), ";".join(families), str(ambiguous).upper()])

    unresolved_path = tables / "UNMATCHED_OR_AMBIGUOUS_MASTER_LOCI.tsv"
    with unresolved_path.open("w", encoding="utf-8", newline="") as handle:
        fields = ["locus_id", "chrom", "start", "end", "strand", "salmonte_target", "status", "coordinate_match_count"]
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(unresolved_rows)


    if duplicate_meta_ids:raise ValueError('TE catalog locus identifiers must be unique.')
    summary=dict(master_rows=meta_rows,unique_master_locus_ids=len(master_ids),duplicate_coordinate_keys=duplicate_coordinate_keys,
      raw_rmsk_rows_scanned=raw_rows_scanned,mapping_status_counts=dict(status_counts),target_mapping_rule_counts=dict(mapping_rule_counts),
      salmonte_reference_targets=len(salmonte_targets),salmonte_targets_with_catalog_loci=len(target_meta),unmatched_or_ambiguous_master_loci=len(unresolved_rows))
    (args.output/'annotation_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    print(f'Annotated {meta_rows:,} TE loci; {len(unresolved_rows):,} unmatched or ambiguous loci.')

if __name__=='__main__':main()
