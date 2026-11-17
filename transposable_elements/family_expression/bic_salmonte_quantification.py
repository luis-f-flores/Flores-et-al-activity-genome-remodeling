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


import argparse
import csv
import datetime as dt
import fcntl
import hashlib
import json
import os
import shutil
import subprocess
import sys
from pathlib import Path


PACKAGE = GB_WORK_ROOT / 'rna_package'
INPUTS = PACKAGE / "inputs" / "P2609-01_fastq"
RUN = PACKAGE / "runs" / "2026-09-19_p2609_salmonte"
INDEX = Path(gb_path('references/salmonte_mm'))
SALMON = Path(gb_path('@tool/salmon'))
KCL_META = Path(gb_path('rna/SalmonTE_out_6h/A6K/aux_info/meta_info.json'))
EXPECTED_CONDITIONS = {"TTX": 3, "Bic": 3, "Bic_Nimo": 3}
THREADS = 4


def now() -> str:
    return dt.datetime.now(dt.timezone.utc).isoformat(timespec="seconds")


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def atomic_json(path: Path, payload: dict) -> None:
    temporary = path.with_name(path.name + ".tmp")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    os.replace(temporary, path)


def input_metadata() -> tuple[list[dict], dict]:
    sheet = INPUTS / "samples.tsv"
    checksum_list = INPUTS / "expected_fastq.md5"
    with sheet.open(newline="", encoding="utf-8") as handle:
        samples = list(csv.DictReader(handle, delimiter="\t"))
    if len(samples) != 9 or len({row["sample_id"] for row in samples}) != 9:
        raise ValueError("Expected exactly nine unique samples")
    by_condition = {condition: sum(row["condition"] == condition for row in samples)
                    for condition in EXPECTED_CONDITIONS}
    if by_condition != EXPECTED_CONDITIONS:
        raise ValueError(f"Unexpected condition counts: {by_condition}")
    checksums = {}
    for line in checksum_list.read_text(encoding="utf-8").splitlines():
        digest, filename = line.split(maxsplit=1)
        checksums[filename.lstrip("*")] = digest.lower()
    if len(checksums) != 18:
        raise ValueError("Expected exactly 18 vendor-derived MD5 entries")
    seen = set()
    for sample in samples:
        sid = sample["sample_id"]
        for mate in (1, 2):
            filename = Path(sample[f"read{mate}"]).name
            if filename != f"{sid}_{mate}.fq.gz":
                raise ValueError(f"Unexpected mate filename for {sid}: {filename}")
            path = INPUTS / filename
            if not path.is_file() or path.stat().st_size == 0 or filename not in checksums:
                raise ValueError(f"Missing, empty, or unregistered FASTQ: {path}")
            seen.add(filename)
    actual = {path.name for path in INPUTS.glob("*.fq.gz")}
    if seen != set(checksums) or seen != actual:
        raise ValueError("FASTQ names differ between samples.tsv, MD5 list and intake folder")
    return samples, {
        "sample_sheet": str(sheet),
        "sample_sheet_sha256": sha256(sheet),
        "md5_list": str(checksum_list),
        "md5_list_sha256": sha256(checksum_list),
        "fastq_count": len(seen),
        "fastq_bytes": sum((INPUTS / filename).stat().st_size for filename in seen),
        "expected_md5": {name: checksums[name] for name in sorted(seen)},
    }


def reference_metadata() -> dict:
    index_info = json.loads((INDEX / "info.json").read_text(encoding="utf-8"))
    old = json.loads(KCL_META.read_text(encoding="utf-8"))
    if index_info["SeqHash"] != old["index_seq_hash"] or index_info["NameHash"] != old["index_name_hash"]:
        raise ValueError("SalmonTE index differs from saved KCl A6K reference")
    version = subprocess.check_output([str(SALMON), "--version"], text=True).strip()
    if version != "salmon 1.10.3" or old["salmon_version"] != "1.10.3":
        raise ValueError(f"Salmon version differs from KCl run: {version}")
    if old["num_valid_targets"] != 301:
        raise ValueError('KCl target table must contain 301 targets')
    return {
        "index": str(INDEX),
        "index_info_sha256": sha256(INDEX / "info.json"),
        "index_seq_hash": index_info["SeqHash"],
        "index_name_hash": index_info["NameHash"],
        "kcl_meta": str(KCL_META),
        "kcl_meta_sha256": sha256(KCL_META),
        "kcl_target_count": old["num_valid_targets"],
        "kcl_inferred_library_type": old["library_types"],
        "salmon_executable": str(SALMON),
        "salmon_version": version,
        "library_type_option": "A",
        "threads_per_sample": THREADS,
    }


def validate_quant(directory: Path) -> dict:
    quant = directory / "quant.sf"
    metadata = directory / "aux_info" / "meta_info.json"
    if not quant.is_file() or not metadata.is_file():
        raise ValueError(f"Missing Salmon output in {directory}")
    with quant.open(encoding="utf-8") as handle:
        header = handle.readline().strip()
        rows = sum(1 for _ in handle)
    if header != "Name\tLength\tEffectiveLength\tTPM\tNumReads" or rows != 301:
        raise ValueError(f"Unexpected quant.sf schema or target count in {directory}: {rows}")
    meta = json.loads(metadata.read_text(encoding="utf-8"))
    if meta.get("num_valid_targets") != 301 or meta.get("num_processed", 0) <= 0:
        raise ValueError(f"Unexpected Salmon metadata in {directory}")
    return {
        "num_processed": meta["num_processed"],
        "num_mapped": meta.get("num_mapped"),
        "percent_mapped": meta.get("percent_mapped"),
        "library_types": meta.get("library_types"),
        "target_count": rows,
        "quant_sha256": sha256(quant),
        "meta_sha256": sha256(metadata),
    }


def write_summary(manifest: dict) -> None:
    summary = RUN / "QC_SUMMARY.tsv"
    fields = ["sample_id", "condition", "status", "num_processed", "num_mapped",
              "percent_mapped", "library_types", "target_count", "quant_sha256"]
    with summary.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t")
        writer.writeheader()
        for sid, record in manifest["samples"].items():
            qc = record.get("qc", {})
            writer.writerow({
                "sample_id": sid, "condition": record["condition"], "status": record["status"],
                "num_processed": qc.get("num_processed", ""),
                "num_mapped": qc.get("num_mapped", ""),
                "percent_mapped": qc.get("percent_mapped", ""),
                "library_types": ",".join(qc.get("library_types") or []),
                "target_count": qc.get("target_count", ""),
                "quant_sha256": qc.get("quant_sha256", ""),
            })


def run() -> None:
    RUN.mkdir(parents=True, exist_ok=True)
    lock_handle = (RUN / "run.lock").open("w")
    try:
        fcntl.flock(lock_handle, fcntl.LOCK_EX | fcntl.LOCK_NB)
    except BlockingIOError as exc:
        raise RuntimeError("Another bicuculline SalmonTE controller is already running") from exc
    samples, inputs = input_metadata()
    reference = reference_metadata()
    free_bytes = shutil.disk_usage(RUN).free
    if free_bytes < 20 * 1024**3:
        raise RuntimeError(f"Less than 20 GiB free at run destination: {free_bytes}")
    manifest_path = RUN / "RUN_MANIFEST.json"
    if manifest_path.exists():
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        if manifest["inputs"]["md5_list_sha256"] != inputs["md5_list_sha256"]:
            raise ValueError("Input MD5 manifest changed since previous launch")
        if manifest["reference"]["index_seq_hash"] != reference["index_seq_hash"]:
            raise ValueError("Reference changed since previous launch")
    else:
        manifest = {
            "run_id": RUN.name,
            "analysis": "P2609-01 SalmonTE family-level quantification",
            "source_script": str(Path(__file__).resolve()),
            "source_script_sha256": sha256(Path(__file__).resolve()),
            "started_utc": now(),
            "inputs": inputs,
            "reference": reference,
            "status": "preflight",
            "samples": {
                row["sample_id"]: {
                    "condition": row["condition"],
                    "read1": str(INPUTS / Path(row["read1"]).name),
                    "read2": str(INPUTS / Path(row["read2"]).name),
                    "status": "pending",
                } for row in samples
            },
        }
    
    checksum = subprocess.run(["md5sum", "-c", str(INPUTS / "expected_fastq.md5")],
                              cwd=INPUTS, text=True, capture_output=True)
    (RUN / "INPUT_MD5_CHECK.log").write_text(checksum.stdout + checksum.stderr, encoding="utf-8")
    if checksum.returncode != 0 or checksum.stdout.count(": OK") != 18:
        manifest["status"] = "input_checksum_failed"
        atomic_json(manifest_path, manifest)
        raise RuntimeError("FASTQ MD5 verification failed; see INPUT_MD5_CHECK.log")
    manifest["inputs"]["verified_md5_count"] = 18
    manifest["inputs"]["verified_utc"] = now()
    manifest["status"] = "running"
    atomic_json(manifest_path, manifest)
    quant_root = RUN / "quant"
    log_root = RUN / "logs"
    quant_root.mkdir(exist_ok=True)
    log_root.mkdir(exist_ok=True)
    try:
        for sample in samples:
            sid = sample["sample_id"]
            record = manifest["samples"][sid]
            final = quant_root / sid
            working = quant_root / (sid + ".work")
            if final.exists():
                record["qc"] = validate_quant(final)
                record["status"] = "complete"
                atomic_json(manifest_path, manifest)
                write_summary(manifest)
                continue
            if working.exists():
                if working.parent != quant_root or working.name != sid + ".work":
                    raise ValueError("Unsafe temporary output path")
                shutil.rmtree(working)
            command = [str(SALMON), "quant", "-i", str(INDEX), "-l", "A",
                       "-1", record["read1"], "-2", record["read2"],
                       "-p", str(THREADS), "-o", str(working)]
            record.update(status="running", started_utc=now(), command=command)
            atomic_json(manifest_path, manifest)
            print(f"{now()} START {sid}", flush=True)
            with (log_root / f"{sid}.log").open("w", encoding="utf-8") as log:
                log.write("Command: " + " ".join(command) + "\n")
                log.flush()
                result = subprocess.run(command, stdout=log, stderr=subprocess.STDOUT)
            if result.returncode != 0:
                raise RuntimeError(f"Salmon failed for {sid} (exit {result.returncode}); see logs/{sid}.log")
            qc = validate_quant(working)
            working.rename(final)
            record.update(status="complete", completed_utc=now(), qc=qc)
            atomic_json(manifest_path, manifest)
            write_summary(manifest)
            print(f"{now()} COMPLETE {sid}: {qc['num_mapped']}/{qc['num_processed']} mapped", flush=True)
        manifest["status"] = "quantification_complete"
        manifest["completed_utc"] = now()
        atomic_json(manifest_path, manifest)
        write_summary(manifest)
        print(f"{now()} Quantification complete for nine libraries; QC_SUMMARY.tsv written", flush=True)
    except Exception as exc:
        manifest["status"] = "failed"
        manifest["failure_utc"] = now()
        manifest["error"] = str(exc)
        atomic_json(manifest_path, manifest)
        write_summary(manifest)
        raise


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("mode", choices=("preflight", "run", "status"))
    args = parser.parse_args()
    if args.mode == "preflight":
        samples, inputs = input_metadata()
        reference = reference_metadata()
        print(json.dumps({"samples": [(row["sample_id"], row["condition"]) for row in samples],
                          "fastq_count": inputs["fastq_count"],
                          "fastq_bytes": inputs["fastq_bytes"],
                          "salmon_version": reference["salmon_version"],
                          "index_seq_hash": reference["index_seq_hash"],
                          "free_gib": round(shutil.disk_usage(PACKAGE).free / 1024**3, 1)}, indent=2))
    elif args.mode == "run":
        run()
    else:
        manifest = RUN / "RUN_MANIFEST.json"
        if not manifest.exists():
            print("No run manifest yet")
            return
        data = json.loads(manifest.read_text(encoding="utf-8"))
        print(json.dumps({"status": data["status"],
                          "samples": {sid: row["status"] for sid, row in data["samples"].items()}}, indent=2))


if __name__ == "__main__":
    try:
        main()
    except Exception as error:
        print(f"ERROR: {error}", file=sys.stderr)
        sys.exit(1)
