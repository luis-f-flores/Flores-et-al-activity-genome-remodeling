"""Generate genomic matches and calculate IAP hypo-DMR enrichment."""
from pathlib import Path
import argparse
import json
import subprocess
import sys


def run(settings, output):
    settings = Path(settings).resolve()
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = json.loads(settings.read_text())
    input_keys = ["context", "rna", "dmrs", "tested_windows", "genes", "exons",
                  "promoters", "kcl_bw", "mann_bw", "bismap", "raw_reports"]
    for key in input_keys:
        path = Path(config[key])
        config[key] = str(path if path.is_absolute() else (settings.parent / path).resolve())
    for key in ["reuse_geometry", "reuse_diagnostics", "reuse_final_sizes"]:
        config.pop(key, None)
    initial = output / "candidate_pools"
    initial_settings = output / "candidate_pool_settings.json"
    initial_settings.write_text(json.dumps(config, indent=2) + "\n")
    script = Path(__file__).with_name("iap_matching.py")
    subprocess.run([sys.executable, str(script), "--settings", str(initial_settings),
                    "--output", str(initial), "--prepare-pools"], check=True)
    config.update(reuse_geometry=str(initial / "candidate_geometry_before_depth.tsv.gz"),
                  reuse_diagnostics=str(initial / "pool_search_diagnostics.tsv"),
                  reuse_final_sizes=str(initial / "final_pool_sizes.tsv"))
    completed_settings = output / "matching_settings.json"
    completed_settings.write_text(json.dumps(config, indent=2) + "\n")
    subprocess.run([sys.executable, str(script), "--settings", str(completed_settings),
                    "--output", str(output / "results")], check=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--settings", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.settings, args.output)
