"""List, inspect and run the paper's numerical analyses and figure components.

Only the Python standard library is needed for list/plan. Scientific dependencies
are checked in the selected runtime before any target is run.
"""
from __future__ import annotations

import argparse
import csv
from dataclasses import dataclass
import os
from pathlib import Path, PureWindowsPath
import shlex
import shutil
import subprocess
import sys
import tomllib

ROOT = Path(__file__).resolve().parent
COMMANDS = {"list", "plan", "check", "run", "prepare"}
ALIASES = {"tpm": "continuous_rna", "heatmap": "kcl_tpm", "b1": "lamin_b1",
           "genic": "genic_overlap", "lad": "lad_overlap", "limma": "regional_limma",
           "damid": "damid_continuous"}
DATASET_ALIASES = {"b1": "chipseq", "chip": "chipseq", "gene_rna": "rna", "public-lads": "public_lads"}


def split(value):
    return [v.strip() for v in value.split(";") if v.strip()]


@dataclass(frozen=True)
class Target:
    id: str
    dataset: str
    kind: str
    starting_point: str
    capability: str
    runtime: str
    script: str
    arguments: str
    cwd: str
    inputs: str
    outputs: str
    depends_on: str
    source_dependencies: str
    python_packages: str
    r_packages: str
    setup: str
    description: str
    note: str


def registry(path=ROOT / "workflows.tsv"):
    with path.open(encoding="utf-8", newline="") as handle:
        targets = [Target(**row) for row in csv.DictReader(handle, delimiter="\t")]
    by_id = {t.id: t for t in targets}
    if len(by_id) != len(targets):
        raise ValueError("Duplicate target ID in workflows.tsv")
    for target in targets:
        for dependency in split(target.depends_on):
            if dependency not in by_id:
                raise ValueError(f"{target.id}: unknown dependency {dependency}")
    ordered(targets, by_id)
    return by_id


def ordered(targets, by_id):
    result, done, active = [], set(), set()
    def visit(target):
        if target.id in active:
            raise ValueError("Dependency cycle at " + target.id)
        if target.id in done:
            return
        active.add(target.id)
        for name in split(target.depends_on):
            visit(by_id[name])
        active.remove(target.id)
        done.add(target.id)
        result.append(target)
    for target in targets:
        visit(target)
    return result


def load_settings(override=None):
    with (ROOT / "config.toml").open("rb") as handle:
        settings = tomllib.load(handle)
    if override:
        with Path(override).expanduser().open("rb") as handle:
            extra = tomllib.load(handle)
        for section, values in extra.items():
            if not isinstance(values, dict):
                raise ValueError(f"Config section {section} must be a table")
            settings.setdefault(section, {}).update(values)
    for key, value in settings["runtime"].items():
        if not isinstance(value, list) or not all(isinstance(x, str) for x in value):
            raise ValueError(f"runtime.{key} must be an array of command arguments")
    return settings


def resolve_path(value):
    path = Path(value).expanduser()
    return path.resolve() if path.is_absolute() else (ROOT / path).resolve()


def context(settings, output=None, work=None, fonts=None, data=None):
    paths = settings["paths"]
    font = fonts or paths.get("font_dir")
    if not font:
        font = "C:/Windows/Fonts" if Path("C:/Windows/Fonts").is_dir() else "/usr/share/fonts/truetype/msttcorefonts"
    values = dict(root=str(ROOT), output=str(resolve_path(output or paths["output_root"])),
                  work=str(resolve_path(work or paths["work_root"])), fonts=str(resolve_path(font)),
                  data=str(resolve_path(data or os.environ.get("GB_DATA_ROOT") or paths["data_root"])))
    for key, value in settings.get("inputs", {}).items():
        values["input_" + key] = str(resolve_path(value)) if value else ""
    return values


def render(value, values):
    try:
        return value.format_map(values)
    except KeyError as error:
        raise ValueError("Unknown configuration placeholder " + str(error)) from error


def target_context(target, values):
    return dict(values, out=str(Path(values["output"]) / target.id))


def normalize_figure(value):
    value = value.lower().replace("figure", "").replace("fig", "")
    value = value.replace("_", "").replace(" ", "").replace(".", "")
    return "Fig_S" + str(int(value[1:])).zfill(2) if value.startswith("s") else "Fig_" + str(int(value))


def figure_routes():
    with (ROOT / "figure_workflows.tsv").open(encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if "targets" not in reader.fieldnames:
            raise ValueError("figure_workflows.tsv requires a targets column for figure selection")
        return {row["figure"]: row for row in reader}


def select_targets(args, by_id):
    selected = []
    if args.target:
        for name in args.target:
            if name in {"all", "analyses", "tables", "plots", "figures", "processing"}:
                kinds = {"analysis", "plot"} if name == "all" else {"analysis"} if name in {"analyses", "tables"} else {"plot"} if name in {"plots", "figures"} else {"processing"}
                selected += [t for t in by_id.values() if t.kind in kinds and
                             t.capability == "execute"]
            else:
                name = ALIASES.get(name, name)
                if name not in by_id:
                    raise ValueError("Unknown target: " + name)
                selected.append(by_id[name])
    else:
        selected = list(by_id.values())
        if args.command == "run":
            selected = [t for t in selected if t.kind in {"analysis", "plot"} and t.capability == "execute"]
    if args.dataset:
        datasets = {DATASET_ALIASES.get(x.lower(), x.lower()) for x in args.dataset}
        selected = [t for t in selected if datasets.intersection(split(t.dataset))]
    if args.figure:
        figures = {normalize_figure(x) for x in args.figure}
        routes = figure_routes()
        unknown = figures - routes.keys()
        if unknown:
            raise ValueError("Unknown figure: " + "; ".join(sorted(unknown)))
        allowed = {name for figure in figures for name in split(routes[figure]["targets"])}
        unknown_targets = allowed - by_id.keys()
        if unknown_targets:
            raise ValueError("Unknown target in figure_map.tsv: " + "; ".join(sorted(unknown_targets)))
        selected = [t for t in selected if t.id in allowed]
        if not selected and args.command != "run" and not args.target and not args.dataset:
            return []  # A known manual figure can still be inspected.
    if not selected:
        raise ValueError("No matching targets. Use 'python run.py list' for target IDs and datasets.")
    return ordered(selected, by_id)


def runtime(target, settings):
    commands = settings["runtime"]
    command = commands.get(target.runtime, [])
    if not command and target.runtime in {"plot_python", "pdf_python"}:
        command = commands.get("python", [])
    if not command and target.runtime in {"python", "plot_python", "pdf_python"}:
        command = [sys.executable]
    command = list(command)
    if command and ("/" in command[0] or "\\" in command[0]) and not Path(command[0]).is_absolute():
        command[0] = str(resolve_path(command[0]))
    return command


def is_wsl(command):
    return bool(command) and Path(command[0]).name.lower() in {"wsl", "wsl.exe"}


def linux_path(value):
    """Translate an individual path argument, never a shell command string."""
    path = PureWindowsPath(value)
    if path.drive and len(path.drive) == 2 and path.drive[1] == ":":
        return "/mnt/" + path.drive[0].lower() + "/" + "/".join(path.parts[1:])
    return value.replace("\\", "/") if value.startswith("/") else value


def runtime_command(prefix, arguments, cwd, values):
    if not is_wsl(prefix):
        return [*prefix, *arguments], {}
    if "--" not in prefix:
        raise ValueError("A WSL runtime must include '--' before its executable")
    index = prefix.index("--")
    command = [*prefix[:index], "--cd", linux_path(str(cwd)), *prefix[index:],
               *[linux_path(str(a)) for a in arguments]]
    env = {"GB_WORK_ROOT": values["work"], "GB_FONT_DIR": values["fonts"], "GB_DATA_ROOT": values["data"]}
    forward = ["GB_WORK_ROOT/p", "GB_FONT_DIR/p", "GB_DATA_ROOT/p", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS",
               "MKL_NUM_THREADS", "MPLBACKEND", "PYTHONIOENCODING", "PYTHONDONTWRITEBYTECODE"]
    env["WSLENV"] = ":".join(filter(None, [os.environ.get("WSLENV", ""), *forward]))
    return command, env


def command_for(target, values, settings):
    local = target_context(target, values)
    script = str(ROOT / target.script) if target.script else ""
    arguments = [render(word, local) for word in shlex.split(target.arguments)]
    arguments = ["--vanilla", script, *arguments] if target.runtime == "rscript" else [script, *arguments]
    return runtime_command(runtime(target, settings), arguments,
                           Path(render(target.cwd or "{root}", local)), values)


def files_for(specification, values):
    import glob
    result = []
    for item in split(specification):
        expanded = render(item, values)
        if "{input_" in item and not expanded:
            result.append((item, False))
        elif any(char in expanded for char in "*?["):
            matches = glob.glob(expanded)
            result += [(str(p), True) for p in matches] or [(expanded, False)]
        else:
            result.append((expanded, bool(expanded) and Path(expanded).exists()))
    return result


def path_key(path):
    """Compare dependency paths consistently, including Windows separators."""
    return os.path.normcase(os.path.normpath(path))


def generated_input(path, produced):
    """Recognize declared files, their containing directories and glob outputs."""
    import fnmatch
    key = path_key(path)
    return key in produced or any(
        item.startswith(key.rstrip(os.sep) + os.sep)
        or fnmatch.fnmatchcase(key, item)
        or fnmatch.fnmatchcase(item, key)
        for item in produced)


def output_specification(target, values):
    return target.outputs


def preflight(targets, values, settings, probe=False):
    problems, probes, produced = {}, {}, set()
    for target in targets:
        local, missing = target_context(target, values), []
        if target.capability != "execute":
            missing.append(target.note or 'Source route; supply the prerequisites listed in the target plan.')
        for source in [target.script, *split(target.source_dependencies)]:
            source = source.split("::", 1)[0]
            if source and not (ROOT / source).is_file():
                missing.append("Missing source: " + source)
        for path, exists in files_for(target.inputs, local):
            if not exists and not generated_input(path, produced):
                missing.append("Missing input: " + path)
        if target.kind == "processing" and target.capability == "execute":
            for path, exists in files_for(target.outputs, local):
                if exists:
                    missing.append("Existing processing output; use a fresh work directory: " + path)
        prefix = runtime(target, settings)
        if target.capability == "execute":
            if not prefix:
                missing.append("Set runtime." + target.runtime + " in your configuration")
            elif shutil.which(prefix[0]) is None and not Path(prefix[0]).is_file():
                missing.append("Runtime executable not found: " + prefix[0])
            elif probe:
                key = (tuple(prefix), target.python_packages, target.r_packages)
                if key not in probes:
                    probes[key] = probe_runtime(target, prefix, values)
                if probes[key]:
                    missing.append(probes[key])
        problems[target.id] = missing
        if not missing:
            produced.update(path_key(render(x, local)) for x in split(output_specification(target, values)))
    return problems


def probe_runtime(target, prefix, values):
    if target.runtime in {"python", "plot_python", "pdf_python", "legacy_hmm_python"}:
        code = "import importlib; " + "; ".join("importlib.import_module(" + repr(p) + ")" for p in split(target.python_packages))
        arguments = ["-c", code]
    elif target.runtime == "rscript":
        names = ",".join('"' + p + '"' for p in split(target.r_packages))
        arguments = ["--vanilla", "-e", f"for(p in c({names})) if(!requireNamespace(p,quietly=TRUE)) stop(p)"]
    else:
        arguments = ["--version"]
    command, extra = runtime_command(prefix, arguments, ROOT, values)
    environment = dict(os.environ, **extra)
    environment.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       MPLBACKEND="Agg", PYTHONIOENCODING="utf8", PYTHONDONTWRITEBYTECODE="1")
    try:
        result = subprocess.run(command, cwd=ROOT, env=environment, text=True,
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT, timeout=45)
    except (OSError, subprocess.TimeoutExpired) as error:
        return "Runtime check failed: " + str(error)
    return "Runtime/package check failed: " + result.stdout.strip()[-1200:] if result.returncode else ""


def show(targets, values, settings, action):
    if action == "list":
        print("target\tdataset\tkind\tstarting_point\tcapability\tdescription")
        for target in targets:
            print("\t".join([target.id, target.dataset, target.kind, target.starting_point,
                             target.capability, target.description]))
        return
    for target in targets:
        local = target_context(target, values)
        print(f"{target.id}: {target.description}")
        print(f"  Dataset: {target.dataset}; start: {target.starting_point}; kind: {target.kind}; capability: {target.capability}")
        if target.depends_on:
            print("  Run after: " + target.depends_on)
        if target.source_dependencies:
            print("  Code dependencies: " + target.source_dependencies)
        if target.script:
            print("  Source: " + target.script)
        for path, _ in files_for(target.inputs, local):
            print("  Input: " + path)
        for item in split(target.outputs):
            print("  Output: " + render(item, local))
        if target.script and runtime(target, settings):
            command, _ = command_for(target, values, settings)
            print("  " + ("Command: " if target.capability == "execute" else "Source command: ") + shlex.join(command))
        print("  Working directory: " + render(target.cwd or "{root}", local))
        if target.note:
            print("  Note: " + target.note)


def protected_output(output):
    output = output.resolve()
    if output == ROOT or output in ROOT.parents:
        raise ValueError("Choose an output directory outside the repository or under results/")
    if ROOT in output.parents and output.relative_to(ROOT).parts[0] != "results":
        raise ValueError("Outputs inside the repository must be under results/")
    if output.exists():
        raise FileExistsError("Choose a new output directory: " + str(output))


def prepare(work):
    work = work.resolve()
    if work == ROOT or work in ROOT.parents or ROOT in work.parents:
        raise ValueError("Choose a working directory outside the code repository")
    with (ROOT / "common/staging.tsv").open(newline="", encoding="utf-8") as handle:
        copies = [(ROOT / r["source"], work / r["destination"]) for r in csv.DictReader(handle, delimiter="\t")]
    for source, destination in copies:
        if destination.exists() and destination.read_bytes() != source.read_bytes():
            raise FileExistsError("Existing modified working file: " + str(destination))
    for source, destination in copies:
        destination.parent.mkdir(parents=True, exist_ok=True)
        if not destination.exists():
            shutil.copyfile(source, destination)
    for folder in ["wgbs/Reports", "wgbs_package/runs", "rna_package/runs", "chipseq/runs", "references", "figures", "rna_clusters", "rna/SalmonTE_DE_1h", "rna/SalmonTE_DE_6h", "rna/SalmonTE_DE_24h"]:
        (work / folder).mkdir(parents=True, exist_ok=True)
    print("Prepared code and configuration in " + str(work))


def setup_target(target, values):
    local = target_context(target, values)
    out = Path(local["out"])
    if target.setup == "limma":
        (out / "qc").mkdir(parents=True)
        (out / "tables").mkdir()
        for item in split(target.inputs):
            source = Path(render(item, local))
            shutil.copyfile(source, out / "tables" / source.name)
    elif target.setup == "public_lad_summary":
        out.mkdir(parents=True)
        for path in split(target.inputs):
            source = Path(render(path, local))
            shutil.copyfile(source, out / source.name)
    elif target.setup == "rna_cluster_GO":
        source = Path(values["work"]) / "rna/K_clustering_Neuron_RNA_seq_11-13-24/filtered_results/Clustered_Genes_Hierarchical_Euclidean_Z.csv"
        destination = Path(values["work"]) / "rna_clusters" / source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        if destination.exists() and destination.read_bytes() != source.read_bytes():
            raise FileExistsError("Existing different GO input: " + str(destination))
        if not destination.exists():
            shutil.copyfile(source, destination)
    elif target.setup == "output_dir":
        out.mkdir(parents=True)


def run_targets(targets, values, settings):
    output = Path(values["output"])
    protected_output(output)
    problems = preflight(targets, values, settings, probe=True)
    if any(problems.values()):
        print_checks(targets, problems)
        raise ValueError("Prerequisite checks failed; no target was run and no output directory was created.")
    output.mkdir(parents=True)
    (output / "logs").mkdir()
    environment = os.environ.copy()
    environment.update(OPENBLAS_NUM_THREADS="1", OMP_NUM_THREADS="1", MKL_NUM_THREADS="1",
                       MPLBACKEND="Agg", PYTHONIOENCODING="utf8", PYTHONDONTWRITEBYTECODE="1",
                       GB_WORK_ROOT=values["work"], GB_FONT_DIR=values["fonts"], GB_DATA_ROOT=values["data"])
    rows = []
    try:
        for target in targets:
            local = target_context(target, values)
            log = output / "logs" / (target.id + ".log")
            row = dict(target=target.id, status="FAIL", returncode="", output=local["out"], log=str(log))
            rows.append(row)
            setup_target(target, values)
            cwd = Path(render(target.cwd or "{root}", local))
            cwd.mkdir(parents=True, exist_ok=True)
            command, extra = command_for(target, values, settings)
            env = dict(environment, **extra)
            if target.id == "lad_overlap":
                env["OUTPUT_DIR"] = local["out"]
                if is_wsl(runtime(target, settings)):
                    env["WSLENV"] += ":OUTPUT_DIR/p"
            with log.open("w", encoding="utf-8") as handle:
                handle.write("Command: " + shlex.join(command) + "\nWorking directory: " + str(cwd) + "\n\n")
                handle.flush()
                result = subprocess.run(command, cwd=cwd, env=env, stdout=handle, stderr=subprocess.STDOUT)
            row["returncode"] = result.returncode
            if result.returncode:
                raise RuntimeError(f"{target.id} failed; see {log}")
            missing = [path for path, exists in files_for(output_specification(target, values), local) if not exists]
            if missing:
                raise RuntimeError(target.id + " did not create expected outputs: " + "; ".join(missing))
            row["status"] = "PASS"
            print(target.id + " PASS", flush=True)
    finally:
        with (output / "run_summary.tsv").open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=["target", "status", "returncode", "output", "log"], delimiter="\t")
            writer.writeheader()
            writer.writerows(rows)


def print_checks(targets, problems):
    for target in targets:
        issues = problems[target.id]
        label = "SOURCE" if target.capability != "execute" else "MISSING" if issues else "READY"
        print(target.id + " " + label)
        for issue in issues:
            print("  " + issue)


def parser():
    p = argparse.ArgumentParser(description=__doc__, epilog="Examples: python run.py list; python run.py run --target kcl_tpm --output results/example")
    sub = p.add_subparsers(dest="command", required=True)
    for name in ["list", "plan", "check", "run", "prepare"]:
        child = sub.add_parser(name)
        child.add_argument("--config", type=Path, help="Explicit local TOML override; relative paths resolve from the repository.")
        child.add_argument("--output", type=Path)
        child.add_argument("--work-root", type=Path)
        child.add_argument("--data-root", type=Path, help="Directory containing external input files.")
        child.add_argument("--fonts", type=Path)
        child.add_argument("--rscript", help="Override the Rscript executable.")
        child.add_argument("--r-via-wsl", action="store_true", help="Use wsl.exe -d Ubuntu -- Rscript.")
        if name != "prepare":
            child.add_argument("--target", action="append", help="Target ID or analyses/plots/all/processing; repeat to combine.")
            child.add_argument("--dataset", action="append", help="wgbs/rna/te/chipseq/public_lads/damid/public_chromatin; repeat to combine.")
            child.add_argument("--figure", action="append", help="Figure ID, for example Fig_5 or Fig_S06.")
        if name == "check":
            child.add_argument("--static", action="store_true", help="Check files and runtime paths without importing packages.")
    return p


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    if not argv:
        argv = ["list"]
    elif argv[0] not in COMMANDS and argv[0] not in {"-h", "--help"}:
        argv.insert(0, "run")
    args = parser().parse_args(argv)
    try:
        settings = load_settings(args.config)
        if args.rscript:
            settings["runtime"]["rscript"] = [args.rscript]
        if args.r_via_wsl:
            settings["runtime"]["rscript"] = ["wsl.exe", "-d", "Ubuntu", "--", "Rscript"]
        values = context(settings, args.output, args.work_root, args.fonts, args.data_root)
        values["selected_figures"] = [normalize_figure(x) for x in (getattr(args, "figure", None) or [])]
        if args.command == "prepare":
            prepare(Path(values["work"]))
            return 0
        targets = select_targets(args, registry())
        if args.figure and args.command in {"list", "plan", "check"}:
            routes = figure_routes()
            for name in args.figure:
                row = routes[normalize_figure(name)]
                print(row["figure"] + ": " + row["title"])
                print("  Analysis code: " + row["analysis_code"])
                print("  Figure code: " + row["plot_code"])
        if args.command in {"list", "plan"}:
            show(targets, values, settings, args.command)
        elif args.command == "check":
            problems = preflight(targets, values, settings, probe=not args.static)
            print_checks(targets, problems)
            return 2 if any(problems.values()) else 0
        else:
            run_targets(targets, values, settings)
        return 0
    except (ValueError, OSError, RuntimeError) as error:
        print("Error: " + str(error), file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
