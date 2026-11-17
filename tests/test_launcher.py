"""Test command selection, dependency ordering and output protection."""
import argparse
from dataclasses import replace
import importlib.util
import io
from pathlib import Path
import sys
import tempfile
import unittest
from contextlib import redirect_stdout
from unittest.mock import patch

ROOT = Path(__file__).resolve().parents[1]
spec = importlib.util.spec_from_file_location("paper_launcher", ROOT / "run.py")
launcher = importlib.util.module_from_spec(spec)
sys.modules[spec.name] = launcher
spec.loader.exec_module(launcher)


class LauncherTests(unittest.TestCase):
    def setUp(self):
        self.targets = launcher.registry()
        self.settings = launcher.load_settings()
        self.settings["runtime"]["python"] = [sys.executable]

    def args(self, **values):
        defaults = dict(command="run", target=None, dataset=None, figure=None)
        defaults.update(values)
        return argparse.Namespace(**defaults)

    def stub(self, folder, exit_code=0):
        source = folder / "stub.py"
        source.write_text("from pathlib import Path\nimport sys\np=Path(sys.argv[1]);p.mkdir(parents=True);(p/'result.tsv').write_text('value\\n7\\n');raise SystemExit(" + str(exit_code) + ")\n")
        return replace(self.targets["damid_continuous"], id="stub", script=str(source),
                       arguments="{out}", inputs="", outputs="{out}/result.tsv",
                       python_packages="", r_packages="")

    def test_external_inputs_and_source_dependencies(self):
        selected = launcher.select_targets(self.args(target=["all"]), self.targets)
        self.assertNotIn("reference_figures", self.targets)
        values = launcher.context(self.settings, data="../example_inputs")
        self.assertNotEqual(Path(values["data"]), ROOT / "data")
        for target in self.targets.values():
            for source in [target.script, *launcher.split(target.source_dependencies)]:
                if source:
                    self.assertTrue((ROOT / source.split("::")[0]).is_file(), source)
        self.assertIn("rna/expression/kcl_tpm.py::calculate", self.targets["kcl_heatmaps"].source_dependencies)
        self.assertEqual(self.targets["promoter_methylation"].depends_on, "continuous_panels")

    def test_dependency_order_cycle_and_authoritative_figure_map(self):
        selected = launcher.ordered([self.targets["rna_cluster_GO"]], self.targets)
        self.assertEqual([t.id for t in selected], ["rna_trajectories", "rna_clusters", "rna_cluster_GO"])
        cyclic = dict(self.targets)
        cyclic["rna_trajectories"] = replace(cyclic["rna_trajectories"], depends_on="rna_cluster_GO")
        with self.assertRaisesRegex(ValueError, "cycle"):
            launcher.ordered([cyclic["rna_cluster_GO"]], cyclic)
        with patch.object(launcher, "figure_routes", return_value={"Fig_2": {"targets": "kcl_heatmaps"}, "Fig_8": {"targets": ""}}):
            selected = launcher.select_targets(self.args(figure=["Figure 2"]), self.targets)
            self.assertEqual(selected[-1].id, "kcl_heatmaps")
            self.assertIn("kcl_tpm_inputs", [t.id for t in selected])
            self.assertEqual(launcher.select_targets(self.args(command="plan", figure=["8"]), self.targets), [])
            with self.assertRaisesRegex(ValueError, "No matching"):
                launcher.select_targets(self.args(figure=["8"]), self.targets)

    def test_generated_inputs_have_upstream_producers(self):
        for target in self.targets.values():
            if target.capability != "execute":
                continue
            upstream = launcher.ordered([target], self.targets)[:-1]
            produced = [value.replace("{out}", "{output}/" + parent.id)
                        for parent in upstream for value in launcher.split(parent.outputs)]
            for value in launcher.split(target.inputs):
                if value.startswith("{output}/"):
                    self.assertTrue(launcher.generated_input(value, {launcher.path_key(p) for p in produced}),
                                    target.id + ": no upstream producer for " + value)

    def test_wsl_translation_keeps_spaced_arguments_and_environment_paths(self):
        command, env = launcher.runtime_command(
            ["wsl.exe", "-d", "Ubuntu", "--", "Rscript"],
            ["--vanilla", r"C:\a b\model.R", r"D:\new run"],
            r"C:\paper code", {"work": r"D:\work root", "fonts": r"C:\Windows\Fonts", "data": r"D:\inputs"})
        self.assertEqual(command[-3:], ["--vanilla", "/mnt/c/a b/model.R", "/mnt/d/new run"])
        self.assertIn("/mnt/c/paper code", command)
        self.assertEqual(env["GB_WORK_ROOT"], r"D:\work root")
        self.assertIn("GB_WORK_ROOT/p", env["WSLENV"])
        self.assertIn("OMP_NUM_THREADS", env["WSLENV"])

    def test_wsl_probe_populates_every_forwarded_runtime_variable(self):
        values = launcher.context(self.settings)
        with patch.object(launcher.subprocess, "run") as execute:
            execute.return_value.returncode = 0
            self.assertEqual(launcher.probe_runtime(self.targets["te_families"],
                             ["wsl.exe", "-d", "Ubuntu", "--", "python3"], values), "")
            environment = execute.call_args.kwargs["env"]
            for name in ["OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"]:
                self.assertEqual(environment[name], "1")
            for name in environment["WSLENV"].split(":"):
                key = name.split("/", 1)[0]
                if key in {"GB_WORK_ROOT", "GB_FONT_DIR", "OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS", "MPLBACKEND", "PYTHONIOENCODING", "PYTHONDONTWRITEBYTECODE"}:
                    self.assertTrue(environment[key])

    def test_missing_prerequisite_blocks_all_commands_before_output_creation(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            target = replace(self.stub(folder), inputs=str(folder / "missing.tsv"))
            values = launcher.context(self.settings, output=folder / "new output")
            with redirect_stdout(io.StringIO()), self.assertRaisesRegex(ValueError, "Prerequisite"):
                launcher.run_targets([target], values, self.settings)
            self.assertFalse(Path(values["output"]).exists())
            self.assertFalse((folder / "new output/stub/result.tsv").exists())

    def test_success_failure_logs_and_existing_outputs_are_preserved(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            target = self.stub(folder)
            values = launcher.context(self.settings, output=folder / "success")
            with redirect_stdout(io.StringIO()):
                launcher.run_targets([target], values, self.settings)
            self.assertIn("PASS", (folder / "success/run_summary.tsv").read_text())
            self.assertIn("Command:", (folder / "success/logs/stub.log").read_text())
            with self.assertRaises(FileExistsError):
                launcher.run_targets([target], values, self.settings)
            target = self.stub(folder, exit_code=3)
            values = launcher.context(self.settings, output=folder / "failure")
            with self.assertRaisesRegex(RuntimeError, "stub failed"):
                launcher.run_targets([target], values, self.settings)
            self.assertIn("FAIL\t3", (folder / "failure/run_summary.tsv").read_text())
            self.assertTrue((folder / "failure/stub/result.tsv").is_file())

    def test_go_input_handoff_preserves_modified_destination(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            values = launcher.context(self.settings, output=folder / "out", work=folder / "work")
            source = folder / "work/rna/K_clustering_Neuron_RNA_seq_11-13-24/filtered_results/Clustered_Genes_Hierarchical_Euclidean_Z.csv"
            source.parent.mkdir(parents=True)
            source.write_bytes(b"gene_id,Cluster\ngene1,2\n")
            launcher.setup_target(self.targets["rna_cluster_GO"], values)
            destination = folder / "work/rna_clusters" / source.name
            self.assertEqual(source.read_bytes(), destination.read_bytes())
            destination.write_bytes(b"author changed\n")
            with self.assertRaises(FileExistsError):
                launcher.setup_target(self.targets["rna_cluster_GO"], values)
            self.assertEqual(destination.read_bytes(), b"author changed\n")

    def test_current_figure_numbers_select_current_generators(self):
        cases = {"Fig_1": {"b1_main1", "b1_signal_spread", "b1_chr19_tracks"},
                 "Fig_3": {"b1_main3", "matched_b1_plot", "b1_locus_tracks"},
                 "Fig_5": {"kcl_heatmaps"}}
        for figure, expected in cases.items():
            selected = launcher.select_targets(self.args(target=["plots"], figure=[figure]), self.targets)
            self.assertTrue(expected <= {target.id for target in selected})

    def test_generated_table_dependencies_are_recognized_before_execution(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            first = self.stub(folder)
            second = replace(first, id="second", depends_on="stub",
                             inputs="{output}/stub/result.tsv", outputs="{out}/result.tsv")
            values = launcher.context(self.settings, output=folder / "new run")
            issues = launcher.preflight([first, second], values, self.settings)
            self.assertFalse(any(issues.values()))

    def test_generated_directories_and_patterns_require_declared_outputs(self):
        declared = {launcher.path_key('/results/models/A.model.json'),
                    launcher.path_key('/results/counts/*.tsv')}
        self.assertTrue(launcher.generated_input('/results/models', declared))
        self.assertTrue(launcher.generated_input('/results/models/*.model.json', declared))
        self.assertTrue(launcher.generated_input('/results/counts/sample.tsv', declared))
        self.assertFalse(launcher.generated_input('/results/models/missing.tsv', declared))
        self.assertFalse(launcher.generated_input('/results/unrelated', declared))

    def test_repository_has_no_scientific_data_or_artwork(self):
        extensions = {".bed", ".npz", ".npy", ".bw", ".bam", ".gz", ".csv", ".pdf", ".png", ".ai", ".pzfx"}
        self.assertFalse([str(p.relative_to(ROOT)) for p in ROOT.rglob("*")
                          if p.is_file() and ".git" not in p.parts and p.suffix.lower() in extensions])


if __name__ == "__main__":
    unittest.main()
