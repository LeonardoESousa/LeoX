import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lx import conf_search as conf


def completed_log(number):
    energy = -40.0 - number * 0.001
    return (
        "Standard orientation:\n -----\n Center Atomic Atomic Coordinates\n"
        " Number Number Type X Y Z\n -----\n"
        f"1 6 0 0.0 0.0 0.0\n2 1 0 {1.0 + number * 0.05} 0.0 0.0\n -----\n"
        f"SCF Done: E(RB3LYP) = {energy:.12f} A.U.\n"
        "Optimization completed.\nStationary point found.\nNormal termination\n"
        "Item Value Threshold Converged?\nStationary point found.\n"
        "Frequencies -- 20.0 100.0 200.0\nTemperature 300.000 Kelvin.\n"
        "Thermal correction to Gibbs Free Energy= 0.01000000\n"
        f"Sum of electronic and thermal Free Energies= {energy + 0.01:.12f}\n"
        "Normal termination\n"
    )


class ResumeTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)
        self.input = self.root / "start.com"
        self.input.write_text("%nprocshared=2\n%mem=1GB\n# b3lyp/6-31G opt freq=noraman\n\n"
                              "Title\n\n0 1\nC 0 0 0\nH 1 0 0\n\n")
        self.batch = self.root / "batch.sh"
        self.batch.write_text("#!/bin/bash\n#SBATCH --time=1-0\nbash \"$1\"\n")
        self.folder = self.root / "Conformational"
        self.calls = []
        for attribute, callback in (("run_command", self.command), ("run_slurm", self.crest),
                                    ("run_gaussian_jobs", self.gaussian)):
            mocker = patch.object(conf, attribute, side_effect=callback)
            mocker.start()
            self.addCleanup(mocker.stop)

    def command(self, command, folder, logfile):
        folder = Path(folder)
        if command[0] == "xtb":
            self.calls.append("xtb")
            (folder / logfile).write_text("GEOMETRY OPTIMIZATION CONVERGED\nnormal termination of xtb\n")
            conf.write_xyz(folder / "xtbopt.xyz", conf.read_xyz(folder / "start.xyz"))
        else:
            self.calls.append("classification")
            (folder / logfile).write_text("CREST terminated normally.\n")
            conf.write_xyz(folder / "reoptimized.xyz.sorted", conf.read_xyz(folder / "reoptimized.xyz"))

    def crest(self, batch, script, folder, nproc, mem):
        self.calls.append("crest")
        folder = Path(folder)
        (folder / "crest.out").write_text("CREST terminated normally.\n")
        structure = conf.read_xyz(folder / "start.xyz")[0]
        conf.write_xyz(folder / "crest_clustered.xyz", [structure, structure])

    def gaussian(self, template, batch, files, folder, max_jobs):
        self.calls.append(("gaussian", tuple(files)))
        for name in files:
            number = int(Path(name).stem.split("-")[1])
            (Path(folder) / Path(name).with_suffix(".log")).write_text(completed_log(number))
            (Path(folder) / f"conformer_{number}.chk").write_bytes(b"checkpoint")
        return []

    def run_search(self, **kwargs):
        return conf.run_workflow(self.input, self.batch, self.batch, workdir=self.folder, **kwargs)

    def test_completed_search_does_no_external_work_on_resume(self):
        self.assertEqual(len(self.run_search()), 2)
        self.calls.clear()
        self.assertEqual(len(self.run_search()), 2)
        self.assertEqual(self.calls, [])

    def test_only_missing_gaussian_job_is_resubmitted(self):
        self.run_search()
        (self.folder / "Geometries/Geometry-2-.log").unlink()
        self.calls.clear()
        self.run_search()
        self.assertEqual(self.calls[0], ("gaussian", ("Geometry-2-.com",)))
        self.assertNotIn("xtb", self.calls)
        self.assertNotIn("crest", self.calls)

    def test_missing_or_modified_report_repeats_classification(self):
        self.run_search()
        (self.folder / "conformation.csv").write_text("interrupted report")
        self.calls.clear()
        self.run_search()
        self.assertEqual(self.calls, ["classification"])

    def test_incomplete_xtb_restarts_in_clean_folder(self):
        self.run_search()
        (self.folder / "xTB/xtb.log").write_text("interrupted")
        (self.folder / "xTB/partial_restart").write_text("stale")
        self.calls.clear()
        self.run_search()
        self.assertEqual(self.calls, ["xtb"])
        self.assertFalse((self.folder / "xTB/partial_restart").exists())
        self.assertFalse((self.folder / "Inputs/Interrupted").exists())

    def test_incomplete_crest_restarts_and_replaces_downstream_inputs(self):
        self.run_search()
        (self.folder / "CREST/crest.out").write_text("interrupted")
        (self.folder / "CREST/partial_restart").write_text("stale")
        self.calls.clear()
        self.run_search()
        self.assertIn("crest", self.calls)
        self.assertIn(("gaussian", ("Geometry-1-.com", "Geometry-2-.com")), self.calls)
        self.assertNotIn("xtb", self.calls)
        self.assertFalse((self.folder / "CREST/partial_restart").exists())
        self.assertFalse((self.folder / "Inputs/Interrupted").exists())

    def test_old_interrupted_archives_are_deleted_on_resume(self):
        self.run_search()
        for name in ("Inputs/Interrupted", "Geometries/Interrupted"):
            directory = self.folder / name
            directory.mkdir(parents=True)
            (directory / "partial.log").write_text("interrupted")
        self.calls.clear()
        self.run_search()
        self.assertEqual(self.calls, [])
        self.assertFalse((self.folder / "Inputs/Interrupted").exists())
        self.assertFalse((self.folder / "Geometries/Interrupted").exists())

    def test_changed_input_is_rejected_before_external_work(self):
        self.run_search()
        self.input.write_text(self.input.read_text().replace("b3lyp", "pbe0"))
        self.calls.clear()
        with self.assertRaisesRegex(ValueError, "different Gaussian input"):
            self.run_search()
        self.assertEqual(self.calls, [])

    def test_changed_threshold_repeats_only_classification(self):
        self.run_search()
        self.calls.clear()
        self.run_search(rthr=0.2)
        self.assertEqual(self.calls, ["classification"])

    def test_old_inputs_are_not_overwritten_on_resume(self):
        self.run_search()
        job = self.folder / "Geometries/Geometry-1-.com"
        original = job.read_text()
        job.write_text(original + "\n! retry input sentinel\n")
        self.calls.clear()
        self.run_search()
        self.assertEqual(self.calls, [])
        self.assertTrue(job.read_text().endswith("! retry input sentinel\n"))


class ResubmissionWatcherTests(unittest.TestCase):
    def test_watcher_always_preserves_frequency_checkpoints(self):
        import os
        with tempfile.TemporaryDirectory() as temporary:
            previous = Path.cwd()
            try:
                os.chdir(temporary)
                Path("Geometry-1-.log").write_text("Normal termination\nNormal termination\n")
                with patch.object(conf.lx.tools, "delchk") as cleanup:
                    watcher = conf.lx.tools.Watcher(".", files=["Geometry-1-.com"], counter=2)
                    watcher.check()
                self.assertEqual(watcher.done, ["Geometry-1-"])
                cleanup.assert_not_called()
            finally:
                os.chdir(previous)

    def test_stale_log_and_status_are_removed_before_watching(self):
        with tempfile.TemporaryDirectory() as temporary:
            folder = Path(temporary)
            (folder / "Geometry-1-.log").write_text("Normal termination\nNormal termination\n")
            (folder / "cmd_0_.sh.status").write_text("0")
            seen = []
            class Watcher:
                error = []
                def __init__(self, files):
                    seen.extend(files)
                    assert not (folder / "Geometry-1-.log").exists()
                    assert not (folder / "cmd_0_.sh.status").exists()
                def run(self, *args):
                    pass
                def hold_watch(self):
                    pass
            with patch.object(conf, "ConformerWatcher", Watcher):
                conf.run_gaussian_jobs(dict(nproc=2, mem="1GB", gaussian="g16"),
                                       "batch.sh", ["Geometry-1-.com"], folder, 1)
            self.assertEqual(seen, ["Geometry-1-.com"])
            self.assertFalse((folder / "Interrupted").exists())


if __name__ == "__main__":
    unittest.main()
