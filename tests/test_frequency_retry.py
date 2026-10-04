import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from lx import conf_search as conf


def optfreq_log(stationary=True, energy=-40.0, imaginary=False):
    return (
        f"SCF Done: E(RB3LYP) = {energy:.12f} A.U.\n"
        "Optimization completed.\nNormal termination\n"
        "Item Value Threshold Converged?\n"
        "Maximum Force 0.00001 0.00045 YES\n"
        "RMS Force 0.00001 0.00030 YES\n"
        "Maximum Displacement 0.01 0.0018 NO\n"
        + ("Stationary point found.\n" if stationary else "")
        + f"Frequencies -- {'-20.0' if imaginary else '20.0'} 100.0 200.0\n"
        "Temperature 300.000 Kelvin.\n"
        "Thermal correction to Gibbs Free Energy= 0.01000000\n"
        f"Sum of electronic and thermal Free Energies= {energy + 0.01:.12f}\n"
        "Normal termination\n"
    )


class FrequencyRetryTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        self.name = "Geometry-1-.com"
        (self.folder / self.name).write_text("original input")
        self.log = self.folder / "Geometry-1-.log"
        self.log.write_text(optfreq_log(False))
        (self.folder / "conformer_1.chk").write_bytes(b"original checkpoint")
        self.template = dict(route="# b3lyp/gen opt=tight freq=noraman guess=mix iop(3/107=100)",
                             bottom="C H 0\n6-31G(d)\n****", nproc=8, mem="2GB", gaussian="g16")

    def test_final_frequency_check_overrides_optimization(self):
        self.assertFalse(conf.frequency_is_stationary(self.log))
        self.log.write_text(optfreq_log(True))
        self.assertTrue(conf.frequency_is_stationary(self.log))
        self.log.write_text(optfreq_log(True, imaginary=True))
        self.assertFalse(conf.frequency_is_stationary(self.log))

    def test_converged_jobs_are_not_resubmitted(self):
        self.log.write_text(optfreq_log(True))
        with patch.object(conf, "run_gaussian_jobs") as runner:
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
        runner.assert_not_called()

    def exercise_retry(self, stationary):
        original = self.log.read_text()
        def run(template, batch, files, folder, max_jobs):
            self.assertEqual((files, max_jobs), ([self.name], 3))
            retry_input = (folder / self.name).read_text()
            first, second = retry_input.split("--Link1--")
            self.assertEqual(folder, self.folder)
            self.assertNotIn("%oldchk", first)
            self.assertIn("%chk=conformer_1.chk", first)
            self.assertIn("opt=readfc guess=read geom=allcheck", first)
            self.assertNotIn("guess=mix", retry_input)
            self.assertIn("freq=noraman temperature=300", second)
            for job in (first, second):
                self.assertIn("iop(3/107=100)", job)
                self.assertIn("C H 0\n6-31G(d)\n****", job)
            self.assertEqual((folder / "Retry/Originals" / self.log.name).read_text(), original)
            (folder / self.log.name).write_text(optfreq_log(stationary, energy=-40.1))
            (folder / "conformer_1.chk").write_bytes(b"retry checkpoint")
            return []
        with patch.object(conf, "run_gaussian_jobs", side_effect=run) as runner, \
             patch.object(conf.lx.parser, "pega_geom",
                          return_value=(np.array([[0., 0., 0.], [1., 0., 0.]]), ["C", "H"])):
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
            self.assertEqual(runner.call_count, 1)
        self.assertIn("-40.100000000000", self.log.read_text())
        self.assertFalse((self.folder / "Retry").exists())
        self.assertEqual((self.folder / "conformer_1.chk").read_bytes(), b"retry checkpoint")

    def test_successful_retry_is_used_and_backups_removed(self):
        self.exercise_retry(True)

    def test_still_nonstationary_retry_is_used_without_looping(self):
        self.exercise_retry(False)

    def test_failed_retry_keeps_original(self):
        original = self.log.read_text()
        def run(template, batch, files, folder, max_jobs):
            (folder / self.log.name).write_text("Error termination\n")
            return ["Geometry-1-"]
        with patch.object(conf, "run_gaussian_jobs", side_effect=run) as runner:
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
            self.assertEqual(runner.call_count, 1)
        self.assertEqual(self.log.read_text(), original)
        self.assertEqual((self.folder / "conformer_1.chk").read_bytes(), b"original checkpoint")

    def test_interrupted_retry_resumes_same_input_and_original_hessian(self):
        original = self.log.read_text()
        calls = []
        def run(template, batch, files, folder, max_jobs):
            calls.append((folder / self.name).read_text())
            self.assertEqual((folder / "conformer_1.chk").read_bytes(), b"original checkpoint")
            if len(calls) == 1:
                self.assertEqual((folder / "Retry/Originals" / self.log.name).read_text(), original)
                (folder / self.log.name).write_text("Optimization completed.\nNormal termination\n")
                (folder / "conformer_1.chk").write_bytes(b"partial optimization checkpoint")
                return ["Geometry-1-"]
            (folder / self.log.name).write_text(optfreq_log(True))
            return []
        with patch.object(conf, "run_gaussian_jobs", side_effect=run), \
             patch.object(conf.lx.parser, "pega_geom",
                          return_value=(np.array([[0., 0., 0.], [1., 0., 0.]]), ["C", "H"])):
            for _ in range(3):
                conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
        self.assertEqual(len(calls), 2)
        self.assertEqual(calls[0], calls[1])
        self.assertFalse((self.folder / "Retry").exists())

    def test_missing_checkpoint_does_not_submit(self):
        (self.folder / "conformer_1.chk").unlink()
        with patch.object(conf, "run_gaussian_jobs") as runner:
            conf.retry_frequency_checks(self.template, "batch.sh", [self.name], self.folder, 3)
        runner.assert_not_called()


if __name__ == "__main__":
    unittest.main()
