import csv
import os
from pathlib import Path
import tempfile
import unittest
from unittest.mock import patch

import numpy as np
from lx import conf_search as conf
from lx import parser


def gaussian_log(geometry, energy=-40.0, correction=0.01, imaginary=False, failed=False):
    rows = ''.join(f'{i} {z} 0 {xyz[0]:.9f} {xyz[1]:.9f} {xyz[2]:.9f}\n'
                   for i, (z, xyz) in enumerate(zip([6, 1], geometry), 1))
    orientation = ' Standard orientation:\n -----\n Center Atomic Atomic Coordinates\n Number Number Type X Y Z\n -----\n' + rows + ' -----\n'
    text = ' Entering Gaussian System\n' + orientation
    text += f' SCF Done: E(RB3LYP) = {energy:.12f} A.U.\n Optimization completed.\n Normal termination\n'
    if failed:
        return text + ' Error termination\n'
    text += ' Frequencies -- ' + ('-20.0' if imaginary else '20.0') + ' 100.0 200.0\n'
    text += ' Temperature 300.000 Kelvin.\n'
    text += f' Thermal correction to Gibbs Free Energy= {correction:.8f}\n'
    text += f' Sum of electronic and thermal Free Energies= {energy+correction:.8f}\n Normal termination\n'
    return text


class WorkflowTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.root = Path(self.temp.name)
        self.geometry = np.array([[0., 0., 0.], [1., 0., 0.]])
        self.input = self.root / 'start.com'
        self.input.write_text('%nprocshared=8\n%mem=2GB\n%chk=original.chk\n'
                              '#p wb97xd/gen opt=(tight, maxcycles=100) freq scrf=(pcm, read) guess=mix\n'
                              ' integral=ultrafine scf=(xqc, maxcycle=500)\n\nTitle\n\n-1 2\n'
                              '6 0.0 0.0 0.0\n1 1.0 0.0 0.0\n\nC H 0\n6-31G(d)\n****\n\nEps=3.0\n')
        self.batch = self.root / 'batch.sh'
        self.batch.write_text('#!/bin/bash\n#SBATCH --time=1-0\nbash "$1"\n')

    def tearDown(self):
        self.temp.cleanup()

    def test_input_settings_and_two_job_tail(self):
        data = parser.read_gaussian_input(self.input)
        self.assertEqual((data['nproc'], data['mem'], data['charge'], data['multiplicity']), (8, '2GB', -1, 2))
        self.assertEqual(data['atoms'], ['C', 'H'])
        folder = self.root / 'inputs'
        folder.mkdir()
        files = conf.make_gaussian_inputs(data, [dict(atoms=data['atoms'], geometry=data['geometry'])]*2, folder)
        generated = (folder / files[0]).read_text()
        opt, freq = generated.split('--Link1--')
        self.assertIn('opt=(tight, maxcycles=100)', opt)
        self.assertNotIn(' freq ', opt)
        self.assertIn('integral=ultrafine scf=(xqc, maxcycle=500)', opt)
        self.assertIn('geom=allcheck guess=read', freq)
        self.assertNotIn('guess=mix', freq)
        for step in (opt, freq):
            self.assertIn('C H 0\n6-31G(d)\n****\n\nEps=3.0', step)
            self.assertIn('%chk=conformer_1.chk', step)
        self.assertIn('%chk=conformer_2.chk', (folder / files[1]).read_text())
        self.assertEqual(parser.pega_geom(str(self.input))[1], ['C', 'H'])

    def test_reject_unsupported_inputs_before_search(self):
        original = self.input.read_text()
        for bad in (original + '\n--Link1--\n', original.replace('6 0.0 0.0 0.0', 'C'),
                    original.replace('guess=mix', 'guess=read')):
            self.input.write_text(bad)
            with self.assertRaises(ValueError):
                parser.read_gaussian_input(self.input)

    def test_gaussian_minima_and_numerical_frequency_geometry(self):
        log = self.root / 'Geometry-1-.log'
        text = gaussian_log(self.geometry)
        # A numerical Hessian's final displaced coordinates/SCF energy must not be used.
        text = text.replace(' Frequencies --', ' Standard orientation:\n -----\n Number Atomic Type X Y Z\n -----\n1 6 0 9 9 9\n2 1 0 8 8 8\n -----\n SCF Done: E(RB3LYP) = -39.0 A.U.\n Frequencies --')
        log.write_text(text)
        result = conf.gaussian_result(log)
        np.testing.assert_allclose(result['geometry'], self.geometry)
        self.assertEqual(result['atoms'], ['C', 'H'])
        self.assertEqual(result['energy'], -40.0)
        for invalid in (gaussian_log(self.geometry, imaginary=True), gaussian_log(self.geometry, failed=True),
                        gaussian_log(self.geometry).replace('300.000 Kelvin', '298.150 Kelvin')):
            log.write_text(invalid)
            with self.assertRaises(ValueError):
                conf.gaussian_result(log)

    def test_correlated_energy_from_thermochemistry(self):
        log = self.root / 'Geometry-1-.log'
        log.write_text(gaussian_log(self.geometry).replace('SCF Done: E(RB3LYP) = -40.000000000000',
                                                          'SCF Done: E(RHF) = -39.000000000000'))
        self.assertAlmostEqual(conf.gaussian_result(log)['energy'], -40.0)

    def test_populations_and_memory(self):
        self.assertEqual([parser.element_symbol(str(z)) for z in range(1, 119)], parser.ELEMENTS[1:])
        _, population = conf.populations([-40., -40. + parser.BOLTZ_EV * 300 / conf.HARTREE_EV])
        self.assertAlmostEqual(population[0]/population[1], np.e)
        self.assertAlmostEqual(sum(population), 100)
        for mem, expected in [('2GB', '2G'), ('100MW', '800M'), ('131072', '1M')]:
            command = conf.sbatch_command(self.batch, self.input, self.root, 8, mem)
            self.assertIn('--mem='+expected, command)
            self.assertIn('--cpus-per-task=8', command)

    def test_frequency_options_and_connectivity_preserved_correctly(self):
        data = parser.read_gaussian_input(self.input)
        data['route'] = '# b3lyp/gen opt=tight freq=(numer, noraman) geom=connectivity'
        data['bottom'] = '1 2 1.0\n2\n\nC H 0\n6-31G(d)\n****'
        opt, freq, tail = conf.gaussian_routes(data)
        self.assertIn('geom=connectivity', opt)
        self.assertIn('freq=(numer, noraman)', freq)
        self.assertNotIn('geom=connectivity', freq)
        self.assertEqual(tail, 'C H 0\n6-31G(d)\n****')

    def test_slurm_failure_is_propagated(self):
        with patch.object(conf.subprocess, 'Popen') as popen:
            popen.return_value.poll.return_value = 9
            popen.return_value.returncode = 9
            with self.assertRaisesRegex(RuntimeError, 'SLURM job failed'):
                conf.run_slurm(self.batch, self.input, self.root, 8, '2GB')

    def test_missing_geometry_is_rejected(self):
        log = self.root / 'Geometry-1-.log'
        log.write_text(gaussian_log(self.geometry).replace(' Standard orientation:', ' Unrelated heading:'))
        with self.assertRaisesRegex(ValueError, 'Missing final Gaussian geometry'):
            conf.gaussian_result(log)

    def test_missing_output_is_reported(self):
        (self.root / 'Geometry-1-.com').write_text('input')
        with self.assertRaisesRegex(ValueError, 'No completed, verified minima'):
            conf.classify_only(self.root)
        self.assertIn('Geometry-1-.log', (self.root / 'rejected_conformers.csv').read_text())

    def test_watcher_finishes_killed_jobs_and_waits_for_both_links(self):
        previous = Path.cwd()
        os.chdir(self.root)
        try:
            Path('Geometry-1-.log').write_text('Normal termination\n')
            watcher = conf.ConformerWatcher(['Geometry-1-.com'])
            watcher.check()
            self.assertEqual(watcher.done, [])
            Path('cmd_0_.sh.status').write_text('137')
            watcher.check()
            self.assertEqual(watcher.error, ['Geometry-1-'])
            Path('cmd_0_.sh.status').unlink()
            Path('Geometry-1-.log').write_text('Normal termination\nNormal termination\n')
            watcher = conf.ConformerWatcher(['Geometry-1-.com'])
            watcher.check()
            self.assertEqual(watcher.done, ['Geometry-1-'])
        finally:
            os.chdir(previous)

    def test_full_workflow_with_simulated_external_programs(self):
        commands = []
        def fake_command(command, folder, logfile):
            commands.append(command)
            folder = Path(folder)
            (folder / logfile).write_text('completed')
            if command[0] == 'xtb':
                conf.write_xyz(folder / 'xtbopt.xyz', conf.read_xyz(folder / 'start.xyz'))
            else:
                ensemble = conf.read_xyz(folder / 'reoptimized.xyz')
                # Collapse two opt jobs to the same minimum, retaining the lower-energy copy.
                conf.write_xyz(folder / 'reoptimized.xyz.sorted', [ensemble[1], ensemble[2]])
        def fake_crest(batch, script, folder, nproc, mem):
            self.assertEqual(nproc, 8)
            self.assertIn('--cluster', Path(script).read_text())
            start = conf.read_xyz(Path(folder) / 'start.xyz')[0]
            conf.write_xyz(Path(folder) / 'crest_clustered.xyz', [start]*4)
        def fake_submit(command):
            wrapper = Path(command[1]).read_text()
            self.assertIn('sbatch --wait', wrapper)
            self.assertIn('--cpus-per-task=8', wrapper)
            self.assertIn('--mem=2G', wrapper)
            number = int(Path(command[2]).name.split('_')[1]) + 1
            geometry = self.geometry.copy()
            if number == 3:
                geometry[1, 0] = 1.2
            Path(f'Geometry-{number}-.log').write_text(gaussian_log(geometry, -40.0 + (0.001 if number == 1 else 0),
                                                                 correction=0.01+(0.002 if number == 3 else 0),
                                                                 failed=number==4))
            return 0
        with patch.object(conf, 'run_command', side_effect=fake_command), \
             patch.object(conf, 'run_slurm', side_effect=fake_crest), \
             patch.object(conf.subprocess, 'call', side_effect=fake_submit), \
             patch.object(conf.time, 'sleep'):
            results = conf.run_workflow(self.input, self.batch, self.batch, max_jobs=1,
                                        workdir=self.root/'search', solvent='toluene')
        self.assertEqual(len(results), 2)
        self.assertEqual(Path(results[0]['source']).name, 'Geometry-2-.log')
        report = (self.root/'search/conformation.lx').read_text()
        rows = [line.split() for line in report.splitlines() if not line.startswith('#')]
        self.assertEqual(len(rows), 2)
        self.assertAlmostEqual(sum(float(row[3]) for row in rows), 100, places=4)
        self.assertEqual([float(row[3]) for row in rows], [50., 50.])
        self.assertGreater(float(rows[0][6]), float(rows[1][6]))
        self.assertIn('Geometry-4-.log', (self.root/'search/rejected_conformers.csv').read_text())
        self.assertIn('--chrg', commands[0])
        self.assertIn('-1', commands[0])
        self.assertIn('--alpb', commands[0])
        self.assertNotIn('--cluster', commands[-1])


if __name__ == '__main__':
    unittest.main()
