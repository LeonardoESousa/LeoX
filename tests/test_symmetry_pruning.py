import csv
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

import numpy as np
from lx import conf_search as conf


class SymmetryTests(unittest.TestCase):
    def setUp(self):
        self.pair = conf.read_xyz(Path(__file__).parent / 'data/mirror_pair.xyz')
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.folder = Path(self.temp.name)
        (self.folder / 'CREST').mkdir()
        (self.folder / 'Geometries').mkdir()

    def test_real_pair_and_permuted_rotated_duplicate(self):
        first, second = self.pair
        self.assertGreater(conf.symmetry_rmsd(first, second), 0.125)
        self.assertLess(conf.symmetry_rmsd(first, second, True), 0.0001)
        permutation = np.arange(80)[::-1]
        duplicate = dict(first, atoms=[first['atoms'][i] for i in permutation],
                         geometry=first['geometry'][permutation] @ np.array([[0,-1,0],[1,0,0],[0,0,1]]) + 4)
        groups = conf.prune_symmetry([first, second, duplicate, second], .125, .05)
        self.assertEqual(len(groups), 1)
        self.assertEqual(len(groups[0]['minima']), 2)
        self.assertEqual(len(conf.prune_symmetry([first, second], .125, .05, False)), 2)
        self.assertEqual(len(conf.prune_symmetry([first, duplicate], .125, .05)[0]['minima']), 1)

    def test_pruning_preserves_original_job_numbers(self):
        self.pair[1]['comment'] = '-2823.735677'
        retained = conf.sampling_prune(self.pair, self.folder / 'CREST', .125, .05, True)
        self.assertEqual(len(retained), 1)
        self.assertEqual(retained[0]['crest_index'], 2)
        template = dict(atoms=retained[0]['atoms'], route='# b3lyp/6-31G opt freq', bottom='',
                        link0=[], nproc=2, mem='1GB', title='Test', charge=0, multiplicity=1)
        files = conf.make_gaussian_inputs(template, retained, self.folder / 'Geometries')
        self.assertEqual(files, ['Geometry-2-.com'])

    def classify(self, results):
        for item in results:
            Path(item['source']).write_text('mock completed log')
        def result(filename):
            return next(item for item in results if Path(item['source']) == filename)
        def cregen(command, folder, logfile):
            (folder / logfile).write_text('CREST terminated normally')
            conf.write_xyz(folder / 'reoptimized.xyz.sorted', conf.read_xyz(folder / 'reoptimized.xyz'))
        with patch.object(conf, 'gaussian_result', side_effect=result), \
             patch.object(conf, 'run_command', side_effect=cregen):
            return conf.classify_only(self.folder)

    def minimum(self, structure, number):
        return dict(structure, energy=-40., gibbs=-39.99,
                    comment='-40.000000000000',
                    source=str(self.folder / 'Geometries' / f'Geometry-{number}-.log'))

    def test_skipped_partner_survives_dft_and_populations(self):
        conf.sampling_prune(self.pair, self.folder / 'CREST', .125, .05, True)
        first = self.minimum(self.pair[0], 1)
        other = self.minimum(dict(atoms=['C','H'], geometry=np.array([[0,0,0],[1,0,0]])), 3)
        # Use the same molecule but a clearly different, non-mirror shape.
        other = dict(first, geometry=first['geometry'] * 1.1, source=other['source'])
        results = self.classify([first, other])
        self.assertEqual([x['multiplicity'] for x in results], [2, 1])
        with open(self.folder / 'conformation.csv') as handle:
            rows = list(csv.DictReader(handle))
        self.assertEqual(rows[0]['Enantiomers'], 'yes')
        self.assertEqual(rows[0]['Partner_Gaussian_log'], 'inferred mirror partner')
        self.assertAlmostEqual(float(rows[0]['PopG_300K_percent']), 200/3, places=4)
        self.assertAlmostEqual(sum(float(row['PopE_300K_percent']) for row in rows), 100, places=4)

    def test_both_dft_partners_sum_actual_gibbs_weights(self):
        first, second = [self.minimum(s, i+1) for i,s in enumerate(self.pair)]
        second['gibbs'] += .0001
        other = dict(first, geometry=first['geometry'] * 1.1,
                     source=str(self.folder/'Geometries/Geometry-4-.log'))
        results = self.classify([first, second, dict(first, source=str(self.folder/'Geometries/Geometry-3-.log')), other])
        self.assertEqual(len(results), 2)
        self.assertEqual(results[0]['multiplicity'], 2)
        self.assertEqual(len(results[0]['minima']), 2)
        self.assertEqual(results[0]['partner_log'], 'Geometry-2-.log')
        with open(self.folder / 'conformation.csv') as handle:
            rows = list(csv.DictReader(handle))
        partner_weight = np.exp(-.0001 * conf.HARTREE_EV / (conf.lx.parser.BOLTZ_EV * 300))
        self.assertAlmostEqual(float(rows[0]['PopG_300K_percent']),
                               100 * (1 + partner_weight) / (2 + partner_weight), places=4)

    def test_achiral_dft_collapse_does_not_keep_multiplicity(self):
        conf.sampling_prune(self.pair, self.folder / 'CREST', .125, .05, True)
        achiral = self.minimum(dict(atoms=['C','H'], geometry=np.array([[0,0,0],[1,0,0]])), 1)
        result = self.classify([achiral])[0]
        self.assertEqual(result['multiplicity'], 1)
        self.assertFalse(result['enantiomers'])

    def test_energy_inconsistency_does_not_remove_structure(self):
        self.pair[1]['comment'] = str(float(self.pair[0]['comment']) + .001)
        self.assertEqual(len(conf.prune_symmetry(self.pair, .125, .05)), 2)


if __name__ == '__main__':
    unittest.main()
