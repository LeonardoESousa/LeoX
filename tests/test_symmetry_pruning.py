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

    def classify(self, results, reverse_cregen=False):
        for item in results:
            Path(item['source']).write_text('mock completed log')
        def result(filename):
            return next(item for item in results if Path(item['source']) == filename)
        def cregen(command, folder, logfile):
            (folder / logfile).write_text('CREST terminated normally')
            structures = conf.read_xyz(folder / 'reoptimized.xyz')
            conf.write_xyz(folder / 'reoptimized.xyz.sorted', structures[::-1] if reverse_cregen else structures)
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
        self.assertEqual([x['multiplicity'] for x in results], [2, 2])
        with open(self.folder / 'conformation.csv') as handle:
            rows = list(csv.DictReader(handle))
        self.assertIn('Multiplicity', rows[0])
        self.assertNotIn('Enantiomers', rows[0])
        self.assertNotIn('Partner_Gaussian_log', rows[0])
        with open(self.folder / 'conformers_manifest.csv') as handle:
            manifest = {row['Gaussian_log']: row for row in csv.DictReader(handle)}
        details = manifest[rows[0]['Gaussian_log']]
        self.assertEqual(details['Group'], rows[0]['Group'])
        self.assertEqual(details['Enantiomers'], 'yes')
        self.assertEqual(details['Partner_Gaussian_log'], 'inferred mirror partner')
        self.assertAlmostEqual(float(rows[0]['PopG_300K_percent']), 50., places=4)
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
                               100 * (1 + partner_weight) / (3 + partner_weight), places=4)

    def test_unsampled_partner_has_same_multiplicity_as_observed_pair(self):
        first = self.minimum(self.pair[0], 1)
        result = self.classify([first])[0]
        self.assertEqual(result['multiplicity'], 2)
        self.assertEqual(result['partner_log'], 'inferred mirror partner')
        self.assertFalse((self.folder / 'CREST/symmetry_groups.json').exists())

    def test_observing_only_the_higher_energy_partner_does_not_bias_populations(self):
        low = self.minimum(self.pair[0], 1)
        high = dict(self.minimum(self.pair[0], 2), geometry=low['geometry'] * 1.1,
                    energy=low['energy'] + .0001, gibbs=low['gibbs'] + .0001,
                    comment='-39.999900000000')
        partner = dict(high, geometry=high['geometry'] * np.array([-1., 1., 1.]),
                       source=str(self.folder / 'Geometries/Geometry-3-.log'))
        results = self.classify([low, high, partner])
        self.assertEqual([item['multiplicity'] for item in results], [2, 2])
        with open(self.folder / 'conformation.csv') as handle:
            rows = list(csv.DictReader(handle))
        for column in ('PopE_300K_percent', 'PopG_300K_percent'):
            self.assertGreater(float(rows[0][column]), float(rows[1][column]))

    def test_achiral_dft_collapse_does_not_keep_multiplicity(self):
        conf.sampling_prune(self.pair, self.folder / 'CREST', .125, .05, True)
        achiral = self.minimum(dict(atoms=['C','H'], geometry=np.array([[0,0,0],[1,0,0]])), 1)
        result = self.classify([achiral])[0]
        self.assertEqual(result['multiplicity'], 1)
        self.assertFalse(result['enantiomers'])

    def test_energy_inconsistency_does_not_remove_structure(self):
        self.pair[1]['comment'] = str(float(self.pair[0]['comment']) + .001)
        self.assertEqual(len(conf.prune_symmetry(self.pair, .125, .05)), 2)

    def test_xyz_blocks_match_csv_groups_even_with_energy_ties_and_reversed_cregen(self):
        first = self.minimum(self.pair[0], 9)
        second = dict(self.minimum(self.pair[0], 2), geometry=first['geometry'] * 1.1)
        third = dict(self.minimum(self.pair[0], 4), geometry=first['geometry'] * 1.2,
                     energy=-40.001, comment='-40.001000000000')
        results = self.classify([first, second, third], reverse_cregen=True)
        with open(self.folder / 'conformation.csv') as handle:
            rows = list(csv.DictReader(handle))
        structures = conf.read_xyz(self.folder / 'conformers_unique.xyz')
        self.assertEqual([row['Gaussian_log'] for row in rows],
                         ['Geometry-4-.log', 'Geometry-2-.log', 'Geometry-9-.log'])
        self.assertEqual(len(structures), len(rows))
        for group, (structure, row, result) in enumerate(zip(structures, rows, results), 1):
            self.assertIn(f'Group={group} ', structure['comment'])
            self.assertIn('Gaussian_log=' + row['Gaussian_log'], structure['comment'])
            self.assertEqual(float(structure['comment'].split()[0]), float(row['E_Hartree']))
            np.testing.assert_allclose(structure['geometry'], result['geometry'], atol=1e-9)

    def test_cached_classification_restores_partner_from_manifest(self):
        first, second = [self.minimum(s, i+1) for i,s in enumerate(self.pair)]
        second['gibbs'] += .0001
        self.classify([first, second])
        signature = {'pruning': conf.PRUNING_VERSION}
        outputs = {name: conf.file_digest(self.folder / name) for name in
                   ('conformation.csv', 'conformers_manifest.csv', 'conformers_unique.xyz')}
        conf.save_json(self.folder / 'CREGEN/classification.json',
                       dict(signature=signature, outputs=outputs))
        def result(filename):
            return dict(next(item for item in (first, second) if Path(item['source']) == filename))
        with patch.object(conf, 'gaussian_result', side_effect=result):
            cached = conf.cached_classification(self.folder, signature)
        self.assertIsNotNone(cached)
        self.assertEqual(cached[0]['multiplicity'], 2)
        self.assertTrue(cached[0]['enantiomers'])
        self.assertEqual(cached[0]['partner_log'], 'Geometry-2-.log')
        self.assertAlmostEqual(cached[0]['minima'][1]['gibbs'], second['gibbs'])


if __name__ == '__main__':
    unittest.main()
