import contextlib
import csv
import io
import os
import tempfile
import unittest
from pathlib import Path
from unittest.mock import patch

from lx import __main__ as menu
from lx import conf_search as conf


class ClassificationMenuTests(unittest.TestCase):
    def setUp(self):
        temporary = tempfile.TemporaryDirectory()
        self.addCleanup(temporary.cleanup)
        self.root = Path(temporary.name)
        previous = Path.cwd()
        os.chdir(self.root)
        self.addCleanup(os.chdir, previous)

    def run_menu(self, folder_answer):
        with patch('builtins.input', side_effect=['5', 'y', folder_answer]), \
             patch.object(menu.lx.tools, 'check_for_updates'), \
             patch.object(conf, 'classify_only', return_value=[]) as classify, \
             contextlib.redirect_stdout(io.StringIO()):
            menu.interface()
        return classify.call_args

    def test_parent_folder_and_saved_charge_spin(self):
        (self.root / 'Conformational/Inputs').mkdir(parents=True)
        (self.root / 'Conformational/Inputs/search.json').write_text('{"charge": -1, "uhf": 2}')
        self.assertEqual(self.run_menu('').args, ('Conformational', 'crest', .125, .05, -1, 2, True))

    def test_current_search_folder_and_custom_path(self):
        (self.root / 'Geometries').mkdir()
        (self.root / 'Conformational').mkdir()
        self.assertEqual(self.run_menu('').args[0], '.')
        self.assertEqual(self.run_menu('PreviousRun').args[0], 'PreviousRun')

    def test_actual_error_is_not_hidden(self):
        error = io.StringIO()
        with patch('builtins.input', side_effect=['5', 'y', 'MissingRun']), \
             patch.object(menu.lx.tools, 'check_for_updates'), \
             patch.object(conf, 'classify_only', side_effect=ValueError('Actual classification error')), \
             contextlib.redirect_stdout(io.StringIO()), contextlib.redirect_stderr(error):
            menu.interface()
        self.assertIn('Actual classification error', error.getvalue())


class GroupPopulationTests(unittest.TestCase):
    def test_higher_energy_pair_can_outweigh_single_minimum(self):
        with tempfile.TemporaryDirectory() as temporary:
            low = dict(energy=-40., gibbs=-39.99, source='low.log', multiplicity=1)
            high = dict(energy=-40. + .01/conf.HARTREE_EV,
                        gibbs=-39.99 + .01/conf.HARTREE_EV, source='high.log',
                        multiplicity=2, enantiomers=True)
            high['minima'] = [dict(high), dict(high)]
            report = Path(temporary) / 'report.csv'
            conf.write_report([low, high], report)
            with report.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertGreater(float(rows[1]['PopE_300K_percent']), float(rows[0]['PopE_300K_percent']))
            self.assertAlmostEqual(sum(float(row['PopE_300K_percent']) for row in rows), 100, places=4)
            high['minima'] = [high['minima'][0]]
            high['multiplicity'] = 1
            conf.write_report([low, high], report)
            with report.open() as handle:
                rows = list(csv.DictReader(handle))
            self.assertLess(float(rows[1]['PopE_300K_percent']), float(rows[0]['PopE_300K_percent']))


if __name__ == '__main__':
    unittest.main()
