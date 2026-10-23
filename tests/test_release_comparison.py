"""Execution metadata must not mask a changed scientific result."""
import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'src'))
from compare_results import compare_file

RELATIVE = 'revision_20261022_r2/card003/results/never_warning_identity.json'


class NestedExecutionMetadataTests(unittest.TestCase):
    def setUp(self):
        self.folder = tempfile.TemporaryDirectory()
        self.addCleanup(self.folder.cleanup)
        self.reference = Path(self.folder.name) / 'reference.json'
        self.actual = Path(self.folder.name) / 'actual.json'
        self.base = {'never_events': 408, 'original_total_result_reused': {
            'elapsed_seconds': 28.4, 'rr': 0.774, 'events': 408, 'seed': 20261022,
            'bootstrap_95_ci': [0.674, 0.870], 'coefficients': {'exposure': -0.256},
            'status': 'estimated', 'design_columns': ['intercept', 'exposure']}}
        self.reference.write_text(json.dumps(self.base), encoding='utf-8')

    def write_actual(self, value):
        self.actual.write_text(json.dumps(value), encoding='utf-8')

    def changed_runtime(self):
        value = copy.deepcopy(self.base)
        value['original_total_result_reused']['elapsed_seconds'] = 21.1
        return value

    def test_exact_runtime_path_is_reported(self):
        self.write_actual(self.changed_runtime())
        result = compare_file(self.reference, self.actual, RELATIVE)
        self.assertEqual(result['excluded_execution_metadata'], {
            'original_total_result_reused.elapsed_seconds': {
                'authoritative': 28.4, 'reproduced': 21.1}})

    def test_scientific_neighbors_remain_strict(self):
        mutations = {'rr': 0.9, 'events': 409, 'seed': 20261023,
                     'bootstrap_95_ci': [0.674, 0.88], 'coefficients': {'exposure': -0.25},
                     'status': 'unestimated', 'design_columns': ['exposure', 'intercept'],
                     'unexpected_estimate': 1}
        for key, changed in mutations.items():
            with self.subTest(field=key):
                value = self.changed_runtime()
                value['original_total_result_reused'][key] = changed
                self.write_actual(value)
                with self.assertRaises((AssertionError, KeyError)):
                    compare_file(self.reference, self.actual, RELATIVE)
        value = self.changed_runtime()
        del value['original_total_result_reused']['rr']
        self.write_actual(value)
        with self.assertRaises((AssertionError, KeyError)):
            compare_file(self.reference, self.actual, RELATIVE)

    def test_exclusion_is_not_global(self):
        self.write_actual(self.changed_runtime())
        with self.assertRaises(AssertionError):
            compare_file(self.reference, self.actual, 'different_file.json')
        value = self.changed_runtime()
        value['elapsed_seconds'] = 12
        self.write_actual(value)
        with self.assertRaises(AssertionError):
            compare_file(self.reference, self.actual, RELATIVE)
        value = self.changed_runtime()
        value['never_events'] = 409
        self.write_actual(value)
        with self.assertRaises(AssertionError):
            compare_file(self.reference, self.actual, RELATIVE)


if __name__ == '__main__':
    unittest.main()
