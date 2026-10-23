"""Software fixtures only; no constructed data enter clinical results."""
import unittest

import numpy as np
import pandas as pd
import statsmodels.api as sm

from revision_clinical import (exact_within_sample, design_variation, draw_children,
                              event_frame, verify_never, mundlak_design)


class RevisionClinicalTests(unittest.TestCase):
    def test_exact_sample_keeps_zero_event_intervals(self):
        f = pd.DataFrame({"PtID": ["a", "a", "b", "b"], "unwarned": [1, 0, 0, 0], "pair_id": ["a1", "a2", "b1", "b2"]})
        self.assertEqual(exact_within_sample(f).pair_id.tolist(), ["a1", "a2"])

    def test_visit_adjustment_removes_pure_common_visit_variation(self):
        f = pd.DataFrame({"PtID": ["a"]*3 + ["b"]*3 + ["c"], "target_visit": [26, 39, 52]*2 + [26],
                          "auc_nmol_l": [1, .8, .6, 2, 1.8, 1.6, 1], "lag_days": [60]*7})
        result, child = design_variation(f)
        self.assertEqual(result["raw_varying_children"], 2)
        self.assertEqual(result["after_visit_varying_children"], 0)
        self.assertAlmostEqual(result["after_visit_sum_squares"], 0, places=20)
        self.assertEqual(child.loc[child.PtID.eq("c"), "within_ss"].iloc[0], 0)

    def test_stratified_child_draw_preserves_arm_counts(self):
        f = pd.DataFrame({"PtID": ["a", "a", "b", "c", "d"], "delivery_arm": ["i", "i", "i", "u", "u"]})
        x = draw_children(f, np.random.default_rng(20261023))
        y = draw_children(f, np.random.default_rng(20261023))
        self.assertEqual(x.tolist(), y.tolist())
        self.assertEqual(sum(p in ("a", "b") for p in x), 2)
        self.assertEqual(sum(p in ("c", "d") for p in x), 2)

    def test_event_ledger_count_disagreement_raises(self):
        f = pd.DataFrame({"PtID": ["a"], "pair_id": ["a1"], "events": [2], "matched": [1], "opportunity_events": [2], "unwarned": [1]})
        e = pd.DataFrame({"pair_id": ["a1"], "onset_time": [100], "warned": [True], "legal_opportunity": [True]})
        with self.assertRaises(AssertionError):
            event_frame(f, e)

    def test_never_identity_requires_identical_exposure(self):
        f = pd.DataFrame({"PtID": ["a"], "pair_id": ["a1"], "observed_bins": [288], "grid_bins": [300],
                          "observed_days": [1.], "events": [2], "start_seconds": [0], "end_seconds": [90000]})
        n = f.assign(unwarned=2, matched=0, alarms=0, unconfirmed=0)
        self.assertEqual(verify_never(f, n)["status"], "PASS")
        n["observed_days"] = 2.
        with self.assertRaises(AssertionError):
            verify_never(f, n)

    def test_expanded_and_grouped_binomial_likelihood_coefficients(self):
        x = np.array([[1., a] for a in [-1, 0, 1, 2]])
        counts = np.array([[1, 3], [2, 3], [3, 2], [4, 1]])
        expanded_x = np.repeat(x, counts.sum(axis=1), axis=0)
        y = np.concatenate([np.r_[np.ones(a), np.zeros(b)] for a, b in counts])
        fit1 = sm.GLM(y, expanded_x, family=sm.families.Binomial()).fit()
        fit2 = sm.GLM(counts, x, family=sm.families.Binomial()).fit()
        np.testing.assert_allclose(fit1.params, fit2.params, atol=1e-10)

    def test_mundlak_within_and_between_add_to_exposure(self):
        f = pd.DataFrame({"PtID": ["a", "a", "b", "b"], "auc_nmol_l": [.1, .3, .5, .9],
                          "target_visit": [26, 39, 26, 52], "lag_days": [60, 61, 59, 62], "age": [8, 8, 10, 10],
                          "management": ["usual care"]*2 + ["intensive/Tandem"]*2, "cohort": ["A"]*2 + ["B"]*2})
        x, removed = mundlak_design(f)
        np.testing.assert_allclose(x.auc_within + x.auc_between, f.auc_nmol_l/.20)
        np.testing.assert_allclose(x.auc_within.groupby(f.PtID).mean(), 0, atol=1e-15)
        self.assertIn("intensive_medtronic", removed)
        self.assertIn("child_mean_target_39", x)


if __name__ == "__main__":
    unittest.main()
