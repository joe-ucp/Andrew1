"""Software checks for input masking without manufacturing outcome labels."""
import sys
from pathlib import Path
import unittest
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from card001_stress import TRACE_BINS,SCORE_START,SCORE_END,complete_windows,mask_positions,mask_trace,trace_score


class MaskStressTests(unittest.TestCase):
    def test_window_selection_depends_only_on_observation_mask(self):
        a=np.ones(TRACE_BINS*10)
        a[TRACE_BINS*3]=np.nan
        b=a*40
        self.assertEqual(complete_windows(a),complete_windows(b))
        self.assertEqual(len(complete_windows(a)),7)
        for start,stop in complete_windows(a):
            self.assertEqual(stop-start,TRACE_BINS)
            self.assertTrue(np.isfinite(a[start:stop]).all())

    def test_seed_schedule_repeatable_and_masks_nested(self):
        first=mask_positions(20261022,4)
        np.testing.assert_array_equal(first,mask_positions(20261022,4))
        for starts in first:
            short=mask_trace(np.ones(TRACE_BINS),starts,5)
            long=mask_trace(np.ones(TRACE_BINS),starts,60)
            self.assertEqual(np.isnan(short).sum(),4)
            self.assertEqual(np.isnan(long).sum(),48)
            self.assertTrue((~np.isnan(short)|np.isnan(long)).all())
            self.assertTrue(np.isfinite(long[:SCORE_START]).all())
            self.assertTrue(np.isfinite(long[SCORE_END:]).all())

    def test_input_mask_changes_warning_but_never_truth_episode(self):
        truth=np.full(TRACE_BINS,100.)
        truth[100:104]=40
        probabilities=np.zeros(TRACE_BINS)
        probabilities[96]=1
        policy={"kind":"M15","threshold":.5}
        baseline=trace_score(truth,truth,{"M15":probabilities},policy)
        masked=truth.copy(); masked[90:102]=np.nan
        after=trace_score(truth,masked,{"M15":probabilities},policy)
        self.assertEqual(baseline["events"],1)
        self.assertEqual(after["events"],1)
        self.assertEqual(baseline["matched"],1)
        self.assertEqual(after["matched"],0)
        self.assertEqual(baseline["truth_observed_days"],after["truth_observed_days"])
        self.assertLess(after["input_observed_days"],baseline["input_observed_days"])

    def test_rejects_masked_outcome_truth(self):
        truth=np.ones(TRACE_BINS)*100; truth[40]=np.nan
        with self.assertRaises(ValueError):
            trace_score(truth,truth,{"M15":np.zeros(TRACE_BINS)},{"kind":"M15","threshold":.5})

if __name__=="__main__":
    unittest.main()
