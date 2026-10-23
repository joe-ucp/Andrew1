"""Independent arithmetic and resampling checks for the causal contrast."""
import sys
from pathlib import Path
import unittest

ROOT=Path(__file__).resolve().parents[1]
sys.path.insert(0,str(ROOT/"src"))
from card002_run import estimate, score_expression, bootstrap
import numpy as np
from scipy.stats import t


class Card002Tests(unittest.TestCase):
    def test_four_cell_contrast_known_effects_and_welch_variance(self):
        y=np.array([40,50,60,10,20,30,30,40,50,20,30,40],dtype=float)
        arm=np.tile([True,True,True,False,False,False],2)
        high=np.repeat([False,True],6)
        out=estimate(y,arm,high)
        self.assertAlmostEqual(out["delta_low"]["estimate_weeks"],30)
        self.assertAlmostEqual(out["delta_high"]["estimate_weeks"],10)
        self.assertAlmostEqual(out["theta"]["estimate_weeks"],20)
        expected_se=np.sqrt(4*100/3)
        self.assertAlmostEqual(out["theta"]["se_weeks"],expected_se)
        self.assertAlmostEqual(out["theta"]["df"],8)
        self.assertAlmostEqual(out["theta"]["ci95_weeks"][0],20-t.ppf(.975,8)*expected_se)

    def test_constant_additive_treatment_benefit_has_zero_interaction(self):
        placebo_low=np.arange(5)*2+15
        placebo_high=np.arange(5)*3+70
        y=np.concatenate([placebo_low+25,placebo_low,placebo_high+25,placebo_high])
        arm=np.tile(np.repeat([True,False],5),2)
        high=np.repeat([False,True],10)
        out=estimate(y,arm,high)
        self.assertEqual(out["theta"]["estimate_weeks"],0)
        self.assertNotEqual(out["placebo_low_minus_high"]["estimate_weeks"],0)

    def test_score_invariant_to_positive_gene_affine_representation(self):
        rng=np.random.default_rng(20)
        x=rng.normal(size=(23,7))
        score,high=score_expression(x)
        other,other_high=score_expression(x*np.arange(1,8)+np.arange(7))
        np.testing.assert_allclose(score,other,atol=1e-14)
        np.testing.assert_array_equal(high,other_high)
        self.assertEqual(high.sum(),11)

    def test_bootstrap_reproducible_and_recomputes_scores(self):
        rng=np.random.default_rng(11)
        x=rng.normal(size=(48,8))
        arm=np.arange(48)%2==0
        y=np.linspace(10,156,48)
        first,summary=bootstrap(x,y,arm,n=50,seed=27)
        second,_=bootstrap(x,y,arm,n=50,seed=27)
        np.testing.assert_array_equal(first,second)
        self.assertEqual(summary["valid"]+summary["unusable"],50)
        self.assertIn("pooled median",summary["recomputed_in_each_resample"])

if __name__=="__main__":
    unittest.main()
