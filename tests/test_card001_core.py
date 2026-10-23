"""Implementation damage tests, not patient analyses."""
import sys
import unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from card001_core import *


class EndpointTests(unittest.TestCase):
    def test_empty_boundary(self):
        self.assertEqual(len(ready([],4)),0)
        self.assertEqual(len(episodes_independent([])[0]),0)
        self.assertEqual(len(labels([])),0)
        for kind in ["M15","M120","missing"]:
            matrix,names=features([],kind)
            self.assertEqual(matrix.shape,(0,len(names)))
        self.assertEqual(len(replay([],[],.5)),0)
        with self.assertRaisesRegex(ValueError,"No participants"):
            paired_bootstrap(np.empty((0,6)))

    def test_recovery_rearming_and_censor(self):
        g=np.array([50]*4+[100]*4+[50]*4+[100]*2+[50]*5+[100]*4+[50]*4+[np.nan]+[50]*4)
        events,censored=episodes(g)
        self.assertEqual(events.tolist(),[8,23])
        self.assertEqual(censored.tolist(),[0,28])
        for a,b in zip(episodes(g),episodes_independent(g)): np.testing.assert_array_equal(a,b)

    def test_random_independent_labels(self):
        rng=np.random.default_rng(42)
        for k in [3,4]:
            for _ in range(300):
                # Long runs stress recovery, interruption, event persistence.
                g=np.repeat(rng.choice([40.,53.,54.,100.,np.nan],40),rng.integers(1,10,40))
                for a,b in zip(episodes(g,k),episodes_independent(g,k)): np.testing.assert_array_equal(a,b)

    def test_15_minutes_means_four_timestamps(self):
        self.assertEqual(episodes([100]*4+[50]*3)[0].size,0)
        self.assertEqual(episodes([100]*4+[50]*4)[0].tolist(),[4])

    def test_gap_low_onset_not_invented(self):
        e,c=episodes([100]*30+[np.nan]+[40]*8)
        self.assertEqual(e.size,0); self.assertEqual(c.tolist(),[31])

    def test_match_lead_bounds_and_refractory(self):
        x=np.full(40,100.)
        a=replay(x,np.ones(40),.5)
        self.assertEqual(a.tolist(),[0,6,12,18,24,30,36])
        ledger,matched=match_alarms([0,1,2],[4,5],x)
        self.assertEqual(len(matched),2)
        self.assertEqual(ledger[-1]["category"],"duplicate")
        self.assertEqual(independent_match_count([0,1,2],[4,5]),2)

    def test_shared_burden_and_no_alarm_low(self):
        x=np.array([100]*30+[np.nan]*20+[50]*4,dtype=float)
        p=np.ones(54); a=replay(x,p,.5)
        self.assertTrue(all(i<30 for i in a))
        score=participant_score(x,[],a,np.isfinite(p))
        self.assertEqual(score["observed_days"],34/288)

    def test_grid_latest_no_future_and_conflict(self):
        r=causal_grid([0,299,301,600,600,900],[100,101,99,60,61,80],start=0,end=900)
        np.testing.assert_allclose(r["glucose"],[100,101,np.nan,80],equal_nan=True)
        self.assertTrue(r["conflict"][2])
        self.assertTrue(np.all(r["source_time"]<=r["time"]))

    def test_delayed_reading(self):
        r=causal_grid([0,300],[100,50],start=0,end=600,delay_seconds=300)
        np.testing.assert_allclose(r["glucose"],[np.nan,100,50],equal_nan=True)

    def test_causal_prefix_attack(self):
        rng=np.random.default_rng(1); x=rng.normal(110,30,200); x[20:30]=np.nan
        altered=x.copy(); altered[101:]=rng.normal(40,5,99)
        for kind in ["M120","M15","missing"]:
            a,_=features(x,kind); b,_=features(altered,kind)
            np.testing.assert_allclose(a[:101],b[:101],equal_nan=True)
        np.testing.assert_array_equal(ready(x,25)[:101],ready(altered,25)[:101])
        t=np.arange(200)*300+150
        a=causal_grid(t,x,start=0,end=60000); b=causal_grid(t,altered,start=0,end=60000)
        np.testing.assert_allclose(a["glucose"][:102],b["glucose"][:102],equal_nan=True)

    def test_training_label_censoring(self):
        g=np.array([100.]*20+[50.]*4+[100.]*20)
        g[5]=np.nan
        y=labels(g)
        self.assertTrue(np.isnan(y[0]))
        self.assertEqual(y[14:18].tolist(),[1.,1.,1.,1.])
        self.assertTrue(np.isnan(y[-1]))

    def test_availability_bound(self):
        x=np.r_[np.full(10,100.),np.full(4,50.),np.nan,np.full(30,100.),np.full(4,50.)]
        e=episodes(x)[0]; available=ready(x,4)&(x>=54)
        alarms=replay(x,np.where(available,1,np.nan),.5)
        s=participant_score(x,e,alarms,available)
        self.assertLessEqual(s["matched"],s["opportunities"])
        self.assertLessEqual(s["opportunities"],s["events"])

    def test_split_identity(self):
        self.assertEqual(participant_split("00042"),participant_split(42))
        self.assertEqual(participant_split(109),("train","generation_exposed"))

    def test_bootstrap_undefined_retained(self):
        result=paired_bootstrap([[0,0,0,1,1,1]],n_boot=20)
        self.assertEqual(result["undefined"],20)
        self.assertIsNone(result["sensitivity_difference_ci"])

if __name__=="__main__": unittest.main()
