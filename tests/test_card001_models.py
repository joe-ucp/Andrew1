import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from card001_core import *
from card001_models import participant_weights,ranking,calibrate
from card001_evaluate import fast_grid,score_policy


class ModelTests(unittest.TestCase):
    def test_equal_patient_weights(self):
        group=np.array([1]*30+[2]*2+[3]*100)
        w=participant_weights(group)
        totals=[w[group==p].sum() for p in [1,2,3]]
        np.testing.assert_allclose(totals,np.repeat(totals[0],3))

    def test_no_test_calibration(self):
        with self.assertRaisesRegex(ValueError,"validation"):
            calibrate([dict(split="test")])

    def test_linear_tie_physical_units(self):
        score=dict(sensitivity=.5,unconfirmed_per_day=.2,total_per_day=.4)
        a=dict(kind="linear",threshold=-60,score=score)
        b=dict(kind="linear",threshold=-61,score=score)
        self.assertEqual(min([a,b],key=ranking)["threshold"],-61)

    def test_fast_replay_independent_damage(self):
        rng=np.random.default_rng(SEED)
        for _ in range(10):
            g=np.repeat(rng.choice([45.,80.,100.,np.nan],40),rng.integers(1,15,40))
            onsets=episodes(g)[0]
            scores=rng.random((4,len(g)))
            scores[rng.random(scores.shape)<.3]=np.nan
            kind=rng.integers(0,5,80)
            cut=rng.choice(np.r_[np.arange(0,1.01,.05),np.inf],80)
            short=rng.choice(np.r_[np.arange(0,1.01,.05),np.inf],80)
            # Availability masks model conditions, excluding low/missing bins.
            for j in range(4): scores[j,~ready(g,25 if j==0 else 4)]=np.nan
            p=dict(glucose=g,onsets=onsets,predictions={k:scores[j] for j,k in enumerate(["M120","M15","linear","missing"])})
            result=fast_grid(g,onsets,scores,kind,cut,short)
            for j in range(80):
                row=dict(kind=["M120","M15","linear","missing","switch"][kind[j]],threshold=cut[j],short_threshold=short[j])
                s=score_policy(p,row)
                self.assertEqual(result[j].tolist(),[s["matched"],s["alarms"]])

if __name__=="__main__": unittest.main()
