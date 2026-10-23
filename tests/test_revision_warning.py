"""Pure synthetic validation of new exploratory diagnostics; no patient data."""
import sys,unittest
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
import numpy as np
from revision_warning import availability_masks,event_availability,select_indices,bootstrap_differences
from card001_models import ranking,admissible
from card003_policy import select_primary


class RevisionWarningTests(unittest.TestCase):
    def test_nested_availability_and_event_operator(self):
        rng=np.random.default_rng(175)
        for length in [0,1,3,4,24,25,300]:
            for _ in range(20):
                g=rng.choice([45.,54.,80.,np.nan],length,p=[.05,.1,.8,.05])
                long,short,switch=availability_masks(g)
                self.assertFalse(np.any(long&~short));np.testing.assert_array_equal(switch,short)
                onsets=np.arange(length)
                a,b,c=[event_availability(x,onsets) for x in [long,short,switch]]
                self.assertFalse(np.any(a&~b));np.testing.assert_array_equal(b,c)

    def test_selector_matches_original_two_cards(self):
        rng=np.random.default_rng(191)
        spec=[]
        for kind in ["M120","M15","linear","missing"]:
            for t in ([-40.,-80.,-100.,np.inf] if kind=="linear" else [0.,.01,.04,np.inf]):
                spec.append(dict(kind=kind,threshold=t))
        for a in [0.,.04,np.inf]:
            for b in [0.,.01,.04,np.inf]:spec.append(dict(kind="switch",threshold=a,short_threshold=b))
        for trial in range(100):
            # Repeated integer counts exercise performance and threshold ties.
            matched=rng.integers(0,5,len(spec));alarm=matched+rng.integers(0,20,len(spec))
            for i,s in enumerate(spec):
                if np.isinf(s["threshold"]) and (s["kind"]!="switch" or np.isinf(s["short_threshold"])):matched[i]=alarm[i]=0
            counts=np.column_stack([matched,alarm]);days=float(rng.integers(3,30));rows=[]
            for i,s in enumerate(spec):rows.append(dict(s,idx=i,score=dict(sensitivity=matched[i]/10,unconfirmed_per_day=(alarm[i]-matched[i])/days,total_per_day=alarm[i]/days)))
            families={k:min([r for r in rows if r["kind"]==k and admissible(r["score"])],key=ranking) for k in ["M120","M15","linear","missing"]}
            comp=min(families.values(),key=ranking)
            switch=min([r for r in rows if r["kind"]=="switch" and admissible(r["score"],comp["score"]["unconfirmed_per_day"]+.1)],key=ranking)
            third=select_primary(rows)["primary"]
            actual=select_indices(counts,days,spec)
            self.assertEqual(actual[:3],(comp["idx"],switch["idx"],third["idx"]))

    def test_cluster_bootstrap_keeps_zero_event_people(self):
        data=np.array([[2,1,0,3,2,10],[0,0,0,1,2,5]],float)
        result=bootstrap_differences(data,1000)
        self.assertGreater(result["undefined"],0)
        np.testing.assert_array_equal(result["sensitivity_difference_ci"],[.5,.5])
        self.assertEqual(result,bootstrap_differences(data,1000))


if __name__=="__main__":unittest.main()
