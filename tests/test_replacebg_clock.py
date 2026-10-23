"""Clock transformation fixtures only; no REPLACE-BG clinical data simulated."""
import sys
from pathlib import Path
import unittest
import numpy as np
import pandas as pd

sys.path.insert(0,str(Path(__file__).resolve().parents[1]/"src"))
from card001_replacebg import relative_clock_seconds,clock_barriers,canonical_arrays,cohort_rows


class ReplaceBGClockTests(unittest.TestCase):
    def test_day_rollover_and_twelve_hour_clocks(self):
        actual=relative_clock_seconds([0,1,0,0,-1],["23:59:59","00:00:00","12:00 AM","12:00 PM","23:00"])
        np.testing.assert_array_equal(actual,[86399,86400,0,43200,-3600])

    def test_invalid_clocks_cannot_silently_become_midnight(self):
        actual=relative_clock_seconds([0,0,0,.5,None,0],["24:00","12:60","00:00 PM","01:00","01:00","5400"])
        self.assertTrue(np.isnan(actual).all())

    def test_clock_reset_marks_barrier_without_time_warp(self):
        times,reasons=clock_barriers([0,300,600,900],[1000,1300,500,800])
        np.testing.assert_array_equal(times,[600])
        self.assertEqual(reasons.tolist(),["internal_clock_reset"])

    def test_conflicting_duplicate_retained_as_missing(self):
        grid,raw=canonical_arrays([0,300,300,600],[100,40,90,110],[0,300,300,600],[1,2,3,4],0,900)
        np.testing.assert_array_equal(grid["time"],[0,300,600])
        self.assertTrue(np.isnan(grid["glucose"][1]))
        self.assertTrue(grid["conflict"][1])
        self.assertEqual(len(raw["time"]),4)

    def test_latest_causal_reading_and_exclusive_end(self):
        grid,raw=canonical_arrays([1,250,301,600,900],[120,100,90,80,60],[1,250,301,600,900],range(5),0,900)
        self.assertTrue(np.isnan(grid["glucose"][0]))
        self.assertEqual(grid["glucose"][1],100)
        self.assertEqual(grid["source_time"][1],250)
        self.assertEqual(grid["glucose"][2],80)
        self.assertEqual(raw["time"][-1],600)
        self.assertTrue((grid["source_time"][np.isfinite(grid["source_time"])]<=grid["time"][np.isfinite(grid["source_time"])]).all())

    def test_missing_participant_keeps_followup_grid(self):
        grid,raw=canonical_arrays([],[],[],[],0,900)
        self.assertEqual(len(grid["time"]),3)
        self.assertTrue(np.isnan(grid["glucose"]).all())
        self.assertEqual(len(raw["time"]),0)

    def test_same_time_internal_disagreement_is_order_invariant(self):
        args=([0,300,300,600],[100,100,100,100],[0,300,10000,600],[1,2,3,4])
        a,araw=canonical_arrays(*args,0,900)
        b,braw=canonical_arrays(args[0],args[1],[0,10000,300,600],[1,3,2,4],0,900)
        np.testing.assert_equal(a["glucose"],b["glucose"])
        np.testing.assert_equal(araw["barrier_times"],braw["barrier_times"])
        self.assertTrue(np.isnan(a["glucose"][1]))
        self.assertEqual(araw["barrier_reasons"].tolist(),["same_time_internal_disagreement"])

    def test_cohort_excludes_exposed_and_nonrandomized_without_glucose(self):
        roster=pd.DataFrame(dict(PtID=[1,2,183],RandDtDaysAfterEnroll=[60,np.nan,60],TrtGroup=["CGM Only"]*3,PtStatus=["Dropped","Completed","Completed"]))
        visits=pd.DataFrame(dict(PtID=[1,1,1,183],VisitDtDaysFromEnroll=[80,90,100,250],VisitMiss=[np.nan,np.nan,1,np.nan]))
        final=pd.DataFrame(dict(PtID=[1],DeathDtDaysAfterEnroll=[np.nan]))
        rows,exclusions=cohort_rows(roster,visits,final,Path("fixture"))
        self.assertEqual([r["id"] for r in rows],[1])
        self.assertEqual(rows[0]["start"],0)
        self.assertEqual(rows[0]["end"],91*86400)
        self.assertEqual(rows[0]["randomization_day"],60)
        self.assertEqual({x["reason"] for x in exclusions},{"not_randomized","generation_schema_exposure"})

    def test_fixed_horizon_convention_has_separate_inclusive_sensitivity(self):
        roster=pd.DataFrame(dict(PtID=[1],RandDtDaysAfterEnroll=[30],TrtGroup=["CGM+BGM"],PtStatus=["Completed"]))
        visits=pd.DataFrame(dict(PtID=[1],VisitDtDaysFromEnroll=[250],VisitMiss=[np.nan]))
        final=pd.DataFrame(dict(PtID=[],DeathDtDaysAfterEnroll=[]))
        primary,_=cohort_rows(roster,visits,final,Path("fixture"))
        inclusive,_=cohort_rows(roster,visits,final,Path("fixture"),days=183)
        self.assertEqual(primary[0]["end"],182*86400)
        self.assertEqual(inclusive[0]["end"],183*86400)

if __name__=="__main__":
    unittest.main()
