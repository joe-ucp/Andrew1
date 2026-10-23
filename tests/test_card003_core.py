"""Software fixtures only. These contain no real or simulated study results."""
import sys
sys.dont_write_bytecode = True
import unittest
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import numpy as np
import pandas as pd
from card003_core import (normalized_auc, clinical_remission, derive_mmt_visit,
    later_interval, overlapping_pairs, independent_episode_labeler,
    residual_identity, fixed_design, gate_assessment, replay_interval)
from card003_core import poisson_separation, fit_association
from card001_core import episodes
from card003_policy import prescribed_folds, select_primary


class ScientificKernels(unittest.TestCase):
    def assay(self):
        return pd.DataFrame(dict(minute=[0,15,30,60,90,120], value=[1]*6,
            units=["nmol/L"]*6, date=["2024-01-01"]*6, status=["Complete"]*6,
            fasting=["Fasting"]*6, temperature=["OK"]*6))

    def test_auc_normalized_units(self):
        self.assertEqual(normalized_auc([1]*6), 1)
        self.assertAlmostEqual(normalized_auc(np.array([0,15,30,60,90,120])/120), .5)
        with self.assertRaises(ValueError): normalized_auc([1,1,np.nan,1,1,1])
        with self.assertRaises(ValueError): normalized_auc([1]*6, [0,15,30,60,90,121])

    def test_analytical_duplicates(self):
        a = self.assay()
        self.assertTrue(derive_mmt_visit(pd.concat([a,a.iloc[:1]]))["eligible"])
        b = a.iloc[:1].copy(); b["value"] = 2
        self.assertFalse(derive_mmt_visit(pd.concat([a,b]))["eligible"])

    def test_fasting_units_and_temperature(self):
        a = self.assay(); a.loc[0,"fasting"] = "Unknown"
        self.assertFalse(derive_mmt_visit(a)["eligible"])
        a = self.assay(); a["units"] = "pmol/L"
        self.assertFalse(derive_mmt_visit(a)["eligible"])
        a = self.assay(); a.loc[3,"temperature"] = "WARM"
        self.assertTrue(derive_mmt_visit(a)["eligible"])
        self.assertFalse(derive_mmt_visit(a,require_temperature_ok=True)["eligible"])

    def test_idaa1c_and_weight(self):
        args = ("2024-01-01",)*4
        a = clinical_remission(7,20,40,"kg",*args)
        b = clinical_remission(7,20,40/.45359237,"lb",*args)
        self.assertEqual(a["idaa1c"],9)
        self.assertEqual(a["remission"],1)
        self.assertAlmostEqual(a["idaa1c"],b["idaa1c"])
        self.assertTrue(np.isnan(clinical_remission(7,20,40,"kg",*args,dose_known=False)["remission"]))

    def test_clinical_date_gate(self):
        a = clinical_remission(7,20,40,"kg","2024-01-01","2024-01-08","2024-01-01","2024-01-01")
        b = clinical_remission(7,20,40,"kg","2024-01-01","2024-01-09","2024-01-01","2024-01-01")
        self.assertEqual(a["remission"],1)
        self.assertTrue(np.isnan(b["remission"]))

    def test_temporal_order_and_ambiguity(self):
        self.assertTrue(later_interval("2024-01-01",["2024-04-01"],"2024-04-01")["eligible"])
        self.assertFalse(later_interval("2024-03-15",["2024-04-01"],"2024-04-01")["eligible"])
        self.assertFalse(later_interval("2024-01-01",["2024-04-01","2024-04-02"],"2024-04-01")["eligible"])
        self.assertFalse(later_interval("2024-01-01",[],"2024-04-01")["eligible"])
        self.assertEqual(overlapping_pairs([(0,2),(1,3),(3,5)]),{0,1})

    def test_event_checker_gaps_and_recovery(self):
        x = [50]*4+[100]*4+[50]*4+[np.nan]+[50]*4+[100]*3+[50]*4+[100]*4+[50]*4
        for got, expected in zip(independent_episode_labeler(x), episodes(x)):
            np.testing.assert_array_equal(got,expected)
        self.assertEqual(len(independent_episode_labeler(x)[0]),2)

    def test_event_checker_random_prefixes(self):
        rng=np.random.default_rng(18)
        # Fixtures stress state transitions, not clinical efficacy.
        for n in (0,1,4,35,500):
            x=rng.choice([np.nan,40,53,54,70],n)
            for a,b in zip(independent_episode_labeler(x),episodes(x)):
                np.testing.assert_array_equal(a,b)

    def test_counterexample(self):
        a=residual_identity(10,9,100); b=residual_identity(20,19,100)
        self.assertEqual(a["unwarned_rate"],b["unwarned_rate"])
        a=residual_identity(10,1,100); b=residual_identity(20,18,100)
        self.assertGreater(a["unwarned_rate"],b["unwarned_rate"])

    def test_failed_gate_cannot_support(self):
        f=pd.DataFrame(dict(PtID=np.arange(200),unwarned=np.ones(200)))
        transfer=dict(overall_sensitivity=.29, delivery_arms={"I":{"unconfirmed_per_day":0,"alarms_per_day":1},"S":{"unconfirmed_per_day":0,"alarms_per_day":1}})
        got=gate_assessment(f,transfer,{"status":"estimated","numeric_H1_condition":True},True)
        self.assertEqual(got["status"],"Inconclusive")

    def test_boundary_replay_preserves_event(self):
        t=np.arange(-9000,6000,300); x=np.full(len(t),100.)
        e=np.flatnonzero(t==0)[0]; x[e:e+4]=50
        p=np.full(len(t),np.nan); p[e-4]=1
        out,events,alarms=replay_interval(t,x,p,.5,np.isfinite(p),0,3000)
        self.assertEqual(out["events"],1)
        self.assertEqual(out["matched"],1)
        self.assertEqual(out["alarms"],0)  # Useful preceding alarm is outside burden window.
        self.assertTrue(events.iloc[0].warned)

    def test_delay_changes_alarm_information_not_event_clock(self):
        t=np.arange(-9000,6000,300);x=np.full(len(t),100.)
        e=np.flatnonzero(t==0)[0];x[e:e+4]=50
        delayed=np.r_[np.nan,x[:-1]]
        p=np.full(len(t),np.nan);p[e]=1
        regular,_,_=replay_interval(t,x,p,.5,np.isfinite(p),0,3000)
        late,_,alarms=replay_interval(t,x,p,.5,np.isfinite(p),0,3000,warning_glucose=delayed)
        self.assertEqual(regular["events"],late["events"])
        self.assertEqual(regular["alarms"],0)
        self.assertEqual(late["alarms"],1)
        self.assertEqual(alarms.iloc[0].category,"late")

    def test_fold_assignment_id_permutation(self):
        ids=["01",2,3,4,5,6,7,8,9,10]
        self.assertEqual(prescribed_folds(ids),prescribed_folds(ids[::-1]))
        self.assertEqual(set(prescribed_folds(ids).values()),set(range(5)))

    def policy_rows(self):
        return [dict(kind=k,threshold=-40 if k=="linear" else .5,
                     short_threshold=.5 if k=="switch" else 0,
                     score=dict(sensitivity=.5,unconfirmed_per_day=.1,total_per_day=.2))
                for k in ["M120","M15","missing","switch","linear"]]

    def test_card003_does_not_inherit_switch_margin(self):
        rows=self.policy_rows()
        switch=next(r for r in rows if r["kind"]=="switch")
        switch["score"]=dict(sensitivity=.6,unconfirmed_per_day=.95,total_per_day=1.1)
        self.assertEqual(select_primary(rows)["primary"]["kind"],"switch")

    def test_policy_simplicity_and_physical_threshold_ties(self):
        rows=self.policy_rows()
        high=dict(rows[-1]); high["threshold"]=-80
        rows.append(high)
        selected=select_primary(rows)["primary"]
        self.assertEqual(selected["kind"],"linear")
        self.assertEqual(selected["threshold"],-80)

    def known_mean_fixture(self):
        # Exactly representable Poisson means, solely an estimator fixture.
        rows=[]; rng=np.random.default_rng(40)
        for child in range(60):
            for visit in [26,39,52]:
                e=int(rng.integers(0,3))
                rows.append(dict(PtID=child,delivery_arm="standard" if child%3==0 else "intensive",
                    auc_nmol_l=.2*e,age=8+child%9,target_visit=visit,
                    management=["usual care","intensive/Tandem","intensive/Medtronic"][child%3],
                    cohort=["A","B"][child%2],lag_days=55+(child+visit)%11,
                    unwarned=2**(5-e),observed_days=1))
        return pd.DataFrame(rows)

    def test_known_mean_estimator_and_fixed_intercepts(self):
        f=self.known_mean_fixture()
        for fixed in [False,True]:
            r,b=fit_association(f,n_boot=3,within_child=fixed)
            self.assertAlmostEqual(r["rr"],.5,places=7)
            self.assertEqual(r["status"],"inconclusive")  # 3 is not 5,000.
            self.assertEqual(r["failed_bootstrap"],0)

    def test_only_constant_companion_lag_may_be_omitted(self):
        f=self.known_mean_fixture()
        with self.assertRaises(ValueError):fixed_design(f,omit_constant_lag=True)
        f["lag_days"]=-28
        x=fixed_design(f,omit_constant_lag=True)
        self.assertNotIn("lag_centered_30d",x.columns)
        self.assertEqual(np.linalg.matrix_rank(x),len(x.columns))

    def test_actual_patient_relabeling_preserves_association(self):
        f=self.known_mean_fixture()
        renamed=f.copy()
        ids=list(f.PtID.unique())
        renamed["PtID"]=renamed.PtID.map({p:f"arbitrary_{len(ids)-i}" for i,p in enumerate(ids)})
        a,_=fit_association(f,n_boot=0)
        b,_=fit_association(renamed,n_boot=0)
        self.assertAlmostEqual(a["rr"],b["rr"],places=12)
        self.assertAlmostEqual(a["robust_se"],b["robust_se"],places=12)

    def test_secondary_constant_nuisance_does_not_drop_exposure(self):
        f=self.known_mean_fixture()
        subset=f[f.management.eq("usual care")]
        result,_=fit_association(subset,n_boot=0,drop_constant_nuisance=True)
        self.assertAlmostEqual(result["rr"],.5,places=7)
        self.assertEqual(result["dropped_constant_nuisance_columns"],["intensive_tandem","intensive_medtronic"])

    def test_refractory_context_requires_known_quiet_run(self):
        from card003_boundary import synchronized
        t=np.arange(-10800,3600,300);g=np.ones(len(t))*100;p=np.ones(len(t))
        self.assertFalse(synchronized(t,g,p,.5,0,24))
        p[24:30]=0
        self.assertTrue(synchronized(t,g,p,.5,0,24))
        self.assertFalse(synchronized(t,g,p,.5,0,27))

    def test_empty_and_allzero_are_not_estimates(self):
        f=self.known_mean_fixture()
        for ff in [f.iloc[0:0],f.assign(unwarned=0)]:
            r,b=fit_association(ff,n_boot=0)
            self.assertEqual(r["status"],"inconclusive")
            self.assertNotIn("rr",r)

    def test_invalid_outcomes_and_exposure_are_rejected(self):
        f=self.known_mean_fixture()
        for ff in [f.assign(unwarned=np.nan),f.assign(unwarned=.5),f.assign(observed_days=0)]:
            r,b=fit_association(ff,n_boot=0)
            self.assertEqual(r["status"],"inconclusive")
            self.assertEqual(r["reason"],"invalid_count_or_observed_exposure")

    def test_ppml_separation_despite_full_column_rank(self):
        X=np.array([[1,0],[1,1],[1,1]])
        self.assertTrue(poisson_separation(X,np.array([1,0,0])))
        self.assertFalse(poisson_separation(X,np.array([1,0,1])))

if __name__=="__main__": unittest.main()
