"""Independent artifact accounting and participant-distribution summaries."""
import argparse,json,hashlib
from pathlib import Path
import numpy as np
import pandas as pd


def quantiles(values):
    values=np.asarray(values,dtype=float); values=values[np.isfinite(values)]
    return dict(n=len(values),minimum=float(values.min()) if len(values) else None,
                q25=float(np.quantile(values,.25)) if len(values) else None,
                median=float(np.median(values)) if len(values) else None,
                q75=float(np.quantile(values,.75)) if len(values) else None,
                maximum=float(values.max()) if len(values) else None)


def audit(root):
    out=root/"results"; distributions=[]; categories=[]; audit_rows=[]
    for cohort in ["internal","external"]:
        frame=pd.read_csv(out/f"{cohort}_participant_scores.csv")
        ledger=pd.read_csv(out/f"{cohort}_alarm_ledger.csv")
        events=pd.read_csv(out/f"{cohort}_event_table.csv")
        summary=pd.read_csv(out/f"{cohort}_summary.csv").set_index("policy")
        for policy,part in frame.groupby("policy",sort=False):
            if (part.matched>part.opportunities).any() or (part.opportunities>part.events).any(): raise AssertionError("Participant coverage bound violated")
            if not np.array_equal(part.alarms,part.matched+part.unconfirmed): raise AssertionError("Alarm accounting differs")
            p=ledger.loc[ledger.policy==policy]
            if len(p)!=part.alarms.sum() or (p.category=="matched").sum()!=part.matched.sum(): raise AssertionError("Full alarm ledger disagrees with totals")
            if events[policy+"_warning"].sum()!=part.matched.sum() or len(events)!=part.events.sum(): raise AssertionError("Independent event table count mismatch")
            np.testing.assert_allclose(part[["events","matched","alarms","unconfirmed","opportunities","observed_days"]].sum(),summary.loc[policy,["events","matched","alarms","unconfirmed","opportunities","observed_days"]],rtol=0,atol=1e-8)
            for measure,values in [("sensitivity",part.matched/part.events.replace(0,np.nan)),("opportunity_coverage",part.opportunities/part.events.replace(0,np.nan)),("unconfirmed_per_day",part.unconfirmed/part.observed_days.replace(0,np.nan)),("total_per_day",part.alarms/part.observed_days.replace(0,np.nan))]:
                distributions.append(dict(cohort=cohort,policy=policy,measure=measure,participants=len(part),no_observed_data=int((part.observed_days==0).sum()),zero_event_participants=int((part.events==0).sum()),**quantiles(values)))
            for category in ["matched","duplicate","late","unobservable_outcome","unmatched_observed"]:
                categories.append(dict(cohort=cohort,policy=policy,category=category,count=int((p.category==category).sum())))
            audit_rows.append(dict(cohort=cohort,policy=policy,participant_coverage_bound="PASS",alarm_accounting="PASS",event_table_identity="PASS",sum_to_summary="PASS"))
        if frame.groupby("id").events.nunique().max()!=1 or frame.groupby("id").observed_days.nunique().max()!=1: raise AssertionError("Policies changed the common denominator")
        sensitivity=out/f"{cohort}_sensitivity_participants.csv"
        if sensitivity.exists():
            attacks=pd.read_csv(sensitivity)
            base=attacks.loc[attacks.variant=="primary"]
            columns=["id","policy","events","matched","opportunities","alarms","unconfirmed","observed_days"]
            pd.testing.assert_frame_equal(frame[columns].sort_values(["id","policy"]).reset_index(drop=True),base[columns].sort_values(["id","policy"]).reset_index(drop=True),check_dtype=False,check_exact=False,rtol=0,atol=1e-8)
            if len(attacks)!=len(frame)*9 or attacks.groupby(["id","policy"]).variant.nunique().min()!=9: raise AssertionError("Missing sensitivity cell")
            audit_rows.append(dict(cohort=cohort,independent_primary_replay="PASS",all_participants_all_nine_variants="PASS"))
    pd.DataFrame(distributions).to_csv(out/"patient_level_distributions.csv",index=False)
    pd.DataFrame(categories).to_csv(out/"alarm_category_counts.csv",index=False)
    payload=dict(status="PASS",checks=audit_rows,source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),patient_interpretation="Event-free participants have undefined individual sensitivity and remain in burden and participant bootstrap denominators")
    (out/"independent_result_accounting_audit.json").write_text(json.dumps(payload,indent=2))
    print(json.dumps(payload))


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data-root",required=True,type=Path)
    audit(p.parse_args().data_root)
