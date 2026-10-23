"""Full threshold-grid and locked-test replay, retaining original denominators."""
from __future__ import annotations
import json,pickle,hashlib,time,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from numba import njit,prange
from card001_core import *
from card001_models import ranking,admissible,predict


@njit(parallel=True,cache=False)
def fast_grid(g,onsets,scores,kinds,cuts,short_cuts):
    """Independent compiled counter checked against reference ledger scorer."""
    result=np.zeros((len(kinds),2),dtype=np.int64)
    for c in prange(len(kinds)):
        kind=kinds[c]; i=0; j=0; matched=0; alarms=0
        while i<len(g):
            if not np.isfinite(g[i]) or g[i]<54:
                i+=1; continue
            if kind==4:
                if np.isfinite(scores[0,i]): p=scores[0,i]; cut=cuts[c]
                else: p=scores[1,i]; cut=short_cuts[c]
            else: p=scores[kind,i]; cut=cuts[c]
            if np.isfinite(p) and p>=cut:
                alarms+=1
                while j<len(onsets) and onsets[j]<i+3: j+=1
                if j<len(onsets) and onsets[j]<=i+6:
                    matched+=1; j+=1
                i+=6
            else: i+=1
        result[c,0]=matched; result[c,1]=alarms
    return result


def specifications():
    rows=[]; cuts=np.r_[np.arange(101)/100,np.inf]
    for kind in ["M120","M15","linear","missing"]:
        for threshold in (np.r_[-np.arange(40,101),np.inf] if kind=="linear" else cuts):
            rows.append(dict(kind=kind,threshold=float(threshold)))
    for long in cuts:
        for short in cuts: rows.append(dict(kind="switch",threshold=float(long),short_threshold=float(short)))
    return rows


def burden_intervals(rows,n_boot=10000):
    data=np.asarray(rows,float); rng=np.random.default_rng(SEED)
    estimates=[]; undefined=0
    for _ in range(n_boot):
        s=data[rng.integers(0,len(data),len(data))].sum(axis=0)
        if s[5]<=0: undefined+=1; continue
        estimates.append([s[3]/s[5],s[4]/s[5],(s[3]+s[1])/s[5],(s[4]+s[2])/s[5]])
    arr=np.asarray(estimates)
    names=["switch_unconfirmed","comparator_unconfirmed","switch_total","comparator_total"]
    return {"undefined_replicates":undefined,**{name:np.quantile(arr[:,j],[.025,.975]).tolist() if len(arr) else None for j,name in enumerate(names)}}


def score_policy(p,policy):
    probs=p["predictions"]; g=p["glucose"]; kind=policy["kind"]
    if kind=="switch":
        long=np.isfinite(probs["M120"])
        score=np.where(long,probs["M120"],probs["M15"])
        cut=np.where(long,policy["threshold"],policy["short_threshold"])
    else: score=probs[kind]; cut=policy["threshold"]
    availability=np.isfinite(score)&np.isfinite(g)&(g>=54)
    alarm=replay(g,score,cut)
    return participant_score(p.get("scoring_glucose",g),p["onsets"],alarm,availability)


def calibrate_run(root):
    from card001_run import load_participants
    frozen_path=root/"calibration.json"
    if frozen_path.exists():
        previous=json.loads(frozen_path.read_text())
        if previous.get("status")!="FROZEN_BEFORE_TEST": raise ValueError("Existing calibration is not a valid frozen artifact")
        for kind,digest in previous["model_sha256"].items():
            if hashlib.sha256((root/"models"/f"{kind}.pkl").read_bytes()).hexdigest()!=digest: raise ValueError("Frozen calibration model identity changed")
        if hashlib.sha256((root/"validation_prediction_manifest.json").read_bytes()).hexdigest()!=previous["prediction_manifest_sha256"]: raise ValueError("Frozen validation prediction identity changed")
        print("Existing immutable calibration retained; no re-selection performed",flush=True)
        return previous
    if list((root/"results").glob("*_comparison.json")):
        raise ValueError("Cannot first calibrate after a test comparison exists")
    provenance=json.loads((root/"validation_prediction_manifest.json").read_text())
    for kind,digest in provenance["model_sha256"].items():
        if hashlib.sha256((root/"models"/f"{kind}.pkl").read_bytes()).hexdigest()!=digest: raise ValueError("Validation predictions belong to different fitted models")
    for section in ["inputs","predictions"]:
        for name,digest in provenance[section].items():
            if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest: raise ValueError(f"Validation {section} bytes changed: {name}")
    for name,digest in provenance["source_sha256"].items():
        if hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()!=digest: raise ValueError("Prediction source changed after validation replay")
    start=time.perf_counter(); rows=specifications()
    kinds=np.array([["M120","M15","linear","missing","switch"].index(r["kind"]) for r in rows])
    cuts=np.array([r["threshold"] for r in rows]); short=np.array([r.get("short_threshold",np.nan) for r in rows])
    counts=np.zeros((len(rows),2),dtype=np.int64); events=0; days=0; patients=[]
    for p in load_participants(root,["validation"]):
        z=np.load(root/"validation_predictions"/f"DCLP3_{p['id']}.npz")
        p["predictions"]={k:z[k] for k in ["M120","M15","linear","missing"]}; p["onsets"]=z["onsets"]
        stack=np.vstack([p["predictions"][k] for k in ["M120","M15","linear","missing"]])
        tallies=fast_grid(p["glucose"],p["onsets"],stack,kinds,cuts,short)
        # Cross-check threshold corners and interior points against independent,
        # complete alarm-ledger replay on every real validation participant.
        for idx in [0,10,50,100,101,120,230,350,368,1000,5000,len(rows)-1]:
            q=score_policy(p,rows[idx])
            if [q["matched"],q["alarms"]]!=tallies[idx].tolist(): raise AssertionError("Compiled and reference replay disagree")
        counts+=tallies; events+=len(p["onsets"]); days+=np.isfinite(p["glucose"]).sum()/288
        patients.append(p)
        print(f"calibration_participant={p['id']} elapsed={time.perf_counter()-start:.1f}s",flush=True)
    if events==0 or days==0: raise ValueError("No identifiable validation event/day denominator")
    for row,tally in zip(rows,counts):
        matched,alarms=map(int,tally)
        row["score"]=dict(events=events,matched=matched,alarms=alarms,observed_days=days,sensitivity=matched/events,unconfirmed_per_day=(alarms-matched)/days,total_per_day=alarms/days)
    best={kind:min([r for r in rows if r["kind"]==kind and admissible(r["score"])],key=ranking) for kind in ["M120","M15","linear","missing"]}
    comparator=min(best.values(),key=ranking)
    margin=comparator["score"]["unconfirmed_per_day"]+.10
    switch=min([r for r in rows if r["kind"]=="switch" and admissible(r["score"],margin)],key=ranking)
    result=dict(comparator=comparator,comparators=best,switch=switch,validation_grid=rows,status="FROZEN_BEFORE_TEST",model_sha256={k:hashlib.sha256((root/"models"/f"{k}.pkl").read_bytes()).hexdigest() for k in ["M120","M15","missing"]},manifest_sha256=hashlib.sha256((root/"participant_manifest.json").read_bytes()).hexdigest(),prediction_manifest_sha256=hashlib.sha256((root/"validation_prediction_manifest.json").read_bytes()).hexdigest(),seconds=time.perf_counter()-start)
    # Development-only paired precision approximation. The absolute gain is not
    # assumed true; report dispersion and the 80%-power normal MDE as descriptive.
    bootstrap=[]
    for p in patients:
        a=score_policy(p,switch); b=score_policy(p,comparator)
        bootstrap.append([a["events"],a["matched"],b["matched"],a["unconfirmed"],b["unconfirmed"],a["observed_days"]])
    precision=paired_bootstrap(bootstrap)
    ci=precision["sensitivity_difference_ci"]
    precision["approximate_80pct_MDE_validation_N"]=(ci[1]-ci[0])/3.92*2.8016 if ci else None
    source_manifest=json.loads((root/"participant_manifest.json").read_text())
    internal_n=sum(p["split"]=="test" for p in source_manifest["participants"])
    precision["validation_participants"]=len(patients)
    precision["internal_allocated_participants"]=internal_n
    precision["approximate_80pct_MDE_internal_N"]=precision["approximate_80pct_MDE_validation_N"]*np.sqrt(len(patients)/internal_n) if ci and internal_n else None
    precision["limitation"]="Development-only normal approximation from bootstrap interval width; not a power guarantee; population shift and alternative effect shape are unidentified"
    result["development_precision"]=precision
    result["frozen_at_utc"]=datetime.now(timezone.utc).isoformat()
    snapshot=root/"calibration_source_snapshot"; snapshot.mkdir(exist_ok=True)
    sources=[Path(__file__).parent/name for name in ["card001_core.py","card001_models.py","card001_evaluate.py","card001_run.py"]]
    result["calibration_source_sha256"]={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    for source in sources: shutil.copy2(source,snapshot/source.name)
    with frozen_path.open("x") as f: json.dump(result,f,indent=2)
    print(json.dumps({k:result[k] for k in ["comparator","switch","development_precision","seconds"]}),flush=True)
    return result


def test_run(root,external_root=None,split_override=None,prefix_override=None):
    from card001_run import load_participants
    executing_source={p.name:p.read_bytes() for p in [Path(__file__).parent/name for name in ["card001_evaluate.py","card001_core.py","card001_models.py","card001_run.py"]]}
    frozen=json.loads((root/"calibration.json").read_text())
    for name in ["card001_core.py","card001_models.py"]:
        if name in frozen.get("calibration_source_sha256",{}) and hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest()!=frozen["calibration_source_sha256"][name]: raise ValueError("Frozen prediction or scoring source changed")
    models={}
    for kind,digest in frozen["model_sha256"].items():
        path=root/"models"/f"{kind}.pkl"
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest: raise ValueError("Frozen model hash mismatch")
        models[kind]=pickle.load(path.open("rb"))
    if hashlib.sha256((root/"participant_manifest.json").read_bytes()).hexdigest()!=frozen["manifest_sha256"]: raise ValueError("Participant manifest changed after freeze")
    output=root/"results"; output.mkdir(exist_ok=True)
    cohort_root=Path(external_root) if external_root is not None else root
    split="external_test" if external_root is not None else "test"
    prefix="external" if external_root is not None else "internal"
    if split_override is not None: split=split_override
    if prefix_override is not None:
        if not prefix_override.replace("_","").isalnum(): raise ValueError("Result prefix must be alphanumeric with underscores")
        prefix=prefix_override
    participant_rows=[]; alarms=[]; event_rows=[]; censored_rows=[]; comparisons=[]; input_hashes={}
    policies={**frozen["comparators"],"switch":frozen["switch"]}
    comparator=frozen["comparator"]["kind"]
    for p in load_participants(cohort_root,[split]):
        input_hashes[p["file"]]=hashlib.sha256(Path(p["file"]).read_bytes()).hexdigest()
        g=p["glucose"]; p["onsets"],unknown=episodes(g)
        independent,unknown2=episodes_independent(g)
        if not np.array_equal(p["onsets"],independent) or not np.array_equal(unknown,unknown2): raise AssertionError("Independent detector disagreement")
        for j,e in enumerate(unknown):
            censored_rows.append(dict(id=p["id"],censored_run_index=j,first_observed_low_time=float(p["time"][e]),onset_unknown=True,reason="sustained_low_seen_before_recovery_established_after_gap_or_start"))
        p["predictions"]={kind:predict(m,g,kind) for kind,m in models.items()}; p["predictions"]["linear"]=-linear_forecast(g)
        per={}
        for name,policy in policies.items():
            s=score_policy(p,policy); per[name]=s
            history_ready=ready(g,25) if name=="M120" else ready(g,4)
            ready_mask=history_ready&(g>=54)
            participant_rows.append(dict(id=p["id"],arm=p["arm"],policy=name,events=s["events"],opportunities=s["opportunities"],matched=s["matched"],alarms=s["alarms"],unconfirmed=s["unconfirmed"],observed_days=s["observed_days"],calendar_days=len(g)/288,model_ready_hours=float(history_ready.sum()/12),alarm_eligible_prediction_hours=float(ready_mask.sum()/12),censored_low_runs=len(unknown)))
            for alarm in s["ledger"]:
                alarms.append(dict(id=p["id"],policy=name,time=float(p["time"][alarm["alarm_index"]]),**alarm))
        a=per["switch"]; b=per[comparator]
        comparisons.append([a["events"],a["matched"],b["matched"],a["unconfirmed"],b["unconfirmed"],a["observed_days"]])
        r120=ready(g,25)&(g>=54)
        for j,e in enumerate(p["onsets"]):
            event_rows.append(dict(id=p["id"],arm=p["arm"],event_index=j,onset_time=float(p["time"][e]),M120_opportunity=bool(r120[max(0,e-6):max(0,e-2)].any()),switch_warning=j in a["matched_event_indices"],comparator_warning=j in b["matched_event_indices"],**{name+"_warning":j in scores["matched_event_indices"] for name,scores in per.items() if name!="switch"}))
    if not participant_rows: raise ValueError("No prespecified test participants")
    frame=pd.DataFrame(participant_rows); frame.to_csv(output/f"{prefix}_participant_scores.csv",index=False)
    pd.DataFrame(alarms).to_csv(output/f"{prefix}_alarm_ledger.csv",index=False)
    pd.DataFrame(censored_rows,columns=["id","censored_run_index","first_observed_low_time","onset_unknown","reason"]).to_csv(output/f"{prefix}_censored_event_table.csv",index=False)
    event_frame=pd.DataFrame(event_rows); event_frame.to_csv(output/f"{prefix}_event_table.csv",index=False)
    if len(event_frame):
        strata=[]
        for group,part in event_frame.groupby("M120_opportunity"):
            strata.append(dict(M120_opportunity=bool(group),events=len(part),**{name:float(part[name].mean()) for name in part.columns if name.endswith("_warning")}))
        pd.DataFrame(strata).to_csv(output/f"{prefix}_gap_exposed_ablation.csv",index=False)
    summary=[]
    for policy,part in frame.groupby("policy"):
        total=part[["events","opportunities","matched","alarms","unconfirmed","observed_days"]].sum()
        summary.append(dict(policy=policy,participants=len(part),participants_with_observed_data=int((part.observed_days>0).sum()),**total.to_dict(),sensitivity=total.matched/total.events if total.events else None,opportunity_coverage=total.opportunities/total.events if total.events else None,conditional_sensitivity=total.matched/total.opportunities if total.opportunities else None,unconfirmed_per_day=total.unconfirmed/total.observed_days if total.observed_days else None,total_per_day=total.alarms/total.observed_days if total.observed_days else None))
    pd.DataFrame(summary).to_csv(output/f"{prefix}_summary.csv",index=False)
    arm_summary=[]
    for (arm,policy),part in frame.groupby(["arm","policy"]):
        s=part[["events","opportunities","matched","alarms","unconfirmed","observed_days"]].sum()
        arm_summary.append(dict(arm=arm,policy=policy,participants=len(part),**s.to_dict(),sensitivity=s.matched/s.events if s.events else None,unconfirmed_per_day=s.unconfirmed/s.observed_days if s.observed_days else None))
    pd.DataFrame(arm_summary).to_csv(output/f"{prefix}_arm_strata.csv",index=False)
    snapshot=output/f"{prefix}_source_snapshot"; snapshot.mkdir(exist_ok=True)
    for name,data in executing_source.items(): (snapshot/name).write_bytes(data)
    provenance=dict(evaluated_at_utc=datetime.now(timezone.utc).isoformat(),calibration_sha256=hashlib.sha256((root/"calibration.json").read_bytes()).hexdigest(),model_sha256=frozen["model_sha256"],inputs=input_hashes,source_sha256={name:hashlib.sha256(data).hexdigest() for name,data in executing_source.items()},source_snapshot=str(snapshot))
    (output/f"{prefix}_evaluation_provenance.json").write_text(json.dumps(provenance,indent=2))
    return finalize_saved_results(root,prefix,cohort_root)


def finalize_saved_results(root,prefix,cohort_root):
    """Summarize retained complete replay ledgers without refitting or re-selection."""
    output=root/"results"; frozen=json.loads((root/"calibration.json").read_text())
    provenance=json.loads((output/f"{prefix}_evaluation_provenance.json").read_text())
    if provenance["calibration_sha256"]!=hashlib.sha256((root/"calibration.json").read_bytes()).hexdigest(): raise ValueError("Saved replay calibration changed")
    frame=pd.read_csv(output/f"{prefix}_participant_scores.csv")
    summary=pd.read_csv(output/f"{prefix}_summary.csv").to_dict("records")
    comparator=frozen["comparator"]["kind"]
    a=frame.loc[frame.policy=="switch"].set_index("id")
    b=frame.loc[frame.policy==comparator].set_index("id").reindex(a.index)
    if not np.array_equal(a.events,b.events) or not np.array_equal(a.observed_days,b.observed_days): raise AssertionError("Saved paired denominators disagree")
    comparisons=np.column_stack([a.events,a.matched,b.matched,a.unconfirmed,b.unconfirmed,a.observed_days])
    ci=paired_bootstrap(comparisons)
    burden_ci=burden_intervals(comparisons)
    by={r["policy"]:r for r in summary}; a=by["switch"]; b=by[comparator]
    difference=a["sensitivity"]-b["sensitivity"] if a["sensitivity"] is not None else None
    margin=a["unconfirmed_per_day"]-b["unconfirmed_per_day"] if a["unconfirmed_per_day"] is not None else None
    floor=a["participants_with_observed_data"]>=20 and a["events"]>=100
    caps=all(r["unconfirmed_per_day"] is not None and r["unconfirmed_per_day"]<=1 and r["total_per_day"]<=2 for r in [a,b])
    dev_mde=frozen["development_precision"]["approximate_80pct_MDE_validation_N"]
    projected_mde=dev_mde*np.sqrt(frozen["development_precision"]["validation_participants"]/a["participants"]) if dev_mde is not None else None
    precision_resolvable=projected_mde is not None and projected_mde<=.05
    margin_refuted=floor and ci["unconfirmed_difference_ci"] is not None and ci["unconfirmed_difference_ci"][0]>.10
    cap_refuted=floor and any(burden_ci[name] is not None and burden_ci[name][0]>(1 if "unconfirmed" in name else 2) for name in ["switch_unconfirmed","comparator_unconfirmed","switch_total","comparator_total"])
    target_refuted=floor and ci["sensitivity_difference_ci"] is not None and ci["sensitivity_difference_ci"][1]<.05
    if not floor or ci["sensitivity_difference_ci"] is None: decision="INCONCLUSIVE_INSUFFICIENT_FEASIBILITY"
    elif not precision_resolvable: decision="INCONCLUSIVE_DEVELOPMENT_PRECISION"
    elif ci["sensitivity_difference_ci"][1]<.05: decision="REFUTED_TARGET_GAIN"
    elif margin_refuted: decision="REFUTED_BURDEN_MARGIN"
    elif cap_refuted: decision="REFUTED_BURDEN_CAP"
    elif not caps or margin>.10: decision="REJECTED_OPERATING_POINT"
    elif difference>=.05 and ci["sensitivity_difference_ci"][0]>0: decision="SUPPORTED_IN_THIS_COHORT"
    else: decision="INCONCLUSIVE_IN_THIS_COHORT"
    conclusion=dict(status=decision,cohort=prefix,analysis_role="primary" if prefix in ["internal","external"] else "sensitivity",cohort_manifest_sha256=hashlib.sha256((cohort_root/"participant_manifest.json").read_bytes()).hexdigest(),comparator=comparator,sensitivity_difference=difference,unconfirmed_difference=margin,feasibility_floor=floor,burden_caps=caps,burden_rate_intervals=burden_ci,statistically_refuted_target_gain=target_refuted,statistically_refuted_margin=margin_refuted,statistically_refuted_cap=cap_refuted,development_projected_MDE=projected_mde,development_precision_resolvable=precision_resolvable,operating_point_rejected=not caps or (margin is not None and margin>.10),**ci,caveat="Retrospective recorded sensor episodes only; not a clinical efficacy or dosing result")
    encode=lambda value:value.item() if isinstance(value,np.generic) else (_ for _ in ()).throw(TypeError(type(value).__name__))
    serialized=json.dumps(conclusion,indent=2,default=encode)
    (output/f"{prefix}_comparison.json").write_text(serialized); print(serialized,flush=True)
    finalization=dict(finalized_at_utc=datetime.now(timezone.utc).isoformat(),source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),replay_table_sha256={name:hashlib.sha256((output/f"{prefix}_{name}.csv").read_bytes()).hexdigest() for name in ["participant_scores","summary","alarm_ledger","event_table"]},description="Summarization of retained frozen replay; no model fitting or threshold selection")
    (output/f"{prefix}_finalization_provenance.json").write_text(json.dumps(finalization,indent=2))
    return conclusion
