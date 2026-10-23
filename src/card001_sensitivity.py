"""Frozen-policy grid/clock/label/reporting-delay attacks on observed timelines."""
from __future__ import annotations
import json,pickle,hashlib,time
from concurrent.futures import ProcessPoolExecutor,as_completed
from pathlib import Path
import numpy as np
import pandas as pd
from card001_core import *
from card001_models import predict
from card001_evaluate import score_policy


def variants(root,p):
    z=np.load(p["file"])
    g=z["glucose"]; t=z["time"]
    yield "primary",g,g,t,4
    yield "three_readings",g,g,t,3
    if "glucose_without_clock_barriers" in z:
        raw=z["glucose_without_clock_barriers"]
        yield "no_clock_barriers",raw,raw,t,4
    # Reporting delay changes available INPUT only. Shifted event truth would
    # conceal the lost lead time and is expressly forbidden here.
    delayed=np.r_[np.nan,g[:-1]] if len(g) else g
    yield "reporting_delay_5min",delayed,g,t,4
    source=root/"clocks"/f"{p['cohort']}_{p['id']}.npz"
    if source.exists():
        a=np.load(source)
        for phase in [60,120,180,240,300]:
            rebuilt=causal_grid(a["time"],a["glucose"],p["start"],p["end"]-1e-6,phase_seconds=phase)
            if "barrier_times" in a:
                bins=np.ceil((a["barrier_times"]-phase)/300)*300+phase
                rebuilt["glucose"][np.isin(rebuilt["time"],bins)]=np.nan
            yield f"grid_phase_{phase}sec",rebuilt["glucose"],rebuilt["glucose"],rebuilt["time"],4
    else:
        # Empty raw stream remains empty under every phase; no low is invented.
        if np.isfinite(g).any(): raise ValueError(f"Raw source required for phase attack: {source}")
        for phase in [60,120,180,240,300]:
            tg=np.arange(np.ceil((p["start"]-phase)/300)*300+phase,p["end"],300)
            empty=np.full(len(tg),np.nan)
            yield f"grid_phase_{phase}sec",empty,empty,tg,4


def initialize_worker(root,cohort_root):
    global _CONTEXT
    root=Path(root); frozen=json.loads((root/"calibration.json").read_text())
    models={}
    for kind,digest in frozen["model_sha256"].items():
        path=root/"models"/f"{kind}.pkl"
        if hashlib.sha256(path.read_bytes()).hexdigest()!=digest: raise ValueError("Frozen model changed")
        models[kind]=pickle.load(path.open("rb"))
    policies={**frozen["comparators"],"switch":frozen["switch"]}
    comparator=frozen["comparator"]["kind"]
    _CONTEXT=(Path(cohort_root),models,policies,comparator)


def participant_replay(p):
    cohort_root,models,policies,comparator=_CONTEXT
    rows=[]; paired={}; causality=[]; prediction_cache={}
    for name,inputs,truth,t,run_length in variants(cohort_root,p):
        event,censored=episodes(truth,run_length)
        event2,censored2=episodes_independent(truth,run_length)
        if not np.array_equal(event,event2) or not np.array_equal(censored,censored2): raise AssertionError("Sensitivity detector mismatch")
        key=hashlib.sha256(np.asarray(inputs,dtype=np.float64).tobytes()).hexdigest()
        if key not in prediction_cache:
            pred={kind:predict(m,inputs,kind) for kind,m in models.items()}
            pred["linear"]=-linear_forecast(inputs)
            prediction_cache[key]=pred
        pred=prediction_cache[key]
        q=dict(glucose=inputs,scoring_glucose=truth,onsets=event,predictions=pred)
        results={}
        for kind,policy in policies.items():
            result=score_policy(q,policy); results[kind]=result
            rows.append(dict(id=p["id"],arm=p["arm"],variant=name,policy=kind,events=result["events"],opportunities=result["opportunities"],matched=result["matched"],alarms=result["alarms"],unconfirmed=result["unconfirmed"],observed_days=result["observed_days"],censored_low_runs=len(censored)))
        a=results["switch"]; b=results[comparator]
        paired[name]=[a["events"],a["matched"],b["matched"],a["unconfirmed"],b["unconfirmed"],a["observed_days"]]
        if name=="primary":
            cut=len(inputs)//2; attack=inputs.copy(); attack[cut:]=45.
            for kind,model in models.items():
                changed=predict(model,attack,kind)
                np.testing.assert_allclose(changed[:cut],pred[kind][:cut],equal_nan=True,rtol=0,atol=0)
            causality.append(dict(id=p["id"],prefix_bins=cut,frozen_model_future_change_test="PASS"))
    return dict(rows=rows,paired=paired,causality=causality)


def run_sensitivity(root,external_root=None,workers=1):
    start=time.perf_counter(); frozen=json.loads((root/"calibration.json").read_text())
    for kind,digest in frozen["model_sha256"].items():
        if hashlib.sha256((root/"models"/f"{kind}.pkl").read_bytes()).hexdigest()!=digest: raise ValueError("Frozen model changed before sensitivity aggregation")
    cohort_root=Path(external_root) if external_root else root
    split="external_test" if external_root else "test"
    prefix="external" if external_root else "internal"
    if not (root/"results"/f"{prefix}_comparison.json").exists(): raise ValueError("Primary locked comparison must precede sensitivity replay")
    manifest=json.loads((cohort_root/"participant_manifest.json").read_text())
    participants=[p for p in manifest["participants"] if p["split"]==split]
    if not participants or len({str(p["id"]) for p in participants})!=len(participants): raise ValueError("Sensitivity requires unique prespecified participants")
    expected_variants={"primary","three_readings","no_clock_barriers","reporting_delay_5min",*[f"grid_phase_{phase}sec" for phase in [60,120,180,240,300]]}
    expected_policies={"M120","M15","linear","missing","switch"}
    def validate_result(p,result):
        if set(result["paired"])!=expected_variants or len(result["rows"])!=45: raise AssertionError("Sensitivity requires every prespecified variant")
        for name in expected_variants:
            rows=[row for row in result["rows"] if row["variant"]==name]
            if len(rows)!=5 or {row["policy"] for row in rows}!=expected_policies or any(str(row["id"])!=str(p["id"]) for row in rows): raise AssertionError("Sensitivity requires all five policies exactly once per participant/variant")
    signature=dict(calibration=hashlib.sha256((root/"calibration.json").read_bytes()).hexdigest(),manifest=hashlib.sha256((cohort_root/"participant_manifest.json").read_bytes()).hexdigest(),source={name:hashlib.sha256((Path(__file__).parent/name).read_bytes()).hexdigest() for name in ["card001_sensitivity.py","card001_evaluate.py","card001_models.py","card001_core.py"]})
    key=hashlib.sha256(json.dumps(signature,sort_keys=True).encode()).hexdigest()
    checkpoint=root/"results"/"sensitivity_checkpoints"/prefix/key[:16]; checkpoint.mkdir(parents=True,exist_ok=True)
    (checkpoint/"commitment.json").write_text(json.dumps(signature,indent=2))
    for name in signature["source"]: (checkpoint/name).write_bytes((Path(__file__).parent/name).read_bytes())
    completed={}; pending=[]; hashes={}
    for p in participants:
        pid=str(p["id"]); source=cohort_root/"clocks"/f"{p['cohort']}_{p['id']}.npz"
        inputs=dict(grid=hashlib.sha256(Path(p["file"]).read_bytes()).hexdigest(),raw_clock=hashlib.sha256(source.read_bytes()).hexdigest() if source.exists() else None)
        hashes[pid]=inputs; cache=checkpoint/f"participant_{pid}.json"
        if cache.exists():
            data=json.loads(cache.read_text())
            if data["inputs"]!=inputs: raise ValueError("Sensitivity checkpoint input changed")
            validate_result(p,data["result"])
            completed[pid]=data["result"]
        else: pending.append(p)
    def retain(p,result):
        validate_result(p,result)
        pid=str(p["id"]); temporary=checkpoint/f"participant_{pid}.tmp"
        temporary.write_text(json.dumps(dict(inputs=hashes[pid],result=result)))
        temporary.replace(checkpoint/f"participant_{pid}.json"); completed[pid]=result
        print(f"{prefix}_sensitivity_participant={pid} completed={len(completed)}/{len(participants)} sec={time.perf_counter()-start:.1f}",flush=True)
    if workers==1:
        initialize_worker(root,cohort_root)
        for p in pending: retain(p,participant_replay(p))
    else:
        with ProcessPoolExecutor(max_workers=workers,initializer=initialize_worker,initargs=(root,cohort_root)) as pool:
            tasks={pool.submit(participant_replay,p):p for p in pending}
            for task in as_completed(tasks): retain(tasks[task],task.result())
    rows=[]; paired={}; causality=[]
    for p in participants: # Original manifest order preserves bootstrap draw identity.
        result=completed[str(p["id"])]; rows.extend(result["rows"]); causality.extend(result["causality"])
        for name,tallies in result["paired"].items(): paired.setdefault(name,[]).append(tallies)
    if len(completed)!=len(participants) or any(len(x)!=len(participants) for x in paired.values()): raise AssertionError("Incomplete sensitivity participants")
    comparator=frozen["comparator"]["kind"]
    frame=pd.DataFrame(rows); out=root/"results"
    frame.to_csv(out/f"{prefix}_sensitivity_participants.csv",index=False)
    summary=[]
    for (variant,policy),part in frame.groupby(["variant","policy"]):
        s=part[["events","matched","opportunities","alarms","unconfirmed","observed_days"]].sum()
        summary.append(dict(variant=variant,policy=policy,**s.to_dict(),sensitivity=s.matched/s.events if s.events else None,unconfirmed_per_day=s.unconfirmed/s.observed_days if s.observed_days else None,total_per_day=s.alarms/s.observed_days if s.observed_days else None))
    pd.DataFrame(summary).to_csv(out/f"{prefix}_sensitivity_summary.csv",index=False)
    result=dict(cohort=prefix,comparator=comparator,paired_intervals={k:paired_bootstrap(v) for k,v in paired.items()},causality=causality,seconds=time.perf_counter()-start,workers=workers,checkpoint=str(checkpoint),execution_signature=signature,delay_semantics="Inputs delayed one bin, outcome truth and observed-day denominator unchanged",outcome_scope="Observed-only sensor glucose; phase and recovery sensitivity may alter identified event count")
    (out/f"{prefix}_sensitivity_audit.json").write_text(json.dumps(result,indent=2)); return result
