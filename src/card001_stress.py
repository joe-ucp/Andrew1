"""Frozen development-only artificial-input-mask stress experiment.

These masks alter model inputs, never original glucose outcome truth. This is
an artificial mechanism check, not held-out or external clinical validation.
"""
from __future__ import annotations
import argparse
import hashlib
import json
import os
from pathlib import Path
import pickle
import sys
import time

import numpy as np
import pandas as pd
from card001_core import episodes, episodes_independent, linear_forecast, match_alarms, replay

SCORE_START=30  # 150 minutes of intact lookback
SCORE_BINS=288  # 24 hours of retained outcome truth, not a midnight-aligned day
SCORE_END=SCORE_START+SCORE_BINS
TRACE_BINS=SCORE_END+10  # through 45 minutes beyond the scored boundary
MAX_WINDOWS_PER_PARTICIPANT=7
DURATIONS=[5,15,30,60]
SEEDS=list(range(20261022,20261042))


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def complete_windows(glucose,max_windows=MAX_WINDOWS_PER_PARTICIPANT):
    """Earliest nonoverlapping complete traces, selected without glucose level."""
    observed=np.isfinite(np.asarray(glucose,float))
    edges=np.r_[0,np.flatnonzero(np.diff(observed.astype(int)))+1,len(observed)]
    windows=[]
    for a,b in zip(edges[:-1],edges[1:]):
        if a==b or not observed[a]:
            continue
        for start in range(int(a),int(b)-TRACE_BINS+1,TRACE_BINS):
            windows.append((start,start+TRACE_BINS))
            if len(windows)==max_windows:
                return windows
    return windows


def mask_positions(seed,n_traces):
    """One random start per six-hour quarter, valid for all four durations."""
    rng=np.random.default_rng(seed)
    return np.array([[SCORE_START+q*72+int(rng.integers(0,61)) for q in range(4)] for _ in range(n_traces)],dtype=int)


def mask_trace(truth,starts,minutes):
    if minutes not in DURATIONS:
        raise ValueError("Only the four prespecified durations are allowed")
    masked=np.asarray(truth,float).copy()
    for start in starts:
        if start<SCORE_START or start+minutes//5>SCORE_END:
            raise ValueError("Mask outside prespecified scoring interval")
        masked[start:start+minutes//5]=np.nan
    return masked


def trace_score(truth,inputs,predictions,policy):
    """Use unchanged outcome truth and its events; replay inputs independently."""
    truth=np.asarray(truth,float); inputs=np.asarray(inputs,float)
    if len(truth)!=TRACE_BINS or len(inputs)!=TRACE_BINS or not np.isfinite(truth).all():
        raise ValueError("Stress traces require complete untouched outcome truth")
    onsets,censored=episodes(truth)
    check,check_censored=episodes_independent(truth)
    if not np.array_equal(onsets,check) or not np.array_equal(censored,check_censored):
        raise AssertionError("Independent truth labelers disagree")
    kind=policy["kind"]
    if kind=="switch":
        long=np.isfinite(predictions["M120"])
        score=np.where(long,predictions["M120"],predictions["M15"])
        threshold=np.where(long,policy["threshold"],policy["short_threshold"])
    else:
        score=predictions[kind]; threshold=policy["threshold"]
    legal=np.isfinite(score)&np.isfinite(inputs)&(inputs>=54)
    alarms=replay(inputs,np.where(legal,score,np.nan),threshold)
    alarms=alarms[alarms<SCORE_END]
    onsets=onsets[onsets<=SCORE_END+6]
    ledger,matched=match_alarms(alarms,onsets,truth)
    matchset=set(matched.tolist())
    inside_events=[j for j,e in enumerate(onsets) if SCORE_START<=e<SCORE_END]
    count=len(inside_events)
    warned=sum(j in matchset for j in inside_events)
    opportunities=sum(bool(legal[max(0,onsets[j]-6):max(0,onsets[j]-2)].any()) for j in inside_events)
    if warned>opportunities or opportunities>count:
        raise AssertionError("Warning sensitivity cannot exceed opportunity coverage")
    scored_alarms=[r for r in ledger if SCORE_START<=r["alarm_index"]<SCORE_END]
    return {"events":count,"matched":warned,"opportunities":opportunities,
            "alarms":len(scored_alarms),"unconfirmed":sum(r["category"]!="matched" for r in scored_alarms),
            "truth_observed_days":SCORE_BINS/288,"input_observed_days":float(np.isfinite(inputs[SCORE_START:SCORE_END]).sum()/288)}


def prepare_schedule(root):
    root=Path(root); folder=root/"results"; folder.mkdir(exist_ok=True)
    path=folder/"mask_stress_schedule.json"
    manifest=json.loads((root/"participant_manifest.json").read_text())
    traces=[]; exclusions=[]
    # Validation belongs to development; no test/external array is opened here.
    for row in sorted(manifest["participants"],key=lambda r:int(r["id"])):
        if row["split"]!="validation":
            continue
        with np.load(row["file"]) as z:
            windows=complete_windows(z["glucose"])
        if not windows:
            exclusions.append({"id":row["id"],"reason":"no_complete_328_bin_trace"})
        source_hash=digest(row["file"])
        for j,(a,b) in enumerate(windows):
            traces.append({"id":row["id"],"window":j,"start_index":a,"stop_index":b,"source_file":row["file"],"source_sha256":source_hash})
    schedule={"experiment":"development-only artificial input masking; original outcome truth retained",
        "participant_manifest_sha256":digest(root/"participant_manifest.json"),
        "selection":"validation participants; earliest up to seven nonoverlapping complete traces per person, independent of glucose values/event counts",
        "trace_bins":TRACE_BINS,"score_start":SCORE_START,"score_end_exclusive":SCORE_END,
        "block_minutes":DURATIONS,"seeds":SEEDS,"blocks_per_scored_day":4,
        "mask_start_rule":"one uniform start in each six-hour quarter with sufficient room for a 60-minute block; paired start times across durations",
        "burden_denominator":"original fully observed truth day; input-observed-time rate also retained",
        "traces":traces,"exclusions":exclusions,
        "positions":{str(seed):mask_positions(seed,len(traces)).tolist() for seed in SEEDS},
        "labels_derived_from_masked_inputs":False,"used_for_model_or_threshold_tuning":False,
        "boundary_convention":"150-minute unmasked initialization; score24h;45-minute true followup; pre-window warning matches retained"}
    if path.exists():
        old=json.loads(path.read_text())
        if old!=schedule:
            raise ValueError("Frozen stress schedule/input identity changed; preserve original and explicitly amend")
    else:
        path.write_text(json.dumps(schedule,indent=2)+"\n")
    return schedule


def predict_traces(models,traces):
    from card001_models import predict
    # Twenty-five NaNs prevent any history from crossing a trace boundary.
    chunks=[]; spans=[]; cursor=0
    for trace in traces:
        chunks.extend([trace,np.full(25,np.nan)])
        spans.append((cursor,cursor+len(trace)))
        cursor+=len(trace)+25
    joined=np.concatenate(chunks)
    whole={kind:predict(model,joined,kind) for kind,model in models.items()}
    whole["linear"]=-linear_forecast(joined)
    return [{kind:values[a:b] for kind,values in whole.items()} for a,b in spans]


def run_stress(root):
    root=Path(root); output=root/"results"; output.mkdir(exist_ok=True)
    started=time.perf_counter()
    schedule=prepare_schedule(root)
    frozen_path=root/"calibration.json"
    frozen=json.loads(frozen_path.read_text())
    if frozen.get("status")!="FROZEN_BEFORE_TEST":
        raise ValueError("Model and policy calibration must be frozen before stress outcomes")
    models={}
    for kind,expected in frozen["model_sha256"].items():
        path=root/"models"/f"{kind}.pkl"
        if digest(path)!=expected:
            raise ValueError("Frozen model checksum mismatch")
        with path.open("rb") as f:
            models[kind]=pickle.load(f)
    traces=[]
    for record in schedule["traces"]:
        if digest(record["source_file"])!=record["source_sha256"]:
            raise ValueError("Frozen stress source trace changed")
        with np.load(record["source_file"]) as z:
            traces.append(z["glucose"][record["start_index"]:record["stop_index"]].copy())
    if not traces:
        raise ValueError("No eligible fully observed development traces; stress experiment unexecuted")
    policies={**frozen["comparators"],"switch":frozen["switch"]}
    rows=[]
    def score_setting(masked,predictions,minutes,seed):
        for record,truth,inputs,probs in zip(schedule["traces"],traces,masked,predictions):
            for name,policy in policies.items():
                row={"id":record["id"],"window":record["window"],"block_minutes":minutes,"seed":seed,"policy":name,**trace_score(truth,inputs,probs,policy)}
                rows.append(row)
    score_setting(traces,predict_traces(models,traces),0,-1)
    for seed in SEEDS:
        starts=schedule["positions"][str(seed)]
        for minutes in DURATIONS:
            masked=[mask_trace(t,s,minutes) for t,s in zip(traces,starts)]
            score_setting(masked,predict_traces(models,masked),minutes,seed)
        print(f"mask_stress_seed={seed}; seconds={time.perf_counter()-started:.1f}",flush=True)
    frame=pd.DataFrame(rows)
    # An invariant proves that masked predictor inputs never redefined truth.
    invariant=frame.groupby(["id","window"]).events.nunique()
    if not invariant.eq(1).all():
        raise AssertionError("Artificial masks changed the outcome event set")
    frame.to_csv(output/"mask_stress_trace_scores.csv",index=False)
    keys=["block_minutes","seed","policy"]
    summary=frame.groupby(keys,sort=True)[["events","matched","opportunities","alarms","unconfirmed","truth_observed_days","input_observed_days"]].sum().reset_index()
    summary["sensitivity"]=summary.matched/summary.events.replace(0,np.nan)
    summary["opportunity_coverage"]=summary.opportunities/summary.events.replace(0,np.nan)
    summary["unconfirmed_per_truth_day"]=summary.unconfirmed/summary.truth_observed_days
    summary["alarms_per_truth_day"]=summary.alarms/summary.truth_observed_days
    summary["unconfirmed_per_input_observed_day"]=summary.unconfirmed/summary.input_observed_days
    truth_events=int(summary.loc[(summary.block_minutes==0)&summary.policy.eq("switch"),"events"].iloc[0])
    summary.to_csv(output/"mask_stress_summary.csv",index=False)
    aggregate=summary.groupby(["block_minutes","policy"])[["sensitivity","opportunity_coverage","unconfirmed_per_truth_day","alarms_per_truth_day"]].agg(["median","min","max"])
    aggregate.to_csv(output/"mask_stress_seed_range.csv")
    comparator=frozen["comparator"]["kind"]
    comparisons=[]
    for (minutes,seed),part in summary.groupby(["block_minutes","seed"]):
        by=part.set_index("policy")
        comparisons.append({"block_minutes":int(minutes),"seed":int(seed),
            "switch_sensitivity":float(by.loc["switch","sensitivity"]),
            "comparator_sensitivity":float(by.loc[comparator,"sensitivity"]),
            "sensitivity_difference":float(by.loc["switch","sensitivity"]-by.loc[comparator,"sensitivity"]),
            "unconfirmed_rate_difference":float(by.loc["switch","unconfirmed_per_truth_day"]-by.loc[comparator,"unconfirmed_per_truth_day"])})
    comparison=pd.DataFrame(comparisons)
    comparison.to_csv(output/"mask_stress_contrasts.csv",index=False)
    contrast_ranges=comparison.groupby("block_minutes")[["switch_sensitivity","comparator_sensitivity","sensitivity_difference","unconfirmed_rate_difference"]].agg(["median","min","max"])
    contrast_ranges.to_csv(output/"mask_stress_contrast_ranges.csv")
    os.environ["MPLCONFIGDIR"]=str(root/"cache"/"matplotlib")
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(10.5,4.4),layout="constrained")
    for name,color in [("M120","#607d8b"),("M15","#ab7445"),("linear","#6d65a0"),("missing","#42875b"),("switch","#036c8c")]:
        part=summary[summary.policy.eq(name)].groupby("block_minutes")
        for ax,metric in zip(axes,["sensitivity","unconfirmed_per_truth_day"]):
            stats=part[metric].agg(["median","min","max"])
            ax.plot(stats.index,stats["median"],marker="o",label=name,color=color,linewidth=1.5)
            ax.fill_between(stats.index,stats["min"],stats["max"],color=color,alpha=.10)
            ax.set_xlabel("Minutes hidden in each of four daily input blocks")
            ax.set_xticks([0,5,15,30,60])
            ax.spines[["top","right"]].set_visible(False)
    axes[0].set_ylabel("Useful-warning sensitivity on unchanged outcomes")
    axes[1].set_ylabel("Unconfirmed alarms per retained truth day")
    axes[1].legend(frameon=False,fontsize=9)
    fig.suptitle("Artificial masks on development traces: medians and ranges over 20 seeds",fontsize=12)
    fig.savefig(output/"mask_stress_figure.png",dpi=180)
    plt.close(fig)
    result={"status":"executed_artificial_development_experiment","hypothesis_status":"not_a_confirmatory_H1_test",
        "n_participants":len({r['id'] for r in schedule['traces']}),"n_complete_windows":len(traces),"truth_events":truth_events,
        "mask_seeds":len(SEEDS),"durations_minutes":DURATIONS,"truth_events_unchanged":True,
        "model_sha256":frozen['model_sha256'],"calibration_sha256":digest(frozen_path),
        "schedule_sha256":digest(output/"mask_stress_schedule.json"),"code_sha256":digest(__file__),
        "comparator":frozen["comparator"]["kind"],"runtime_seconds":time.perf_counter()-started,
        "interpretation":"Artificial masking isolates input-history loss on a selected fully observed development subset. Seeds are mask realizations, not independent patients. No efficacy, transfer or native-gap hypothesis is established by this experiment."}
    (output/"mask_stress_results.json").write_text(json.dumps(result,indent=2)+"\n")
    text="# Artificial input-mask stress experiment\n\n"
    text+=f"Executed four block durations (5, 15, 30, 60 minutes) across 20 fixed mask seeds on {len(traces)} fully observed development windows from {result['n_participants']} validation participants. Each selected 24-hour scoring interval had 150-minute lookback and 45-minute outcome followup. Original outcome events and truth-time denominators remained identical for every mask.\n\n"
    text+=f"The preselected subset contains **{truth_events} qualifying truth events**. This small event denominator limits interpretation of every sensitivity estimate; rows and seeds are not independent patient replication.\n\n"
    text+="No model, threshold or segment selection was tuned using stress results. The first up to seven eligible complete windows per person were selected chronologically, independent of low-event counts. This is a bounded fully observed subset, not a native-gap test or independent clinical evaluation. Four input blocks occur per day, one per six-hour quarter; mask start positions are paired across durations.\n\n"
    text+="Seed-level results and median/minimum/maximum summaries are in the accompanying CSV files. Ranges describe artificial-mask variability; they are not participant confidence intervals. Alarm rates use retained truth days to isolate the artificial availability intervention; rates per remaining input-observed day are also saved. Native-gap primary results retain their original common observed-sensor-day denominator.\n\n"
    text+="Sensitivity columns show separate policy medians. The contrast column is the median of paired seed differences, which need not equal the difference of those separate medians.\n\n"
    text+=f"The frozen comparator is **{comparator}**. The primary two-policy descriptive comparison is:\n\n"
    text+="| Block duration | Switch sensitivity, median | Comparator sensitivity, median | Difference, median [seed range] | Unconfirmed-rate difference, median [seed range] |\n|---|---:|---:|---:|---:|\n"
    def fmt(x):
        return f"{x:.3f}" if np.isfinite(x) else "NA"
    for minutes,values in contrast_ranges.iterrows():
        text+=f"| {minutes} min | {fmt(values[('switch_sensitivity','median')])} | {fmt(values[('comparator_sensitivity','median')])} | {fmt(values[('sensitivity_difference','median')])} [{fmt(values[('sensitivity_difference','min')])}, {fmt(values[('sensitivity_difference','max')])}] | {fmt(values[('unconfirmed_rate_difference','median')])} [{fmt(values[('unconfirmed_rate_difference','min')])}, {fmt(values[('unconfirmed_rate_difference','max')])}] |\n"
    text+="\n![Artificial-mask stress](mask_stress_figure.png)\n\n"
    text+="The unmasked baseline is explicitly reported alongside all durations, all policies and all 20 seeds. No favorable synthetic result can replace failure of the original two-cohort hypothesis.\n"
    (output/"mask_stress_report.md").write_text(text,encoding="utf-8")
    print(json.dumps(result,indent=2))
    return result


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True)
    parser.add_argument("--freeze-only",action="store_true")
    parser.add_argument("--threads",type=int,default=2)
    args=parser.parse_args()
    root=args.data_root.resolve()
    if args.threads<1:
        raise ValueError("Thread count must be positive")
    os.environ["OMP_NUM_THREADS"]=str(args.threads)
    if args.freeze_only:
        s=prepare_schedule(root)
        print(json.dumps({"status":"schedule_frozen_no_stress_outcomes_computed","traces":len(s['traces']),"seeds":len(SEEDS)}))
    else:
        run_stress(root)


if __name__=="__main__":
    main()
