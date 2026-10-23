"""Complete CLVer replay and clinical analysis under the independently frozen policy.

Run after card003_data.py --prepare and card003_policy.freeze_from_card001.
Raw, processed and result artifacts use the explicitly supplied study root.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import argparse
import json
import math
import pickle
import platform
import time
from datetime import datetime, timezone
from pathlib import Path
import numpy as np
import pandas as pd
from card003_core import (fit_association, fixed_design, gate_assessment,
                          replay_interval, sha256_file, SEED)
from card003_data import window_grid
from card001_core import linear_forecast
from card001_models import predict


def clean(value):
    if isinstance(value,dict): return {str(k):clean(v) for k,v in value.items()}
    if isinstance(value,(list,tuple,np.ndarray)): return [clean(v) for v in value]
    if isinstance(value,(np.integer,)): return int(value)
    if isinstance(value,(np.bool_,)): return bool(value)
    if isinstance(value,(float,np.floating)):
        return float(value) if math.isfinite(value) else None
    return value


def dump(path,value):
    path.parent.mkdir(parents=True,exist_ok=True)
    path.write_text(json.dumps(clean(value),indent=2,allow_nan=False),encoding="utf-8")


def load_frozen(root):
    path=root/"models"/"frozen_policy.json"
    frozen=json.loads(path.read_text())
    if frozen["status"]!="FROZEN_BEFORE_CLVER_REPLAY" or frozen["CLVer_outcomes_used"]:
        raise ValueError("No independently frozen Card003 policy")
    models={}
    for kind,record in frozen["models"].items():
        modelpath=Path(record["path"])
        if sha256_file(modelpath)!=record["sha256"]:raise ValueError("Frozen model hash changed")
        with modelpath.open("rb") as f: models[kind]=pickle.load(f)
    return frozen,models,sha256_file(path)


def policy_arrays(predictions,choice):
    kind=choice["kind"]
    if kind=="never":
        n=len(predictions["M15"])
        return np.zeros(n),np.inf,np.ones(n,bool)
    if kind=="switch":
        long=np.isfinite(predictions["M120"])
        p=np.where(long,predictions["M120"],predictions["M15"])
        cut=np.where(long,choice["threshold"],choice["short_threshold"])
    else:p=predictions[kind];cut=choice["threshold"]
    return p,cut,np.isfinite(p)


def predict_grid(models,g):
    return {**{kind:predict(model,g,kind) for kind,model in models.items()},"linear":-linear_forecast(g)}


def transfer_table(frame):
    def total(f):
        events=int(f.events.sum());days=float(f.observed_days.sum());matched=int(f.matched.sum())
        return dict(intervals=len(f),children=int(f.PtID.nunique()),events=events,matched=matched,
                    unwarned=int(f.unwarned.sum()),observed_days=days,alarms=int(f.alarms.sum()),
                    unconfirmed=int(f.unconfirmed.sum()),sensitivity=matched/events if events else np.nan,
                    unconfirmed_per_day=f.unconfirmed.sum()/days if days else np.nan,
                    alarms_per_day=f.alarms.sum()/days if days else np.nan,
                    total_per_100days=events/days*100 if days else np.nan,
                    unwarned_per_100days=f.unwarned.sum()/days*100 if days else np.nan,
                    opportunity_coverage=f.opportunity_events.sum()/events if events else np.nan,
                    residual_children=int(f.loc[f.unwarned>0,"PtID"].nunique()))
    result=total(frame)
    result["overall_sensitivity"]=result["sensitivity"]
    result["delivery_arms"]={str(k):total(g) for k,g in frame.groupby("delivery_arm")}
    result["management_groups"]={str(k):total(g) for k,g in frame.groupby("management")}
    return result


def replay_frames(root,frame,models,frozen,*,variant="primary",phase=0,delay=0,
                  exclude_blank=False,run_length=4,all_policies=False,include_date_only=False):
    out=root/"results";predfolder=root/"processed"/"predictions"/variant
    predfolder.mkdir(parents=True,exist_ok=True)
    choices=dict(frozen["choices"]["frozen_family"])
    choices["never"]={"kind":"never","threshold":np.inf}
    primary_kind=frozen["choices"]["primary"]["kind"]
    if not all_policies:choices={primary_kind:choices[primary_kind]}
    rows={k:[] for k in choices};events={k:[] for k in choices};alarms={k:[] for k in choices}
    for j,row in enumerate(frame.to_dict("records")):
        start,end=row["start_seconds"],row["end_seconds"]
        if phase==0 and not exclude_blank and not include_date_only:
            with np.load(row["grid_path"]) as data:grid={k:data[k] for k in data.files}
        else:grid=window_grid(root,str(row["PtID"]),start,end,phase_minutes=phase,exclude_blank=exclude_blank,include_date_only=include_date_only)
        target_g=grid["glucose"];t=grid["time"]
        if delay:
            # Delay input availability only; physiological/event clock stays fixed.
            delayed=window_grid(root,str(row["PtID"]),start,end,phase_minutes=phase,delay_minutes=delay,exclude_blank=exclude_blank,include_date_only=include_date_only)
            if not np.array_equal(t,delayed["time"]):raise AssertionError("Delay must preserve event grid")
            warning_g=delayed["glucose"]
        else:warning_g=target_g
        path=predfolder/f"{row['pair_id']}.npz"
        cache_valid=False
        if path.exists():
            with np.load(path) as cached:
                if str(cached["policy_hash"] )!=frozen["artifact_hash"]:raise ValueError("Stale policy prediction cache")
                cache_valid=np.array_equal(cached["glucose"],warning_g,equal_nan=True) and np.array_equal(cached["time"],t)
                if cache_valid:predictions={k:cached[k] for k in ("M120","M15","missing","linear")}
            if not cache_valid:
                archived=root/"history"/"obsolete_predictions"/variant/(path.stem+"_"+sha256_file(path)+".npz")
                archived.parent.mkdir(parents=True,exist_ok=True)
                if not archived.exists():archived.write_bytes(path.read_bytes())
        if not cache_valid:
            predictions=predict_grid(models,warning_g)
            np.savez_compressed(path,**predictions,glucose=warning_g,time=t,policy_hash=frozen["artifact_hash"])
        for kind,choice in choices.items():
            p,cut,available=policy_arrays(predictions,choice)
            metrics,ee,aa=replay_interval(t,target_g,p,cut,available,start,end,run_length,warning_g)
            r={**row,**metrics,"policy":kind,"variant":variant}
            if metrics["unwarned"]+metrics["matched"]!=metrics["events"]:raise AssertionError("Broken residual identity")
            rows[kind].append(r)
            events[kind].extend(ee.assign(PtID=row["PtID"],pair_id=row["pair_id"],variant=variant,policy=kind).to_dict("records"))
            alarms[kind].extend(aa.assign(PtID=row["PtID"],pair_id=row["pair_id"],variant=variant,policy=kind).to_dict("records"))
        if (j+1)%40==0:print(f"replay {variant}: {j+1}/{len(frame)} intervals",flush=True)
    for kind in choices:
        pd.DataFrame(events[kind]).to_csv(out/f"events_{variant}_{kind}.csv",index=False)
        pd.DataFrame(alarms[kind]).to_csv(out/f"alarms_{variant}_{kind}.csv",index=False)
    pd.DataFrame(events[primary_kind]).to_csv(out/f"events_{variant}.csv",index=False)
    pd.DataFrame(alarms[primary_kind]).to_csv(out/f"alarms_{variant}.csv",index=False)
    answer={k:pd.DataFrame(v) for k,v in rows.items()}
    for kind,table in answer.items():table.to_csv(out/f"intervals_{variant}_{kind}.csv",index=False)
    return answer


def fit_and_save(root,name,frame,*,n_boot=0,**kwargs):
    t=time.perf_counter()
    result,boot=fit_association(frame,n_boot=n_boot,**kwargs)
    result["elapsed_seconds"]=time.perf_counter()-t
    result["analysis_name"]=name
    result["inferential_role"]="primary_bootstrap" if name=="primary" else "companion_bootstrap" if n_boot==5000 else "sensitivity_cluster_sandwich_only"
    if n_boot==0 and "rr" in result:
        result["status"]="secondary_estimate_robust_interval"
        result["reason"]="not_a_primary_5000_bootstrap_inference"
    if len(boot):boot.to_csv(root/"results"/f"bootstrap_{name}.csv.gz",index=False,compression="gzip")
    dump(root/"results"/f"model_{name}.json",result)
    print(f"fit {name}: {result.get('rr')} CI={result.get('bootstrap_95_ci',result.get('robust_95_ci'))} status={result['status']} sec={result['elapsed_seconds']:.2f}",flush=True)
    return result


def descriptive_tables(root,primary):
    out=root/"results"
    clinical=pd.read_csv(root/"processed"/"clinical_derivations.csv",dtype={"PtID":str})
    clinical["clinical_state"]=clinical.remission.map({0:"not remission",1:"remission"}).fillna("unknown")
    pd.crosstab(clinical.index_visit,clinical.clinical_state,dropna=False).to_csv(out/"remission_visit_frequency.csv")
    transitions=[]
    for pid,g in clinical.sort_values("index_visit").groupby("PtID"):
        for a,b in zip(g.to_dict("records"),g.to_dict("records")[1:]):
            transitions.append(dict(PtID=pid,index_visit=a["index_visit"],target_visit=b["index_visit"],
                from_state=a["clinical_state"],to_state=b["clinical_state"]))
    tf=pd.DataFrame(transitions);tf.to_csv(out/"remission_transitions.csv",index=False)
    tf.groupby(["index_visit","target_visit","from_state","to_state"]).size().rename("n").reset_index().to_csv(out/"remission_transition_counts.csv",index=False)
    state=primary.copy()
    state["clinical_state"]=state.remission.map({0:"not remission",1:"remission"}).fillna("unknown")
    state["peak_above_0_30"]=np.where(state.peak_nmol_l.isna(),"unknown",np.where(state.peak_nmol_l>.30,"yes","no"))
    state["auc_band"]=pd.cut(state.auc_nmol_l,[-np.inf,.2,.4,.6,np.inf],right=False,labels=["<0.20","0.20–<0.40","0.40–<0.60","≥0.60"])
    summaries={}
    for name,grouping in [("remission",["clinical_state"]),("remission_peak",["clinical_state","peak_above_0_30"]),("auc_bands",["auc_band"]),("management",["management"])]:
        records=[]
        for keys,g in state.groupby(grouping,observed=False,dropna=False):
            if not isinstance(keys,tuple):keys=(keys,)
            summary=transfer_table(g)
            summary.pop("delivery_arms");summary.pop("management_groups")
            records.append({**dict(zip(grouping,keys)),**summary})
        table=pd.DataFrame(records);table.to_csv(out/f"descriptive_{name}.csv",index=False);summaries[name]=records
    primary.auc_nmol_l.describe(percentiles=[.05,.25,.5,.75,.95]).to_csv(out/"auc_distribution.csv")
    # Inclusion/missingness is summarized without imputing unobserved outcomes.
    allpairs=pd.read_csv(root/"processed"/"interval_derivations.csv",dtype={"PtID":str})
    included=set(primary.pair_id)
    allpairs["primary_included"]=allpairs.pair_id.isin(included)
    characteristics=allpairs.groupby("primary_included").agg(intervals=("pair_id","size"),children=("PtID","nunique"),auc_mean=("auc_nmol_l","mean"),age_mean=("age","mean"),known_remission=("remission","count"),remission_mean=("remission","mean"))
    characteristics.to_csv(out/"inclusion_characteristics.csv")
    return summaries


def predicted_curve(root,frame):
    import os
    config=Path(root)/"tools"/"matplotlib_cache";config.mkdir(parents=True,exist_ok=True)
    os.environ["MPLCONFIGDIR"]=str(config)
    import statsmodels.api as sm
    from statsmodels.stats.sandwich_covariance import cov_cluster
    X=fixed_design(frame)
    fit=sm.GLM(frame.unwarned.to_numpy(),X.to_numpy(),family=sm.families.Poisson(),offset=np.log(frame.observed_days)).fit(maxiter=100)
    covariance=cov_cluster(fit,frame.PtID.to_numpy())
    rows=[]
    for value in np.linspace(frame.auc_nmol_l.min(),frame.auc_nmol_l.max(),101):
        xx=X.copy();xx["exposure"]=value/.20
        rates=np.exp(xx.to_numpy()@fit.params)
        mean=float(rates.mean());gradient=(rates[:,None]*xx.to_numpy()).mean(axis=0)/mean
        se=float(np.sqrt(gradient@covariance@gradient));limits=100*np.exp(np.log(mean)+np.array([-1,1])*1.95996398454*se)
        rows.append(dict(auc_nmol_l=value,standardized_unwarned_per_100days=mean*100,robust_lower=limits[0],robust_upper=limits[1]))
    curve=pd.DataFrame(rows);curve.to_csv(root/"results"/"predicted_rate_curve.csv",index=False)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig,axes=plt.subplots(1,2,figsize=(9,3.7))
    axes[0].hist(frame.auc_nmol_l,bins=20,color="#32658a");axes[0].set(xlabel="Index normalized C-peptide AUC (nmol/L)",ylabel="Intervals")
    axes[1].plot(curve.auc_nmol_l,curve.standardized_unwarned_per_100days,color="#32658a")
    axes[1].fill_between(curve.auc_nmol_l,curve.robust_lower,curve.robust_upper,color="#32658a",alpha=.18)
    axes[1].set(xlabel="Index normalized C-peptide AUC (nmol/L)",ylabel="Predicted unwarned events\nper 100 observed days")
    fig.suptitle("CLVer observational association; warning-transfer burden gate failed",fontsize=10)
    fig.tight_layout();fig.savefig(root/"results"/"exposure_and_predicted_rate.png",dpi=180);plt.close(fig)


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root",type=Path,required=True)
    ap.add_argument("--stage",choices=["primary","complete"],default="complete")
    args=ap.parse_args();root=args.data_root.resolve();start=time.perf_counter()
    frozen,models,digest=load_frozen(root);frozen["artifact_hash"]=digest
    frame=pd.read_csv(root/"processed"/"intervals_all_coverage.csv",dtype={"PtID":str})
    commitment=dict(started_utc=datetime.now(timezone.utc).isoformat(),frozen_policy_sha256=digest,
        source_sha256=sha256_file(root/"raw"/"CLVerPublicDataset.zip"),
        input_intervals_sha256=sha256_file(root/"processed"/"intervals_all_coverage.csv"),
        code_sha256={p.name:sha256_file(p) for p in Path(__file__).parent.glob("card003*.py")},
        version=dict(python=platform.python_version(),numpy=np.__version__,pandas=pd.__version__),
        bootstrap_primary=5000,bootstrap_seed=SEED,sensitivity_uncertainty="child-cluster sandwich intervals; not primary bootstrap claims")
    dump(root/"results"/"analysis_commitment.json",commitment)
    # Preserve the exact executing scientific sources even when the report
    # wrapper or methods prose is edited later in the same research session.
    snapshot=root/"results"/"analysis_source_snapshots"/sha256_file(root/"results"/"analysis_commitment.json")
    snapshot.mkdir(parents=True,exist_ok=True)
    for name,source_digest in commitment["code_sha256"].items():
        original=Path(__file__).parent/name
        if sha256_file(original)!=source_digest:raise ValueError("Source changed during analysis commitment")
        (snapshot/name).write_bytes(original.read_bytes())
    tables=replay_frames(root,frame,models,frozen,all_policies=True)
    kind=frozen["choices"]["primary"]["kind"]
    all_primary=tables[kind];primary=all_primary[all_primary.coverage>=.7].copy()
    primary.to_csv(root/"results"/"primary_intervals.csv",index=False)
    transfer=transfer_table(primary);dump(root/"results"/"transfer.json",transfer)
    # Benchmark the first committed 100 draws; final run repeats identical seed
    # and retains every one of the required 5000 draws.
    benchmark=fit_and_save(root,"benchmark_100",primary,n_boot=100)
    result=fit_and_save(root,"primary",primary,n_boot=5000)
    gates=gate_assessment(primary,transfer,result,True)
    dump(root/"results"/"primary_gates.json",gates)
    descriptions=descriptive_tables(root,primary)
    if "rr" in result:predicted_curve(root,primary)
    companions={}
    companions["total_events"]=fit_and_save(root,"total_events",primary,outcome="events",n_boot=5000)
    known=primary[primary.remission.notna()]
    companions["remission"]=fit_and_save(root,"remission",known,exposure="remission",n_boot=5000)
    if args.stage=="primary":
        dump(root/"results"/"analysis_summary.json",dict(complete=False,status="primary_done_sensitivities_pending",primary=result,gates=gates,transfer=transfer,companions=companions))
        return
    sensitivities={}
    def assess(name,f,**kwargs):
        if not kwargs.get("within_child",False):kwargs["drop_constant_nuisance"]=True
        sensitivities[name]=fit_and_save(root,name,f,**kwargs)
        sensitivities[name]["transfer"]=transfer_table(f)
        return sensitivities[name]
    for policy,table in tables.items():assess("policy_"+policy,table[table.coverage>=.7])
    assess("coverage_120hours",all_primary[all_primary.observed_days>=5])
    assess("coverage_90percent",all_primary[all_primary.coverage>=.9])
    assess("DeviceCGM_only",primary[primary.source.eq("DeviceCGM")])
    assess("exclude_any_clock_barrier_interval",primary[primary.clock_barrier_bins.eq(0)])
    assess("temperature_OK",primary[primary.temperature_all_ok.eq(True)])
    assess("augmented_HbA1c_insulin",primary.dropna(subset=["hba1c","insulin_u_kg_day"]),augmented=True)
    within=assess("within_child",primary,within_child=True)
    pos=primary[primary.groupby("PtID").unwarned.transform("sum")>0]
    within["within_child_variation"]=pos.groupby("PtID").auc_nmol_l.agg(["min","max","count"]).assign(range=lambda f:f["max"]-f["min"]).to_dict("index")
    contemporaneous=primary[primary.contemporaneous_assay_eligible.eq(True)&primary.contemporaneous_auc_nmol_l.notna()].copy()
    contemporaneous["auc_nmol_l"]=contemporaneous.contemporaneous_auc_nmol_l;contemporaneous["lag_days"]=-28
    assess("contemporaneous",contemporaneous,omit_constant_lag=True)
    assess("quadratic_AUC",primary,quadratic=True)
    lowchild=primary.groupby("PtID").auc_nmol_l.mean().sort_values(kind="stable").index[0]
    highchild=primary.groupby("PtID").events.sum().sort_values(ascending=False,kind="stable").index[0]
    assess("drop_lowest_mean_secretion_child",primary[primary.PtID.ne(lowchild)])["excluded_child"]=lowchild
    assess("drop_highest_total_event_child",primary[primary.PtID.ne(highchild)])["excluded_child"]=highchild
    # Constant nuisance columns are algebraically redundant with the intercept;
    # removing them preserves the model column space and is secondary-only.
    for name,g in primary.groupby("management"):
        r=assess("management_"+name.replace("/","_").replace(" ","_"),g)
        r["exploratory_limitation"]="Within-management exploratory fit; constant nuisance columns removed without changing model column space."
    for variant,kwargs in [("phase1",dict(phase=1)),("phase2",dict(phase=2)),("phase3",dict(phase=3)),("phase4",dict(phase=4)),("delay5",dict(delay=5)),("run3",dict(run_length=3)),("blank_units_excluded",dict(exclude_blank=True)),("date_only_midnight_assumption",dict(include_date_only=True))]:
        table=replay_frames(root,frame,models,frozen,variant=variant,**kwargs)[kind]
        assess(variant,table[table.coverage>=.7])
    # Phase0 is the primary pipeline, explicitly included in sensitivity matrix.
    sensitivities["phase0"]={**result,"transfer":transfer}
    dump(root/"results"/"sensitivity_results.json",sensitivities)
    sensitivity_rows=[]
    for name,r in sensitivities.items():
        ci=r.get("robust_95_ci",[None,None])
        sensitivity_rows.append(dict(analysis=name,children=r.get("n_children"),intervals=r.get("n_intervals"),
            events=r.get("events"),rr=r.get("rr"),robust_lower=ci[0],robust_upper=ci[1],status=r["status"],reason=r.get("reason")))
    pd.DataFrame(sensitivity_rows).to_csv(root/"results"/"sensitivity_matrix.csv",index=False)
    try:
        import psutil
        memory=psutil.Process().memory_info()._asdict()
    except ImportError:memory={"peak_memory_unavailable":True}
    summary=dict(complete=True,status="completed_empirical_run",primary=result,gates=gates,transfer=transfer,
        companions=companions,sensitivity_names=list(sensitivities),description=descriptions,
        primary_policy=kind,primary_policy_artifact_sha256=digest,
        elapsed_seconds=time.perf_counter()-start,memory=memory,
        limitations=["Observational C-peptide/remission associations, not treatment effects.",
        "VisitInfo date is a protocol-based proxy for undated physical-exam weight.",
        "Blank DeviceCGM units assumed mg/dL; explicit exclusion sensitivity retained.",
        "Clock discontinuities receive an unknown separating bin; later source segments remain observed.",
        "Secondary sensitivity intervals use clustered sandwich uncertainty, not primary 5000-draw bootstrap.",
        "Replay does not identify event prevention, actual caregiver awareness or clinical deployment effects."])
    dump(root/"results"/"analysis_summary.json",summary)
    dump(root/"results"/"result_hashes.json",{str(p.relative_to(root)):sha256_file(p) for p in (root/"results").glob("*") if p.is_file() and p.name!="result_hashes.json"})
    print(json.dumps(clean(dict(status=summary["status"],gates=gates,primary=result,elapsed_seconds=summary["elapsed_seconds"]))),flush=True)


if __name__=="__main__":main()
