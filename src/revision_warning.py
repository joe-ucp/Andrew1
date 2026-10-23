"""Post-primary, fixed-model warning diagnostics. Original outputs are read-only."""
from __future__ import annotations
import argparse,hashlib,io,json,os,pickle,shutil,time,zipfile
from collections import Counter
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from card001_core import ready,episodes,linear_forecast,paired_bootstrap
from card001_models import predict
from card001_evaluate import score_policy,fast_grid,specifications

ORIGINAL=None
OUT=None
CLINICAL=None
SEED=20261023
POLICIES={"M15_004":dict(kind="M15",threshold=.04),"M15_001":dict(kind="M15",threshold=.01),"switch_004_004":dict(kind="switch",threshold=.04,short_threshold=.04),"switch_004_001":dict(kind="switch",threshold=.04,short_threshold=.01),"M120_004":dict(kind="M120",threshold=.04)}
CONTRASTS=[("switch_004_001","M15_004","Frozen primary contrast"),("switch_004_004","M15_004","Model/history choice at common numeric threshold"),("switch_004_001","switch_004_004","Fallback threshold change with fixed model choice"),("M15_001","M15_004","Short-model threshold change"),("M120_004","M15_004","Long versus short standalone at common numeric threshold")]

def sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def write_json(path,value):
    def safe(v):
        if isinstance(v,dict):return {str(k):safe(x) for k,x in v.items()}
        if isinstance(v,(list,tuple)):return [safe(x) for x in v]
        if isinstance(v,np.generic):return safe(v.item())
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    Path(path).write_text(json.dumps(safe(value),indent=2,allow_nan=False),encoding="utf-8")

def protocol():
    return dict(status="EXPLORATORY_AFTER_PRIMARY_RESULTS",seed=SEED,selection_bootstrap_draws=1000,selection_population="All 29 original DCLP3 validation participants, including zero-observation participants",selection_target="Conditional stability given all originally fitted models; no model refitting, no external choice",availability="Strict 25-bin readiness implies strict 4-bin readiness. With identical current>=54 gating, L is a subset of S, so L union S=S. Applying any over a common pre-event window preserves containment and equality.",paired_table="Original frozen M15 versus switch flags: both, M15_only, switch_only, neither; all and M120-opportunity strata",ablations=POLICIES,contrasts=CONTRASTS,ablation_bootstrap_draws=10000,ablation_scope="Original held-out DCLP3 and REPLACE-BG; fixed models and thresholds, no favorable selection",demographics="Original enrollment age and recorded gender/race/ethnicity; duration=age minus reported diagnosis age; local screening HbA1c from exactly named Screening/Screening Visit; no later visits substituted",interpretation="Exploratory computational policy contrasts. Shared refractory state couples branches. Equal numeric probability thresholds do not equalize calibration or burden. No clinical intervention effect, independent training uncertainty, or corrected confirmatory claim is identified.")

def freeze():
    OUT.mkdir(parents=True,exist_ok=True); path=OUT/"exploratory_protocol.json"
    specification=protocol()
    if path.exists():
        old=json.loads(path.read_text());old.pop("frozen_utc",None)
        if old!=json.loads(json.dumps(specification)):raise ValueError("Exploratory protocol changed; create a further version")
    else:write_json(path,dict(specification,frozen_utc=datetime.now(timezone.utc).isoformat()))
    protocol_dir=OUT/"protocol";protocol_dir.mkdir(exist_ok=True)
    (protocol_dir/"revision_warning.md").write_text("# Exploratory warning revision protocol\n\nFrozen before executing these new diagnostics, after the original results were known. Original models, thresholds and outputs remain controlling.\n\n"+"\n\n".join(f"**{k}**: {json.dumps(v,ensure_ascii=False) if isinstance(v,(dict,list)) else v}" for k,v in specification.items())+"\n",encoding="utf-8")
    src=Path(__file__); dest=OUT/"source_snapshots"/sha(src);dest.mkdir(parents=True,exist_ok=True)
    snapshot=dest/src.name
    if snapshot.exists():
        if sha(snapshot)!=sha(src):raise ValueError("Revision source snapshot changed")
    else:shutil.copy2(src,snapshot)
    frozen=json.loads((ORIGINAL/"calibration.json").read_text())
    for k,h in frozen["model_sha256"].items():
        if sha(ORIGINAL/"models"/f"{k}.pkl")!=h:raise ValueError("Original fitted model changed")
    if sha(ORIGINAL/"participant_manifest.json")!=frozen["manifest_sha256"]:raise ValueError("Original participant manifest changed")
    dependencies=[ORIGINAL/"validation_prediction_manifest.json",ORIGINAL/"participant_manifest.json",ORIGINAL/"external"/"replacebg"/"participant_manifest.json"]
    if CLINICAL is not None:
        dependencies += [CLINICAL/"models"/"frozen_policy.json",CLINICAL/"results"/"primary_intervals.csv",CLINICAL/"results"/"events_primary.csv"]
    dependencies += [ORIGINAL/"results"/f"{c}_{suffix}" for c in ["internal","external"] for suffix in ["comparison.json","summary.csv","participant_scores.csv","alarm_ledger.csv","event_table.csv","censored_event_table.csv"]]
    write_json(OUT/"source_provenance.json",dict(protocol_sha256=sha(path),revision_source_sha256=sha(src),original_calibration_sha256=sha(ORIGINAL/"calibration.json"),model_sha256=frozen["model_sha256"],original_source_sha256={n:sha(src.parent/n) for n in ["card001_core.py","card001_models.py","card001_evaluate.py","card003_policy.py"]},original_input_sha256={str(p):sha(p) for p in dependencies}))

def cohort_manifest(cohort):
    root=ORIGINAL if cohort=="internal" else ORIGINAL/"external"/"replacebg"
    rows=json.loads((root/"participant_manifest.json").read_text())["participants"]
    return [p for p in rows if p["split"]==("test" if cohort=="internal" else "external_test")]

def availability_masks(g):
    allowed=np.isfinite(g)&(g>=54);long=ready(g,25)&allowed;short=ready(g,4)&allowed
    return long,short,np.where(long,True,short)

def event_availability(mask,onsets):return np.array([mask[max(0,int(e)-6):max(0,int(e)-2)].any() for e in onsets],bool)

def availability():
    rows=[]
    def check(cohort,pid,interval,g,onsets,scored=None,source_path=None):
        long,short,switch=availability_masks(g)
        if np.any(long&~short) or not np.array_equal(switch,short):raise AssertionError("Bin availability containment failed")
        a,b,c=[event_availability(m,onsets) for m in [long,short,switch]]
        if np.any(a&~b) or not np.array_equal(b,c):raise AssertionError("Event availability containment failed")
        if scored is None:scored=np.ones(len(g),bool)
        rows.append(dict(cohort=cohort,participant_id=pid,interval=interval,input_path=str(source_path),input_sha256=sha(source_path),scored_bins=int(scored.sum()),observed_bins=int((scored&np.isfinite(g)).sum()),long_prediction_bins=int((long&scored).sum()),short_prediction_bins=int((short&scored).sum()),switch_prediction_bins=int((switch&scored).sum()),events=len(onsets),long_opportunity_events=int(a.sum()),short_opportunity_events=int(b.sum()),switch_opportunity_events=int(c.sum()),subset_violations=0,equality_violations=0))
    all_dclp=json.loads((ORIGINAL/"participant_manifest.json").read_text())["participants"]
    for p in all_dclp:
        z=np.load(p["file"]);g=z["glucose"];check("DCLP3_"+p["split"],p["id"],"full",g,episodes(g)[0],source_path=p["file"])
    for p in cohort_manifest("external"):
        z=np.load(p["file"]);g=z["glucose"];check("REPLACE-BG",p["id"],"full",g,episodes(g)[0],source_path=p["file"])
    if CLINICAL is not None:
        clver=CLINICAL/"results";intervals=pd.read_csv(clver/"primary_intervals.csv");events=pd.read_csv(clver/"events_primary.csv")
        for p in intervals.itertuples():
            z=np.load(p.grid_path);g=z["glucose"];t=z["time"];e=events.loc[(events.pair_id==p.pair_id)&(events.policy=="M15")]
            onsets=e.onset_index.to_numpy(dtype=int)
            np.testing.assert_array_equal(t[onsets],e.onset_time.to_numpy())
            if len(onsets)!=p.events:raise AssertionError("CLVer preserved event identity differs")
            np.testing.assert_array_equal(event_availability(availability_masks(g)[1],onsets),e.legal_opportunity.to_numpy(bool))
            check("CLVer_primary",p.PtID,p.pair_id,g,onsets,(t>=p.start_seconds)&(t<p.end_seconds),source_path=p.grid_path)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/"availability_containment_by_interval.csv",index=False)
    summary=frame.groupby("cohort").agg(intervals=("interval","size"),participants=("participant_id","nunique"),events=("events","sum"),long_opportunity_events=("long_opportunity_events","sum"),short_opportunity_events=("short_opportunity_events","sum"),switch_opportunity_events=("switch_opportunity_events","sum"),subset_violations=("subset_violations","sum"),equality_violations=("equality_violations","sum")).reset_index()
    summary.to_csv(OUT/"availability_containment_summary.csv",index=False)
    write_json(OUT/"availability_proof.json",dict(status="PASS",cohorts=summary.cohort.tolist(),CLVer_evaluated=CLINICAL is not None,proof=protocol()["availability"],interpretation="Switching adds no prediction opportunity relative to the always-short comparator. Any warning difference must arise from model scores, threshold choice and alarm-state interactions."))

def paired():
    rows=[];ledger=[]
    for cohort in ["internal","external"]:
        p=pd.read_csv(ORIGINAL/"results"/f"{cohort}_event_table.csv")
        p["pair_category"]=np.select([p.M15_warning&p.switch_warning,p.M15_warning&~p.switch_warning,~p.M15_warning&p.switch_warning],["both","M15_only","switch_only"],default="neither")
        p["cohort"]=cohort;ledger.append(p)
        for label,sub in [("all",p),("no_long_opportunity",p.loc[~p.M120_opportunity]),("long_opportunity",p.loc[p.M120_opportunity])]:
            c=sub.pair_category.value_counts();d={n:int(c.get(n,0)) for n in ["both","M15_only","switch_only","neither"]}
            rows.append(dict(cohort=cohort,stratum=label,events=len(sub),**d,net_switch_minus_M15=d["switch_only"]-d["M15_only"],difference=(d["switch_only"]-d["M15_only"])/len(sub) if len(sub) else None))
    pd.concat(ledger).to_csv(OUT/"preserved_event_pairs.csv",index=False);pd.DataFrame(rows).to_csv(OUT/"paired_event_discordance.csv",index=False)

def ablation():
    models={k:pickle.load((ORIGINAL/"models"/f"{k}.pkl").open("rb")) for k in ["M120","M15"]};rows=[];identities={}
    model_hashes={k:sha(ORIGINAL/"models"/f"{k}.pkl") for k in models}
    prediction_source_hashes={n:sha(Path(__file__).parent/n) for n in ["card001_core.py","card001_models.py"]}
    pred_dir=OUT/"fixed_predictions";pred_dir.mkdir(exist_ok=True)
    for cohort in ["internal","external"]:
        original=pd.read_csv(ORIGINAL/"results"/f"{cohort}_participant_scores.csv")
        for p in cohort_manifest(cohort):
            z=np.load(p["file"]);g=z["glucose"];onsets=episodes(g)[0]
            path=pred_dir/f"{cohort}_{p['id']}.npz"
            sidecar=path.with_suffix(".json"); expected=dict(input_path=p["file"],input_sha256=sha(p["file"]),models=model_hashes,source_sha256=prediction_source_hashes)
            if path.exists():
                if not sidecar.exists():raise ValueError("Unverified prediction cache cannot be reused")
                saved=json.loads(sidecar.read_text());digest=saved.pop("prediction_sha256")
                if saved!=expected or digest!=sha(path):raise ValueError("Prediction cache provenance changed")
                q=np.load(path);pred={k:q[k] for k in models}
            else:
                pred={k:predict(m,g,k) for k,m in models.items()};np.savez_compressed(path,**pred)
                write_json(sidecar,dict(expected,prediction_sha256=sha(path)))
            long,short,_=availability_masks(g)
            np.testing.assert_array_equal(np.isfinite(pred["M120"]),long);np.testing.assert_array_equal(np.isfinite(pred["M15"]),short)
            identities[str(path)]=dict(prediction_sha256=sha(path),provenance_sha256=sha(sidecar),**expected)
            patient=dict(glucose=g,onsets=onsets,predictions=pred)
            for name,policy in POLICIES.items():
                s=score_policy(patient,policy)
                row=dict(cohort=cohort,id=p["id"],arm=p["arm"],policy=name,**{k:s[k] for k in ["events","opportunities","matched","alarms","unconfirmed","observed_days"]});rows.append(row)
                old_name={"M15_004":"M15","M120_004":"M120","switch_004_001":"switch"}.get(name)
                if old_name:
                    old=original.loc[(original.id==p["id"])&(original.policy==old_name)].iloc[0]
                    for k in ["events","opportunities","matched","alarms","unconfirmed","observed_days"]:
                        if not np.isclose(row[k],old[k],rtol=0,atol=1e-9):raise AssertionError("Frozen original replay changed")
            print(f"fixed_ablation {cohort} participant={p['id']}",flush=True)
    frame=pd.DataFrame(rows);frame.to_csv(OUT/"fixed_ablation_participants.csv",index=False)
    summary=[];contrasts=[]
    for (cohort,name),sub in frame.groupby(["cohort","policy"]):
        s=sub[["events","matched","opportunities","alarms","unconfirmed","observed_days"]].sum()
        summary.append(dict(cohort=cohort,policy=name,**s.to_dict(),sensitivity=s.matched/s.events,unconfirmed_per_day=s.unconfirmed/s.observed_days,total_per_day=s.alarms/s.observed_days))
    pd.DataFrame(summary).to_csv(OUT/"fixed_ablation_summary.csv",index=False)
    for cohort in ["internal","external"]:
        for a_name,b_name,label in CONTRASTS:
            a=frame.loc[(frame.cohort==cohort)&(frame.policy==a_name)].set_index("id");b=frame.loc[(frame.cohort==cohort)&(frame.policy==b_name)].set_index("id").reindex(a.index)
            data=np.column_stack([a.events,a.matched,b.matched,a.unconfirmed,b.unconfirmed,a.observed_days]);ci=bootstrap_differences(data,10000)
            contrasts.append(dict(cohort=cohort,a=a_name,b=b_name,contrast=label,sensitivity_difference=(a.matched.sum()-b.matched.sum())/a.events.sum(),unconfirmed_difference=(a.unconfirmed.sum()-b.unconfirmed.sum())/a.observed_days.sum(),**ci))
    write_json(OUT/"fixed_ablation_contrasts.json",contrasts)
    write_json(OUT/"fixed_prediction_provenance.json",dict(models={k:sha(ORIGINAL/"models"/f"{k}.pkl") for k in models},predictions=identities,original_replay_reconciled=True))

def bootstrap_differences(data,n):
    rng=np.random.default_rng(SEED);draws=rng.integers(0,len(data),size=(n,len(data)));s=data[draws].sum(axis=1);valid=(s[:,0]>0)&(s[:,5]>0)
    d=(s[valid,1]-s[valid,2])/s[valid,0];u=(s[valid,3]-s[valid,4])/s[valid,5]
    return dict(replicates=n,undefined=int((~valid).sum()),sensitivity_difference_ci=np.quantile(d,[.025,.975]).tolist(),unconfirmed_difference_ci=np.quantile(u,[.025,.975]).tolist())

def select_indices(counts,days,spec):
    matched=counts[:,0];alarms=counts[:,1];u=alarms-matched
    kinds=np.array([r["kind"] for r in spec]);cut=np.array([r["threshold"] for r in spec]);short=np.array([r.get("short_threshold",0) for r in spec]);physical=np.where((kinds=="linear")&np.isfinite(cut),-cut,cut)
    feasible=(u/days<=1)&(alarms/days<=2)
    def best(mask):
        ix=np.flatnonzero(mask)
        if not len(ix):raise AssertionError("Never-alarm feasibility missing")
        order=np.lexsort((ix,-short[ix],-physical[ix],alarms[ix],u[ix],-matched[ix]));return int(ix[order[0]])
    family={k:best(feasible&(kinds==k)) for k in ["M120","M15","linear","missing","switch"]}
    comp=best(feasible&(kinds!="switch"))
    switch=best(feasible&(kinds=="switch")&(u/days<=u[comp]/days+.10+1e-12))
    priority={"linear":0,"M15":1,"M120":2,"missing":3,"switch":4}
    primary003=min(family.values(),key=lambda i:(-matched[i],u[i],alarms[i],priority[spec[i]["kind"]]))
    return comp,switch,primary003,family

def stability():
    rows=json.loads((ORIGINAL/"participant_manifest.json").read_text())["participants"];rows=[p for p in rows if p["split"]=="validation"]
    provenance=json.loads((ORIGINAL/"validation_prediction_manifest.json").read_text())
    for category in ["inputs","predictions"]:
        for p,h in provenance[category].items():
            if sha(p)!=h:raise AssertionError("Original validation artifact changed")
    spec=specifications();kinds=np.array([["M120","M15","linear","missing","switch"].index(s["kind"]) for s in spec]);cuts=np.array([s["threshold"] for s in spec]);short=np.array([s.get("short_threshold",np.nan) for s in spec]);tallies=[];days=[];events=[]
    for p in rows:
        z=np.load(p["file"]);q=np.load(ORIGINAL/"validation_predictions"/f"DCLP3_{p['id']}.npz");stack=np.vstack([q[k] for k in ["M120","M15","linear","missing"]]);tallies.append(fast_grid(z["glucose"],q["onsets"],stack,kinds,cuts,short));days.append(np.isfinite(z["glucose"]).sum()/288);events.append(len(q["onsets"]))
    tallies=np.asarray(tallies);days=np.asarray(days);events=np.asarray(events);np.savez_compressed(OUT/"validation_person_full_threshold_tallies.npz",counts=tallies,days=days,events=events,ids=[p["id"] for p in rows])
    totals=tallies.sum(axis=0);frozen=json.loads((ORIGINAL/"calibration.json").read_text());np.testing.assert_array_equal(totals,np.array([[r["score"]["matched"],r["score"]["alarms"]] for r in frozen["validation_grid"]]))
    original=select_indices(totals,days.sum(),spec)
    from card003_policy import select_primary
    third_frozen=select_primary(frozen["validation_grid"])["primary"]
    if CLINICAL is not None:
        clinical_frozen=json.loads((CLINICAL/"models"/"frozen_policy.json").read_text())["choices"]["primary"]
        for key in ["kind","threshold","short_threshold"]:
            if clinical_frozen.get(key)!=third_frozen.get(key):raise AssertionError("Clinical frozen policy differs from independently reconstructed validation selection")
    for idx,old in zip(original[:3],[frozen["comparator"],frozen["switch"],third_frozen]):
        if any(old.get(k)!=v for k,v in spec[idx].items()):raise AssertionError("Resampling selector differs from original frozen choice")
    rng=np.random.default_rng(SEED);draws=rng.integers(0,len(rows),size=(1000,len(rows)));weights=np.array([np.bincount(d,minlength=len(rows)) for d in draws]);np.savez_compressed(OUT/"validation_selection_draws.npz",indices=draws,weights=weights)
    records=[];families=[];undefined=0
    def identity(i):
        r=spec[i];label=lambda v:f"{v:g}" if np.isfinite(v) else "never";return f"{r['kind']}:{label(r['threshold'])}"+(f":{label(r['short_threshold'])}" if "short_threshold" in r else "")
    for j,w in enumerate(weights):
        e=int(w@events);d=float(w@days)
        if e<=0 or d<=0:undefined+=1;continue
        count=np.tensordot(w,tallies,axes=(0,0));comp,sw,third,fam=select_indices(count,d,spec)
        records.append(dict(draw=j,events=e,observed_days=d,comparator=identity(comp),comparator_family=spec[comp]["kind"],switch=identity(sw),card003_primary=identity(third),card003_family=spec[third]["kind"],same_comparator=comp==original[0],same_switch=sw==original[1],same_joint=(comp,sw)==original[:2],same_card003=third==original[2]))
        families.extend(dict(draw=j,family=k,selected=identity(i)) for k,i in fam.items())
    frame=pd.DataFrame(records);frame.to_csv(OUT/"selection_stability_draws.csv",index=False);pd.DataFrame(families).to_csv(OUT/"selection_stability_family_thresholds.csv",index=False)
    freq=[]
    for key in ["comparator","comparator_family","switch","card003_primary","card003_family"]:
        for choice,n in frame[key].value_counts().items():freq.append(dict(selection=key,choice=choice,count=int(n),fraction=n/len(frame)))
    pd.DataFrame(freq).to_csv(OUT/"selection_stability_frequencies.csv",index=False)
    joint=frame.groupby(["comparator","switch"]).size().rename("count").reset_index();joint["fraction"]=joint["count"]/len(frame);joint.to_csv(OUT/"selection_stability_joint_choices.csv",index=False)
    write_json(OUT/"selection_stability_summary.json",dict(seed=SEED,requested_draws=1000,valid_draws=len(frame),undefined=undefined,participants=len(rows),operating_points=len(spec),original_comparator=identity(original[0]),original_switch=identity(original[1]),original_card003=identity(original[2]),same_comparator=float(frame.same_comparator.mean()),same_switch=float(frame.same_switch.mean()),same_joint=float(frame.same_joint.mean()),same_card003=float(frame.same_card003.mean()),scope="Conditional validation-person selection stability given originally trained models. No external choice, no refitting, no corrected confirmatory inference."))

def demographics():
    def table(archive,name):
        data=zipfile.ZipFile(archive).read(name);enc="utf-16" if data[:2] in [b"\xff\xfe",b"\xfe\xff"] else "utf-8-sig"
        return pd.read_csv(io.BytesIO(data),sep="|",encoding=enc)
    persons=[];summary=[];sources=[]
    for cohort in ["DCLP3_all","internal","external"]:
        dclp=cohort!="external";source_root=ORIGINAL if dclp else ORIGINAL/"external"/"replacebg"
        archive=Path(json.loads((source_root/"participant_manifest.json").read_text())["archive"])
        screen_name="Data Files/DiabScreening_a.txt" if dclp else "Data Tables/HScreening.txt";roster_name="Data Files/PtRoster_a.txt" if dclp else "Data Tables/HPtRoster.txt";a1c_name="Data Files/DiabLocalHbA1c_a.txt" if dclp else "Data Tables/HLocalHbA1c.txt"
        screen=table(archive,screen_name);roster=table(archive,roster_name);a1c=table(archive,a1c_name);a1c=a1c.loc[a1c.Visit.eq("Screening" if dclp else "Screening Visit")]
        if screen.PtID.duplicated().any() or roster.PtID.duplicated().any() or a1c.PtID.duplicated().any():raise ValueError("Demographic source not unique per person")
        screen=screen.set_index("PtID");roster=roster.set_index("PtID");a1c=a1c.set_index("PtID")
        if cohort=="DCLP3_all":
            population=json.loads((ORIGINAL/"participant_manifest.json").read_text())["participants"];source_scores=[]
            for p in population:
                g=np.load(p["file"])["glucose"];onsets,censored=episodes(g)
                source_scores.append(dict(id=p["id"],observed_days=np.isfinite(g).sum()/288,calendar_days=len(g)/288,events=len(onsets),censored_low_runs=len(censored),matched=np.nan))
            scores=pd.DataFrame(source_scores).set_index("id")
        else:
            population=cohort_manifest(cohort);scores=pd.read_csv(ORIGINAL/"results"/f"{cohort}_participant_scores.csv");scores=scores.loc[scores.policy.eq("M15")].set_index("id")
        for p in population:
            pid=p["id"];s=screen.loc[pid];r=roster.loc[pid];v=scores.loc[pid];age=pd.to_numeric(s["AgeAtEnrollment"] if dclp else r["AgeAsOfEnrollDt"],errors="coerce");diag=pd.to_numeric(s.DiagAge,errors="coerce");duration=age-diag
            persons.append(dict(cohort=cohort,id=pid,split=p["split"],arm=p["arm"],age_at_enrollment=age,recorded_gender=s.Gender,recorded_race=s.Race,recorded_ethnicity=s.Ethnicity,diabetes_duration_years=duration if duration>=0 else np.nan,screening_local_hba1c=pd.to_numeric(a1c.loc[pid,"HbA1cTestRes"],errors="coerce") if pid in a1c.index else np.nan,observed_days=v.observed_days,calendar_days=v.calendar_days,observed_fraction=v.observed_days/v.calendar_days,events=v.events,event_rate_per_observed_day=v.events/v.observed_days if v.observed_days>0 else np.nan,censored_low_runs=v.censored_low_runs,patient_sensitivity=v.matched/v.events if v.events>0 else np.nan))
        sources.append(dict(cohort=cohort,archive=str(archive),sha256=sha(archive),screening_table=screen_name,roster_table=roster_name,hba1c_table=a1c_name,hba1c_visit="Screening" if dclp else "Screening Visit",duration_definition="enrollment age minus reported diagnosis age; negative values unavailable"))
    frame=pd.DataFrame(persons);frame.to_csv(OUT/"cohort_characteristics_by_person.csv",index=False)
    for cohort,part in frame.groupby("cohort"):
        for name in ["age_at_enrollment","diabetes_duration_years","screening_local_hba1c","observed_days","observed_fraction","events","event_rate_per_observed_day","patient_sensitivity"]:
            if cohort=="DCLP3_all" and name=="patient_sensitivity":continue
            x=part[name].dropna();summary.append(dict(cohort=cohort,measure=name,n=len(x),missing=len(part)-len(x),mean=x.mean(),sd=x.std(),median=x.median(),q25=x.quantile(.25),q75=x.quantile(.75),minimum=x.min(),maximum=x.max()))
    pd.DataFrame(summary).to_csv(OUT/"cohort_characteristics_continuous.csv",index=False)
    categorical=[]
    for cohort,part in frame.groupby("cohort"):
        for field in ["recorded_gender","recorded_race","recorded_ethnicity","arm"]:
            for category,n in part[field].fillna("MISSING").value_counts().items():categorical.append(dict(cohort=cohort,field=field,category=category,n=int(n),denominator=len(part),fraction=n/len(part)))
    pd.DataFrame(categorical).to_csv(OUT/"cohort_characteristics_categorical.csv",index=False);write_json(OUT/"demographic_source_provenance.json",sources)
    full=frame.loc[frame.cohort=="DCLP3_all"];split_rows=[]
    for label,p in [("all",full)]+list(full.groupby("split")):
        split_rows.append(dict(split=label,allocated_people=len(p),observed_people=int((p.observed_days>0).sum()),event_free_people=int((p.events==0).sum()),events=int(p.events.sum()),censored_low_runs=int(p.censored_low_runs.sum()),observed_days=p.observed_days.sum(),calendar_days=p.calendar_days.sum(),aggregate_observed_fraction=p.observed_days.sum()/p.calendar_days.sum(),events_per_100_observed_days=100*p.events.sum()/p.observed_days.sum()))
    pd.DataFrame(split_rows).to_csv(OUT/"DCLP3_full_and_split_coverage.csv",index=False)
    summaries=[];alarm=[];person_policy=[]
    for cohort in ["internal","external"]:
        p=frame.loc[frame.cohort==cohort];score=pd.read_csv(ORIGINAL/"results"/f"{cohort}_summary.csv");ledger=pd.read_csv(ORIGINAL/"results"/f"{cohort}_alarm_ledger.csv");individual=pd.read_csv(ORIGINAL/"results"/f"{cohort}_participant_scores.csv")
        for policy in ["M120","M15","switch"]:
            s=score.loc[score.policy==policy].iloc[0]
            people=individual.loc[individual.policy==policy];positive=people.loc[people.events>0];sensitivity=positive.matched/positive.events
            person_policy.append(dict(cohort=cohort,policy=policy,people=len(people),people_with_events=len(positive),event_free_people=int((people.events==0).sum()),episode_weighted_sensitivity=s.sensitivity,person_mean_sensitivity=sensitivity.mean(),person_median_sensitivity=sensitivity.median(),person_q25_sensitivity=sensitivity.quantile(.25),person_q75_sensitivity=sensitivity.quantile(.75)))
            summaries.append(dict(cohort=cohort,policy=policy,allocated_people=len(p),observed_people=int((p.observed_days>0).sum()),event_free_people=int((p.events==0).sum()),events=int(p.events.sum()),censored_low_runs=int(p.censored_low_runs.sum()),observed_days=p.observed_days.sum(),calendar_days=p.calendar_days.sum(),aggregate_observed_fraction=p.observed_days.sum()/p.calendar_days.sum(),events_per_observed_day=p.events.sum()/p.observed_days.sum(),episode_weighted_sensitivity=s.sensitivity))
            for category in ["matched","duplicate","late","unobservable_outcome","unmatched_observed"]:alarm.append(dict(cohort=cohort,policy=policy,category=category,n=int(((ledger.policy==policy)&(ledger.category==category)).sum())))
    pd.DataFrame(summaries).to_csv(OUT/"source_coverage_event_summary.csv",index=False);pd.DataFrame(alarm).to_csv(OUT/"alarm_category_summary.csv",index=False);pd.DataFrame(person_policy).to_csv(OUT/"event_weighted_and_person_sensitivity.csv",index=False)

def figure_and_report():
    required=["availability_proof.json","paired_event_discordance.csv","fixed_ablation_contrasts.json","fixed_prediction_provenance.json","selection_stability_summary.json","cohort_characteristics_continuous.csv","DCLP3_full_and_split_coverage.csv","event_weighted_and_person_sensitivity.csv"]
    for name in required:
        if not (OUT/name).is_file():raise FileNotFoundError(f"Required revision stage missing: {name}")
    os.environ.setdefault("MPLCONFIGDIR",str(OUT/"matplotlib_cache"));import matplotlib
    matplotlib.use("Agg");import matplotlib.pyplot as plt
    pair=pd.read_csv(OUT/"paired_event_discordance.csv");plt.rcParams.update({"font.size":11.5,"axes.titlesize":12})
    fig=plt.figure(figsize=(10,8.4),layout="constrained");layout=fig.add_gridspec(2,2,height_ratios=[1.5,1]);axes=[fig.add_subplot(layout[0,:]),fig.add_subplot(layout[1,0]),fig.add_subplot(layout[1,1])]
    rows=pair.loc[pair.stratum!="all"].reset_index(drop=True);y=np.arange(len(rows));left=-100*rows.M15_only/rows.events;right=100*rows.switch_only/rows.events
    axes[0].barh(y,left,color="#8d4b35",label="M15 only");axes[0].barh(y,right,color="#28728e",label="Switch only");axes[0].axvline(0,color="black",lw=.7)
    axes[0].set(yticks=y,yticklabels=[("DCLP3" if r.cohort=="internal" else "REPLACE-BG")+" / "+("no long opportunity" if r.stratum=="no_long_opportunity" else "long opportunity")+f"\n(n={r.events})" for r in rows.itertuples()],xlabel="Discordant warnings (% of stratum events)",title="Paired event discordance")
    axes[0].legend(frameon=False,loc="lower right",ncol=2);axes[0].invert_yaxis();axes[0].set_xlim(-23,42)
    for i,r in rows.iterrows():axes[0].text(left.iloc[i]-.4,i,str(r.M15_only),ha="right",va="center",fontsize=11.5);axes[0].text(right.iloc[i]+.4,i,str(r.switch_only),ha="left",va="center",fontsize=11.5)
    for j,cohort in enumerate(["internal","external"]):
        c=json.loads((ORIGINAL/"results"/f"{cohort}_comparison.json").read_text())
        for ax,key,scale in [(axes[1],"sensitivity_difference",100),(axes[2],"unconfirmed_difference",1)]:
            point=c[key]*scale;ci=np.array(c[key+"_ci"])*scale;ax.errorbar(point,j,xerr=[[point-ci[0]],[ci[1]-point]],fmt="o",color="#183b56",capsize=4)
    for ax in axes[1:]:ax.axvline(0,color="black",lw=.7);ax.set(yticks=[0,1],yticklabels=["DCLP3","REPLACE-BG"],ylim=(-.7,1.7));ax.invert_yaxis()
    axes[1].axvline(5,color="#777",ls="--",lw=.8);axes[1].set(title="Frozen sensitivity contrast",xlabel="Switch − M15\n(percentage points)")
    axes[2].axvline(.1,color="#777",ls="--",lw=.8);axes[2].set(title="Frozen alarm-burden contrast",xlabel="Switch − M15\nUnconfirmed alarms / observed day",xticks=[-.05,0,.05,.10])
    for ax in axes:ax.spines[["top","right"]].set_visible(False)
    fig.suptitle("History switching: no added availability over M15,\nfewer useful warnings overall",fontsize=12.5)
    fig.savefig(OUT/"principal_warning_diagnostics.png",dpi=180);plt.close(fig)
    stability=json.loads((OUT/"selection_stability_summary.json").read_text());abl=pd.read_csv(OUT/"fixed_ablation_summary.csv");availability=pd.read_csv(OUT/"availability_containment_summary.csv")
    family=pd.read_csv(OUT/"selection_stability_family_thresholds.csv");frequency=family.groupby(["family","selected"]).size().rename("count").reset_index();frequency["fraction"]=frequency["count"]/stability["valid_draws"];frequency.to_csv(OUT/"selection_stability_within_family_frequencies.csv",index=False)
    lines=["# Reconstructed exploratory warning diagnostics", "",
        "All numeric results below were reconstructed from the caller-supplied source archives and the newly fitted, frozen warning models. These diagnostics use the published post-primary exploratory protocol; they do not retune models or operating points.", "",
        "## Availability and paired outcomes", "", protocol()["availability"], "",
        "|Cohort|Intervals|People|Events|Long opportunity|Short/switch opportunity|",
        "|---|---:|---:|---:|---:|---:|"]
    for r in availability.itertuples():lines.append(f"|{r.cohort}|{r.intervals}|{r.participants}|{r.events}|{r.long_opportunity_events}|{r.short_opportunity_events}|")
    if CLINICAL is None:lines += ["", "CLVer was not supplied to this warning-only run. Pass --clinical-root after its independent reconstruction to include its containment check."]
    lines += ["", "|Cohort|Stratum|Events|Both|M15 only|Switch only|Neither|Net switch−M15|", "|---|---|---:|---:|---:|---:|---:|---:|"]
    for r in pair.itertuples():lines.append(f"|{r.cohort}|{r.stratum}|{r.events}|{r.both}|{r.M15_only}|{r.switch_only}|{r.neither}|{r.net_switch_minus_M15:+d}|")
    lines += ["", "## Fixed-threshold replay", "", "|Cohort|Policy|Sensitivity|Unconfirmed/day|All alarms/day|", "|---|---|---:|---:|---:|"]
    for r in abl.itertuples():lines.append(f"|{r.cohort}|{r.policy}|{r.sensitivity:.6f}|{r.unconfirmed_per_day:.6f}|{r.total_per_day:.6f}|")
    lines += ["", "The equal-threshold switch versus M15 contrast changes model choice at common numeric thresholds. Changing the fallback threshold preserves model choice. Refractory-state interactions remain in each replay; equal numeric thresholds do not equalize calibration or burden. Full 10,000-resample contrast intervals are in fixed_ablation_contrasts.json.", "", "## Conditional selection stability", "",
        f"Of {stability['requested_draws']} validation-person resamples (seed {SEED}), {stability['valid_draws']} were defined. The original comparator setting recurred in {stability['same_comparator']:.1%}, switch pair in {stability['same_switch']:.1%}, their joint choice in {stability['same_joint']:.1%}, and the separately specified clinical primary setting in {stability['same_card003']:.1%}.", "", stability["scope"], "",
        f"All {stability['participants']} allocated validation participants were resampled, including zero-observation participants. Tallies for all {stability['operating_points']} operating points reproduce the complete validation ledger before resampling. No external outcome enters selection.", "", "## Cohort summaries and interpretation", "",
        "Enrollment demographics, screening HbA1c, derived duration, coverage, event frequency, alarm categories, censored-low counts and person-level summaries are saved with source-table provenance. Demographic denominators use allocated participants. Individual sensitivity is undefined for event-free participants; their burden remains included.", "",
        "The principal figure uses paired discordant-event counts and the frozen primary participant-bootstrap intervals. Dashed lines mark the +5 percentage-point sensitivity target and +0.10/day unconfirmed-alarm margin. The linear comparator's selected never-alarm sentinel is omitted from this figure. No finite linear-grid point was feasible in the original study; reconstructed selections are retained in calibration.json.", "",
        "These intervals are conditional on fitted models. Actual receipt latency and physiology during unobserved intervals remain unknown. Cross-cohort differences are descriptive and do not identify sensitivity deterioration caused by transport.", ""]
    (OUT/"revision_warning_report.md").write_text("\n".join(lines),encoding="utf-8")
    write_json(OUT/"revision_completion.json",dict(status="COMPLETE_EXPLORATORY_DIAGNOSTICS",completed_utc=datetime.now(timezone.utc).isoformat(),original_results_modified=False,source_sha256=sha(Path(__file__)),artifacts={p.name:sha(p) for p in OUT.glob("*") if p.is_file() and p.name!="revision_completion.json"}))

def main():
    global ORIGINAL,OUT,CLINICAL
    p=argparse.ArgumentParser();p.add_argument("--warning-root",type=Path,required=True);p.add_argument("--output-root",type=Path,required=True);p.add_argument("--clinical-root",type=Path);p.add_argument("--stage",choices=["freeze","availability","paired","ablation","stability","demographics","report","all"],default="all");a=p.parse_args()
    ORIGINAL=a.warning_root.resolve();OUT=a.output_root.resolve();CLINICAL=a.clinical_root.resolve() if a.clinical_root else None
    freeze()
    jobs={"availability":availability,"paired":paired,"ablation":ablation,"stability":stability,"demographics":demographics,"report":figure_and_report}
    if a.stage=="all":
        for job in jobs.values():job()
    elif a.stage in jobs:jobs[a.stage]()

if __name__=="__main__":main()
