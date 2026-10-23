"""Verify refractory synchronization; extend secondary policy warmup as needed."""
from pathlib import Path
import json,argparse
import numpy as np
import pandas as pd
from card003_analyze import load_frozen,predict_grid,policy_arrays,fit_and_save,transfer_table,dump
from card003_data import window_grid
from card003_core import replay_interval,sha256_file


def synchronized(t,g,p,cut,start,known_from):
    """Six fully determined no-candidate bins erase any prior 30-minute timer."""
    candidate=np.isfinite(g)&(g>=54)&np.isfinite(p)&(p>=cut)
    before=np.flatnonzero(t<=start-1800)
    end=before[-1] if len(before) else -1
    quiet=~candidate[known_from:end+1]
    return bool(len(quiet)>=6 and (np.convolve(quiet.astype(int),np.ones(6,dtype=int),"valid")>=6).any())


def audit_primary(root,frame,frozen):
    variants=["primary","phase1","phase2","phase3","phase4","delay5","run3","blank_units_excluded","date_only_midnight_assumption"]
    choice=frozen["choices"]["primary"];kind=choice["kind"];rows=[]
    for variant in variants:
        for row in frame.to_dict("records"):
            path=root/"processed"/"predictions"/variant/f"{row['pair_id']}.npz"
            if not path.exists():raise ValueError("Required variant prediction is incomplete: "+str(path))
            with np.load(path) as cache:
                p,cut,_=policy_arrays(cache,choice)
                known=24 if kind in ("M120","missing","switch") else 3
                if variant=="delay5":known+=1
                ok=synchronized(cache["time"],cache["glucose"],p,cut,row["start_seconds"],known)
            rows.append(dict(pair_id=row["pair_id"],policy=kind,variant=variant,synchronized=ok))
    table=pd.DataFrame(rows);table.to_csv(root/"results"/"primary_refractory_synchronization.csv",index=False)
    if not table.synchronized.all():raise ValueError("A primary-policy variant requires more refractory context; extend before completion")
    return dict(interval_variants=len(table),all_synchronized=True)


def extend_families(root):
    root=Path(root);frozen,models,digest=load_frozen(root)
    frame=pd.read_csv(root/"processed"/"intervals_all_coverage.csv",dtype={"PtID":str})
    folder=root/"results"/"boundary_extension";folder.mkdir(exist_ok=True)
    choices={k:v for k,v in frozen["choices"]["frozen_family"].items() if k in ("M120","missing","switch")}
    rows={k:[] for k in choices};events={k:[] for k in choices};alarms={k:[] for k in choices};audit=[]
    cachefolder=root/"processed"/"predictions"/"extended_boundary";cachefolder.mkdir(parents=True,exist_ok=True)
    for j,row in enumerate(frame.to_dict("records")):
        for lookback in (360,1440,4320,20160):
            grid=window_grid(root,row["PtID"],row["start_seconds"],row["end_seconds"],lookback_minutes=lookback)
            path=cachefolder/f"{row['pair_id']}_{lookback}.npz"
            if path.exists():
                with np.load(path) as old:
                    if str(old["policy_hash"])!=digest or not np.array_equal(old["glucose"],grid["glucose"],equal_nan=True):raise ValueError("Stale extended-boundary cache")
                    predictions={k:old[k] for k in ("M120","M15","missing","linear")}
            else:
                predictions=predict_grid(models,grid["glucose"])
                np.savez_compressed(path,**predictions,glucose=grid["glucose"],time=grid["time"],policy_hash=digest)
            proofs={}
            for kind,choice in choices.items():
                p,cut,available=policy_arrays(predictions,choice)
                proofs[kind]=synchronized(grid["time"],grid["glucose"],p,cut,row["start_seconds"],24)
            if all(proofs.values()):break
        if not all(proofs.values()):raise ValueError("Even14-day warmup did not synchronize; use full-source replay")
        for kind,choice in choices.items():
            p,cut,available=policy_arrays(predictions,choice)
            metrics,ee,aa=replay_interval(grid["time"],grid["glucose"],p,cut,available,row["start_seconds"],row["end_seconds"])
            if metrics["events"]!=row["events"]:raise AssertionError("Extended history changes event denominator; audit primary event boundaries")
            rows[kind].append({**row,**metrics,"policy":kind,"variant":"extended_boundary"})
            events[kind].extend(ee.assign(PtID=row["PtID"],pair_id=row["pair_id"],policy=kind).to_dict("records"))
            alarms[kind].extend(aa.assign(PtID=row["PtID"],pair_id=row["pair_id"],policy=kind).to_dict("records"))
            audit.append(dict(pair_id=row["pair_id"],policy=kind,lookback_minutes=lookback,synchronized=proofs[kind]))
        if (j+1)%50==0:print(f"extended policy boundary: {j+1}/{len(frame)}",flush=True)
    result={}
    for kind,data in rows.items():
        table=pd.DataFrame(data);table.to_csv(folder/f"intervals_{kind}.csv",index=False)
        pd.DataFrame(events[kind]).to_csv(folder/f"events_{kind}.csv",index=False)
        pd.DataFrame(alarms[kind]).to_csv(folder/f"alarms_{kind}.csv",index=False)
        primary=table[table.coverage>=.7]
        r=fit_and_save(root,"extended_boundary_"+kind,primary,drop_constant_nuisance=True)
        r["transfer"]=transfer_table(primary)
        r["boundary_provenance"]="At least six no-candidate bins after full feature history and before earliest scored warning; initial timer cannot affect scored alarms"
        result["policy_"+kind]=r
    pd.DataFrame(audit).to_csv(folder/"synchronization_audit.csv",index=False)
    dump(folder/"corrected_policy_results.json",result)
    return result


def finalize(root):
    root=Path(root)
    frozen,_,_=load_frozen(root)
    frame=pd.read_csv(root/"processed"/"intervals_all_coverage.csv",dtype={"PtID":str})
    primaryproof=audit_primary(root,frame,frozen)
    corrected=json.loads((root/"results"/"boundary_extension"/"corrected_policy_results.json").read_text())
    for key in corrected:
        kind=key.removeprefix("policy_")
        extended=pd.read_csv(root/"results"/"boundary_extension"/f"intervals_{kind}.csv")
        comparison=frame[["pair_id","events"]].merge(extended[["pair_id","events"]],on="pair_id",suffixes=("_original","_extended"),validate="one_to_one")
        if len(comparison)!=len(frame) or not comparison.events_original.eq(comparison.events_extended).all():
            raise AssertionError("Extended history changed original event denominator")
    path=root/"results"/"sensitivity_results.json"
    old=json.loads(path.read_text());archive=root/"history"/"before_secondary_boundary_correction.json"
    if not archive.exists():archive.write_bytes(path.read_bytes())
    differences=[]
    for key,result in corrected.items():
        differences.append(dict(policy=key,old_rr=old[key].get("rr"),corrected_rr=result.get("rr"),old_unwarned=old[key]["transfer"]["unwarned"],corrected_unwarned=result["transfer"]["unwarned"]))
        old[key]=result
    dump(path,old)
    rows=[]
    for name,r in old.items():
        ci=r.get("robust_95_ci",[None,None]);rows.append(dict(analysis=name,children=r.get("n_children"),intervals=r.get("n_intervals"),events=r.get("events"),rr=r.get("rr"),robust_lower=ci[0],robust_upper=ci[1],status=r["status"],reason=r.get("reason")))
    pd.DataFrame(rows).to_csv(root/"results"/"sensitivity_matrix.csv",index=False)
    proof=dict(primary=primaryproof,secondary_changes=differences,status="PASS")
    dump(root/"results"/"boundary_integrity_closure.json",proof)
    return proof


if __name__=="__main__":
    ap=argparse.ArgumentParser();ap.add_argument("--data-root",type=Path,required=True);ap.add_argument("--finalize-only",action="store_true");a=ap.parse_args()
    if not a.finalize_only:extend_families(a.data_root)
    if (a.data_root/"results"/"analysis_summary.json").exists():print(json.dumps(finalize(a.data_root)))
