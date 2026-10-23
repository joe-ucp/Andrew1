"""Staged warning-study execution in a caller-selected research directory."""
from __future__ import annotations
import argparse,hashlib,json,os,pickle,time,sys,shutil
from pathlib import Path
import numpy as np
from card001_core import *
from card001_models import fit_model,predict


def load_participants(root,splits):
    manifest=json.loads((root/"participant_manifest.json").read_text())
    result=[]
    for row in manifest["participants"]:
        if row["split"] not in splits: continue
        a=np.load(row["file"])
        result.append(dict(row,glucose=a["glucose"],time=a["time"]))
    return result


def train(root):
    started=time.perf_counter(); out=root/"models"; out.mkdir(exist_ok=True)
    participants=load_participants(root,["train"])
    import sklearn,scipy,pandas
    sources=[Path(__file__).parent/name for name in ["card001_core.py","card001_models.py","card001_data.py"]]
    signature=dict(manifest_sha256=hashlib.sha256((root/"participant_manifest.json").read_bytes()).hexdigest(),code_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in sources},versions={"numpy":np.__version__,"sklearn":sklearn.__version__,"scipy":scipy.__version__,"pandas":pandas.__version__},seed=SEED,learner_grid={"max_leaf_nodes":[7,15],"l2_regularization":[0,1],"learning_rate":.05,"max_iter":200,"min_samples_leaf":50,"early_stopping":False},folds="Card003 sorted sha256 hash cyclic 5-folds, all data from person together",status="TRAINING_NO_TEST_OUTCOMES_OPENED")
    commitment=out/"training_commitment.json"
    if commitment.exists():
        old=json.loads(commitment.read_text())
        if old!=signature:
            raise ValueError("Training source/code commitment changed; refuse to reuse or relabel models. Use a fresh versioned output directory.")
    else:
        if list(out.glob("*.pkl")): raise ValueError("Models exist without original commitment; cannot establish training provenance")
        commitment.write_text(json.dumps(signature,indent=2))
        snapshot=out/"source_snapshot"; snapshot.mkdir(exist_ok=True)
        for p in sources: shutil.copy2(p,snapshot/p.name)
    for kind in ["M120","M15","missing"]:
        if (out/f"{kind}.pkl").exists():
            print(f"existing_model={kind}; use separate output directory for protocol changes",flush=True); continue
        model,report=fit_model(participants,kind)
        with (out/f"{kind}.pkl").open("wb") as f: pickle.dump(model,f)
        (out/f"{kind}_fit.json").write_text(json.dumps(report,indent=2))
        print(f"completed_model={kind} cumulative_seconds={time.perf_counter()-started:.1f}",flush=True)
    predict_validation(root)


def predict_validation(root):
    if (root/"calibration.json").exists():
        frozen=json.loads((root/"calibration.json").read_text())
        path=root/"validation_prediction_manifest.json"
        if hashlib.sha256(path.read_bytes()).hexdigest()!=frozen["prediction_manifest_sha256"]: raise ValueError("Frozen validation prediction manifest changed")
        provenance=json.loads(path.read_text())
        for section in ["inputs","predictions"]:
            for name,digest in provenance[section].items():
                if hashlib.sha256(Path(name).read_bytes()).hexdigest()!=digest: raise ValueError("Frozen validation artifact changed")
        for kind,digest in frozen["model_sha256"].items():
            if hashlib.sha256((root/"models"/f"{kind}.pkl").read_bytes()).hexdigest()!=digest: raise ValueError("Frozen model changed")
        print("Existing frozen validation predictions retained",flush=True)
        return
    folder=root/"validation_predictions"; folder.mkdir(exist_ok=True)
    models={kind:pickle.load((root/"models"/f"{kind}.pkl").open("rb")) for kind in ["M120","M15","missing"]}
    provenance=dict(model_sha256={kind:hashlib.sha256((root/"models"/f"{kind}.pkl").read_bytes()).hexdigest() for kind in models},inputs={},predictions={},source_sha256={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in [Path(__file__).parent/"card001_core.py",Path(__file__).parent/"card001_models.py"]})
    for p in load_participants(root,["validation"]):
        predictions={kind:predict(m,p["glucose"],kind) for kind,m in models.items()}
        predictions["linear"]=-linear_forecast(p["glucose"])
        onsets,censored=episodes(p["glucose"])
        onsets2,censored2=episodes_independent(p["glucose"])
        if not np.array_equal(onsets,onsets2) or not np.array_equal(censored,censored2): raise AssertionError("Independent event detector failure")
        destination=folder/f"DCLP3_{p['id']}.npz"
        np.savez_compressed(destination,**predictions,onsets=onsets,censored=censored)
        provenance["inputs"][p["file"]]=hashlib.sha256(Path(p["file"]).read_bytes()).hexdigest()
        provenance["predictions"][str(destination)]=hashlib.sha256(destination.read_bytes()).hexdigest()
    (root/"validation_prediction_manifest.json").write_text(json.dumps(provenance,indent=2))
    print("validation_predictions_ready",flush=True)


def status(root):
    archives=list((root/"raw").glob("*.zip"))
    model_count=len(list((root/"models").glob("*.pkl")))
    state="SOURCE_AVAILABLE" if archives else "ACCESS_BLOCKED"
    if (root/"models"/"training_commitment.json").exists(): state="MODEL_FITTING_OR_INTERRUPTED"
    if model_count==3: state="MODELS_FITTED_CALIBRATION_PENDING"
    if (root/"calibration.json").exists(): state="POLICIES_FROZEN_TEST_PENDING"
    comparisons={}
    for cohort in ["internal","external"]:
        path=root/"results"/f"{cohort}_comparison.json"
        if path.exists(): comparisons[cohort]=json.loads(path.read_text())
    if comparisons: state="INTERNAL_EVALUATED_EXTERNAL_PENDING" if len(comparisons)==1 else "BOTH_PRIMARY_COHORTS_EVALUATED"
    overall="NOT_EVALUATED"
    if comparisons:
        decisions=[r["status"] for r in comparisons.values()]
        if any(d.startswith("REFUTED_") for d in decisions): overall="REFUTED"
        elif len(comparisons)==2 and all(x=="SUPPORTED_IN_THIS_COHORT" for x in decisions): overall="SUPPORTED_PRIMARY_ENDPOINTS_REQUIRES_ROBUSTNESS_REVIEW"
        elif len(comparisons)==2: overall="INCONCLUSIVE_OR_REJECTED_OPERATING_POINT"
        else: overall="INTERNAL_ONLY_EXTERNAL_NOT_EVALUATED"
    payload=dict(card="001",scientific_status=overall,execution_status=state,model_count=model_count,comparisons=comparisons,archives=[str(p) for p in archives],reason="Only archived empirical comparisons determine scientific status; code and synthetic validation do not",updated=time.strftime("%Y-%m-%dT%H:%M:%S%z"))
    (root/"execution_status.json").write_text(json.dumps(payload,indent=2)); print(json.dumps(payload))


def main():
    parser=argparse.ArgumentParser(); parser.add_argument("--data-root",type=Path,required=True); parser.add_argument("--external-root",type=Path); parser.add_argument("--cohort-split"); parser.add_argument("--result-prefix"); parser.add_argument("--stage",choices=["status","train","predict-validation","calibrate","test","sensitivity"],default="status")
    parser.add_argument("--workers",type=int,choices=range(1,9),default=1)
    args=parser.parse_args(); root=args.data_root
    root.mkdir(parents=True,exist_ok=True)
    if args.stage=="train": train(root)
    elif args.stage=="predict-validation": predict_validation(root)
    elif args.stage=="status": status(root)
    else:
        from card001_evaluate import calibrate_run,test_run
        if args.stage=="calibrate": calibrate_run(root)
        elif args.stage=="test": test_run(root,args.external_root,args.cohort_split,args.result_prefix)
        else:
            from card001_sensitivity import run_sensitivity
            run_sensitivity(root,args.external_root,args.workers)

if __name__=="__main__": main()
