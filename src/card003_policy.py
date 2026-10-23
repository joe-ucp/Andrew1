"""Independent Card003 policy choice from the DCLP3-only validation ledger.

Card001's switch adds a comparator difference constraint. Card003 must select
from the entire capped family without that extra constraint, never from CLVer.
"""
from __future__ import annotations
import hashlib
import math
import json
from datetime import datetime, timezone
from pathlib import Path

PRIORITY = {"linear":0,"M15":1,"M120":2,"missing":3,"switch":4}


def prescribed_folds(participant_ids):
    ids = {str(int(p)) for p in participant_ids}
    ordered = sorted(ids, key=lambda p:hashlib.sha256(("t1d-card-003-fold:"+p).encode()).hexdigest())
    return {p:k%5 for k,p in enumerate(ordered)}


def _performance(row):
    s = row["score"]
    return (-s["sensitivity"],s["unconfirmed_per_day"],s["total_per_day"])


def _threshold_order(row):
    # The shared core represents linear forecast alarms as -forecast >= -cut.
    physical = -row["threshold"] if row["kind"] == "linear" else row["threshold"]
    if not math.isfinite(row["threshold"]):
        physical = math.inf  # explicit never-alarm at an exact zero-event tie
    return (-physical,-row.get("short_threshold",0))


def select_primary(validation_grid):
    """All policies face ≤1 unconfirmed and ≤2 issued alarms/observed day."""
    by_kind = {}
    for kind in PRIORITY:
        candidates = []
        for r in validation_grid:
            if r["kind"] != kind:
                continue
            s = r["score"]
            if s["sensitivity"] is not None and s["unconfirmed_per_day"] is not None and s["unconfirmed_per_day"] <= 1 and s["total_per_day"] <= 2:
                candidates.append(r)
        if not candidates:
            raise ValueError(f"No feasible DCLP3 validation operating point for {kind}")
        by_kind[kind] = min(candidates,key=lambda r:(_performance(r),_threshold_order(r)))
    # Numeric threshold scales differ across model classes; each class already
    # used its threshold tie-break. Exact cross-family performance ties use the
    # explicitly prespecified simplicity priority, not unlike-unit comparison.
    primary = min(by_kind.values(),key=lambda r:(_performance(r),PRIORITY[r["kind"]]))
    return {"primary":primary,"frozen_family":by_kind,
            "dataset":"DCLP3","selection_data":"validation participants only",
            "CLVer_used_for_selection":False,
            "family_tie_priority":list(PRIORITY),
            "switch_additional_margin_constraint":False,
            "threshold_ties":"higher physical threshold within policy; switch M120 then M15",
            "cross_family_threshold_note":"Different score units are not compared; performance ties use the prespecified family priority."}


def freeze_from_card001(source_root, destination_root):
    """Validate provenance, then persist the separate Card003 operating points.

    This function never opens CLVer data or DCLP3 held-out participant arrays.
    Existing differing freeze files cause failure instead of silent replacement.
    """
    from card003_core import sha256_file
    from card001_core import participant_split
    source, dest = Path(source_root), Path(destination_root)
    frozen_path = source / "calibration.json"
    calibration = json.loads(frozen_path.read_text())
    if calibration.get("status") != "FROZEN_BEFORE_TEST":
        raise ValueError("DCLP3 validation artifact is not declared frozen")
    manifest_path = source / "participant_manifest.json"
    if sha256_file(manifest_path) != calibration["manifest_sha256"]:
        raise ValueError("Participant manifest changed after DCLP3 freeze")
    manifest = json.loads(manifest_path.read_text())
    participants = manifest["participants"]
    for p in participants:
        if p["cohort"] != "DCLP3" or p["split"] != participant_split(p["id"])[0]:
            raise ValueError("Policy manifest contains wrong cohort or split")
    training_ids = [p["id"] for p in participants if p["split"] == "train"]
    folds = prescribed_folds(training_ids)
    models, fits = {}, {}
    for kind in ("M120", "M15", "missing"):
        path = source / "models" / f"{kind}.pkl"
        digest = sha256_file(path)
        if digest != calibration["model_sha256"][kind]:
            raise ValueError("Fitted model differs from frozen artifact")
        fit_path = source / "models" / f"{kind}_fit.json"
        fit = json.loads(fit_path.read_text())
        if fit["folds"] != folds:
            raise ValueError("Folds do not match all assigned training participants, including zero-observation people")
        candidates = fit["candidates"]
        if {(p["max_leaf_nodes"],p["l2_regularization"]) for p in candidates} != {(7,0),(7,1),(15,0),(15,1)}:
            raise ValueError("Incomplete prespecified model grid")
        if any(len(p["fold_losses"]) != 5 for p in candidates):
            raise ValueError("Incomplete prespecified grouped cross-validation")
        selected = min(candidates,key=lambda r:(r["weighted_log_loss"],r["max_leaf_nodes"],-r["l2_regularization"]))
        if selected != fit["selected"]:
            raise ValueError("Model hyperparameter selection mismatch")
        models[kind] = {"path":str(path),"sha256":digest}
        fits[kind] = {"path":str(fit_path),"sha256":sha256_file(fit_path)}
    grid = calibration["validation_grid"]
    expected_counts = {"M120":102,"M15":102,"missing":102,"linear":62,"switch":10404}
    if {k:sum(r["kind"]==k for r in grid) for k in expected_counts} != expected_counts:
        raise ValueError("Validation ledger does not contain the complete threshold family")
    probabilities = [x / 100 for x in range(101)] + [math.inf]
    expected_settings = {(k,t,None) for k in ("M120","M15","missing") for t in probabilities}
    expected_settings |= {("linear",float(-t),None) for t in range(40,101)} | {("linear",math.inf,None)}
    expected_settings |= {("switch",a,b) for a in probabilities for b in probabilities}
    actual_settings = {(r["kind"],r["threshold"],r.get("short_threshold")) for r in grid}
    if actual_settings != expected_settings:
        raise ValueError("Threshold setting identities differ from the prespecified grid")
    choices = select_primary(grid)
    commitment_path = source / "models" / "training_commitment.json"
    commitment = json.loads(commitment_path.read_text())
    expected = dict(max_leaf_nodes=[7,15],l2_regularization=[0,1],learning_rate=.05,
                    max_iter=200,min_samples_leaf=50,early_stopping=False)
    if commitment["learner_grid"] != expected or commitment["seed"] != 20261022:
        raise ValueError("Learner settings differ from Card003 commitment")
    if commitment["manifest_sha256"] != calibration["manifest_sha256"]:
        raise ValueError("Training and calibration used different participant manifests")
    for name, digest in commitment["code_sha256"].items():
        if sha256_file(source/"models"/"source_snapshot"/name) != digest:
            raise ValueError("Training source snapshot changed")
    provenance = dict(
        choices=choices,models=models,fit_reports=fits,
        calibration_sha256=sha256_file(frozen_path),calibration_path=str(frozen_path),
        participant_manifest_sha256=calibration["manifest_sha256"],
        training_commitment_sha256=sha256_file(commitment_path),
        versions=commitment["versions"],folds=folds,
        selection_code_sha256=sha256_file(Path(__file__)),
        status="FROZEN_BEFORE_CLVER_REPLAY",CLVer_outcomes_used=False)
    dest.mkdir(parents=True,exist_ok=True)
    out = dest / "frozen_policy.json"
    if out.exists():
        old = json.loads(out.read_text())
        old.pop("frozen_utc",None)
        if old != provenance:
            raise ValueError("Card003 policy already frozen differently; preserve it and explicitly version any reformulation")
    else:
        provenance["frozen_utc"] = datetime.now(timezone.utc).isoformat()
        out.write_text(json.dumps(provenance,indent=2),encoding="utf-8")
    return out
