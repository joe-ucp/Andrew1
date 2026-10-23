"""Scientific kernels for the CLVer secretion and hypoglycemia analysis. No data acquisition or warning fitting.

Small constructed arrays in tests are software fixtures, never patient substitutes.
Real analysis inputs and outputs reside under explicitly configured data roots.
"""
from __future__ import annotations

import hashlib
import json
import math
import warnings
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 20261022
MMTT_TIMES = (0, 15, 30, 60, 90, 120)
EXPECTED_CLVER_SHA256 = "8421a9a9ff7dd1485ce0e4bdb984b5ea037300a0e7216f9f445334e56d32d385"


def sha256_file(path: Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1024 * 1024), b""):
            h.update(b)
    return h.hexdigest()


def normalized_auc(values, times=MMTT_TIMES) -> float:
    """Trapezoidal two-hour mean concentration, in input concentration units."""
    v = np.asarray(values, dtype=float)
    t = np.asarray(times, dtype=float)
    if tuple(t) != MMTT_TIMES or v.shape != (6,):
        raise ValueError("Exactly the six prespecified nominal MMTT times are required")
    if not np.isfinite(v).all() or (v < 0).any():
        raise ValueError("Missing, nonnumeric, infinite or negative MMTT value")
    return float(np.sum((v[:-1] + v[1:]) * np.diff(t) / 2.0) / 120.0)


def clinical_remission(hba1c, daily_units, weight, weight_units,
                       index_date, hba1c_date, insulin_date, exam_date,
                       dose_known=True, weight_known=True):
    """IDAA1c and companion state; missing is never treated as nonremission."""
    result = {"idaa1c": np.nan, "remission": np.nan,
              "insulin_u_kg_day": np.nan, "weight_kg": np.nan, "reason": ""}
    try:
        h, d, w = map(float, (hba1c, daily_units, weight))
    except (ValueError, TypeError):
        result["reason"] = "nonnumeric_clinical_component"
        return result
    if not dose_known or not weight_known:
        result["reason"] = "explicit_unknown_component"
        return result
    unit = str(weight_units).strip().casefold()
    if unit in ("lb", "lbs", "pound", "pounds"):
        w *= 0.45359237
    elif unit not in ("kg", "kilogram", "kilograms"):
        result["reason"] = "unsupported_weight_unit"
        return result
    if not all(np.isfinite([h, d, w])) or h <= 0 or d < 0 or w <= 0:
        result["reason"] = "invalid_clinical_component"
        return result
    dates = [pd.to_datetime(x, errors="coerce") for x in
             (index_date, hba1c_date, insulin_date, exam_date)]
    if any(pd.isna(x) for x in dates):
        result["reason"] = "missing_clinical_date"
        return result
    if any(abs((x.normalize() - dates[0].normalize()).days) > 7 for x in dates[1:]):
        result["reason"] = "clinical_date_outside_seven_days"
        return result
    idaa = h + 4.0 * d / w
    result.update(idaa1c=idaa, remission=int(idaa <= 9), insulin_u_kg_day=d / w,
                  weight_kg=w, reason="known")
    return result


def derive_mmt_visit(samples: pd.DataFrame, *, fasting_codes=("Fasting",),
                     require_temperature_ok=False):
    """Derive one visit from canonical rows; codes supplied from source glossary.

    Columns: minute,value,units,date,status,fasting,temperature. A fasting MMTT
    requires its baseline sample to have the documented fasting code. Challenge
    samples cannot independently establish fasting at challenge initiation.
    Identical analytical duplicates are collapsed; conflicting samples exclude.
    """
    cols = ["minute", "value", "units", "date", "status", "fasting", "temperature"]
    missing = set(cols) - set(samples.columns)
    if missing:
        raise ValueError(f"Missing canonical assay columns: {sorted(missing)}")
    s = samples[cols].drop_duplicates().copy()
    reasons = []
    if len(s) != 6 or sorted(s.minute.tolist()) != list(MMTT_TIMES):
        reasons.append("incomplete_or_conflicting_nominal_samples")
    values = pd.to_numeric(s.value, errors="coerce")
    if not np.isfinite(values).all() or (values < 0).any():
        reasons.append("nonnumeric_or_negative_assay")
    if set(s.units.astype(str).str.strip()) != {"nmol/L"}:
        reasons.append("unsupported_or_inconsistent_assay_unit")
    dates = pd.to_datetime(s.date, errors="coerce")
    if dates.isna().any() or dates.dt.normalize().nunique() != 1:
        reasons.append("missing_or_inconsistent_collection_date")
    if not s.status.astype(str).str.strip().str.casefold().eq("complete").all():
        reasons.append("incomplete_result_status")
    baseline = s.loc[s.minute.eq(0)]
    if len(baseline) != 1 or baseline.iloc[0].fasting not in fasting_codes:
        reasons.append("fasting_not_established")
    if require_temperature_ok and not s.temperature.eq("OK").all():
        reasons.append("temperature_not_ok")
    out = {"eligible": not reasons, "reason": ";".join(reasons) or "eligible",
           "auc_nmol_l": np.nan, "peak_nmol_l": np.nan,
           "collection_date": None, "temperature_all_ok": bool(s.temperature.eq("OK").all())}
    if not reasons:
        ordered = s.assign(value=values).sort_values("minute")
        out.update(auc_nmol_l=normalized_auc(ordered.value),
                   peak_nmol_l=float(ordered.value.max()),
                   collection_date=dates.iloc[0].normalize().isoformat())
        for m, v in zip(ordered.minute, ordered.value):
            out[f"cpep_{m}"] = float(v)
    return out


def later_interval(index_date, target_anchor_dates, target_visit_date):
    """Visit anchor derives only from unique target C-PEP-0 collection dates."""
    target_dates = pd.to_datetime(pd.Series(list(target_anchor_dates)), errors="coerce")
    if target_dates.isna().any() or target_dates.dt.normalize().nunique() != 1:
        return {"eligible": False, "reason": "absent_or_ambiguous_target_anchor"}
    target = target_dates.iloc[0].normalize()
    visit = pd.to_datetime(target_visit_date, errors="coerce")
    index = pd.to_datetime(index_date, errors="coerce")
    if pd.isna(index) or pd.isna(visit) or abs((target - visit.normalize()).days) > 7:
        return {"eligible": False, "reason": "target_visit_date_mismatch"}
    start = target - pd.Timedelta(days=28)
    if start <= index.normalize():
        return {"eligible": False, "reason": "outcome_does_not_follow_index"}
    return {"eligible": True, "reason": "eligible", "start": start,
            "end": target, "lag_days": (start - index.normalize()).days}


def overlapping_pairs(intervals):
    """Both members of overlapping half-open intervals are excluded."""
    bad = set()
    ordered = sorted(enumerate(intervals), key=lambda x: x[1][0])
    for k, (i, (a, b)) in enumerate(ordered):
        for j, (c, d) in ordered[k + 1:]:
            if c >= b:
                break
            if a < d and c < b:
                bad.update((i, j))
    return bad


def independent_episode_labeler(glucose, run_length=4):
    """Run-length based cross-check, independent from streaming card001 labeler.

    Recovery requires a complete >=54 run before a complete <54 run can begin
    an episode. Gaps erase recovery information. The first qualifying low run
    without recovery is censored, never assigned an artificial onset.
    """
    x = np.asarray(glucose, dtype=float)
    state = np.where(~np.isfinite(x), 0, np.where(x < 54, -1, 1))
    boundaries = np.r_[0, np.flatnonzero(np.diff(state) != 0) + 1, len(x)]
    events, censored = [], []
    armed = False
    episode_active = False
    for a, b in zip(boundaries[:-1], boundaries[1:]):
        if a == b:
            continue
        if state[a] == 0:
            armed = False
            episode_active = False
        elif state[a] == 1 and b - a >= run_length:
            armed = True
            episode_active = False
        elif state[a] == -1 and b - a >= run_length:
            if armed:
                events.append(a)
                armed = False
                episode_active = True
            elif not episode_active:
                censored.append(a)
                episode_active = True
    return np.asarray(events, dtype=int), np.asarray(censored, dtype=int)


def residual_identity(total, matched, observed_days):
    if total < matched or min(total, matched) < 0 or observed_days <= 0:
        raise ValueError("Invalid event or exposure accounting")
    residual = total - matched
    return {"total": total, "matched": matched, "unwarned": residual,
            "total_rate": total / observed_days,
            "unwarned_rate": residual / observed_days,
            "sensitivity": matched / total if total else np.nan}


def replay_interval(times, glucose, probabilities, threshold, availability, start, end,
                    run_length=4, warning_glucose=None):
    """Replay an already frozen policy, retaining boundaries and all event onsets.

    Arrays must include at least 150 minutes of lookback and 45 minutes of
    follow-up grid positions. Missing positions must be represented by NaN.
    Neither this function nor its event ledger selects or calibrates a policy.
    """
    from card001_core import episodes, match_alarms, replay
    t, x, p = np.asarray(times), np.asarray(glucose), np.asarray(probabilities)
    available = np.asarray(availability, dtype=bool)
    if any(a.ndim != 1 for a in (t,x,p,available)) or not np.isfinite([start,end]).all() or end <= start:
        raise ValueError("Replay requires one-dimensional arrays and an ordered finite interval")
    if not (len(t) == len(x) == len(p) == len(available)):
        raise ValueError("Misaligned replay arrays")
    if len(t) < 2 or not np.all(np.diff(t) == 300):
        raise ValueError("Expected complete five-minute grid, with NaN for missing data")
    if t[0] > start - 9000 or t[-1] < end + 2700:
        raise ValueError("Boundary grid must retain 150-minute lookback and 45-minute follow-up")
    wx = x if warning_glucose is None else np.asarray(warning_glucose,dtype=float)
    if wx.shape != x.shape:
        raise ValueError("Misaligned delayed warning glucose")
    onsets, censored = episodes(x,run_length)
    check, check_censored = independent_episode_labeler(x,run_length)
    if not np.array_equal(onsets, check) or not np.array_equal(censored, check_censored):
        raise AssertionError("Independent episode labelers disagree")
    legal = available & np.isfinite(wx) & (wx >= 54) & np.isfinite(p)
    alarms = replay(wx, np.where(legal, p, np.nan), threshold)
    # No scored alarm after interval end, but outcome context extends 30 minutes.
    alarms = alarms[t[alarms] < end]
    context_onsets = onsets[t[onsets] <= end + 1800]
    alarm_ledger, matched_indices = match_alarms(alarms, context_onsets, x)
    matched = set(matched_indices.tolist())
    event_rows = []
    for j, onset in enumerate(context_onsets):
        if not start <= t[onset] < end:
            continue
        opportunity = bool(legal[max(0, onset - 6):max(0, onset - 2)].any())
        event_rows.append({"onset_time": float(t[onset]), "onset_index": int(onset),
                           "warned": j in matched, "legal_opportunity": opportunity,
                           "unwarned_reason": "warned" if j in matched else
                           "no_legal_opportunity" if not opportunity else "no_qualifying_alarm"})
    alarm_rows = []
    for row in alarm_ledger:
        ai = row["alarm_index"]
        if start <= t[ai] < end:
            alarm_rows.append({**row, "alarm_time": float(t[ai])})
    inside = (t >= start) & (t < end)
    observed = int(np.isfinite(x[inside]).sum())
    n = len(event_rows)
    m = sum(e["warned"] for e in event_rows)
    n_alarms = len(alarm_rows)
    unconfirmed = sum(a["category"] != "matched" for a in alarm_rows)
    days = observed / 288
    result = {"events": n, "matched": m, "unwarned": n - m,
              "observed_bins": observed, "grid_bins": int(inside.sum()),
              "observed_days": days, "coverage": observed / inside.sum(),
              "alarms": n_alarms, "unconfirmed": unconfirmed,
              "opportunity_events": sum(e["legal_opportunity"] for e in event_rows),
              "censored_runs": int(np.sum((t[censored] >= start) & (t[censored] < end))),
              "sensitivity": m / n if n else np.nan,
              "alarms_per_day": n_alarms / days if days else np.nan,
              "unconfirmed_per_day": unconfirmed / days if days else np.nan}
    return result, pd.DataFrame(event_rows), pd.DataFrame(alarm_rows)


def fixed_design(frame, exposure="auc_nmol_l", within_child=False, augmented=False,
                 quadratic=False, omit_constant_lag=False):
    """Prespecified design; no outcome-based selection and no fitted penalty."""
    f = frame.copy()
    X = pd.DataFrame(index=f.index)
    if within_child:
        X = pd.get_dummies(f.PtID.astype(str), prefix="child", dtype=float)
    else:
        X["intercept"] = 1.0
    X["exposure"] = f[exposure].astype(float) / (0.20 if exposure == "auc_nmol_l" else 1.0)
    if quadratic:
        X["exposure_squared"] = X.exposure ** 2
    X["target_39"] = f.target_visit.eq(39).astype(float)
    X["target_52"] = f.target_visit.eq(52).astype(float)
    X["lag_centered_30d"] = (f.lag_days.astype(float) - 60.0) / 30.0
    if omit_constant_lag:
        if f.lag_days.nunique() != 1:
            raise ValueError("Only an algebraically constant companion lag can be omitted")
        X = X.drop(columns="lag_centered_30d")
    if not within_child:
        X["age_centered_12"] = f.age.astype(float) - 12.0
        X["intensive_tandem"] = f.management.eq("intensive/Tandem").astype(float)
        X["intensive_medtronic"] = f.management.eq("intensive/Medtronic").astype(float)
        X["cohort_B"] = f.cohort.eq("B").astype(float)
        if not f.management.isin(["usual care", "intensive/Tandem", "intensive/Medtronic"]).all():
            raise ValueError("Unknown management category")
        if not f.cohort.isin(["A", "B"]).all():
            raise ValueError("Unknown cohort")
    if not f.target_visit.isin([26, 39, 52]).all():
        raise ValueError("Unknown target visit")
    if augmented:
        X["hba1c"] = f.hba1c.astype(float)
        X["insulin_u_kg_day"] = f.insulin_u_kg_day.astype(float)
    if not np.isfinite(X.to_numpy()).all():
        raise ValueError("Missing essential covariate")
    return X


def poisson_separation(design, counts):
    """Detect a PPML nonexistence direction; rank alone is insufficient.

    A direction zero on positive-count rows and nonpositive on zero-count rows,
    with a strict decrease for at least one zero, improves likelihood without a
    finite maximizer. Bounded coordinates fix the otherwise arbitrary scale.
    """
    from scipy.optimize import linprog
    a = np.asarray(design, dtype=float)
    y = np.asarray(counts, dtype=float)
    zeros, positives = a[y == 0], a[y > 0]
    if not len(zeros):
        return False
    lp = linprog(zeros.sum(axis=0), A_ub=zeros, b_ub=np.zeros(len(zeros)),
                 A_eq=positives if len(positives) else None,
                 b_eq=np.zeros(len(positives)) if len(positives) else None,
                 bounds=[(-1,1)] * a.shape[1], method="highs")
    if not lp.success:
        raise ValueError("PPML existence check failed: " + lp.message)
    return bool(lp.fun < -1e-8)


def _checked_poisson_fit(sm, y, design, offset):
    if poisson_separation(design, y):
        raise ValueError("separated_Poisson_mean_no_finite_MLE")
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        fit = sm.GLM(y, design, family=sm.families.Poisson(), offset=offset).fit(maxiter=100)
    if not fit.converged or not np.isfinite(fit.params).all():
        raise ValueError("nonconverged_or_nonfinite")
    # Numerical identifiability is also checked at the fitted likelihood.
    information = design.T @ (fit.mu[:,None] * design)
    if np.linalg.matrix_rank(information) < information.shape[1]:
        raise ValueError("singular_observed_information")
    notes = sorted(set(type(w.message).__name__ + ": " + str(w.message) for w in caught))
    return fit, notes


def fit_association(frame, outcome="unwarned", exposure="auc_nmol_l", n_boot=5000,
                    seed=SEED, within_child=False, augmented=False, quadratic=False,
                    omit_constant_lag=False, drop_constant_nuisance=False):
    """PPML and delivery-arm-stratified child bootstrap; failed fits retained."""
    import statsmodels.api as sm
    from statsmodels.stats.sandwich_covariance import cov_cluster
    f = frame.copy().reset_index(drop=True)
    raw_y = pd.to_numeric(f[outcome],errors="coerce").to_numpy(dtype=float)
    raw_days = pd.to_numeric(f.observed_days,errors="coerce").to_numpy(dtype=float)
    if not np.isfinite(raw_y).all() or (raw_y < 0).any() or (raw_y != np.floor(raw_y)).any() or not np.isfinite(raw_days).all() or (raw_days <= 0).any():
        return {"outcome":outcome,"exposure":exposure,"n_intervals":len(f),
                "status":"inconclusive","reason":"invalid_count_or_observed_exposure"},pd.DataFrame()
    f[outcome] = raw_y
    f["observed_days"] = raw_days
    if within_child:
        f = f[f.groupby("PtID")[outcome].transform("sum") > 0].copy()
    result = {"outcome": outcome, "exposure": exposure, "n_intervals": len(f),
              "n_children": int(f.PtID.nunique()), "events": int(f[outcome].sum()),
              "observed_days": float(f.observed_days.sum()), "n_boot": n_boot,
              "seed": seed, "status": "inconclusive"}
    if f.empty or f[outcome].sum() <= 0:
        result["reason"] = "no_intervals_or_no_positive_events"
        return result, pd.DataFrame()
    try:
        X = fixed_design(f, exposure, within_child, augmented, quadratic, omit_constant_lag)
    except ValueError as exc:
        result["reason"] = "invalid_design:" + str(exc)
        return result, pd.DataFrame()
    arr = X.to_numpy(dtype=float)
    dropped=[]
    if drop_constant_nuisance:
        if within_child:
            raise ValueError("Do not drop fixed child intercepts")
        dropped=[c for c in X.columns if c not in ("intercept","exposure","exposure_squared") and X[c].nunique()==1]
        X=X.drop(columns=dropped)
        arr=X.to_numpy(dtype=float)
    result["dropped_constant_nuisance_columns"] = dropped
    result["design_columns"] = list(X.columns)
    result["omitted_algebraically_constant_lag"] = bool(omit_constant_lag)
    if np.linalg.matrix_rank(arr) < arr.shape[1]:
        result["reason"] = "rank_deficient_primary_design"
        return result, pd.DataFrame()
    y = f[outcome].to_numpy(dtype=float)
    offset = np.log(f.observed_days.to_numpy(dtype=float))
    if not np.isfinite(offset).all() or (y < 0).any():
        result["reason"] = "invalid_outcome_or_exposure"
        return result, pd.DataFrame()
    try:
        fit, fit_notes = _checked_poisson_fit(sm, y, arr, offset)
        k = list(X.columns).index("exposure")
        beta = float(fit.params[k])
        cov = cov_cluster(fit, f.PtID.to_numpy())
        se = float(np.sqrt(cov[k, k]))
        if not np.isfinite(cov).all() or not np.isfinite(se):
            raise ValueError("nonfinite_cluster_covariance")
        with np.errstate(over="raise", invalid="raise"):
            rr = float(np.exp(beta))
            robust_ci = np.exp([beta - 1.95996398454 * se, beta + 1.95996398454 * se]).tolist()
        if not np.isfinite([rr, *robust_ci]).all() or min([rr,*robust_ci]) <= 0:
            raise ValueError("nonfinite_transformed_estimate_or_interval")
    except Exception as exc:
        result["reason"] = "primary_fit_failure:" + str(exc)
        return result, pd.DataFrame()
    result.update(beta=beta, rr=rr, robust_se=se, robust_95_ci=robust_ci,
                  exposure_min=float(f[exposure].min()), exposure_max=float(f[exposure].max()))
    result["coefficients"]={name:float(value) for name,value in zip(X.columns,fit.params)}
    if quadratic:
        result["rr_interpretation"]="Linear coefficient at exposure zero; not a constant per-0.20 RR when squared term is present"
        b2=float(fit.params[list(X.columns).index("exposure_squared")])
        result["local_per_0_20_rr"]={str(q):float(np.exp(beta+b2*(2*(float(f[exposure].quantile(q))/.20)+1))) for q in (.05,.25,.5,.75,.95)}
    result["fit_warnings"] = fit_notes
    rng = np.random.default_rng(seed)
    rows = {p: np.flatnonzero(f.PtID.to_numpy() == p) for p in f.PtID.unique()}
    per_child = f[["PtID", "delivery_arm"]].drop_duplicates()
    if per_child.PtID.duplicated().any():
        raise ValueError("Child has conflicting randomized delivery arm")
    strata = [g.PtID.to_numpy() for _, g in per_child.groupby("delivery_arm", sort=True)]
    boots = []
    for b in range(n_boot):
        picks = np.concatenate([rng.choice(ids, len(ids), replace=True) for ids in strata])
        draw = json.dumps([str(p) for p in picks],separators=(",",":"))
        ix = np.concatenate([rows[p] for p in picks])
        try:
            if within_child:
                # Fixed intercepts must distinguish repeated copies of a sampled child.
                bf = pd.concat([f.iloc[rows[p]].assign(PtID=f"boot_{j}") for j, p in enumerate(picks)], ignore_index=True)
                bdesign = fixed_design(bf, exposure, True, augmented, quadratic, omit_constant_lag)
                bx = bdesign.to_numpy(dtype=float)
                by = bf[outcome].to_numpy(dtype=float)
                bo = np.log(bf.observed_days.to_numpy(dtype=float))
                bk = list(bdesign.columns).index("exposure")
            else:
                bx, by, bo, bk = arr[ix], y[ix], offset[ix], k
            if np.linalg.matrix_rank(bx) < bx.shape[1]:
                raise ValueError("rank_deficient")
            bfit, bnotes = _checked_poisson_fit(sm, by, bx, bo)
            if not bfit.converged or not np.isfinite(bfit.params[bk]):
                raise ValueError("nonconverged")
            with np.errstate(over="raise", invalid="raise"):
                boot_rr = float(np.exp(bfit.params[bk]))
            if not np.isfinite(boot_rr) or boot_rr <= 0:
                raise ValueError("nonfinite_bootstrap_RR")
            boots.append({"replicate": b, "rr": boot_rr, "failure": "", "warnings": ";".join(bnotes),"sampled_child_ids":draw})
        except Exception as exc:
            boots.append({"replicate": b, "rr": np.nan, "failure": str(exc), "warnings": "","sampled_child_ids":draw})
    boot = pd.DataFrame(boots, columns=["replicate", "rr", "failure", "warnings","sampled_child_ids"])
    failed = int(boot.rr.isna().sum()) if len(boot) else 0
    result["failed_bootstrap"] = failed
    if n_boot and (n_boot - failed) > 0:
        ci = np.nanpercentile(boot.rr, [2.5, 97.5]).tolist()
        result["bootstrap_95_ci"] = ci
        result["numeric_H1_condition"] = result["rr"] <= .8 and ci[0] <= .8 and ci[1] < 1
        result["full_20pct_magnitude_95pct"] = ci[1] <= .8
        result["target_magnitude_excluded"] = ci[0] > .8
    result["status"] = "estimated" if n_boot == 5000 and failed / n_boot <= .05 else "inconclusive"
    result["reason"] = "complete" if result["status"] == "estimated" else "bootstrap_incomplete_or_unstable"
    return result, boot


def gate_assessment(intervals, transfer, model_result, integrity_pass):
    """The qualified claim is a conjunction; never rescue a failed transfer."""
    floors = len(intervals) >= 200 and intervals.PtID.nunique() >= 60
    residual = intervals.unwarned.sum() >= 100 and intervals.loc[intervals.unwarned > 0, "PtID"].nunique() >= 30
    transfer_ok = transfer.get("overall_sensitivity", -1) >= .3
    arms = transfer.get("delivery_arms", {})
    transfer_ok &= len(arms) == 2 and all(v["unconfirmed_per_day"] <= 1 and v["alarms_per_day"] <= 2 for v in arms.values())
    gates = {"sample_floor": bool(floors), "residual_floor": bool(residual),
             "transfer": bool(transfer_ok), "integrity": bool(integrity_pass),
             "identified_model": model_result.get("status") == "estimated"}
    status = "Inconclusive"
    if all(gates.values()):
        if model_result.get("numeric_H1_condition"):
            status = "Supported in this cohort and policy"
        elif model_result.get("target_magnitude_excluded"):
            status = "Targeted association refuted"
    return {"status": status, "gates": gates}
