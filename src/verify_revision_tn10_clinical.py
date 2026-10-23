"""Independent arithmetic review; imports neither TN10 production module."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from scipy.stats import t as student_t


def digest(p):
    return hashlib.sha256(Path(p).read_bytes()).hexdigest()


def main(base, revised, out, raw_root=None):
    base, revised, out = map(lambda p: Path(p).resolve(), (base, revised, out))
    raw_root = Path(raw_root).resolve() if raw_root is not None else base / "raw"
    out.mkdir(parents=True, exist_ok=True)
    repo = Path(__file__).resolve().parents[1]
    sources = [raw_root/"TN10_metadata.csv", base/"results/analysis_dataset.tsv", base/"results/primary_results.json",
               base/"results/continuous_models.json", base/"results/bootstrap_draws.tsv", base/"results/bootstrap_results.json",
               revised/"followup_identification_range.json", revised/"tn10_curve_steps.csv", revised/"tn10_curve_risk_counts.csv",
               repo/"src/revision_tn10.py"]
    before = {str(p): digest(p) for p in sources}
    data = pd.read_csv(base/"results/analysis_dataset.tsv", sep="\t")
    raw = pd.read_csv(raw_root/"TN10_metadata.csv")
    meta = raw.loc[raw["Time.Point"].eq("Baseline")].set_index("libraryID").loc[data.libraryID]
    assert len(data) == 73 and len(meta) == 73 and meta.index.is_unique
    event = pd.to_numeric(meta.DIAG_T1D).to_numpy(int)
    duration = pd.to_numeric(meta.Survival).to_numpy(float)
    arm = meta.TREATMENTDESC.eq("Teplizumab").to_numpy(int)
    np.testing.assert_array_equal(data.event, event)
    np.testing.assert_array_equal(data.treatment, arm)
    np.testing.assert_allclose(data.survival_weeks, duration, atol=1e-10, rtol=0)
    start = pd.to_datetime(meta.TREATMENTDATE, format="%m/%d/%Y")
    date = pd.to_datetime(meta.IDDMDATE.loc[event == 1], format="%m/%d/%Y")
    elapsed = (date-start.loc[event == 1]).dt.total_seconds().to_numpy()/604800
    err = float(np.max(np.abs(elapsed-duration[event == 1])))
    assert err < 1e-9 and sum(event == 0) == 19
    assert meta.IDDMDATE.loc[event == 0].eq(".").all()
    assert len(np.unique(duration[event == 0])) == 1
    assert np.isclose(duration[event == 0][0], np.round(max(duration[event == 1]), 4), atol=1e-10)
    assert not any("last" in str(c).lower() or "follow" in str(c).lower() for c in meta.columns)
    y = np.minimum(duration, 156)
    np.testing.assert_allclose(data.restricted_weeks, y, atol=1e-10, rtol=0)
    high = data.high.to_numpy(bool)
    groups = [(False, 1, "low_treated", 1), (False, 0, "low_placebo", -1),
              (True, 1, "high_treated", -1), (True, 0, "high_placebo", 1)]
    stored = json.loads((base/"results/primary_results.json").read_text())
    bounds = json.loads((revised/"followup_identification_range.json").read_text())
    steps = pd.read_csv(revised/"tn10_curve_steps.csv")
    risks = pd.read_csv(revised/"tn10_curve_risk_counts.csv")
    means, parts, weights, summaries = [], [], np.zeros(len(data)), []
    for h, a, key, sign in groups:
        ix = (high == h) & (arm == a)
        yy = y[ix]
        n = int(sum(ix))
        weights[ix] = sign / n
        means.append(yy.mean())
        parts.append(float(yy.var(ddof=1)/n))
        assert abs(yy.mean()-stored["cells"][key]["mean_weeks"]) < 1e-10
        s = steps[steps.group.eq(key)].sort_values("weeks")
        for row in s.itertuples():
            expected = 1 - np.sum(ix & (event == 1) & (duration <= row.weeks)) / n
            assert abs(expected-row.stage3_free) < 1e-12
        area = float(np.sum(np.diff(s.weeks.to_numpy())*s.stage3_free.to_numpy()[:-1]))
        assert abs(area-yy.mean()) < 1e-10
        r = risks[risks.group.eq(key)]
        for row in r.itertuples():
            assert row.encoded_at_risk == np.sum(ix & (duration >= row.week))
        summaries.append(dict(group=key, n=n, diagnoses_by_156=int(sum(ix & (event == 1) & (duration <= 156))),
                              unknown_followup=int(sum(ix & (event == 0))), mean_weeks=float(yy.mean()), curve_area=area))
    theta = float(weights @ y)
    var = sum(parts)
    dfs = [s["n"]-1 for s in summaries]
    df = var**2/sum(p*p/d for p,d in zip(parts, dfs))
    ci = theta + np.array([-1, 1])*student_t.ppf(.975, df)*np.sqrt(var)
    np.testing.assert_allclose([theta, df, np.sqrt(var)], [stored["theta"]["estimate_weeks"], stored["theta"]["df"], stored["theta"]["se_weeks"]], atol=1e-10)
    np.testing.assert_allclose(ci, stored["theta"]["ci95_weeks"], atol=1e-10)
    uncertain = event == 0
    known = float(weights[~uncertain] @ y[~uncertain])
    lower = known + float(np.minimum(0, 156*weights[uncertain]).sum())
    upper = known + float(np.maximum(0, 156*weights[uncertain]).sum())
    np.testing.assert_allclose([lower, upper], bounds["theta_range_weeks"], atol=1e-10)
    z = ((data.score-data.score.mean())/data.score.std(ddof=1)).to_numpy()
    age = pd.to_numeric(meta.AGE_AT_RAND, errors="coerce").to_numpy(float)
    hb = pd.to_numeric(meta.HBA1C, errors="coerce").to_numpy(float)
    complete = np.isfinite(age) & np.isfinite(hb)
    continuous = json.loads((base/"results/continuous_models.json").read_text())
    continuous_checks = []
    for key, mask, adjusted in [("unadjusted_full", np.ones(len(y),bool),False), ("unadjusted_complete",complete,False), ("adjusted_complete",complete,True)]:
        X = np.column_stack([np.ones(len(y)), arm, z, arm*z])
        if adjusted:
            X = np.column_stack([X, age, hb])
        X, yy = X[mask], y[mask]
        beta = np.linalg.lstsq(X, yy, rcond=None)[0]
        inverse = np.linalg.inv(X.T @ X)
        leverage = np.sum((X @ inverse) * X, axis=1)
        residual = yy-X @ beta
        hc3 = inverse @ (X.T @ ((residual/(1-leverage))[:,None]**2*X)) @ inverse
        se = float(np.sqrt(hc3[3,3]))
        interval = beta[3] + np.array([-1,1])*student_t.ppf(.975, len(yy)-X.shape[1])*se
        expected = continuous[key]["coefficients"]["treatment_x_score_z"]
        np.testing.assert_allclose([beta[3],se], [expected["estimate"],expected["se"]], atol=1e-9)
        np.testing.assert_allclose(interval, expected["ci95"], atol=1e-9)
        continuous_checks.append(dict(model=key, n=len(yy), beta=float(beta[3]), se=se, ci95=interval.tolist()))
    draws = pd.read_csv(base/"results/bootstrap_draws.tsv", sep="\t")
    bootstrap = json.loads((base/"results/bootstrap_results.json").read_text())
    assert len(draws) == 10000 and np.isfinite(draws.to_numpy()).all()
    bootci = np.percentile(draws.theta, [2.5,97.5])
    np.testing.assert_allclose(bootci, bootstrap["percentile_ci95_weeks"]["theta"], atol=1e-10)
    after = {str(p): digest(p) for p in sources}
    assert after == before
    result = dict(status="PASS", imported_production_analysis=False, input_sha256=before, original_files_unchanged=True,
        event_date_max_error_weeks=err, non_event_n=19, common_non_event_duration_weeks=float(duration[event==0][0]),
        actual_last_contact_available=False, cells=summaries, theta_weeks=theta, welch_ci95=ci.tolist(),
        conservative_sample_identification_range=[lower,upper], continuous_manual_hc3=continuous_checks,
        saved_bootstrap_percentiles=bootci.tolist(), bootstrap_scope="Percentiles of all stored draws verified; gene-level bootstrap itself not rerun in this review",
        score_scope="Uses original audited score labels; does not independently rederive gene list or count normalization",
        risk_scope="Encoded immediately-before risk counts only; does not verify missing clinical last-contact information")
    (out/"tn10_independent_clinical_review.json").write_text(json.dumps(result,indent=2,allow_nan=False),encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in ("input_sha256", "cells")},indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--original-root', type=Path, required=True)
    parser.add_argument('--revision-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    parser.add_argument('--raw-root', type=Path)
    args = parser.parse_args()
    main(args.original_root, args.revision_root, args.output_root, args.raw_root)
