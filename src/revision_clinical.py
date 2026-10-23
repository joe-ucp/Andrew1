"""Read-only follow-up to completed CLVer analysis with explicit input/output roots."""
from __future__ import annotations

import argparse
import datetime as dt
import json
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.special import expit
from statsmodels.stats.sandwich_covariance import cov_cluster

from card003_core import fixed_design, _checked_poisson_fit, sha256_file

SEED = 20261023
N_BOOT = 5000
Z = 1.959963984540054
BAND_NAMES = ["<0.20", "0.20–<0.40", "0.40–<0.60", "≥0.60"]


def save_json(path, obj):
    Path(path).write_text(json.dumps(obj, indent=2, allow_nan=False), encoding="utf-8")


def read_frame(path):
    return pd.read_csv(path, dtype={"PtID": str, "pair_id": str})


def exact_within_sample(frame):
    """Outcome-positive children, retaining all their zero-event intervals too."""
    return frame.loc[frame.groupby("PtID").unwarned.transform("sum").gt(0)].copy().reset_index(drop=True)


def design_variation(frame):
    f = frame.reset_index(drop=True).copy()
    ids = f.PtID.to_numpy()
    a = f.auc_nmol_l.to_numpy(float)
    x = a - f.groupby("PtID").auc_nmol_l.transform("mean").to_numpy()
    visit = np.column_stack([f.target_visit.eq(v) for v in (39, 52)]).astype(float)
    lag = ((f.lag_days - 60) / 30).to_numpy(float)[:, None]
    def residualize(nuisance):
        centered = nuisance.copy()
        for child in pd.unique(ids):
            ix = ids == child
            centered[ix] -= centered[ix].mean(axis=0)
        return x - centered @ np.linalg.lstsq(centered, x, rcond=None)[0]
    visit_residual = residualize(visit)
    full_residual = residualize(np.column_stack([visit, lag]))
    rows = []
    for child in pd.unique(ids):
        ix = ids == child
        rows.append(dict(PtID=child, intervals=int(ix.sum()), auc_min=float(a[ix].min()),
                         auc_max=float(a[ix].max()), auc_range=float(np.ptp(a[ix])),
                         within_ss=float(x[ix] @ x[ix]),
                         after_visit_ss=float(visit_residual[ix] @ visit_residual[ix]),
                         after_visit_lag_ss=float(full_residual[ix] @ full_residual[ix])))
    bychild = pd.DataFrame(rows)
    multiple = bychild.intervals.gt(1)
    ss = float(x @ x)
    summary = dict(children=len(bychild), intervals=len(f),
                   intervals_per_child={str(int(k)): int(v) for k, v in bychild.intervals.value_counts().sort_index().items()},
                   multiple_interval_children=int(multiple.sum()),
                   raw_varying_children=int((multiple & bychild.auc_range.gt(1e-10)).sum()),
                   after_visit_varying_children=int((multiple & bychild.after_visit_ss.gt(1e-20)).sum()),
                   after_visit_lag_varying_children=int((multiple & bychild.after_visit_lag_ss.gt(1e-20)).sum()),
                   within_auc_sum_squares=ss,
                   after_visit_sum_squares=float(visit_residual @ visit_residual),
                   after_visit_lag_sum_squares=float(full_residual @ full_residual),
                   after_visit_retained_fraction=float(visit_residual @ visit_residual / ss) if ss > 0 else None,
                   after_visit_lag_retained_fraction=float(full_residual @ full_residual / ss) if ss > 0 else None,
                   repeated_child_auc_range_quantiles={str(q): float(bychild.loc[multiple, "auc_range"].quantile(q)) for q in (.25, .5, .75)},
                   interpretation="Unweighted design geometry, not power or causal identification; tolerance 1e-10 nmol/L")
    return summary, bychild


def mundlak_design(f):
    X = fixed_design(f)
    mean = f.groupby("PtID").auc_nmol_l.transform("mean") / .20
    X["exposure"] = X.exposure - mean
    X = X.rename(columns={"exposure": "auc_within"})
    X["auc_between"] = mean
    for col in ["target_39", "target_52", "lag_centered_30d"]:
        X["child_mean_" + col] = X[col].groupby(f.PtID).transform("mean")
    dropped = [c for c in X if c not in ("intercept", "auc_within", "auc_between") and X[c].nunique() == 1]
    return X.drop(columns=dropped), dropped


def checked_estimate(f, X, outcome="unwarned", family="poisson", contrast_columns=("exposure",)):
    arr = X.to_numpy(float)
    result = dict(status="not_estimable", inferential_role="exploratory_child_cluster_sandwich", outcome=outcome,
                  n_rows=len(f), n_children=int(f.PtID.nunique()), design_columns=list(X),
                  design_rank=int(np.linalg.matrix_rank(arr)), events=int(f[outcome].sum()))
    if result["design_rank"] != arr.shape[1]:
        result["reason"] = "rank_deficient_design"
        return result, None, None
    try:
        if family == "poisson":
            fit, notes = _checked_poisson_fit(sm, f[outcome].to_numpy(float), arr, np.log(f.observed_days.to_numpy(float)))
            result["observed_days"] = float(f.observed_days.sum())
            kind = "rate_ratio"
        else:
            with warnings.catch_warnings(record=True) as caught:
                warnings.simplefilter("always")
                fit = sm.GLM(f[outcome].to_numpy(float), arr, family=sm.families.Binomial()).fit(maxiter=100)
            notes = sorted(set(type(w.message).__name__ + ": " + str(w.message) for w in caught))
            if not fit.converged or not np.isfinite(fit.params).all() or any("separat" in w.lower() for w in notes):
                raise ValueError("nonconvergence_nonfinite_or_separation")
            info = arr.T @ ((fit.mu * (1 - fit.mu))[:, None] * arr)
            if np.linalg.matrix_rank(info) != arr.shape[1]:
                raise ValueError("singular_observed_information")
            kind = "conditional_odds_ratio"
        cov = cov_cluster(fit, f.PtID.to_numpy())
        if not np.isfinite(cov).all() or (np.diag(cov) < 0).any():
            raise ValueError("invalid_cluster_covariance")
        contrasts = {}
        for c in contrast_columns:
            k = list(X).index(c)
            b, se = float(fit.params[k]), float(np.sqrt(cov[k, k]))
            transformed = np.exp([b, b-Z*se, b+Z*se])
            if not np.isfinite(transformed).all() or (transformed <= 0).any():
                raise ValueError("nonfinite_transformed_estimate")
            contrasts[c] = dict(beta=b, robust_se=se, estimate=float(transformed[0]),
                                robust_95_ci=transformed[1:].tolist(), measure=kind)
        result.update(status="estimated_exploratory", contrasts=contrasts,
                      coefficients={k: float(v) for k, v in zip(X, fit.params)},
                      cluster_covariance=cov.tolist(), fit_warnings=notes,
                      finite_sample_correction="G/(G-1) times (N-1)/(N-p); N is rows of fitted likelihood")
        return result, fit, cov
    except Exception as exc:
        result["reason"] = str(exc)
        return result, None, None


def event_frame(f, ledger):
    e = ledger.loc[ledger.pair_id.isin(f.pair_id)].copy()
    if e.duplicated(["pair_id", "onset_time"]).any():
        raise AssertionError("Duplicate event ledger key")
    for col in ("warned", "legal_opportunity"):
        if e[col].dtype != bool:
            e[col] = e[col].astype(str).map({"True": True, "False": False})
        if e[col].isna().any():
            raise AssertionError("Unknown event truth value")
    counts = e.groupby("pair_id").agg(events=("warned", "size"), matched=("warned", "sum"), opportunity_events=("legal_opportunity", "sum"))
    aligned = counts.reindex(f.pair_id, fill_value=0)
    for col in aligned:
        np.testing.assert_array_equal(aligned[col].to_numpy(), f[col].to_numpy())
    e = e[["pair_id", "onset_time", "warned", "legal_opportunity"]].merge(f, on="pair_id", validate="many_to_one")
    e["unwarned_binary"] = (~e.warned).astype(int)
    assert int(e.unwarned_binary.sum()) == int(f.unwarned.sum())
    return e


def verify_never(f, never):
    n = never.set_index("pair_id").loc[f.pair_id].reset_index()
    for col in ["PtID", "observed_bins", "grid_bins", "observed_days", "events", "start_seconds", "end_seconds"]:
        np.testing.assert_array_equal(f[col].to_numpy(), n[col].to_numpy())
    np.testing.assert_array_equal(f.events.to_numpy(), n.unwarned.to_numpy())
    assert (n.matched == 0).all() and (n.alarms == 0).all() and (n.unconfirmed == 0).all()
    return dict(status="PASS", intervals=len(f), children=int(f.PtID.nunique()),
                exact_identity="never-warning unwarned count equals all-event count in every identical interval and exposure",
                never_events=int(n.unwarned.sum()), never_alarms=int(n.alarms.sum()),
                consequence="Identical design, response and offset imply identical PPML and original bootstrap interval; reuse is not independent evidence")


def draw_children(f, rng):
    per = f[["PtID", "delivery_arm"]].drop_duplicates()
    if per.PtID.duplicated().any():
        raise ValueError("Child appears in multiple randomized arms")
    return np.concatenate([rng.choice(g.PtID.to_numpy(), len(g), replace=True) for _, g in per.groupby("delivery_arm", sort=True)])


def bands_for(f):
    return np.searchsorted([.20, .40, .60], f.auc_nmol_l.to_numpy(float), side="right")


def band_table(f):
    rows = []
    for k, name in enumerate(BAND_NAMES):
        b = f.loc[bands_for(f) == k]
        totals = {c: int(b[c].sum()) for c in ("events", "matched", "unwarned", "alarms", "unconfirmed", "opportunity_events")}
        days = float(b.observed_days.sum())
        rows.append(dict(band=name, band_index=k, children=int(b.PtID.nunique()), intervals=len(b), observed_days=days,
                         auc_median=float(b.auc_nmol_l.median()), auc_min=float(b.auc_nmol_l.min()), auc_max=float(b.auc_nmol_l.max()),
                         **totals, total_per_100days=100*totals["events"]/days, unwarned_per_100days=100*totals["unwarned"]/days,
                         sensitivity=totals["matched"]/totals["events"], opportunity_coverage=totals["opportunity_events"]/totals["events"],
                         alarms_per_day=totals["alarms"]/days, unconfirmed_per_day=totals["unconfirmed"]/days))
    return pd.DataFrame(rows)


def paired_bootstrap(f, n_boot=N_BOOT, seed=SEED):
    X = fixed_design(f).to_numpy(float)
    k = list(fixed_design(f)).index("exposure")
    y = {name: f[name].to_numpy(float) for name in ("unwarned", "events")}
    offset = np.log(f.observed_days.to_numpy(float))
    rows = {p: np.flatnonzero(f.PtID.to_numpy() == p) for p in f.PtID.unique()}
    band, days = bands_for(f), f.observed_days.to_numpy(float)
    rng = np.random.default_rng(seed)
    out, bandout = [], []
    for b in range(n_boot):
        picks = draw_children(f, rng)
        ix = np.concatenate([rows[p] for p in picks])
        row = dict(replicate=b, sampled_child_ids=json.dumps(picks.tolist(), separators=(",", ":")))
        bx = X[ix]
        for name in ("unwarned", "events"):
            try:
                if np.linalg.matrix_rank(bx) != bx.shape[1]:
                    raise ValueError("rank_deficient")
                fit, notes = _checked_poisson_fit(sm, y[name][ix], bx, offset[ix])
                beta = float(fit.params[k])
                if not np.isfinite(beta) or not np.isfinite(np.exp(beta)):
                    raise ValueError("nonfinite_transformed_estimate")
                row[name + "_beta"] = beta
                row[name + "_failure"] = ""
                row[name + "_warnings"] = ";".join(notes)
            except Exception as exc:
                row[name + "_beta"] = np.nan
                row[name + "_failure"] = str(exc)
                row[name + "_warnings"] = ""
        row["difference_log_rr"] = row["unwarned_beta"] - row["events_beta"]
        out.append(row)
        brow = {"replicate": b}
        for j in range(4):
            ids = ix[band[ix] == j]
            denom = days[ids].sum()
            brow[f"band{j}_days"] = denom
            for name in ("unwarned", "events"):
                brow[f"band{j}_{name}_per_100days"] = 100 * y[name][ids].sum() / denom if denom > 0 else np.nan
        bandout.append(brow)
        if (b+1) % 1000 == 0:
            print(f"Paired bootstrap {b+1}/{n_boot}", flush=True)
    return pd.DataFrame(out), pd.DataFrame(bandout)


def cohort_context(f, roster, temporal, allpairs):
    included = set(f.PtID)
    eligible = temporal
    child = f.drop_duplicates("PtID")
    subset = exact_within_sample(f).drop_duplicates("PtID")
    context = dict(flow=[
        dict(stage="randomized", children=int(roster.PtID.nunique()), pairs=None),
        dict(stage="possible forward pairs", children=int(allpairs.PtID.nunique()), pairs=len(allpairs)),
        dict(stage="assay and temporal gates", children=int(eligible.PtID.nunique()), pairs=len(eligible)),
        dict(stage="at least 70% observed CGM", children=len(included), pairs=len(f)),
        dict(stage="original fixed-child sample", children=len(subset), pairs=len(exact_within_sample(f)))],
        sex_variable="Not present in released randomized roster or canonical analysis input; not inferred",
        treatment_context="Recorded management and existing treatment/device alarms influenced observed trajectories; replay alarms were not deployed",
        populations={})
    for name, c in [("randomized", roster), ("primary", child), ("within_child_sample", subset)]:
        c = c.drop_duplicates("PtID")
        context["populations"][name] = dict(children=len(c), age_median=float(c.age.median()),
            age_q25=float(c.age.quantile(.25)), age_q75=float(c.age.quantile(.75)), age_min=float(c.age.min()), age_max=float(c.age.max()),
            management=c.management.value_counts().to_dict(), source=c.source.value_counts().to_dict(), cohort=c.cohort.value_counts().to_dict(),
            delivery_arm=c.delivery_arm.value_counts().to_dict())
    context["primary_intervals"] = dict(target_visit={str(k): int(v) for k, v in f.target_visit.value_counts().sort_index().items()},
        observed_days=float(f.observed_days.sum()), auc_min=float(f.auc_nmol_l.min()), auc_max=float(f.auc_nmol_l.max()),
        auc_q25=float(f.auc_nmol_l.quantile(.25)), auc_median=float(f.auc_nmol_l.median()), auc_q75=float(f.auc_nmol_l.quantile(.75)),
        lag_median=float(f.lag_days.median()), lag_q25=float(f.lag_days.quantile(.25)), lag_q75=float(f.lag_days.quantile(.75)))
    return context


def standard_curve(f, fit, cov):
    values = np.linspace(f.auc_nmol_l.min(), f.auc_nmol_l.max(), 250)
    rows = []
    for a in values:
        X = fixed_design(f.assign(auc_nmol_l=a)).to_numpy(float)
        rates = np.exp(X @ fit.params)
        mean = rates.mean()
        gradient = (rates[:, None] * X).mean(axis=0) / mean
        se = float(np.sqrt(gradient @ cov @ gradient))
        rows.append(dict(auc_nmol_l=a, rate_per_100days=100*mean,
                         lower=100*mean*np.exp(-Z*se), upper=100*mean*np.exp(Z*se)))
    return pd.DataFrame(rows)


def plot_results(f, bands, curves, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10})
    fig, axes = plt.subplots(2, 2, figsize=(11.2, 7.2), gridspec_kw={"height_ratios": [3, 1]}, sharex="col")
    colors = {"events": "#247BA0", "unwarned": "#B6465F"}
    for j, name in enumerate(("events", "unwarned")):
        ax, hx, curve, color = axes[0, j], axes[1, j], curves[name], colors[name]
        rate = "total_per_100days" if name == "events" else "unwarned_per_100days"
        ax.plot(curve.auc_nmol_l, curve.rate_per_100days, color=color, label="Adjusted fitted mean")
        ax.fill_between(curve.auc_nmol_l, curve.lower, curve.upper, color=color, alpha=.15, label="95% cluster CI")
        ax.errorbar(bands.auc_median, bands[rate], yerr=np.vstack([bands[rate]-bands[rate+"_lower"], bands[rate+"_upper"]-bands[rate]]),
                    fmt="o", color="#222222", capsize=4, label="Empirical band rate (95% bootstrap CI)")
        ax.set_title("All observed events / never-warning" if name == "events" else "Events unwarned by frozen M15")
        ax.set_ylabel("Episodes per 100 observed-sensor days")
        ax.set_ylim(bottom=0)
        ax.grid(alpha=.18)
        ax.legend(fontsize=8, loc="upper right", frameon=False)
        hx.hist(f.auc_nmol_l, bins=np.linspace(f.auc_nmol_l.min(), f.auc_nmol_l.max(), 26), color=color, alpha=.65)
        for cutoff in [.2, .4, .6]:
            hx.axvline(cutoff, color="#777777", lw=.8, linestyle="--")
        hx.set_ylabel("Intervals")
        hx.set_xlabel("Earlier normalized MMTT C-peptide AUC (nmol/L)")
        hx.set_xlim(f.auc_nmol_l.min(), f.auc_nmol_l.max())
        for a in (ax, hx):
            a.spines[["top", "right"]].set_visible(False)
    fig.suptitle("Exploratory clinical follow-up: event frequency and replay residual risk", fontsize=14, y=.98)
    fig.text(.5, .015, "Points: fixed AUC bands at their observed medians. Curves: equal-interval standardized means, not causal effects.\nOverall observed support only; support can be sparse or absent within a covariate stratum. Warning-burden transfer gate remains failed.", ha="center", fontsize=8.5)
    fig.tight_layout(rect=[0, .065, 1, .945])
    fig.savefig(out / "joint_auc_rates.png", dpi=180)
    plt.close(fig)


def build_report(summary, bands, out, repo):
    def fmt(c):
        return f"{c['estimate']:.3f} ({c['robust_95_ci'][0]:.3f}–{c['robust_95_ci'][1]:.3f})"
    pooled = summary["same_sample_pooled"]["contrasts"]["exposure"]
    event = summary["conditional_warning"]["contrasts"]["exposure"]
    paired = summary["paired_comparison"]
    v = summary["variation"]["within_child_sample"]
    ctx = summary["cohort_context"]
    primary = ctx["populations"]["primary"]
    btxt = "\n".join(f"| {r.band} | {r.children} / {r.intervals} | {r.events} / {r.unwarned} | {r.total_per_100days:.3f} / {r.unwarned_per_100days:.3f} | {r.sensitivity*100:.1f}% | {r.opportunity_coverage*100:.1f}% | {r.alarms_per_day:.3f} / {r.unconfirmed_per_day:.3f} |" for r in bands.itertuples())
    m = summary["mundlak"]
    mtext = (f"Full correlated-effects decomposition gives within-child RR {fmt(m['contrasts']['auc_within'])} and between-child RR {fmt(m['contrasts']['auc_between'])}." if m["status"] == "estimated_exploratory" else "The correlated-effects decomposition is not estimable: " + m.get("reason", "unknown"))
    text = f"""# Exploratory clinical reanalysis, revision r2

These follow-up analyses were declared in `protocol/revision_clinical.md` before new estimation, after review of the original results. The original source, frozen policy and all original scientific artifacts remain unchanged. The qualified warning-transfer claim remains **inconclusive because its alarm-burden gate failed**. None of these post-review analyses changes that decision.

## Pooled versus within-child association on identical participants

The exact original fixed-child sample contains 58 children and 164 intervals, retaining all their zero-event intervals, 180 unwarned events, and {summary['same_sample_pooled']['observed_days']:.3f} observed-sensor days. Original-form pooled PPML on these identical rows gives **RR {fmt(pooled)} per 0.20 nmol/L**. The preserved original fixed-child estimate is **RR 1.035 (0.631–1.697)**. Both displayed intervals here are child-cluster sandwich intervals. The original full-sample pooled RR is 0.764 (original 5,000-draw bootstrap interval 0.626–0.896). Restricting to outcome-positive children is itself outcome-based sample selection; it does not create a causal comparison.

Of the 58 fixed-child-sample children, {v['multiple_interval_children']} have multiple included intervals and {v['raw_varying_children']} have nonzero within-child AUC variation. Counts by number of intervals: {json.dumps(v['intervals_per_child'])}. After removing shared target-visit variation, {v['after_visit_varying_children']} children retain nonzero residual variation and {100*v['after_visit_retained_fraction']:.1f}% of the raw within-child AUC sum of squares remains. Further adjusting lag leaves {100*v['after_visit_lag_retained_fraction']:.1f}%. These quantify numerical design information, not power; nonzero variation can still be small. The repeated-child AUC range median is {v['repeated_child_auc_range_quantiles']['0.5']:.3f} nmol/L (IQR {v['repeated_child_auc_range_quantiles']['0.25']:.3f}–{v['repeated_child_auc_range_quantiles']['0.75']:.3f}). Full-sample and child-level diagnostic tables are retained.

{mtext} This model augments the original adjustment with each child's mean AUC, mean visit indicators and mean lag, separating AUC deviation from child mean. It is exploratory, does not eliminate time-varying confounding, and is not asserted to equal fixed-child PPML algebraically.

## Event frequency versus conditional warning probability

The audited primary ledger contains {summary['conditional_warning']['n_rows']} qualifying events in {summary['conditional_warning']['event_bearing_intervals']} intervals from {summary['conditional_warning']['n_children']} children: 180 unwarned and 228 usefully warned. Event-level logistic regression with the original pooled covariates gives **conditional odds ratio {fmt(event)} per 0.20 nmol/L** for being unwarned given an observed qualifying event. This is an odds ratio, not a rate ratio. The event-level and grouped-binomial coefficients agree to {summary['conditional_warning']['grouped_binomial_max_beta_difference']:.2g}. The event-level child-cluster correction uses 408 event rows. Zero-event intervals contribute no trials to this conditional likelihood but remain in every rate analysis.

Conditioning on an observed event selects a subset of trajectories; event occurrence and observation may depend on secretion, management and unmeasured factors. This conditional association is not causal warnability, event prevention, awareness or an intervention effect. The warning policy and its information timing were unchanged.

Using the same 5,000 child draws stratified by delivery arm (new seed 20261023), the paired difference **log RR(unwarned) − log RR(total)** is **{paired['difference_log_rr']:.4f}**, percentile 95% CI **{paired['difference_log_rr_95_ci'][0]:.4f}–{paired['difference_log_rr_95_ci'][1]:.4f}**. The **ratio of RRs is {paired['ratio_of_rrs']:.3f} ({paired['ratio_of_rrs_95_ci'][0]:.3f}–{paired['ratio_of_rrs_95_ci'][1]:.3f})**. Joint failures: {paired['joint_failures']}/5000; unwarned fit failures: {paired['unwarned_failures']}; total fit failures: {paired['events_failures']}. An interval containing one does not establish equivalence. Separate marginal intervals from these new draws are in the JSON; they do not replace the original primary interval.

The never-warning control has exactly the same outcome as all-event frequency in every primary interval and identical offsets: 408 residual events, zero alarms and zero matches. Its adjusted estimate is therefore the original all-event RR **0.774 (0.674–0.870)**, with the same original seed-20261022 bootstrap interval reused only after identity assertions. This mathematical identity is not independent corroboration.

## Joint descriptive patterns

| Earlier AUC (nmol/L) | Children / intervals | Total / unwarned events | Total / unwarned per 100 days | Useful sensitivity | Legal opportunity | All / unconfirmed alarms per day |
|---|---:|---:|---:|---:|---:|---:|
{btxt}

Child counts overlap across bands. These rates use observed-sensor days, and sensitivity/opportunity use all qualifying observed events. The table preserves the joint pattern of event frequency, warning opportunity and alarm burden rather than interpreting residual rate alone. Full tables include exposure, matched counts and bootstrap rate intervals.

![Observed and modeled total and unwarned event rates]({str(out / 'joint_auc_rates.png')})

Points are empirical fixed-band rates at band medians, with 95% percentile intervals from the paired child draws. Curves set AUC to each value across observed overall support and average fitted rates equally over all 276 intervals' original covariates. Shading is the child-cluster delta-method 95% interval for each standardized fitted mean, not an individual prediction interval. Overall AUC support is {ctx['primary_intervals']['auc_min']:.3f}–{ctx['primary_intervals']['auc_max']:.3f} nmol/L; particular covariate strata can lack support at those extremes. Neither empirical nor fitted contrasts establish intervention effects.

## Cohort and source context

The flow is 113 randomized children → 339 possible forward pairs → 301 assay/temporal pairs from 108 children → 276 intervals from 102 children meeting at least 70% observed coverage. The primary children's enrollment age median is {primary['age_median']:.1f} years (IQR {primary['age_q25']:.1f}–{primary['age_q75']:.1f}; range {primary['age_min']:.0f}–{primary['age_max']:.0f}). Child-level management counts: {json.dumps(primary['management'])}; CGM sources: {json.dumps(primary['source'])}; cohorts: {json.dumps(primary['cohort'])}. Sex is absent from the released roster/canonical analysis file used here and was not inferred. Detailed randomized-versus-included context and target-visit counts are saved.

The final strict-time primary excludes unresolved date-only measurements. Blank DeviceCGM units retain the original explicit mg/dL assumption; clock discontinuities remain masked separators with following consistent local-time segments retained. Absolute clock validity is not proven. The original source audits, pre-audit sensitivity history, boundary synchronization proof and policy hashes remain in force. This revision neither replays a different policy nor repairs the transfer failure by retuning.

## Reproduction and preserved provenance

Run `python src/revision_clinical.py --source-root STUDY_ROOT --output-root NEW_REVISION_ROOT` with the pinned dependencies. Protocol, revision script and all input hashes were captured before estimation in `analysis_commitment.json`; the completion manifest verifies original protected artifacts unchanged. New numerical data and figures are in `{out}` and the report is under the revision output root. All analyses here are exploratory, unadjusted for multiplicity, and should be interpreted together with the original null/imprecise fixed-child and clinical-remission findings.
"""
    report_dir=out.parent/"reports"
    report_dir.mkdir(parents=True,exist_ok=True)
    (report_dir / "revision_clinical_r2.md").write_text(text, encoding="utf-8")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--source-root", type=Path, required=True)
    ap.add_argument("--output-root", type=Path, required=True)
    ap.add_argument("--repo", type=Path, default=Path(__file__).resolve().parents[1])
    args = ap.parse_args()
    source, outroot, repo = args.source_root.resolve(), args.output_root.resolve(), args.repo.resolve()
    if source == outroot or source in outroot.parents or outroot in source.parents:
        raise ValueError("New numerical outputs must be in a separate revision directory")
    out = outroot / "results"
    out.mkdir(parents=True, exist_ok=True)
    if (outroot / "analysis_commitment.json").exists():
        raise FileExistsError("This revision already has a commitment; preserve completed/partial artifacts and choose a new directory")
    protected = [p for p in (source / "results").rglob("*") if p.is_file()]
    protected += [p for p in (source / "models").rglob("*") if p.is_file()]
    protected += [source / "processed" / name for name in ("randomized_roster.csv", "interval_derivations.csv", "intervals_all_coverage.csv")]
    protected += sorted((repo / "src").glob("card003*.py")) + sorted((repo / "reports").glob("card003*.md")) + [repo / "protocol" / "card003_methods.md"]
    before = {str(p): sha256_file(p) for p in sorted(set(protected))}
    commitment = dict(phase="DECLARED_BEFORE_REVISION_ESTIMATION", exploratory=True, utc=dt.datetime.now(dt.timezone.utc).isoformat(),
                      seed=SEED, paired_bootstrap_replicates=N_BOOT, protocol_sha256=sha256_file(repo / "protocol" / "revision_clinical.md"),
                      revision_script_sha256=sha256_file(Path(__file__)), source_root=str(source), output_root=str(outroot), protected_sha256=before)
    save_json(outroot / "analysis_commitment.json", commitment)
    start = time.perf_counter()
    f = read_frame(source / "results" / "primary_intervals.csv")
    assert len(f) == 276 and f.PtID.nunique() == 102 and not f.pair_id.duplicated().any()
    assert f.events.sum() == 408 and f.unwarned.sum() == 180
    original_primary = json.loads((source / "results" / "model_primary.json").read_text())
    original_total = json.loads((source / "results" / "model_total_events.json").read_text())
    original_within = json.loads((source / "results" / "model_within_child.json").read_text())
    sub = exact_within_sample(f)
    assert len(sub) == original_within["n_intervals"] == 164 and sub.PtID.nunique() == original_within["n_children"] == 58
    assert int(sub.unwarned.sum()) == original_within["events"]
    assert abs(sub.observed_days.sum() - original_within["observed_days"]) < 1e-8
    sub.to_csv(out / "exact_within_child_sample.csv", index=False)
    summary = dict(exploratory=True, original_primary_result=original_primary, original_fixed_child_result=original_within,
                   original_total_result=original_total, original_primary_decision="INCONCLUSIVE: warning-burden transfer gate failed", variation={})
    same, _, _ = checked_estimate(sub, fixed_design(sub))
    summary["same_sample_pooled"] = same
    save_json(out / "same_sample_pooled.json", same)
    for name, sample in (("full_primary", f), ("within_child_sample", sub)):
        summary["variation"][name], bychild = design_variation(sample)
        bychild.to_csv(out / f"auc_variation_{name}.csv", index=False)
    save_json(out / "auc_variation.json", summary["variation"])
    mx, removed = mundlak_design(f)
    mundlak, _, _ = checked_estimate(f, mx, contrast_columns=("auc_within", "auc_between"))
    mundlak["dropped_algebraically_constant_columns"] = removed
    summary["mundlak"] = mundlak
    save_json(out / "mundlak.json", mundlak)
    e = event_frame(f, read_frame(source / "results" / "events_primary.csv"))
    logistic, eventfit, eventcov = checked_estimate(e, fixed_design(e), outcome="unwarned_binary", family="binomial")
    if eventfit is None:
        raise ValueError("Conditional logistic model was not estimable: " + logistic.get("reason", "unknown"))
    positive = f.loc[f.events.gt(0)]
    grouped = sm.GLM(positive[["unwarned", "matched"]].to_numpy(), fixed_design(positive).to_numpy(), family=sm.families.Binomial()).fit(maxiter=100)
    diff = float(np.max(np.abs(grouped.params-eventfit.params)))
    np.testing.assert_allclose(grouped.params, eventfit.params, atol=1e-8, rtol=1e-8)
    logistic.update(event_bearing_intervals=int(e.pair_id.nunique()), grouped_binomial_max_beta_difference=diff,
                    conditioning="Observed qualifying event; not a causal effect or rate ratio", zero_event_intervals_retained_in_rate_models=int(f.events.eq(0).sum()))
    summary["conditional_warning"] = logistic
    e.to_csv(out / "conditional_event_rows.csv", index=False)
    save_json(out / "conditional_warning_logistic.json", logistic)
    ident = verify_never(f, read_frame(source / "results" / "intervals_primary_never.csv"))
    ident["original_total_result_reused"] = original_total
    summary["never_warning_identity"] = ident
    save_json(out / "never_warning_identity.json", ident)
    fits, covs = {}, {}
    for name, original in (("unwarned", original_primary), ("events", original_total)):
        model, fits[name], covs[name] = checked_estimate(f, fixed_design(f), outcome=name)
        assert abs(model["contrasts"]["exposure"]["beta"] - original["beta"]) < 1e-9
        save_json(out / ("reproduced_" + name + ".json"), model)
    boot, bboot = paired_bootstrap(f)
    boot.to_csv(out / "paired_child_bootstrap.csv.gz", index=False, compression="gzip")
    bboot.to_csv(out / "paired_band_rates.csv.gz", index=False, compression="gzip")
    good = boot.difference_log_rr.notna()
    difference = original_primary["beta"] - original_total["beta"]
    diff_ci = np.percentile(boot.loc[good, "difference_log_rr"], [2.5, 97.5])
    comparison = dict(status="estimated_exploratory" if (~good).mean() <= .05 else "unstable_bootstrap", seed=SEED, n_boot=N_BOOT,
        joint_failures=int((~good).sum()), unwarned_failures=int(boot.unwarned_beta.isna().sum()), events_failures=int(boot.events_beta.isna().sum()),
        difference_log_rr=float(difference), difference_log_rr_95_ci=diff_ci.tolist(), ratio_of_rrs=float(np.exp(difference)),
        ratio_of_rrs_95_ci=np.exp(diff_ci).tolist(), inference="Paired whole-child bootstrap; interval containing one does not establish equivalence",
        marginal_rr_95_ci={name: np.exp(np.percentile(boot.loc[boot[name+"_beta"].notna(), name+"_beta"], [2.5, 97.5])).tolist() for name in ("unwarned", "events")})
    summary["paired_comparison"] = comparison
    save_json(out / "paired_comparison.json", comparison)
    bands = band_table(f)
    for j in range(4):
        for name, rate in (("events", "total_per_100days"), ("unwarned", "unwarned_per_100days")):
            vals = bboot[f"band{j}_{name}_per_100days"]
            lo, hi = np.nanpercentile(vals, [2.5, 97.5])
            bands.loc[j, rate+"_lower"] = lo
            bands.loc[j, rate+"_upper"] = hi
            bands.loc[j, rate+"_undefined_draws"] = int(vals.isna().sum())
    bands.to_csv(out / "joint_auc_bands.csv", index=False)
    curves = {name: standard_curve(f, fits[name], covs[name]) for name in ("events", "unwarned")}
    for name, curve in curves.items():
        curve.to_csv(out / ("fitted_curve_" + name + ".csv"), index=False)
    plot_results(f, bands, curves, out)
    summary["cohort_context"] = cohort_context(f, read_frame(source / "processed" / "randomized_roster.csv"),
                                             read_frame(source / "processed" / "intervals_all_coverage.csv"),
                                             read_frame(source / "processed" / "interval_derivations.csv"))
    save_json(out / "cohort_context.json", summary["cohort_context"])
    build_report(summary, bands, out, repo)
    changed = [p for p, h in before.items() if sha256_file(Path(p)) != h]
    if changed:
        raise AssertionError("Original protected artifact changed: " + repr(changed))
    summary["original_artifact_preservation"] = dict(status="PASS", protected_files=len(before), changed=[])
    summary["elapsed_seconds"] = time.perf_counter() - start
    summary["complete"] = True
    save_json(out / "revision_clinical_summary.json", summary)
    save_json(outroot / "completion_manifest.json", dict(status="COMPLETE_EXPLORATORY", utc=dt.datetime.now(dt.timezone.utc).isoformat(),
        original_preservation=summary["original_artifact_preservation"], protocol_sha256=sha256_file(repo / "protocol" / "revision_clinical.md"),
        script_sha256=sha256_file(Path(__file__)), commitment_sha256=sha256_file(outroot / "analysis_commitment.json"),
        outputs_sha256={str(p.relative_to(outroot)): sha256_file(p) for p in out.rglob("*") if p.is_file()}))
    print(json.dumps(dict(complete=True, same_sample=summary["same_sample_pooled"]["contrasts"], logistic=logistic["contrasts"], paired=comparison,
                         elapsed_seconds=summary["elapsed_seconds"], original_preservation=summary["original_artifact_preservation"]), indent=2), flush=True)


if __name__ == "__main__":
    main()
