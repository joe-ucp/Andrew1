"""Extract event/opportunity/warning accounting and saved nonlinear contrasts.

Reads completed study ledgers and fit summaries. No model fitting, prediction,
policy selection, participant exclusion or modification of original evidence.
Only aggregate outputs are written, beneath the revision-3 paper data folder.
"""
from pathlib import Path
from datetime import datetime, timezone
import argparse
import hashlib
import json
import math

import numpy as np
import pandas as pd


def sha(path):
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def truth(series):
    if series.dtype == bool:
        return series.to_numpy()
    mapped = series.map({"True": True, "False": False})
    assert mapped.notna().all(), "Unrecognized Boolean ledger value"
    return mapped.to_numpy(bool)


def close(a, b, atol=1e-9):
    assert math.isclose(float(a), float(b), rel_tol=1e-12, abs_tol=atol), (a, b)


def decomposition(cohort, policy, people, observed_people, intervals, E, O, W, bins):
    E, O, W, bins = map(int, (E, O, W, bins))
    assert 0 <= W <= O <= E and bins > 0 and E > 0 and O > 0
    days = bins / 288
    without_opportunity, opportunity_unwarned = E - O, O - W
    assert E - W == without_opportunity + opportunity_unwarned
    rate = (E / days) * (1 - (O / E) * (W / O))
    close(rate, (E - W) / days)
    return dict(cohort=cohort, policy=policy, allocated_people=int(people),
        observed_people=int(observed_people), intervals=int(intervals), events=E,
        opportunity_events=O, warned_events=W, unwarned_events=E-W,
        without_opportunity=without_opportunity,
        opportunity_without_warning=opportunity_unwarned,
        opportunity_fraction=O/E, conversion_given_opportunity=W/O,
        useful_warning_sensitivity=W/E,
        fraction_unwarned_with_opportunity=opportunity_unwarned/(E-W),
        observed_bins=bins, observed_days=days,
        total_event_rate_per_100days=100*E/days,
        unwarned_rate_per_100days=100*(E-W)/days,
        no_opportunity_rate_per_100days=100*without_opportunity/days,
        opportunity_unwarned_rate_per_100days=100*opportunity_unwarned/days)


def alarm_row(cohort, policy, ledger, days, expected_alarms, expected_unconfirmed, W):
    categories = ["matched", "duplicate", "late", "unobservable_outcome", "unmatched_observed"]
    assert set(ledger.category.unique()).issubset(categories)
    counts = {c: int(ledger.category.eq(c).sum()) for c in categories}
    assert sum(counts.values()) == len(ledger) == expected_alarms
    assert counts["matched"] == W
    assert expected_unconfirmed == len(ledger) - W
    return dict(cohort=cohort, policy=policy, observed_days=days,
        **counts, alarms=len(ledger), unconfirmed=int(expected_unconfirmed),
        unmatched_observed_per_day=counts["unmatched_observed"]/days,
        late_per_day=counts["late"]/days,
        unobservable_outcome_per_day=counts["unobservable_outcome"]/days,
        unconfirmed_per_day=expected_unconfirmed/days,
        total_alarms_per_day=len(ledger)/days)


def run(root, output_root=None):
    root = Path(root).resolve()
    output = (Path(output_root).resolve() if output_root is not None
              else root / "revision_20261022_r3/paper/endpoints")
    output.mkdir(parents=True, exist_ok=True)
    inputs = {}

    def source(relative):
        p = root / relative
        if relative not in inputs:
            inputs[relative] = dict(bytes=p.stat().st_size, sha256=sha(p))
        return p

    def csv(relative, **kwargs):
        return pd.read_csv(source(relative), **kwargs)

    def record(relative):
        return json.loads(source(relative).read_text(encoding="utf-8"))

    availability = csv("revision_20261022_r2/card001/availability_containment_by_interval.csv")
    assert (availability[["subset_violations", "equality_violations"]] == 0).all().all()
    rows, alarm_rows = [], []
    for cohort, label, allocated in [("internal", "DCLP3_test", 26), ("external", "REPLACE-BG", 225)]:
        events = csv(f"card001/results/{cohort}_event_table.csv")
        ledger = csv(f"card001/results/{cohort}_alarm_ledger.csv")
        people = csv(f"card001/results/{cohort}_participant_scores.csv")
        summary = csv(f"card001/results/{cohort}_summary.csv").set_index("policy")
        ready = availability.loc[availability.cohort.eq(label)].set_index("participant_id")
        assert len(ready) == allocated and ready.index.is_unique
        assert not events.duplicated(["id", "event_index"]).any()
        previous = None
        for policy in ["M15", "switch"]:
            p = people.loc[people.policy.eq(policy)].set_index("id")
            assert len(p) == allocated and p.index.is_unique and set(p.index) == set(ready.index)
            alarms = ledger.loc[ledger.policy.eq(policy)]
            assert not alarms.duplicated(["id", "alarm_index"]).any()
            for pid, person in p.iterrows():
                e = events.loc[events.id.eq(pid)]
                a = alarms.loc[alarms.id.eq(pid)]
                warned = truth(e[f"{policy}_warning"])
                assert len(e) == person.events == ready.loc[pid, "events"]
                assert int(warned.sum()) == person.matched
                assert 0 <= person.matched <= person.opportunities <= person.events
                assert person.opportunities == ready.loc[pid, "short_opportunity_events"]
                assert person.opportunities == ready.loc[pid, "switch_opportunity_events"]
                assert int(truth(e.M120_opportunity).sum()) == ready.loc[pid, "long_opportunity_events"]
                assert len(a) == person.alarms
                matched = a.loc[a.category.eq("matched"), "event_index"]
                assert matched.notna().all() and matched.is_unique
                assert set(matched.astype(int)) == set(e.loc[warned, "event_index"].astype(int))
                assert len(a) - len(matched) == person.unconfirmed
                close(person.observed_days*288, ready.loc[pid, "observed_bins"], atol=1e-6)
            totals = p[["events", "opportunities", "matched", "alarms", "unconfirmed", "observed_days"]].sum()
            for field, value in totals.items():
                close(value, summary.loc[policy, field])
            bins = int(ready.observed_bins.sum())
            row = decomposition(label, policy, allocated, int((p.observed_days > 0).sum()),
                allocated, totals.events, totals.opportunities, totals.matched, bins)
            close(row["observed_days"], totals.observed_days)
            current = (row["events"], row["opportunity_events"], bins)
            if previous is not None:
                assert current == previous, "Policy denominator/opportunity changed"
            previous = current
            rows.append(row)
            alarm_rows.append(alarm_row(label, policy, alarms, row["observed_days"],
                int(totals.alarms), int(totals.unconfirmed), row["warned_events"]))

    intervals = csv("card003/results/primary_intervals.csv")
    events = csv("card003/results/events_primary.csv")
    ledger = csv("card003/results/alarms_primary.csv")
    transfer = record("card003/results/transfer.json")
    assert len(intervals) == 276 and intervals.PtID.nunique() == 102 and intervals.pair_id.is_unique
    # Replay ledgers also contain the original low-coverage windows. The already
    # frozen primary interval list defines the analysis; do not add those rows.
    clinical_ledger_scope = dict(source_event_rows=len(events), source_alarm_rows=len(ledger),
        included_primary_intervals=len(intervals))
    included_pairs = set(intervals.pair_id)
    events = events.loc[events.pair_id.isin(included_pairs)].copy()
    ledger = ledger.loc[ledger.pair_id.isin(included_pairs)].copy()
    clinical_ledger_scope.update(included_event_rows=len(events), included_alarm_rows=len(ledger))
    assert not events.duplicated(["pair_id", "onset_index"]).any()
    assert not ledger.duplicated(["pair_id", "alarm_index"]).any()
    assert set(events.pair_id).issubset(set(intervals.pair_id))
    assert set(ledger.pair_id).issubset(set(intervals.pair_id))
    for table in (intervals, events, ledger):
        assert set(table.policy) == {"M15"} and set(table.variant) == {"primary"}
    assert not np.any(truth(events.warned) & ~truth(events.legal_opportunity))
    for row in intervals.itertuples():
        e = events.loc[events.pair_id.eq(row.pair_id)]
        a = ledger.loc[ledger.pair_id.eq(row.pair_id)]
        E, O, W = len(e), int(truth(e.legal_opportunity).sum()), int(truth(e.warned).sum())
        assert (E, O, W, E-W) == (row.events, row.opportunity_events, row.matched, row.unwarned)
        assert len(a) == row.alarms and a.category.eq("matched").sum() == W
        assert len(a)-W == row.unconfirmed
        close(row.observed_days*288, row.observed_bins, atol=1e-6)
    E, O, W = len(events), int(truth(events.legal_opportunity).sum()), int(truth(events.warned).sum())
    row = decomposition("CLVer", "M15", 102, 102, len(intervals), E, O, W, intervals.observed_bins.sum())
    for field, key in [("events", "events"), ("warned_events", "matched"),
                       ("unwarned_events", "unwarned"), ("observed_days", "observed_days")]:
        close(row[field], transfer[key])
    rows.append(row)
    alarm_rows.append(alarm_row("CLVer", "M15", ledger, row["observed_days"],
        int(intervals.alarms.sum()), int(intervals.unconfirmed.sum()), W))
    assert (E, O, W) == (408, 406, 228)
    assert alarm_rows[3]["unmatched_observed"] == 46744
    endpoints, alarms = pd.DataFrame(rows), pd.DataFrame(alarm_rows)

    # Algebraic contrasts from the existing quadratic coefficients/covariance.
    # The saved diagnostics already contain the five quantile contrasts.
    quadratic = record("card003/results/quadratic_diagnostics.json")
    model = record("card003/results/model_quadratic_AUC.json")
    primary = record("card003/results/model_primary.json")
    for name, value in model["coefficients"].items():
        close(value, quadratic["coefficients"][name])
    columns = quadratic["covariance_columns"]
    covariance = np.array(quadratic["child_cluster_covariance"])
    beta = np.array([quadratic["coefficients"][name] for name in columns])
    k, q = columns.index("exposure"), columns.index("exposure_squared")
    minimum, maximum = intervals.auc_nmol_l.min(), intervals.auc_nmol_l.max()
    close(minimum, model["exposure_min"]); close(maximum, model["exposure_max"])
    contrasts = []

    def contrast(x, role, quantile=None):
        vector = np.zeros(len(columns)); vector[k] = 1; vector[q] = 2*x/.20+1
        logrr = float(vector @ beta); se = float(np.sqrt(vector @ covariance @ vector))
        rr = math.exp(logrr)
        lo, hi = np.exp([logrr-1.95996398454*se, logrr+1.95996398454*se])
        left = intervals.auc_nmol_l.between(x-.10, x+.10)
        right = intervals.auc_nmol_l.between(x+.10, x+.30)
        return dict(selection=role, quantile=quantile, auc_from=x, auc_to=x+.20,
            inside_observed_support=bool(x >= minimum and x+.20 <= maximum),
            quadratic_rr=rr, quadratic_cluster_ci_lower=float(lo), quadratic_cluster_ci_upper=float(hi),
            loglinear_rr=primary["rr"], loglinear_cluster_ci_lower=primary["robust_95_ci"][0],
            loglinear_cluster_ci_upper=primary["robust_95_ci"][1],
            loglinear_primary_bootstrap_ci_lower=primary["bootstrap_95_ci"][0],
            loglinear_primary_bootstrap_ci_upper=primary["bootstrap_95_ci"][1],
            intervals_within_0_10_of_start=int(left.sum()), children_within_0_10_of_start=int(intervals.loc[left,"PtID"].nunique()),
            intervals_within_0_10_of_end=int(right.sum()), children_within_0_10_of_end=int(intervals.loc[right,"PtID"].nunique()))

    for saved in quadratic["local_per_0_20_contrasts"]:
        close(saved["auc_nmol_l"], intervals.auc_nmol_l.quantile(saved["quantile"]))
        r = contrast(saved["auc_nmol_l"], "existing_quantile_contrast", saved["quantile"])
        close(r["quadratic_rr"], saved["rr"])
        close(r["quadratic_cluster_ci_lower"], saved["robust_95_ci"][0])
        close(r["quadratic_cluster_ci_upper"], saved["robust_95_ci"][1])
        close(r["quadratic_rr"], model["local_per_0_20_rr"][str(saved["quantile"])] )
        contrasts.append(r)
    for x in [.2, .4, .6, .8, 1., 1.2]:
        contrasts.append(contrast(x, "round_grid_algebra_from_saved_fit"))
    nonlinear = pd.DataFrame(contrasts)
    assert nonlinear.inside_observed_support.all()

    endpoints.to_csv(output/"event_opportunity_warning_decomposition.csv", index=False)
    alarms.to_csv(output/"alarm_composition.csv", index=False)
    nonlinear.to_csv(output/"saved_quadratic_contrasts.csv", index=False)
    lines = ["# Event, opportunity and useful-warning accounting", "", "Author: UCP Technology LLC.", "",
        "E counts qualifying observed events, O counts those with at least one legal input opportunity, and W counts those with a useful warning. "
        "Unwarned events satisfy E−W=(E−O)+(O−W). Rates use the same observed-sensor-day denominator within each cohort. "
        "This accounting identity does not decompose regression coefficients or identify causal effects.", "",
        "| Cohort / policy | E | O | W | E−W | E−O | O−W | O/E | W/O | Unwarned with opportunity |",
        "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['cohort']} / {r['policy']} | {r['events']} | {r['opportunity_events']} | {r['warned_events']} | {r['unwarned_events']} | {r['without_opportunity']} | {r['opportunity_without_warning']} | {100*r['opportunity_fraction']:.3f}% | {100*r['conversion_given_opportunity']:.3f}% | {100*r['fraction_unwarned_with_opportunity']:.3f}% |")
    lines += ["", "| Cohort / policy | Observed days | Total events / 100 days | Unwarned / 100 days | No opportunity / 100 days | Opportunity but unwarned / 100 days |", "|---|---:|---:|---:|---:|---:|"]
    for r in rows:
        lines.append(f"| {r['cohort']} / {r['policy']} | {r['observed_days']:.9f} | {r['total_event_rate_per_100days']:.6f} | {r['unwarned_rate_per_100days']:.6f} | {r['no_opportunity_rate_per_100days']:.6f} | {r['opportunity_unwarned_rate_per_100days']:.6f} |")
    a = alarm_rows[3]
    lines += ["", f"The external switch issued {a['unmatched_observed']:,} unmatched alarms with observed follow-up over {a['observed_days']:.12f} observed days: **{a['unmatched_observed_per_day']:.9f}/day**. "
        f"The complete unconfirmed count is {a['unconfirmed']:,}, or {a['unconfirmed_per_day']:.9f}/day, including {a['late']:,} late and {a['unobservable_outcome']:,} unobservable-outcome alarms. "
        "The unmatched-observed point estimate alone exceeds the study-defined one-per-day cap. Alarm categories do not establish whether individual alarms were clinically unnecessary.", "",
        "# Existing quadratic AUC sensitivity", "", "The five quantile contrasts below were already saved with the original sensitivity fit. "
        "Each is a finite +0.20 nmol/L contrast with a child-cluster sandwich interval; no model was refit. "
        "The CSV also gives round-grid contrasts calculated from those same coefficients/covariance, and descriptive ±0.10 nmol/L neighborhood counts. "
        "Pooled support does not establish support in every covariate stratum.", "",
        "| Starting quantile | AUC from → to (nmol/L) | Quadratic RR | Cluster 95% interval |",
        "|---|---|---:|---|"]
    for r in contrasts[:5]:
        lines.append(f"| {r['quantile']:.0%} | {r['auc_from']:.6f} → {r['auc_to']:.6f} | {r['quadratic_rr']:.6f} | {r['quadratic_cluster_ci_lower']:.6f}–{r['quadratic_cluster_ci_upper']:.6f} |")
    lines += ["", f"The constant log-linear RR is {primary['rr']:.9f}, with cluster interval {primary['robust_95_ci'][0]:.9f}–{primary['robust_95_ci'][1]:.9f} and the original primary bootstrap interval {primary['bootstrap_95_ci'][0]:.9f}–{primary['bootstrap_95_ci'][1]:.9f}. "
        "The saved quadratic point contrasts decline more steeply at higher AUC; the median-AUC contrast is close to the constant log-linear estimate. "
        "Uncertainty includes one at the 5th, 25th and 95th starting quantiles. These observational sensitivity estimates leave the primary transfer gate and inference unchanged.", "",
        "Verification reconciled each warning-study participant and each clinical interval across event, alarm and aggregate ledgers; M15/switch opportunities additionally matched the preserved availability audit. "
        "All input hashes are checked unchanged after extraction. The warning-study event table contains the long-history opportunity flag; its short/switch opportunity counts are taken from participant ledgers and reconciled to the separate availability audit."]
    (output/"endpoint_decomposition_report.md").write_text("\n".join(lines)+"\n", encoding="utf-8")
    for relative, meta in inputs.items():
        p = root / relative
        assert p.stat().st_size == meta["bytes"] and sha(p) == meta["sha256"], f"Input changed: {p}"
    payload = dict(status="PASS", created_utc=datetime.now(timezone.utc).isoformat(),
        scope="Aggregate extraction and algebra from preserved event/alarm ledgers and model coefficients; no fitting or policy replay",
        original_inputs_unchanged=True, source_sha256=sha(Path(__file__)), inputs=inputs,
        clinical_ledger_scope=clinical_ledger_scope,
        checks=dict(warning_participant_ledger_reconciliations=502, clinical_interval_reconciliations=276,
            accounting_identity_all_rows=True, identical_policy_denominators=True,
            existing_quadratic_contrasts_reproduced=5, new_model_fits=0),
        outputs={name:sha(output/name) for name in ["event_opportunity_warning_decomposition.csv",
            "alarm_composition.csv", "saved_quadratic_contrasts.csv", "endpoint_decomposition_report.md"]})
    (output/"extraction_manifest.json").write_text(json.dumps(payload, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    print(endpoints[["cohort","policy","events","opportunity_events","warned_events","unwarned_events","without_opportunity","opportunity_without_warning","conversion_given_opportunity"]].to_string(index=False))
    print(json.dumps(dict(status="PASS", output=str(output), external_switch_unmatched_observed_per_day=a["unmatched_observed_per_day"]), indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, required=True)
    parser.add_argument("--output-root", type=Path,
                        help="Output directory; defaults to the paper/endpoints folder under the data root")
    args = parser.parse_args()
    run(args.data_root, args.output_root)
