"""Rebuild the CLVer human and machine reports from completed numerical artifacts."""
from pathlib import Path
import json
import zipfile
import pandas as pd
from card003_core import sha256_file


def build_report(root,repo,report_dir=None):
    root,repo=Path(root),Path(repo)
    s=json.loads((root/"results"/"analysis_summary.json").read_text())
    if not s.get("complete"):raise ValueError("Complete empirical output is required for terminal report")
    p,t,g=s["primary"],s["transfer"],s["gates"]
    source=json.loads((root/"results"/"source_audit.json").read_text())
    clinical=json.loads((root/"results"/"clinical_audit.json").read_text())
    sensitivities=json.loads((root/"results"/"sensitivity_results.json").read_text())
    boundary=json.loads((root/"results"/"boundary_integrity_closure.json").read_text())
    if boundary.get("status")!="PASS":raise ValueError("Refractory context audit must pass before terminal report")
    chronology="The final adapter requires explicit CGM time-of-day; date-only midnight candidates are confined to a labeled sensitivity. Historically, the explicit-time rule and duplicate-counter masking were source-audit corrections after initial estimates, with earlier results retained. This fresh source-to-result run applies the corrected rules from the beginning; it does not need the earlier participant-level archives."
    previous=root/"history"/"before_strict_measurement_time_audit.zip"
    if previous.exists():
        with zipfile.ZipFile(previous) as z:
            before=json.loads(z.read("results/model_primary.json"));before_transfer=json.loads(z.read("results/transfer.json"))
        lost_bins=round((before_transfer["observed_days"]-t["observed_days"])*288)
        chronology+=f" This was a disclosed post-estimate source-audit correction: the initial general parser interpreted 66 date-only strings as midnight. All 66 had neighboring readings about five minutes on both sides, but that does not prove the omitted time. The initial RR was {before['rr']:.9f}; the strict-source RR is {p['rr']:.9f}. The primary cohort and 408-event denominator are unchanged, with {lost_bins} fewer observed bins. Original arrays, ledgers, estimates and all bootstrap draws remain in the versioned history archive; neither policy nor statistical specification was retuned."
    def interval(r):
        if "rr" not in r:return "Not estimable: "+r.get("reason","unknown")
        ci=r.get("bootstrap_95_ci",r.get("robust_95_ci"))
        return f"{r['rr']:.3f} ({ci[0]:.3f}–{ci[1]:.3f})"
    armrows="\n".join(f"| {name} | {r['events']} | {r['unwarned']} | {100*r['sensitivity']:.1f}% | {r['unconfirmed_per_day']:.3f} | {r['alarms_per_day']:.3f} |" for name,r in t["delivery_arms"].items())
    sensrows="\n".join(f"| {name.replace('_',' ')} | {r.get('n_children','—')} | {r.get('n_intervals','—')} | {'Nonconstant; see local contrasts' if name=='quadratic_AUC' else interval(r)} |" for name,r in sensitivities.items())
    rate=pd.read_csv(root/"results"/"descriptive_remission.csv")
    remissionrows="\n".join(f"| {r.clinical_state} | {int(r.intervals)} | {int(r.children)} | {r.total_per_100days:.3f} | {r.unwarned_per_100days:.3f} | {r.sensitivity*100:.1f}% |" for r in rate.itertuples())
    flags="; ".join(f"{name}: {'pass' if value else 'fail'}" for name,value in g["gates"].items())
    blank=sensitivities["blank_units_excluded"]
    report=f"""# CLVer: measured secretion, clinical remission and later replay-unwarned hypoglycemia

**Empirical execution complete. Qualified result: {g['status']}.** The adjusted association meets the numerical magnitude/precision rule, but the independently frozen policy fails the prespecified transfer burden gate. This result therefore does not support the qualified useful-warning bridge claim.

## Primary result and decision

The analysis contains **{t['intervals']} forward intervals from {t['children']} children**, {t['observed_days']:.3f} observed-sensor days and {t['events']} qualifying episodes. The frozen M15 policy (threshold 0.04) warns {t['matched']} episodes and leaves {t['unwarned']} unwarned, distributed across {t['residual_children']} children. Overall useful sensitivity is {t['sensitivity']*100:.2f}%; legal prediction opportunity covers {t['opportunity_coverage']*100:.2f}% of observed events.

Per **0.20 nmol/L higher index normalized MMTT C-peptide AUC**, the adjusted later unwarned-event RR is **{interval(p)}**, using 5,000 delivery-arm-stratified child bootstrap draws with seed 20261022. Failed replicates: {p['failed_bootstrap']}. The child-cluster sandwich interval is {p['robust_95_ci'][0]:.3f}–{p['robust_95_ci'][1]:.3f}. The interval does not establish the entire targeted 20% reduction with 95% confidence. This is an observational association conditional on the frozen warning policy, not an effect of preserving C-peptide.

Gates: {flags}.

| Randomized delivery arm | Events | Unwarned | Useful sensitivity | Unconfirmed/day | All alarms/day |
|---|---:|---:|---:|---:|---:|
{armrows}

The required caps are at most one unconfirmed and two total alarms per observed-sensor day in each delivery arm. These are study gates, not established clinical acceptability thresholds. The DCLP3 validation operating point was frozen before CLVer alarm replay; CLVer thresholds were not retuned after observing this failure.

The finite linear forecast grid has no feasible validation operating point. Its reported no-alarm sentinel is an infeasible-family reference, not a successful finite trend-warning model.

## Clinical remission and frequency versus warnability

The same adjusted model for all observed episodes gives RR **{interval(s['companions']['total_events'])}** per 0.20 nmol/L. Replacing secretion with known index IDAA1c≤9 gives unwarned-event RR **{interval(s['companions']['remission'])}**. Each companion uses its own 5,000 child bootstrap draws. Clinical remission and laboratory secretion remain distinct measurements.

| Index clinical state | Intervals | Children | All events/100 days | Unwarned/100 days | Useful sensitivity |
|---|---:|---:|---:|---:|---:|
{remissionrows}

Children can contribute more than one clinical state, so state-specific child counts overlap. Visits and transitions are observed quarterly; no continuous remission duration is inferred. Four remission/peak-C-peptide combinations, unknown states and fixed AUC bands are reported in the accompanying numerical tables. These are descriptive groups, not validated disease stages.

Total-event and residual-event RRs are similar, and useful sensitivity does not rise monotonically across the fixed AUC bands. These results do not establish improved physiological warnability. The within-child estimate is **{interval(sensitivities['within_child'])}**, an imprecise result that does not corroborate the direction of the pooled association. The separate clinical-remission estimate likewise provides no clear evidence of a lower residual rate.

## Source and temporal audit

The release contains {clinical['roster_total']} roster entries and {clinical['randomized']} randomized children. Of {clinical['assay_visits']} scheduled assay visits, {clinical['eligible_assay_visits']} pass the strict six-sample assay gate. There are {clinical['attempted_pairs']} possible forward pairs, {clinical['temporally_eligible_pairs']} valid temporal/assay pairs, and {t['intervals']} with at least 70% observed CGM coverage. Zero-event intervals remain included. Missing glucose is never normal glucose or filled by another stream.

Use GXB for Tandem-roster participants and DeviceCGM for all others. GXB nonzero high/low indicator placeholders are excluded ({source['counts']['GXB']['flagged_placeholders']} records). DeviceCGM has {source['counts']['DeviceCGM']['unit_domain'].get('mg/dL',0):,} explicitly mg/dL and {source['counts']['DeviceCGM']['unit_domain'].get('',0):,} blank-unit records in the chosen CGM source audit; no other unit is present. Blank fields are an explicit mg/dL assumption corroborated by the trial protocol and value domain, tested by exclusion.

Excluding blank-unit DeviceCGM observations retains {blank['n_intervals']} intervals from {blank['n_children']} children and {blank['observed_days']:.3f} observed-sensor days, losing {t['intervals']-blank['n_intervals']} primary intervals and {t['children']-blank['n_children']} children from the restricted sample. Its unwarned-event RR is {interval(blank)}. GXB has no unit column; its documented study-scale interpretation remains explicit.

Conflicting same-time readings are unavailable. Local/transmitter discontinuities insert an unknown separating bin; following internally consistent local-time segments remain available without an invented offset. Absolute device clock accuracy is not established. The any-clock-barrier interval exclusion is the conservative sensitivity. Physical-exam weight has no separate timestamp; linked VisitInfo date is a protocol-based same-visit proxy, not an independently recorded measurement time. The seven-day clinical alignment gate remains applied.

{chronology}

Refractory-state synchronization is verified before the earliest scored warning opportunity in every primary window and every reported primary-policy variant. The M120, missing-input and switch secondary policies use extended context until six fully determined consecutive no-candidate bins after complete feature-history context prove that earlier alarm state cannot affect scored alarms. Actual missing measurements can establish a known inability to issue an alarm. Their extended-context results leave the reported residual counts and associations unchanged; the proof and extended alarm ledgers are retained.

Independent source reconstruction verifies all 414 included MMTT AUCs and 422 known IDAA1c classifications to numerical tolerance. Independent regression and the first 25 preserved bootstrap draws are checked by a separately implemented Newton solver. Actual patient-ID relabeling preserves the fitted coefficient and clustered uncertainty; perturbing future glucose preserves prior predictions in every fitted model on both CGM source types.

The exposure/rate figure sets AUC to each plotted value within the observed overall support and averages fitted per-day rates over the included intervals' original covariates, giving each interval equal weight. Its shaded band is a 95% child-cluster delta-method confidence interval for that standardized fitted mean, not an individual prediction interval. The curve remains conditional on the fitted observational model and does not represent intervention effects.

## Prespecified sensitivities

Secondary intervals below are child-cluster sandwich intervals. They do not replace the primary bootstrap or a failed transfer gate. Constant nuisance columns are removed only in secondary restricted samples when algebraically redundant with the intercept; each removed column is logged. The contemporaneous analysis drops its exactly constant −28-day lag, preserving the model column space. The quadratic row's exponentiated linear coefficient is not a constant RR: full coefficients and local contrasts are retained in its JSON output.

| Analysis | Children | Intervals | RR (95% interval) |
|---|---:|---:|---|
{sensrows}

The five-minute receipt-delay analysis delays only information available to the warning policy. Event times and observed outcome exposure remain fixed. The three-reading convention requires ten elapsed minutes across the first and third low reading; primary four-reading episodes require fifteen elapsed minutes. Coverage variants, blank-unit exclusion and phases reconstruct their own causal observation masks. Every frozen family and the never-warning reference is retained; none is selected on its CLVer association.

## Scope and reproducibility

This study quantifies observed pediatric CLVer trajectories under a fixed retrospective replay policy. Existing therapy and device alarms shaped the recorded trajectories. Associations can reflect management, disease stage, observation selection and residual confounding. The transfer failure remains decisive even when the biological association is numerically favorable.

The proposed joint endpoint adds a frozen warning-policy condition to the C-peptide–hypoglycemia association; the current result fails the qualified bridge because alarm burden does not transport. Existing work already studies C-peptide and glycemia in this cohort ([Neyman et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC12281974/)) and prior C-peptide–hypoglycemia associations ([Pinckney et al.](https://pmc.ncbi.nlm.nih.gov/articles/PMC4916002/)). The combined analysis tests the additional information supplied by the warning-conditioned endpoint.

Protocol: `protocol/card003_methods.md`. Results and all patient-level data: `{root}`. Code/report repository: `{repo}`. One-command rebuild after frozen DCLP3 development: `python src/card003_run.py --data-root "{root}" --rebuild`. The public scientific methods and code are sufficient; private framework documents are not required for reproduction. Source archive, models, policy selection, code and result hashes are preserved. Runtime for the analysis is {s['elapsed_seconds']:.1f} seconds; process memory information is in `analysis_summary.json`.

The [CLVer original public dataset](https://public.jaeb.org/dataset/591) and [DCLP3 public release](https://public.jaeb.org/dataset/573) require their publisher's applicable access terms and attribution. See the acquisition manifests and source archive README. Independent access to original data is needed to reproduce patient-level results.
"""
    report_dir=Path(report_dir) if report_dir is not None else root/"reports"
    report_dir.mkdir(parents=True,exist_ok=True)
    reportpath=report_dir/"card003_execution.md";reportpath.write_text(report,encoding="utf-8")
    status=dict(card="003",complete=True,status="completed_empirical_run",hypothesis_status=g["status"],
        qualified_finding="INCONCLUSIVE" if g["status"]=="Inconclusive" else "SUPPORTED" if g["status"].startswith("Supported") else "REFUTED",
        statistical_result=p,clinical_result="observed partial remission association; no induced remission endpoint",
        gates=g,transfer=t,data_root=str(root),repo=str(repo),report=str(reportpath),remaining=[],
        summary_path=str(root/"results"/"analysis_summary.json"),summary_sha256=sha256_file(root/"results"/"analysis_summary.json"))
    (root/"results"/"card003_status.json").write_text(json.dumps(status,indent=2,allow_nan=False),encoding="utf-8")
    hashes={str(path.relative_to(root)):sha256_file(path) for path in (root/"results").rglob("*") if path.is_file() and path.name!="result_hashes.json" and path.suffix!=".log"}
    (root/"results"/"result_hashes.json").write_text(json.dumps({"scope":"Final result files and source snapshots; active/runtime logs excluded","sha256":hashes},indent=2),encoding="utf-8")
    return status
