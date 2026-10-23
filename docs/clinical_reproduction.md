# CLVer original-source reproduction

This component regenerates the clinical study from the original CLVer archive and newly trained DCLP3 warning models. It does not use released aggregate results, retained participant tables, saved clinical fits, or an author's data directory as analysis inputs. Aggregate release evidence is used only by the separate comparison step.

## Inputs and execution order

Obtain `CLVerPublicDataset.zip` through the Jaeb CLVer data-access page, <https://public.jaeb.org/dataset/591>, completing the provider's registration and terms. The analysis uses the 86,726,379-byte archive with SHA-256:

```
8421a9a9ff7dd1485ce0e4bdb984b5ea037300a0e7216f9f445334e56d32d385
```

The warning component must first regenerate its DCLP3 development models and complete validation grid. Its output directory supplies `calibration.json`, the three fitted model artifacts, and their source/split provenance. CLVer policy selection independently checks all 10,772 validation settings and uses the clinical protocol's selection rule. It does not import the warning study's superiority decision.

Use the release's pinned Python environment and run from the repository root. `INPUTS` and `RUN` below are absolute directories selected by the reproducer; no particular drive, account name, or private framework is required. Keep `RUN` separate from both the code checkout and the original-input directory.

```sh
python -B src/reproduce_clinical.py --run-root RUN --archive INPUTS/CLVerPublicDataset.zip --warning-root RUN/card001 --stage complete
```

For overlapping source preparation with warning-model training:

```sh
python -B src/reproduce_clinical.py --run-root RUN --archive INPUTS/CLVerPublicDataset.zip --stage prepare
python -B src/reproduce_clinical.py --run-root RUN --warning-root RUN/card001 --stage analyze
python -B src/reproduce_clinical.py --run-root RUN --stage revision
python -B src/reproduce_clinical.py --run-root RUN --stage verify
```

The `prepare` stage copies the checksum-matching original ZIP into `RUN/card003/raw`, reconstructs the source tables and grids, and runs independent clinical and clock checks. The original input is checked unchanged. `analyze` verifies the prepared-file and source-code hashes, freezes the independently selected policy, runs all primary/companion/sensitivity analyses and replay-boundary corrections, and generates the clinical report. `revision` runs the declared exploratory diagnostics, independent verification, and publication figure. `verify` repeats independent source, primary-estimate, bootstrap-subset, clock, refractory, and revision checks. Each subprocess has its own log in `RUN/card003/logs`.

Use a fresh run directory for a full rerun. Source preparation refuses an occupied clinical study directory; analysis refuses an existing frozen policy or completed analysis; the exploratory revision refuses to overwrite an existing analysis commitment. Other study components may already exist beneath the shared run root. A failed stage leaves its log and partial artifacts available for diagnosis; elapsed time is not evidence of completion. Only successful final status files establish completion. Runtime caches, temporary files, participant records, predictions, fits, bootstrap draws and figures are all written beneath `RUN`.

## Source rules and historical corrections

`protocol/card003_methods.md` and `protocol/revision_clinical.md` specify the complete estimands, adjustment, gates, selection and sensitivity analyses. The executable source adapter is `src/card003_data.py`.

- Normalized stimulated C-peptide AUC uses the required mixed-meal assay times, fixed trapezoidal quadrature divided by 120 minutes, and the assay/fasting/date eligibility rules. The clinical remission measure is IDAA1c, kept separate from laboratory secretion.
- Weight lacks an independent measurement date in the release. Its linked visit date is an explicit protocol-based proxy, subject to the seven-day clinical linkage rule. This is not a recovered actual weight timestamp.
- Blank `DeviceCGM` units use the stated mg/dL assumption. The exclusion sensitivity removes these measurements and recomputes coverage and outcomes; nonblank units and value-domain evidence are audited.
- The primary excludes all 66 date-only timestamp candidates. Their interpretation as midnight remains a separate sensitivity. Three same-time groups with conflicting internal counters are unavailable under a source-independent grouping rule.
- Clock discontinuities become unknown separating bins; subsequent consistent recorded-local-time segments remain available. This does not assert that absolute clock time was externally verified.
- Causal grid construction never uses a future sensor record. Source conflicts and barriers remain unavailable. Clinical intervals and alarm/event ledgers are regenerated from these grids.
- Refractory state is synchronized before scored alarms using six fully determined consecutive no-candidate bins after complete feature-history context. Missing inputs can make a bin a known no-candidate; the witness need not consist of fully observed glucose. Longer-history companion policies receive the necessary extra context.

The strict date/counter and replay-context rules were source/replay corrections during the original research run. A fresh reproduction starts with the corrected rules. It need not import the earlier participant-level outputs to reconstruct final results. The midnight sensitivity preserves the alternative date-only interpretation; no warning threshold is retuned in CLVer.

## Principal claims and result map

Paths in this table are relative to `RUN`. The main rate contrasts are per **+0.20 nmol/L** earlier normalized AUC. A completed association estimate does not override the warning-transfer burden gate.

| Claim or diagnostic | Inputs and computation | Generated evidence |
|---|---|---|
| Cohort flow, assay/remission linkage and sensor exposure | Original roster, visit, sample, insulin, examination and CGM tables; `card003_data.py`, independent `verify_card003.py` | `card003/results/interval_audit.json`, `independent_clinical_verification.json`; canonical tables in `card003/processed` |
| Independently frozen clinical warning policy | New DCLP3 validation grid/models; `card003_policy.py` | `card003/models/frozen_policy.json`, binding hashes and paths to the newly trained warning models |
| Warning transfer, observed residual burden and qualified hypothesis gates | Frozen policy and causal clinical windows; `card003_analyze.py`, `card003_core.py` | `card003/results/transfer.json`, `primary_gates.json`, event/alarm/interval ledgers |
| Primary unwarned-event association, all-event companion and remission association | Prespecified PPML adjustment and whole-child bootstrap; 5,000 draws each, seed 20261022 | `card003/results/model_primary.json`, `model_total_events.json`, `model_remission.json`; corresponding `bootstrap_*.csv.gz` |
| Source, timing, policy, coverage, nonlinear, management and within-child sensitivities | Same frozen policy or prespecified companion policies, all source rules retained; clustered sandwich uncertainty for secondary fits | `card003/results/sensitivity_results.json`, `sensitivity_matrix.csv`, `invariance_diagnostics.json`, `quadratic_diagnostics.json`, boundary-extension outputs |
| Exact-sample pooled versus within-child estimates, within-child information, Mundlak model | Original freshly regenerated primary interval set; `revision_clinical.py` | `revision_20261022_r2/card003/results/same_sample_pooled.json`, `auc_variation.json`, `mundlak.json` |
| Conditional unwarned odds among observed events | Event-level logistic likelihood and equivalent grouped-binomial check; child-cluster uncertainty | `revision_20261022_r2/card003/results/conditional_warning_logistic.json` |
| Unwarned versus total association comparison | Same 5,000 delivery-arm-stratified child resamples for both models, seed 20261023 | `revision_20261022_r2/card003/results/paired_comparison.json`, `paired_child_bootstrap.csv.gz` |
| Never-warning identity and joint empirical/fitted rates | Assert identical outcomes/offsets before reusing all-event estimate; empirical band rates with paired bootstrap; standardized fitted means with cluster delta-method intervals | `revision_20261022_r2/card003/results/never_warning_identity.json`, `joint_auc_bands.csv`, `fitted_curve_events.csv`, `fitted_curve_unwarned.csv` |
| Publication Figure 2 | Saved newly generated bands, curves and observed AUC support; `revision_clinical_figure.py`, no additional refit | `figures/secretion_rates.png` |
| E/O/W accounting and nonlinear finite contrasts | Newly generated warning/clinical ledgers plus warning availability diagnostics; `paper_endpoint_decomposition.py` after both components finish | `revision_20261022_r3/paper/endpoints/event_opportunity_warning_decomposition.csv`, `alarm_composition.csv`, `saved_quadratic_contrasts.csv` |

The rate models retain zero-event intervals. Conditional logistic inference is explicitly conditional on an observed qualifying event and yields an **odds ratio**, not a rate ratio or causal warning effect. A paired ratio-of-rate-ratios interval containing one does not establish equivalence. Never-warning and all-event estimates are mathematically identical once the asserted endpoints and exposure match. The model curves show uncertainty in standardized fitted means, not individual prediction intervals. TN10 is a separate component and is never numerically combined with these clinical rates into a biological treatment effect.

## Verification

The independent source verifier reconstructs clinical values directly from original archive tables without importing the production derivation. The independent clock verifier compares literal source timestamps with retained flags/epochs and checks canonical-grid causality. The primary verifier uses an independent Newton fit and explicit cluster sandwich and reconstructs a preserved bootstrap subset. The replay verifier checks the state-synchronization witness for each clinical interval/variant. The revision verifier independently reconstructs coefficients and uncertainty.

The new portability tests verify separation of code, original inputs and generated outputs without requiring a specific drive. The existing domain tests cover assay gates, temporal alignment, episode rules, replay, finite-model handling, policy selection and revision estimands.

The completed original-source run on 23 October 2026 regenerated the full clinical analysis, source and replay corrections, revision diagnostics, and Figure 2. The source verifier reconstructed 414 eligible AUC visits and 422 known remission visits. It checked 301 interval grids, 5,935,232 explicit-time rows, 66 preserved date-only candidates and three same-time counter-conflict groups, with zero observed noncausal/conflict/barrier/date-only-source bins. The [scientific verification record](../provenance/scientific_verification.json) retains the completed run and its result comparisons.

The [current repository verification record](../provenance/verification.json) identifies the checkout checks performed during maintenance. Those checks are distinct from rerunning clinical inference; see [verification scope](VERIFICATION.md).
