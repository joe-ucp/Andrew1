# Retained insulin secretion and later hypoglycemia outside warning coverage

## Question and interpretation

This observational study asks whether laboratory-measured endogenous insulin secretion predicts a lower subsequent rate of significant low-CGM episodes left without a timely warning by an independently selected warning policy. It uses CLVer for biological and clinical measurements and DCLP3 for all warning-policy development. The trials' participants are never joined across studies.

The exposure is the mean concentration represented by the two-hour mixed-meal tolerance test (MMTT) C-peptide AUC. The outcome is a subset of episodes identified from recorded CGM, not hypoglycemia prevented by a deployed intervention. A replay-unwarned event does not establish absence of symptoms, caregiver awareness or an actual device alarm.

Clinical partial remission is assessed separately with insulin-dose-adjusted HbA1c (IDAA1c). Observing remission at a visit does not demonstrate that treatment induced or prolonged it. The separate cohorts estimate disease delay and replay-defined warning coverage; they do not identify a combined teplizumab-plus-warning effect or convert disease-delay time into prevented hypoglycemia.

## Primary hypothesis and decision rule

The primary association is the adjusted rate ratio (RR) for later unwarned episodes per **0.20 nmol/L higher index normalized MMTT C-peptide AUC**. The target is RR at or below 0.80.

Qualified support requires all of the following:

1. The point estimate is at most 0.80, the lower endpoint of the child-bootstrap 95% interval is at most 0.80, and the upper endpoint is below 1.00. Separately state whether the entire interval is below 0.80; the primary support rule alone does not establish the full magnitude with 95% confidence.
2. CLVer useful-warning sensitivity is at least 30% overall. Within each randomized delivery arm, unconfirmed alarms are at most one and all alarms at most two per observed-sensor day.
3. At least 60 children and 200 eligible intervals contribute. At least 100 unwarned episodes occur across at least 30 children.
4. Temporal ordering, source identity, event accounting and absence of future-dependent predictions pass verification. Credible preprocessing corrections do not reverse the interpretation.

A lower confidence endpoint above 0.80 excludes the targeted 20% reduction under the declared model, while allowing smaller associations. A failed transfer, precision, identification or integrity requirement gives an inconclusive qualified result. A transfer failure does not show that biological preservation lacks clinical value. These thresholds are design choices, not clinical acceptability standards or a formal power calculation.

## Sources and provenance

Use the original [CLVer public archive](https://public.jaeb.org/dataset/591) and [DCLP3 release 3](https://public.jaeb.org/dataset/573), with their READMEs, glossaries and protocols. Complete the publisher's applicable access and usage requirements. Preserve original bytes, URLs, acquisition dates, byte counts, SHA-256 digests and required study attribution. Source data are kept outside the code repository and excluded from redistributable packages unless allowed by their terms.

The committed CLVer archive digest is `8421a9a9ff7dd1485ce0e4bdb984b5ea037300a0e7216f9f445334e56d32d385`. A different release requires an explicit version audit. Source member digests are recorded in the acquisition code and checked independently.

The analysis needs CLVer `PtRoster`, `SampleResults`, `InsulinUse`, `DiabPhysExam`, `VisitInfo`, `DeviceCGM` and `DevicePumpCGMDATAGXB`. DCLP3 uses original `DexcomClarityCGM_a` and roster/visit records. Original receiver and pump streams must not be merged to manufacture primary predictor availability.

All data, intermediate arrays, model files, ledgers and numerical results reside in the user-configured output tree, separate from the read-only source inputs and release directory. Reproduction uses the original provider files, these scientific methods and the pinned software environment.

## Laboratory reconstruction

Include randomized CLVer participants. The prognostic index visits are weeks 13, 26 and 39. Require one complete numeric nonnegative C-peptide result at each nominal minute 0, 15, 30, 60, 90 and 120. Collapse analytically identical duplicates. Conflicting samples or an incomplete six-sample set make that visit unavailable.

All six measurements must have compatible collection dates and concentration units of nmol/L. Require complete laboratory result status and documented fasting. The original source repeats `STASampleFastStatus` on each sample; the raw adapter conservatively requires literal `Fasting` on all six samples. Explicit nonfasting or unknown status on any required sample excludes the primary exposure. The generic assay kernel checks the baseline field, while this stronger source-specific requirement is enforced by the raw adapter.

Temperature flags alone do not invalidate an otherwise eligible primary visit. A separate sensitivity requires all six temperature indicators to be `OK`.

Compute the trapezoidal integral over the six nominal sampling times and divide by 120 minutes. The resulting exposure is a mean concentration in nmol/L, not an unnormalized nmol·min/L area. Preserve all contributing values and the peak concentration. Missing assays are not imputed. C-peptide is not a direct measurement of beta-cell mass or autoimmune control.

## Clinical partial remission

Calculate `IDAA1c = HbA1c [%] + 4 × insulin dose [U/kg/day]`; a value at or below 9 defines the declared clinical state. Obtain central HbA1c from `GLYHB`. Link insulin reports and examination weight through participant and `ParentLoginVisitID`, then match the scheduled visit and inspect dates.

Insulin is reported as absolute daily units and must be divided by body weight in kilograms. Use kilograms directly or convert pounds with exactly 0.45359237 kg/lb. Require positive known weight, nonnegative known daily insulin and a valid HbA1c. Each component's date must be within seven calendar days of the index MMTT. Unresolved duplicate clinical joins, unknown components or incompatible dates give an unknown remission state, not nonremission.

Unknown clinical remission does not remove an otherwise valid C-peptide interval from the primary model. Report the clinical state's visit frequencies, observed transitions and disagreement with peak stimulated C-peptide above 0.30 nmol/L. The peak threshold is a companion laboratory comparison, not a replacement remission definition. Do not infer daily state or remission duration between quarterly observations.

## Later outcome windows

Pair visits 13→26, 26→39 and 39→52. The target visit anchor is its unique baseline C-peptide collection date. Its numeric assay value is not an exposure or warning feature. An absent or ambiguous anchor, or disagreement with the clinical visit date by more than seven days, excludes the pair with a recorded reason.

The scored interval is the half-open 28-day period ending at midnight at the target date: `[target midnight −28 days, target midnight)`. Its start must follow the index MMTT collection date strictly. Do not reposition intervals according to event counts. If two outcome intervals overlap for a child, exclude both affected pairs in the primary analysis.

The construction conditions on attendance with an identifiable later MMTT date. Report this selection and the laboratory/clinical characteristics of excluded intervals. Include valid zero-event intervals.

## Source, clock and observation rules

Use `DevicePumpCGMDATAGXB` for roster participants assigned `DeviceTrtGroup=Tandem`; use `DeviceCGM` for everyone else. A missing chosen stream stays missing. Report coverage by randomized management arm and source.

Use records marked usable. DeviceCGM additionally requires `RecordType=CGM`. GXB requires numeric glucose and `HighLowIndicator=0`; nonzero-flag zero placeholders are unavailable, not physiological glucose zeroes. Do not clip otherwise valid numerical extremes.

Audit glucose units. Blank unit fields may use the study's implicit mg/dL encoding only with explicit documentation; repeat the analysis excluding blank-unit records. Evidence of another encoding requires a documented conversion or exclusion.

Interpret recorded local timestamps consistently. Examine transmitter information for resets and incompatible clock increments, and exclude unresolved clock segments. Record exactly how each clock exclusion changes source availability. The DCLP3 development adapter uses adjusted timestamps when present and original timestamps otherwise. Differences in device measurement time and actual receipt time remain a limitation of recorded-time replay.

Drop identical duplicate records. Conflicting values at an identical source time are unavailable. On a five-minute grid anchored at local midnight, a decision bin can use only the latest valid source observation in `(t−5 minutes,t]`. Keep original source and observation time. Do not interpolate, carry values forward across bins, fill gaps or use another source to fill a primary gap.

The primary 28-day interval requires at least 70% observed bins. Coverage is determined before warning replay and independently of the number of low events. The burden denominator is observed bins multiplied by 5/1,440 days. Missing CGM is not normal glucose; inference concerns observed intervals and episodes, not unseen lows over all calendar time.

Retain at least 150 minutes of preceding grid positions to initialize prediction, episode and alarm state. Missing lookback remains missing. Retain at least 45 minutes of following positions so alarms near the end can be classified. Record boundary losses and use the same state initialization for separated intervals.

## Episode and alarm accounting

Four consecutive observed five-minute readings below 54 mg/dL establish an episode; the first reading is its onset. Four consecutive readings at or above 54 establish recovery and rearm detection. Begin unarmed. A missing bin resets state to unknown until recovery is established. A sustained low first observed after a gap has an unknown onset and is reported in a separate censoring ledger, not manufactured into a new qualifying event. A confirmed onset survives subsequent censoring of episode duration.

Derive episodes independently of every prediction policy. Every observable qualifying onset in an eligible interval remains in the denominator, including onsets without a legal input opportunity.

A useful alarm is issued 15–30 minutes before onset, including the endpoints. Match alarms in time order to the earliest eligible unmatched episode, with at most one match per alarm and episode. Enforce a 30-minute refractory period after every issued alarm. No predictive alarm may be issued with missing current glucose or current glucose below 54.

For onset coverage, allow matching alarms in the 30 minutes preceding the scored interval. For alarm burden, count alarms inside the scored interval and allow matching to episodes up to 30 minutes afterward with enough additional readings to confirm them. Unobservable outcomes are unconfirmed/unknown, not physiological false positives. Preserve matched, late, duplicate, unmatched-observed and unobservable-outcome alarm categories.

For each interval report total events, warned events, unwarned events, observed-sensor days, opportunity coverage, useful sensitivity and alarm burden. Split unwarned events into no legal prediction opportunity and available prediction without a qualifying alarm. The identity `unwarned rate = total rate × (1 − sensitivity)` holds only with identical event sets and exposure denominators; it does not identify the residual association from either factor alone.

## Warning development entirely outside CLVer

Use DCLP3 CGM from randomization through day 182, truncated by documented withdrawal/final follow-up. Use its original Dexcom export, the causal grid, unchanged event definition and explicit clock audit.

Force previously inspected participant IDs `1,3,4,5,6,7,8,10,12,13,14,15,19,20,22,32,36,39,50,51,53,100,101,102,103,104,105,106,107,108,109` into training. For all others, canonicalize the integer ID, hash `t1d-card-001:` followed by that ID with SHA-256, convert the first eight hex digits to an integer and reduce modulo 10,000. Assign below 6,000 to training, 6,000–7,999 to validation, and at least 8,000 to held-out testing. Preserve entire people in one partition.

The model family consists of a two-hour classifier (M120), a 15-minute classifier (M15), an availability switch using M120 when available and M15 otherwise, a missing-input classifier, and a 15-minute linear trend extrapolated 30 minutes ahead. M120 requires 25 consecutive bins and M15 four. The missing-input classifier may run whenever M15 is available. All policies abstain when their required inputs are unavailable.

Classifiers predict episode onset 15–30 minutes ahead. A negative label requires sufficient observed outcome context through 45 minutes to exclude the target onset. Omit and count uncertain labels. Future readings may define labels but never predictor features.

Features are five-minute glucose lags, current value, ordinary least-squares slope against elapsed minutes, mean, population standard deviation, minimum and maximum. M120 also uses the summaries over trailing 15, 30, 60 and 120 minutes. The missing-input alternative uses the M120 lag positions, observation indicators, time since observation, observed-only summaries and counts; undefined summaries remain missing. Exclude patient identifiers, whole-record normalization, future insulin/carbohydrates and outcome-derived features.

Use histogram gradient boosting with maximum leaves 7 or 15, L2 strength 0 or 1, learning rate 0.05, 200 iterations, at least 50 samples per leaf, seed 20261022 and disabled random-row early stopping. Weight training records so each contributing participant has equal total weight.

Assign grouped cross-validation folds by sorting **all assigned training participant IDs**, including participants without eligible model rows, by SHA-256 of `t1d-card-003-fold:` followed by canonical ID, then cycle through folds 0–4. Keep each person's rows together. Choose hyperparameters by participant-weighted log loss; exact ties prefer fewer leaves and then larger L2. Publish folds, row exclusions, software versions and fitted settings.

On DCLP3 validation only, sweep classifier probabilities from 0 to 1 in 0.01 increments plus an explicit never-alarm choice. Sweep linear forecast cutoffs 40–100 mg/dL in increments of one. Evaluate all switch threshold pairs. Replay actual refractory and matching rules at every candidate point.

For each family member maximize sensitivity subject to at most one unconfirmed and two total alarms per observed-sensor day. Ties prefer fewer unconfirmed alarms, then fewer total alarms, then a higher threshold. For switch threshold ties compare M120 followed by M15. The implementation compares physical forecast thresholds for the linear member rather than the negated internal score. Once each family member's threshold tie is resolved, exact cross-family performance ties use this priority: linear, M15, M120, missing-input, switch. Different score units are not numerically ranked across families.

Select the highest-sensitivity feasible family member as the primary policy. This selection has **no additional switch-versus-comparator burden-difference constraint**. Save every family member for sensitivity analyses. Publish model, threshold, source, environment and code hashes before any CLVer alarm replay. CLVer never supplies fitting, threshold selection, calibration or subgroup choices. DCLP3 held-out outcomes are not used to select this policy.

## Primary association model and uncertainty

For interval `j` in child `i`, let `Uij` be unwarned events, `Dij` observed-sensor days and `Cij` index normalized AUC. Fit the Poisson log-link mean model by maximum pseudolikelihood:

`log E(Uij|Xij) = log(Dij) + intercept + beta × Cij/0.20 + age + target visit + management + cohort + lag`.

Age is continuous, centered at 12 years. Target visit is categorical with levels 26, 39 and 52 weeks. Management is usual care, intensive/Tandem or intensive/Medtronic. Cohort is A or B and does not identify verapamil/placebo. Lag is days from index MMTT to outcome-window start, centered at 60 and scaled by 30. Essential covariate uncertainty excludes the affected interval with a reason. No automated variable selection, outcome-selected transformation, primary interaction, zero inflation or unplanned penalty is used.

Use 5,000 child-cluster bootstrap replicates with seed 20261022, stratified by randomized delivery arm. Resample complete child interval sets and refit. Report the percentile interval for `exp(beta)` and child-cluster sandwich uncertainty. Preserve all failed replicates. More than 5% failures, rank deficiency, nonidentified Poisson means or nonfinite estimates/intervals give an inconclusive inferential result. The implementation additionally checks finite-MLE existence by separation-direction linear programming and information-matrix rank; it does not hide separation with a penalty.

Uncertainty is conditional on the independently frozen warning policy. It does not describe retraining variation, deployment effects or the causal effect of raising C-peptide.

## Mandatory companion analyses

1. Refit the same adjustment model for all observed qualifying events. Compare total and unwarned rates, sensitivity, opportunity coverage and burden together. If residual rates differ while conditional warnability is unchanged, describe a frequency association.
2. Replace continuous secretion with clinical remission at the index visit among records with known IDAA1c. Report its RR and absolute events per 100 observed-sensor days, state frequencies and observed transitions.
3. Cross-tabulate clinical remission against peak C-peptide above 0.30 nmol/L, including missingness. Report residual burden in the four measured combinations; suppress unstable subgroup fits.
4. Add index HbA1c and insulin U/kg/day to the primary model as a different conditional prognostic adjustment, not a purified causal model.
5. Estimate a within-child sensitivity with child intercepts, index AUC, target-visit indicators and lag. Children with all-zero outcomes do not identify the conditional count slope; report informative children and within-child exposure variation.
6. Repeat residual associations for every frozen policy, including the switch, strongest nonswitch member and never-warning reference. Do not choose the most favorable CLVer association. Management-arm analyses are exploratory and do not establish treatment synergy.
7. Repeat using same-visit secretion and its preceding 28-day interval, explicitly contemporaneous rather than prognostic.
8. Show exposure distribution and model-predicted rates only within observed support. Describe empirical rates in fixed AUC bands below 0.20, 0.20–below 0.40, 0.40–below 0.60 and at least 0.60 nmol/L. If the relation is visibly nonmonotonic, add the prespecified quadratic-AUC sensitivity and explain the limitation of one RR.

None of these analyses replaces a failed primary conjunction.

## Robustness and damage tests

Run and report each of these prespecified alternatives without choosing a favorable operating point:

- Primary non-Tandem DeviceCGM source subset.
- Minimum 120 observed hours in 28 days and minimum 90% grid coverage, alongside primary 70% coverage.
- Exclusion of blank-unit CGM records.
- Temperature-OK assays only.
- Five grid phases: 0, 1, 2, 3 and 4 minutes.
- A five-minute reporting delay.
- Three consecutive low readings as the alternate event convention, with its different elapsed duration stated.
- Separate leave-one-child-out exclusions for the lowest-secretion child and highest-event-count child; report whether a few children determine the slope.
- Every independently frozen warning family member and the never-warning reference.
- Conditional HbA1c/insulin adjustment, within-child model, contemporaneous window and quadratic-AUC sensitivity described above.

Implementation validation includes independent episode labelers, one-to-one laboratory joins, AUC dimension checks, insulin/weight conversion fixtures, rejection of flagged zero placeholders, exact exposure accounting, prefix invariance when future glucose changes, and invariance to participant-ID relabeling in outcome calculations. Separately preserve the deliberately ID-based development splits. Verify the arithmetic counterexample showing that fewer total lows need not imply fewer unwarned lows. Software fixtures are never patient evidence.

## Required evidence package and completion boundary

Publish cohort flow, assay/remission derivations, interval/source coverage, exposure distribution, reason-coded clinical/source/time exclusions, event and alarm ledgers, transfer table, primary RR and interval, remission comparisons, sensitivity matrix, bootstrap draws/failures, software/code/source hashes and runtime/memory information. Include a reproducible command and study attribution. Report the primary result as supported in this cohort and policy, targeted association refuted, or inconclusive, with the conditions stated.

The authorized CLVer archive is now supplied and source reconstruction has been executed. Strict assay/time/coverage gates retain 276 intervals from 102 children. Independent original-table reconstruction agrees on all 414 eligible MMTT AUCs and 422 known clinical remission states. DCLP3 development used the complete committed grid and independently selected M15 at threshold 0.04 for the clinical study before CLVer alarm replay. The empirical primary association and 5,000-draw uncertainty have been executed; the terminal numerical report and sensitivity ledger record the resulting evidence. Completion and qualified scientific support remain separate: an alarm-burden transfer failure cannot be repaired by a favorable association.

## Recorded implementation clarifications

The source has no independent physical-exam date. `VisitInfo.VisitDt`, linked by child and `ParentLoginVisitID`, is used as a protocol-based same-visit weight-date proxy. The CRF collects weight during that visit and the protocol specifies an examination at each visit. This supports the proxy but does not establish an actual independent weighing timestamp. The seven-day alignment rule is retained.

Source-date parsing accepts documented calendar strings, rejects numeric timestamps without an explicit origin, and preserves original source timestamps. For DeviceCGM, an adjacent positive local-time increment no larger than 30 minutes triggers an unknown separating bin when the transmitter counter resets or its positive increment differs by more than 60 seconds. No clock offset is invented. Subsequent internally consistent recorded-local-time segments are retained; absolute local clock accuracy is not claimed. A conservative sensitivity excludes any interval with a flagged boundary in its scored 28-day window. Lookback/follow-up boundaries remain visible in the saved grid arrays.

An independent audit after the initial estimates identified 66 CGM source strings with dates but no explicit time-of-day. The initial general parser had interpreted these as midnight. All 66 have neighboring observations about five minutes on both sides, supporting a formatting-omission explanation without proving it. The final primary source adapter conservatively excludes them as unresolved measurement times and requires the literal month/day/year, 12-hour clock, seconds and AM/PM format. Three duplicate source-time groups with distinct finite transmitter counters are also masked, including when their glucose values agree. The pre-audit sources, intervals, estimates and draws are preserved in a versioned history archive; the explicit midnight interpretation remains a separate sensitivity. The fitted warning policy, thresholds, clinical eligibility and statistical design remain unchanged. All dependent final estimates are recomputed.

The complete finite linear forecast grid is 40–100 mg/dL. No finite linear operating point satisfies the DCLP3 validation burden cap. Its reported no-alarm sentinel is an explicit infeasible-family fallback/reference, not a feasible finite trend-warning model. The separate never-warning reference is retained, and neither sentinel affects selection of M15 as the primary policy.

The five-minute reporting-delay sensitivity delays warning inputs and availability only. Its event grid and observed exposure remain the original outcome clock. Legality of the warning uses the delayed visible glucose, so a delayed alarm can be counted as late relative to the unchanged event onset.

Finite lookback is not assumed to erase the refractory timer automatically. Every primary-policy interval variant is checked for six fully determined consecutive no-candidate bins after complete feature-history context and before the earliest scored warning opportunity. Actual missing measurements can establish a known inability to alarm. All 2,709 interval variants pass. M120, missing-input and switch secondary replays extend lookback until the same proof holds; all their event denominators and residual associations remain unchanged. Both initial and extended ledgers are retained.

In restricted secondary samples, a nuisance column constant across all rows is algebraically redundant with the intercept and may be removed without changing the model column space. Each such column is logged. The exposure and squared exposure are never removed. The primary design never uses this option. In the contemporaneous companion, lag is exactly −28 days and its redundant column is omitted explicitly. Remaining rank deficiency or separation is reported, not repaired with a penalty.

The primary, all-event and clinical-remission estimates each use 5,000 child-bootstrap draws. Other prespecified sensitivities report child-cluster sandwich intervals and are labeled as secondary. The quadratic model's linear coefficient is not a constant per-0.20 rate ratio; retain its full coefficients and local exposure contrasts. The lowest-secretion child is determined by the mean of that child's included index AUC values; the highest-event child is determined by the sum of included total qualifying events, with deterministic ID ordering for ties. Both deletions are separate diagnostics.
