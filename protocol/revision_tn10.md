# TN10 post-results diagnostic specification

Revision date: 22 October 2026. This is a post-results presentation and provenance revision. It is not a new preregistration. The frozen scientific configuration, numerical outputs and decision rules remain unchanged.

## Diagnostic scope

- Estimate TN10 as a separate scientific component with no pooled causal synthesis with glucose-warning studies. The authoritative combined manuscript retains this separation.
- Give the fixed-group Welch interval, score-reconstruction bootstrap interval and adjusted continuous interaction comparable prominence.
- Tabulate existing continuous, complete-case, normalization and horizon results without choosing a replacement primary analysis.
- Plot the four released treatment-by-score outcome curves through 156 weeks with explicit risk-set conventions. Recompute their areas and reconcile those areas to the existing restricted means.
- Plot the 26-week margin only in the dichotomized interaction panel, where it applies. Do not impose it on subgroup treatment effects or continuous-score interactions.
- Audit original endpoint provenance. Distinguish verification of encoded durations from verification of individual clinical follow-up.
- Bound the finite-sample interaction by allowing each of the 19 non-event endpoints without verified last-contact to vary in [0,156] weeks, holding all 54 recorded diagnoses fixed. This is a post-results provenance sensitivity, not a confidence interval, plausible censoring model, or new efficacy hypothesis.

## Scientific object retained

The 359-gene AbATE list projects to 312 exact TN10 symbols. Primary scoring remains pooled baseline log2(1+CPM), per-gene sample-SD standardization, an unweighted mean and a strict-above-median high group. The study population remains 73 baseline RNA libraries, interpreted as 73 participants under the original publication's description. The outcome remains the released duration restricted at 156 weeks. The effect is the low-minus-high difference between randomized-arm mean differences; the numerical target remains 26 weeks. Neither score, endpoint horizon nor hypothesis is selected anew.

## Newly explicit provenance limitation

The original checks establish zero early censoring *encoded* in the release. All 19 non-events nevertheless have the same Survival value, 551.4286 weeks, equal to the largest event duration rounded to four decimals. No individual last-contact field is provided. The original authors' public TN10 code consumes Survival from an R object without deriving censor times. Thus clinical completeness to 156 weeks cannot be independently verified for those 19 records. New curves, risk sets and restricted means are labeled conditional on the released durations. The prior `all_73_primary_outcomes_observed` assertion is not independent confirmation of clinical ascertainment and is narrowed by this audit; the original wording is retained in the historical record outside this release.

The conservative finite-sample range is −14.1427489177 to 136.7754329004 weeks. Its lower endpoint assigns zero to the unverifiable low-treated/high-placebo records and 156 to the low-placebo/high-treated records; its upper endpoint reverses those assignments. This linear contrast bound does not include sampling uncertainty. It deliberately does not use potentially informative later RNA sampling to tighten the endpoint box.

## Storage and execution

New numerical outputs and figures are written under `RUN/revision_20261022_r2/card002`. The portable `src/reproduce_tn10.py` reconstructs the original primary arithmetic from the three provider files, then runs `src/revision_tn10.py` and independent numerical checks. It records source identities and checks input invariance after execution. The full source path recomputes the frozen analysis; the separate retained-intermediate path does not refit it. Curves have no newly selected tests or confidence bands. Continuous and bootstrap analyses retain the original specifications.

The three-year stage-3-free endpoint is bounded and horizon-dependent. Untreated prognosis and randomized treatment-benefit modification remain separate. The estimated contrast concerns stage-3-free time through 156 weeks under the released endpoint encoding; it does not identify biological resistance or a treatment-selection rule.
