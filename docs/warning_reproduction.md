# Warning-study original-source reproduction

This pipeline reconstructs the DCLP3 development and held-out analyses, external
REPLACE-BG transport, prespecified robustness checks, and subsequent warning
diagnostics from the two original public ZIP archives. It fits all three models
from scratch. A previously fitted model, processed CGM file, validation tally or
result cache is not an input.

The reported claim concerns the frozen history-switching policy relative to its
DCLP3-selected always-short comparator, with observable event denominators and
paired warning/alarm-burden accounting. Exploratory fixed-threshold contrasts
and selection-stability resampling remain separate from that original test.

## Inputs and environment

Obtain the original public datasets under the source providers' access terms:

| Archive | SHA-256 |
|---|---|
| DCLP3 Public Dataset - Release 3 - 2022-08-04.zip | `6dbe8580fe38612622f0ebf5da58fed04b4ed556d1e06bf3a67813791222f5dc` |
| Replace-BG Dataset.zip | `6d4bb4e47e71fff72cd3eb94c5a5e60c9f32e066606bb3883821debf528dc004` |

The wrapper rejects other archive bytes before creating the work directory.
Source archives are read in place. The output directory must be separate from
both the repository directory and the directory holding the source archives.

Use Python 3.12 and install `requirements-lock.txt` in a clean environment.
The tested interpreter is Python 3.12.14; numerical versions are NumPy 2.5.4,
pandas 3.1.0, SciPy 1.18.1, scikit-learn 1.7.2 and Numba 0.68.0. The wrapper
records the executable, versions, thread settings and every subprocess command.
Eight OpenMP threads and one BLAS thread were used for the reconstruction.

## Commands

Replace the three path values with local locations. The following single command
reconstructs all warning-specific results:

```text
python -B src/reproduce_warning.py --dclp3-archive INPUTS/DCLP3.zip --replacebg-archive INPUTS/REPLACE_BG.zip --work-root RUN/card001 --diagnostics-root RUN/revision_20261022_r2/card001 --threads 8 --workers 2 --stage all
```

The archive filenames may differ; their hashes must match. To reproduce the
combined paper, run the stages below in coordination with the independent
clinical pipeline:

1. `--stage core` prepares DCLP3, fits all 63 cross-validation/final models,
   predicts validation data, freezes all 10,772 threshold comparisons, prepares
   REPLACE-BG, and replays both held-out cohorts.
2. Freeze the clinical policy from `RUN/card001` and complete its original-source
   CLVer reconstruction using the clinical entry point. No clinical outcome is
   used to choose warning models or thresholds.
3. `--stage robustness` runs all nine endpoint/clock/phase/delay variants,
   aggregate-stream and inclusive-endpoint checks, three conservative-clock
   final refits at the original hyperparameters, and the frozen artificial-mask
   experiment. It also reconstructs source flow and development opportunities.
4. `--stage diagnostics --clinical-root RUN/card003` generates the warning
   diagnostics and checks availability containment in the reconstructed CLVer
   windows. Supply `--clinical-root` only after clinical primary results exist.

All stages require the two archive arguments and the same work-root/settings.
Without `--clinical-root`, diagnostics explicitly identify their scope as DCLP3
and REPLACE-BG; they do not silently substitute old CLVer results. The separate
clinical validation choice is still reconstructed directly from the full DCLP3
validation ledger. Repeated completed stages reuse only source/version/command
and output-hash verified checkpoints. An incomplete or mismatched checkpoint
causes an error and requires a fresh work directory.

The retained `development_opportunity_audit.json` field `timing_note` describes
the **original study's** audit during fitting, before validation selection and
test replay. It is historical context, not the timing of this new reproduction.
Here that audit is regenerated during the robustness stage and still reads only
training and validation participants. The new run's execution logs are the
authoritative chronology; no retrospective prespecification claim is made.

## Scientific choices retained

The release retains the exact participant exclusions, SHA-256 allocation and
fivefold groups, random seeds, model features, inverse-participant training
weights, four hyperparameter settings, 200 boosting iterations, threshold grids,
selection ordering, refractory state, event matching, censoring rules and
bootstrap draw counts. There is no row cap, participant subsampling, parameter
retuning or replacement of real outcomes by synthetic data.

The diagnostic ablations retain thresholds M15 0.01/0.04, M120 0.04, and switch
0.04/0.01 versus 0.04/0.04. The 1,000 validation-participant resamples use seed
20261023 and condition on the newly fitted models. Input nesting establishes
that switch availability equals M15 availability; observed warning differences
arise through scores, thresholds and alarm-state interactions.

## Outputs

- `card001/models/`: all fitted models, hyperparameter reports, training
  commitment and the exact training-source snapshot.
- `card001/calibration.json`: complete 10,772-point ledger, frozen policy
  settings and development-precision assessment.
- `card001/results/`: participant scores, complete event/alarm/censored-low
  ledgers, paired bootstrap comparisons, nine-variant robustness results,
  aggregate and inclusive-endpoint results, and participant distributions.
- `card001/conservative_clock/`: modified source-clock arrays, three final
  refits and original-threshold sensitivity results.
- `RUN/revision_20261022_r2/card001/`: paired event discordance, fixed-threshold
  ablations, complete validation-person tallies and selection draws, cohort
  characteristics/coverage, `principal_warning_diagnostics.png`, and a report
  rendered from reconstructed values.
- `card001/execution_logs/`: command, environment, stage timing, source and
  generated-output hash receipts; no original-owner cache is consulted.

## Verification

The completed original-source run on 23 October 2026 regenerated the DCLP3
models and full validation grid, internal and external replays, robustness
analyses, warning diagnostics and publication figure. The historical
[scientific verification record](../provenance/scientific_verification.json)
records the successful full computation and authoritative comparisons. Source
preparation also reproduced all 168 DCLP3 participant allocations, time bounds
and canonical array fields exactly in an independent comparison.

The [current repository record](../provenance/verification.json) reports checkout
and test checks separately. Repository cleanup does not constitute a new model
fit; see [verification scope](VERIFICATION.md).

Dataset attribution: the DCLP3 Study Group and Jaeb Center, and the REPLACE-BG
Study Group, T1D Exchange and Jaeb Center supplied the public source datasets.
The analyses, content and conclusions are the responsibility of UCP Technology
LLC and were not reviewed or approved by those source organizations.
