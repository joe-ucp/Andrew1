# TN10: complete reconstruction from original inputs

The offline entry point reconstructs the TN10 analysis from three original files. It does not read a previous analysis, fitted model, derived patient table, or author cache. The six-input directory used by the complete analysis may contain additional dataset archives; this entry point reads only the three files below.

## Inputs and environment

Obtain the original files from their source, retain the exact basenames, and place them in a directory separate from both this checkout and the new output directory. Source checksums and sizes are checked before output creation.

| Original input | Bytes | SHA-256 | Source |
|---|---:|---|---|
| TN10_metadata.csv | 44,999 | `dabd7beff35f82aef231f9c9e0fd5b91974cdbc9f175edd8831c475824517586` | [Zenodo record 17063348](https://zenodo.org/records/17063348/files/TN10_metadata.csv?download=1) |
| TN10_raw_counts.csv | 8,508,942 | `74085ba79ef79c44d5e25679e36f3615562c0a644700f2332426d71b942e1ee8` | [Zenodo record 17063348](https://zenodo.org/records/17063348/files/TN10_raw_counts.csv?download=1) |
| JCI176403_Table6.pdf | 210,945 | `93b00f6ccb5688eea28f936cfc8e72117e65f70b612c656531cea11c5b51f90c` | [Original supplemental table](https://www.jci.org/articles/view/176403/sd/pdf/render/6) |

Use Python 3.12 and a fresh virtual environment with `requirements-lock.txt`. The validated environment was Python 3.12.14 on 64-bit Windows. The runner checks the installed versions of NumPy 2.5.4, pandas 3.1.0, SciPy 1.18.1, statsmodels 0.15.0, Matplotlib 3.11.2, pdfplumber 0.11.9, and pypdf 6.10.0 against the lock. The lock also pins their dependencies and the other study dependencies.

From the repository root, using that environment's interpreter:

```text
python -m pip install -r requirements-lock.txt
python -B src/reproduce_tn10.py --inputs-root ORIGINAL_INPUT_DIRECTORY --output-root NEW_STUDY_OUTPUT_DIRECTORY
```

Replace the two directory placeholders with actual paths. Outputs may share a root with other studies, but the three TN10 destination subtrees must be empty. Input, output, and repository trees must be disjoint. The runner refuses overlapping paths, an occupied TN10 output subtree, a missing or changed original, or a mismatched direct dependency. It performs no downloads. Child processes use the active interpreter, remove inherited `PYTHONPATH` and `PYTHONHOME`, disable user-site imports and bytecode writes, and keep temporary files and Matplotlib caches within the new output tree.

## Computations and outputs

The five sequential stages are:

| Stage | Source | New outputs beneath the chosen output root |
|---|---|---|
| Primary reconstruction and all frozen sensitivities | `src/card002_run.py` | `card002/results/`, `card002/frozen_config.json`, and source/environment records |
| Separate primary arithmetic verification | `src/verify_card002.py` | `card002/results/independent_verification.json` |
| Outcome curves, continuous-model table, and follow-up identification bounds | `src/revision_tn10.py` | `revision_20261022_r2/card002/` |
| Separate date, cell, curve, risk-count, bound, and HC3 verification | `src/verify_revision_tn10_clinical.py` | `card002/verification/tn10_independent_clinical_review.json` |
| Full bootstrap assignment replay and minimum-follow-up profile | `src/paper_tn10_diagnostics.py` | `revision_20261022_r3/paper/tn10/` |

`card002/reproduction.json` records the actual commands, source hashes, input hashes and sizes, installed versions, interpreter, command logs, elapsed times, output hashes, and completion status. Each stage's combined output is retained in `card002/reproduction_logs/`. Detailed stage manifests retain their own provenance. An interrupted or failed run keeps its failure record and must be restarted in another fresh output tree.

The primary calculation re-extracts the 359 negatively associated genes from the original PDF and uses the 312 exact symbols present in the count matrix. CPM normalization uses all 13,021 supplied gene rows. The analysis retains the 73 baseline participants, including 42 assigned teplizumab and 31 placebo. It computes pooled sample-standardized log2(1 + CPM), the unweighted gene score, and the strict-above-pooled-median group assignment. The outcome is the released endpoint truncated at 156 weeks.

All 10,000 stratified bootstrap draws use seed 20261022. Placebo is sampled before treatment within each draw. Gene means, sample standard deviations, scores, and the median split are recomputed in every draw. The later diagnostic stage replays those same draws, verifies all 30,000 effect values, and records individual group-assignment changes. The normalizations, horizons, continuous models, complete-case comparison, leave-one-out checks, design calculation, and follow-up profiles retain the original scientific definitions.

The first verifier uses a separate PDF parser, pandas score calculation, and an OLS HC2 contrast. The second uses separate source-date arithmetic and manual OLS with HC3 covariance, and independently checks all curve steps and risk counts. These are numerical implementation checks.

## Verification

The completed original-source run on 23 October 2026 executed all five TN10 stages, including the full 10,000-draw score reconstruction and assignment replay. The independent verifiers checked the primary contrast, date arithmetic, outcome curves, risk counts, follow-up bounds and continuous-model uncertainty. The complete research run passed 82 tests, all 54 aggregate evidence comparisons, all four publication-figure pixel comparisons and seven exact intermediate-reconstruction outputs. Its recorded stage time was 78.05 minutes for all three investigations, not a TN10-only runtime.

The [scientific verification record](../provenance/scientific_verification.json) preserves that completed run. The [current repository verification record](../provenance/verification.json) records checkout and test checks separately. Neither a repository hash check nor a unit test is presented as a new TN10 refit. See [verification scope](VERIFICATION.md).

## Interpretation retained by the reproduction

The fixed-group interaction estimate is 69.1754329004 weeks, with Welch 95% interval 20.4923156479–117.8585501529. The full score-reconstruction bootstrap interval is 1.4158912764–104.8659037837 weeks. The adjusted continuous interaction is −16.0165647879 weeks per score standard deviation, with interval −46.6021729670–14.5690433913. These quantify different model specifications and must not be substituted for one another.

Nineteen coded non-events lack an independently documented individual last-contact date in the released inputs. Their common encoded duration is not proof of complete clinical follow-up. The reconstructed endpoints and outcome curves remain conditional on that release encoding. The identification range −14.1427489177 to 136.7754329004 weeks varies those 19 truncated endpoints between 0 and 156 weeks while holding the observed event durations fixed; it is not a confidence interval. The program explicitly records the follow-up limitation and retains the original hypothesis status. Portability changes do not alter the estimands, uncertainty procedures, or scientific conclusions.
