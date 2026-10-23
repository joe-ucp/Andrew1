# Verification scope

Two records describe different checks on the same research project:

| Record | Scope | What it establishes |
|---|---|---|
| [Current repository verification](../provenance/verification.json) | Checks performed on the current checkout, including the maintained file manifest, tests and documentation | The condition of the repository after maintenance, with actual check outcomes recorded separately from scientific execution. |
| [Completed scientific verification](../provenance/scientific_verification.json) | The original-source reproduction completed on 23 October 2026 at 06:11:03 UTC | Recalculation of the three investigations from six provider files in the declared numerical environment, followed by authoritative result and figure comparisons. |

The completed scientific run passed **82 tests**, **54 aggregate evidence comparisons**, and **four exact figure-pixel comparisons**. A separate reconstruction from its newly generated intermediates matched **seven committed output hashes**. The recorded original-source stage time was **78.05 minutes**. These figures describe that completed run; they are not substituted for the current checkout's discovered test count or presented as a new analysis during repository cleanup.

The paper, frozen scientific choices and findings remain the reference. Documentation and repository maintenance do not retune models, change event definitions, replace input datasets or initiate a new research question. Changes to numerical analysis logic or scientific settings require a new, explicitly identified execution and comparison record.

## Project-root verification results

On 23 October 2026, the project root was checked as an ordinary Git repository with local `_UCP/` and `_private/` directories present:

- `git add .` included all **148 content files**, including all **54 scientific evidence files** and the **26 files under nested `results/` directories**. No required file was ignored.
- Content verification passed with root Git metadata and excluded private directories present: **147 recorded file hashes** checked. The Git index contained only the reviewed public files; no private history or remote was configured.
- **104 tests passed; two were skipped** because this Windows account cannot create symbolic links. Actual Windows directory-junction rejection tests passed. The suite also rejected staged private files, deleted private files in history, and private history on another branch. No test failed or errored.
- Comparisons against the retained completed run passed for **54 evidence files and four exact figure pixel arrays**.
- Reconstruction from those retained intermediates passed **seven exact output hashes**; input files remained unchanged.
- All **108 relative links and fragments** across 15 Markdown files resolved after the README editorial update.

The 54 scientific evidence files and eight manuscript/figure files are byte-identical to the completed scientific run's reference files. This cleanup reran implementation tests, result verification and intermediate reconstruction; it did not refit the full participant-level analyses. The exact outcomes and executable command forms are in the current verification record above.

## Checks you can run

Run these commands from the repository root using the pinned environment described in [REPRODUCING.md](REPRODUCING.md):

```text
python -I -B src/release.py verify
python -I -B src/release.py tests --output-root /path/to/new-test-output
```

The first checks the content allowlist and recorded file hashes. Root Git metadata and optional local `_UCP/` and `_private/` directories are outside the content manifest; the Git index and all locally reachable commit trees are checked for private or non-allowlisted paths. The second runs the included tests and records the actual counts, results and imported module locations. Neither command refits participant-level analyses.

For the full scientific path, acquire the [six original provider files](DATASETS.md), then run:

```text
python -I -B src/release.py reproduce --inputs-root /path/to/t1d-inputs --output-root /path/to/new-t1d-run
```

Only a successful `reproduction_result.json`, supported by stage logs and `result_comparison.json`, establishes completion of that new run. To reconstruct the specified calculations from a completed run's intermediates, use the separate command in [REPRODUCING.md](REPRODUCING.md#reconstruct-from-retained-intermediates).

## Interpretation

The numerical comparisons retain exact integer, text, schema and figure-pixel checks and declared floating-point tolerances. Explicitly excluded execution metadata are reported by field; estimates, event counts, confidence intervals and decisions remain compared. A reported mismatch must be investigated without changing scientific commitments to make the check pass.

Reproduction establishes agreement with the computational findings. It preserves the inferior frozen switching policy, the observational nature of the CLVer association and the uncertainty in TN10 treatment-effect modification. It does not turn those results into demonstrated patient benefit or validated treatment selection.
