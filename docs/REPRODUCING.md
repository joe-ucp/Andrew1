# Reproduction instructions

Run commands from the repository root. The three levels below answer different questions. A manifest check is not a refit, and reconstruction from fitted-model summaries is not reproduction from original datasets.

## Environment

The numerical reference environment is **Python 3.12.14**, with every package pinned in [requirements-lock.txt](../requirements-lock.txt). The validation record gives the tested operating system and package versions. Windows is the exercised platform; paths are configurable and no drive letter is required. Other platforms have not been separately tested.

PowerShell:

```powershell
python -m venv ../t1d-env
& ../t1d-env/Scripts/python.exe -m pip install -r requirements-lock.txt
& ../t1d-env/Scripts/python.exe -m pip check
```

POSIX shell equivalent (not separately platform-tested):

```sh
python3.12 -m venv ../t1d-env
../t1d-env/bin/python -m pip install -r requirements-lock.txt
../t1d-env/bin/python -m pip check
```

In subsequent examples, `python` means that environment's interpreter. Use its full path rather than relying on another Python installation. `-I` ignores inherited Python paths and user-site packages. Subprocesses receive only this repository's scientific source directory on `PYTHONPATH`; temporary files and numerical caches are directed into the new output tree. No installed libraries or scripts are taken from an author's data cache.

The full warning-model grid is computationally substantial: its original 63 fits took about 31 minutes with eight threads, before calibration, sensitivity replay, clinical bootstraps and diagnostics. Reserve several GB of memory and ample free disk space for grids, features, predictions and resamples. Actual validation runtimes are in the verification record; no reduced bootstrap or development grid is offered as an equivalent reproduction.

## Verify repository files

```text
python -I -B src/release.py verify
python -I -B src/release.py tests --output-root /path/to/new-test-output
```

The first command uses Python's standard library and, when root Git metadata is present, the Git command-line tool. It checks the explicit content allowlist, recorded SHA-256 values, unexpected content files and forbidden symlink members. It prunes root `.git/` and the optional local `_UCP/` and `_private/` directories before scanning public content, then checks the Git index and all reachable commit trees for private or non-allowlisted paths. Those two private directories must remain ignored and untracked. A checkout must have complete local history for this audit. The `.gitignore` permits exactly the reviewed public files, so ordinary `git add .` includes the public repository content. `PACKAGE_SHA256.json` covers the listed public content except itself, avoiding a circular hash.

The [current verification record](../provenance/verification.json) reports the checks performed on this repository. The [scientific verification record](../provenance/scientific_verification.json) preserves the completed original-source run and intermediate reconstruction. That historical run finished on 23 October 2026 at 06:11:03 UTC: 82 tests, 54 aggregate comparisons, four exact figure matches and seven exact reconstruction outputs passed. Its recorded stage time was 78.05 minutes. Checkout cleanup and integrity checks do not imply a new participant-level analysis; [verification scope](VERIFICATION.md) explains which checks establish which result.

The second command uses the pinned environment and runs all included scientific and release-safety tests. It records discovered test counts and imported module locations; it does not enforce a historical count or include private workflow tests. Preserve its `tests.json` and `tests.log`.

## Reproduce from original datasets

First [acquire the six original files](DATASETS.md). Input, output and repository trees must be separate. `--output-root` must not already exist.

```text
python -I -B src/release.py inputs --inputs-root /path/to/t1d-inputs
python -I -B src/release.py reproduce --inputs-root /path/to/t1d-inputs --output-root /path/to/new-t1d-run
```

This is the complete scientific path. It runs the public tests; rebuilds TN10; prepares and trains DCLP3; calibrates the frozen choices; replays internal and external warning cohorts and prespecified sensitivities; reconstructs CLVer and its bootstrap and diagnostic analyses; runs independent source and numerical checks; regenerates endpoint accounting and four scientific figures; and compares against all retained aggregate evidence.

Generated `commands.json` records actual subprocess arguments and return codes. Stage logs are under `logs/` and study-specific log directories. `reproduction_result.json` is written as PASS only after all stages, source-input invariance checks and result comparisons pass. `result_comparison.json` lists each scientific comparison, any explicitly excluded execution metadata and any discrepancy.

Scientific choices are preserved: participant separation and deterministic splits; the full committed training/calibration grids; fixed score thresholds; event/readiness/refractory and source-clock rules; original model specifications; primary 10,000/5,000-draw resampling counts and seeds; and separate exploratory diagnostic seeds. Details and narrower commands are in [warning reproduction](warning_reproduction.md), [clinical reproduction](clinical_reproduction.md), and [TN10 reproduction](tn10_reproduction.md).

### Inspect output without recomputing

```text
python -I -B src/compare_results.py --run-root /path/to/completed-t1d-run --output /path/to/new-comparison.json
```

Exact byte comparisons are attempted first. For structured numerical outputs, the comparator checks fields and ordered rows with `atol=1e-10`, `rtol=1e-10`, sufficient for binary floating-point/serialization roundoff without changing displayed scientific conclusions. Integral values and strings are exact. Explicit exclusions cover execution times (including the exact nested path `never_warning_identity.json/original_total_result_reused/elapsed_seconds`), location-dependent cohort-manifest hashes, and one corrected TN10 provenance sentence (`original_assertion_narrowed`); their old and new values are reported. All ascertainment counts and scientific estimates remain compared. Participant-specific parameters omitted from the public aggregate projection are not silently exported. Figures require exact pixel equality in the declared environment. Any failure stops the run and remains visible in the report.

## Reconstruct from retained intermediates

After the full original-source run, or when a reviewer already holds its equivalent authorized intermediate files:

```text
python -I -B src/reconstruct_paper.py --data-root /path/to/completed-t1d-run --raw-root /path/to/t1d-inputs --output-root /path/to/new-diagnostics
```

This reruns ledger accounting, contrasts from saved quadratic coefficients, TN10 bootstrap label reconstruction and the follow-up profile. It uses the 25 dependencies listed in [external_dependencies.json](../provenance/external_dependencies.json), two of which are original TN10 CSVs resolved through `--raw-root`. All other inputs are generated by the full source pipeline. It requires no private workflow material or owner cache.

The manifest retains historical original-run checksums for traceability. Rebuilt provenance fields can change file bytes; such differences are reported per input. Every input is hashed before and after execution. Seven committed aggregate outputs must match [reconstruction_expected.json](../provenance/reconstruction_expected.json) byte for byte. This path does **not** establish that upstream models were refitted: use the original-source path for that claim.

## Rebuild manuscript formats

The authoritative [paper.md](../manuscript/paper.md) and four figures are bundled with designated PDF, Word and TeX representations. Format rebuilding is separate from statistical reproduction. It additionally requires Pandoc 3.9, XeLaTeX and the fonts in [document_build_environment.json](../provenance/document_build_environment.json). These external tools are not bundled or installed by the numerical quick start.

```text
python -I -B src/build_paper_formats.py --source-root . --output-root /path/to/new-paper-build --work-root /path/to/new-paper-work --pandoc /path/to/pandoc --xelatex /path/to/xelatex
```

Omit `--xelatex` to build TeX and Word only. Both destination directories must be new and separate from the repository. The source manuscript and designated PDF are never overwritten. Format/build checks retained from the authoritative manuscript are labeled historical; the current repository's checks are in [verification.json](../provenance/verification.json).

## Failure handling and interpretation

Do not retune thresholds, alter seeds, drop observations or broaden tolerances after a failed comparison. Save the log and compare source identity, environment, preprocessing and calculation separately. Use a new run directory after repairing a failed stage, or the documented guarded per-study checkpoints where explicitly supported. Existing frozen output is never a fallback analysis input for the original-source path.

A computational match reproduces the reported evidence. It does not resolve missing individual TN10 last-contact dates, observational confounding in CLVer, clinical warning effects that replay cannot observe, or the failed warning-burden constraints. These are limitations of the completed research questions, not reasons to replace their results.
