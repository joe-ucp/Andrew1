"""Reproduce CLVer from its original archive and newly developed DCLP3 models.

No retained participant records or fitted clinical results are required. All
runtime artifacts are written below --run-root; the supplied archive is read-only.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import time


EXPECTED_ARCHIVE_SHA256 = "8421a9a9ff7dd1485ce0e4bdb984b5ea037300a0e7216f9f445334e56d32d385"
REPO = Path(__file__).resolve().parents[1]


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def write_json(path, value):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def validate_layout(run_root, archive=None, repo=REPO):
    """Keep generated study records outside source code and original-input folders."""
    run_root, repo = Path(run_root).resolve(), Path(repo).resolve()
    if run_root == repo or run_root.is_relative_to(repo) or repo.is_relative_to(run_root):
        raise ValueError("The run root and source repository must be separate directories")
    if archive is not None:
        archive = Path(archive).resolve()
        target = run_root / "card003" / "raw" / "CLVerPublicDataset.zip"
        if archive != target and (archive.parent == run_root
                or archive.parent.is_relative_to(run_root)
                or run_root.is_relative_to(archive.parent)):
            raise ValueError("Keep the original-input and generated-output directories separate")
    return run_root


def environment(run_root):
    env = dict(os.environ)
    env.pop("PYTHONHOME", None)
    env.update(PYTHONDONTWRITEBYTECODE="1", PYTHONUTF8="1", PYTHONNOUSERSITE="1",
               OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
               OMP_NUM_THREADS="4", NUMBA_NUM_THREADS="8")
    for key, name in (("MPLCONFIGDIR", "matplotlib"), ("NUMBA_CACHE_DIR", "numba"),
                      ("TEMP", "tmp"), ("TMP", "tmp")):
        folder = run_root / "runtime_cache" / name
        folder.mkdir(parents=True, exist_ok=True)
        env[key] = str(folder)
    # Dependencies come from the selected pinned environment, never inherited
    # author caches, user-site packages or another checkout.
    env["PYTHONPATH"] = str(REPO / "src")
    return env


def validate_stage_state(run_root, stage):
    """A preparation or fit must not replace a prior scientific commitment."""
    study = Path(run_root) / "card003"
    if stage in ("prepare", "complete") and study.exists() and any(study.iterdir()):
        raise FileExistsError("Use a fresh clinical study directory for source preparation")
    if stage == "analyze":
        for committed in (study / "models" / "frozen_policy.json",
                          study / "results" / "analysis_summary.json",
                          study / "clinical_execution.json"):
            if committed.exists():
                raise FileExistsError("Existing clinical commitment is protected: " + str(committed))


def execute(run_root, label, arguments, env):
    folder = run_root / "card003" / "logs"
    folder.mkdir(parents=True, exist_ok=True)
    logfile = folder / (label + ".log")
    command = [sys.executable, "-B"] + [str(x) for x in arguments]
    start = time.perf_counter()
    print("Starting " + label, flush=True)
    with logfile.open("w", encoding="utf-8") as stream:
        result = subprocess.run(command, cwd=REPO, env=env, stdout=stream,
                                stderr=subprocess.STDOUT)
    record = dict(command=command, elapsed_seconds=time.perf_counter()-start,
                  returncode=result.returncode, log=str(logfile), log_sha256=sha(logfile))
    if result.returncode:
        raise RuntimeError(f"{label} failed; inspect {logfile}")
    print(f"Finished {label} in {record['elapsed_seconds']:.1f} seconds", flush=True)
    return record


def prepare(run_root, archive, env):
    study = run_root / "card003"
    for name in ("raw", "processed", "results", "history", "logs"):
        (study / name).mkdir(parents=True, exist_ok=True)
    target = study / "raw" / "CLVerPublicDataset.zip"
    source = Path(archive).resolve() if archive is not None else target
    if not source.is_file() or sha(source) != EXPECTED_ARCHIVE_SHA256:
        raise ValueError("Supply the checksum-matching original CLVer public archive")
    if target.exists():
        if sha(target) != EXPECTED_ARCHIVE_SHA256:
            raise ValueError("The existing run archive has a different source identity")
    elif source != target:
        shutil.copyfile(source, target)
    else:
        raise FileNotFoundError(target)
    source_before = sha(source)
    records = [execute(run_root, "source_prepare", ["src/card003_data.py",
               "--archive", target, "--output", study / "processed", "--prepare"], env),
               execute(run_root, "source_clinical_verification", ["src/verify_card003.py",
               "--data-root", study, "--clinical-only"], env),
               execute(run_root, "source_clock_verification", ["src/verify_card003_clocks.py",
               "--data-root", study], env)]
    if sha(source) != source_before:
        raise ValueError("The original source archive changed")
    outputs = [p for p in (study / "processed").rglob("*") if p.is_file()
               and "predictions" not in p.relative_to(study).parts]
    source_names = ["card003_data.py", "card003_core.py", "card001_core.py"]
    record = dict(status="PREPARED_FROM_ORIGINAL_ARCHIVE", archive_sha256=source_before,
                  archive_source=str(source), original_archive_unchanged=True,
                  code_sha256={n:sha(REPO / "src" / n) for n in source_names},
                  processed_sha256={p.relative_to(study).as_posix():sha(p) for p in outputs},
                  commands=records,
                  interval_audit=json.loads((study / "results" / "interval_audit.json").read_text()))
    write_json(study / "source_reproduction.json", record)
    return record


def verify_prepared(study):
    path = study / "source_reproduction.json"
    if not path.exists():
        raise ValueError("Run --stage prepare before analysis")
    record = json.loads(path.read_text())
    if sha(study / "raw" / "CLVerPublicDataset.zip") != EXPECTED_ARCHIVE_SHA256:
        raise ValueError("Original archive identity changed")
    for name, digest in record["processed_sha256"].items():
        if sha(study / name) != digest:
            raise ValueError("Prepared source changed: " + name)
    for name, digest in record["code_sha256"].items():
        if sha(REPO / "src" / name) != digest:
            raise ValueError("Source construction code changed after preparation: " + name)
    return record


def analyze(run_root, warning_root, env):
    study = run_root / "card003"
    verify_prepared(study)
    warning_root = Path(warning_root).resolve()
    if not (warning_root / "calibration.json").is_file():
        raise FileNotFoundError("Complete the newly trained DCLP3 validation grid first: " + str(warning_root))
    records = [execute(run_root, "independent_policy_freeze", ["src/card003_run.py",
               "--data-root", study, "--freeze-policy-from", warning_root], env),
               execute(run_root, "full_clinical_analysis", ["src/card003_run.py",
               "--data-root", study], env),
               execute(run_root, "independent_refractory_verification", ["src/verify_card003_boundary.py",
               "--data-root", study], env)]
    status = json.loads((study / "results" / "card003_status.json").read_text())
    if not status.get("complete"):
        raise ValueError("The clinical analysis has remaining work")
    write_json(study / "clinical_execution.json", dict(status="COMPLETE", commands=records,
               warning_root=str(warning_root), frozen_policy_sha256=sha(study / "models" / "frozen_policy.json"),
               scientific_status=status["hypothesis_status"]))
    return records


def revision(run_root, env):
    study = run_root / "card003"
    destination = run_root / "revision_20261022_r2" / "card003"
    records = [execute(run_root, "clinical_revision", ["src/revision_clinical.py",
               "--source-root", study, "--output-root", destination, "--repo", REPO], env),
               execute(run_root, "clinical_revision_independent", ["src/verify_revision_clinical.py",
               "--data-root", run_root], env),
               execute(run_root, "clinical_publication_figure", ["src/revision_clinical_figure.py",
               "--source-root", study, "--revision-root", destination,
               "--output", run_root / "figures" / "secretion_rates.png"], env)]
    write_json(destination / "source_reproduction_execution.json", dict(status="COMPLETE", commands=records))
    return records


def verify(run_root, env):
    study = run_root / "card003"
    verify_prepared(study)
    records = [execute(run_root, "independent_clinical_final", ["src/verify_card003.py",
               "--data-root", study, "--expected-result", study / "results" / "model_primary.json",
               "--bootstrap-draws", study / "results" / "bootstrap_primary.csv.gz"], env),
               execute(run_root, "independent_clock_final", ["src/verify_card003_clocks.py",
               "--data-root", study], env),
               execute(run_root, "independent_refractory_final", ["src/verify_card003_boundary.py",
               "--data-root", study], env),
               execute(run_root, "independent_revision_final", ["src/verify_revision_clinical.py",
               "--data-root", run_root], env)]
    return records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--run-root", type=Path, required=True)
    parser.add_argument("--archive", type=Path, help="Original CLVer ZIP; required for a new prepare/complete run")
    parser.add_argument("--warning-root", type=Path, help="New DCLP3 development output; default RUN_ROOT/card001")
    parser.add_argument("--stage", choices=["prepare", "analyze", "revision", "verify", "complete"], default="complete")
    args = parser.parse_args()
    root = validate_layout(args.run_root, args.archive)
    validate_stage_state(root, args.stage)
    root.mkdir(parents=True, exist_ok=True)
    env = environment(root)
    warning = args.warning_root or root / "card001"
    start = time.perf_counter()
    if args.stage in ("prepare", "complete"):
        prepare(root, args.archive, env)
    if args.stage in ("analyze", "complete"):
        analyze(root, warning, env)
    if args.stage in ("revision", "complete"):
        revision(root, env)
    if args.stage in ("verify", "complete"):
        verify(root, env)
    result = dict(status="COMPLETE" if args.stage=="complete" else "STAGE_COMPLETE", stage=args.stage,
                  completed_utc=datetime.now(timezone.utc).isoformat(), elapsed_seconds=time.perf_counter()-start,
                  run_root=str(root), source_repository=str(REPO), python=sys.version,
                  source_to_result_inputs="Original CLVer ZIP and new DCLP3 training/calibration outputs; no retained clinical records or fits")
    write_json(root / "card003" / ("reproduction_" + args.stage + ".json"), result)
    print(json.dumps(result, indent=2), flush=True)


if __name__ == "__main__":
    main()
