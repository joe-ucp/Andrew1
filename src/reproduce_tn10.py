"""Rebuild the complete TN10 analysis from three checksum-verified originals.

The input directory is read-only. This offline runner uses only this checkout,
the active Python environment, and newly generated artifacts beneath the chosen
output root. It never downloads inputs or reads an author's previous analysis.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys
import time

from card002_run import SOURCES


def sha(path):
    with Path(path).open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def run(inputs_root, output_root):
    started = time.perf_counter()
    inputs_root, output_root = Path(inputs_root).resolve(), Path(output_root).resolve()
    repo = Path(__file__).resolve().parents[1]
    if output_root == repo or repo in output_root.parents or output_root in repo.parents:
        raise ValueError("Output and package trees must be separate")
    if (inputs_root == output_root or inputs_root in output_root.parents
            or output_root in inputs_root.parents):
        raise ValueError("Input and output trees must be separate")
    card = output_root / "card002"
    revision = output_root / "revision_20261022_r2" / "card002"
    diagnostics = output_root / "revision_20261022_r3" / "paper" / "tn10"
    for destination in (card, revision, diagnostics):
        if destination.exists() and any(destination.iterdir()):
            raise FileExistsError(f"A fresh TN10 output tree is required: {destination}")
    inputs = {}
    for name, (url, expected, size) in SOURCES.items():
        path = inputs_root / name
        if not path.is_file() or path.stat().st_size != size or sha(path) != expected:
            raise ValueError(f"Original input identity mismatch: {name}")
        inputs[name] = {"sha256": expected, "bytes": size, "source_url": url}

    # Dependencies are supplied by the active environment, never a data-tree
    # sys.path insertion. Verify the exact declared TN10 dependency versions.
    lock = repo / "requirements-lock.txt"
    pinned = dict(line.strip().split("==", 1) for line in lock.read_text().splitlines()
                  if line.strip() and not line.lstrip().startswith("#"))
    required = ["numpy", "pandas", "scipy", "statsmodels", "matplotlib", "pdfplumber", "pypdf"]
    versions = {name: importlib.metadata.version(name) for name in required}
    for name, version in versions.items():
        if version != pinned.get(name):
            raise RuntimeError(f"Dependency mismatch for {name}: installed {version}, required {pinned.get(name)}")

    logs = card / "reproduction_logs"
    logs.mkdir(parents=True, exist_ok=False)
    cache = card / "cache"
    (cache / "tmp").mkdir(parents=True)
    environment = os.environ.copy()
    environment.pop("PYTHONPATH", None)
    environment.pop("PYTHONHOME", None)
    environment.update(PYTHONDONTWRITEBYTECODE="1", PYTHONNOUSERSITE="1", PYTHONIOENCODING="utf-8",
                       MPLCONFIGDIR=str(cache / "matplotlib"), TMP=str(cache / "tmp"),
                       TEMP=str(cache / "tmp"))
    raw_args = ["--raw-root", str(inputs_root)]
    stages = [
        ("primary", "card002_run.py", ["--data-root", str(card), *raw_args]),
        ("independent_primary", "verify_card002.py", ["--data-root", str(card), *raw_args]),
        ("curves_and_ascertainment", "revision_tn10.py", ["--original-root", str(card),
          "--output-root", str(revision), *raw_args]),
        ("independent_curves_and_hc3", "verify_revision_tn10_clinical.py", ["--original-root", str(card),
          "--revision-root", str(revision), "--output-root", str(card / "verification"), *raw_args]),
        ("full_bootstrap_diagnostics", "paper_tn10_diagnostics.py", ["--data-root", str(output_root),
          "--original-root", str(card), "--revision2-root", str(revision),
          "--output-root", str(diagnostics), *raw_args]),
    ]
    source_names = ["reproduce_tn10.py"] + [script for _, script, _ in stages]
    receipt = {
        "status": "RUNNING", "started_utc": datetime.now(timezone.utc).isoformat(),
        "scope": "Offline, full raw-input TN10 reconstruction; no prior author results used as analysis input",
        "python": sys.version, "executable": sys.executable, "python_prefix": sys.prefix,
        "python_base_prefix": sys.base_prefix, "dependencies": versions, "lock_sha256": sha(lock),
        "source_sha256": {name: sha(repo / "src" / name) for name in source_names},
        "original_inputs": inputs, "raw_input_directory": str(inputs_root),
        "output_root": str(output_root), "bootstrap_repetitions": 10000, "bootstrap_seed": 20261022,
        "commands": [], "scientific_inference": "Conditional on released endpoints; individual last contact is unavailable for 19 non-events",
    }
    manifest = card / "reproduction.json"
    dump(manifest, receipt)
    for label, script, arguments in stages:
        command = [sys.executable, "-B", str(repo / "src" / script), *arguments]
        stage_start = time.perf_counter()
        print(f"START {label}", flush=True)
        log = logs / f"{label}.log"
        with log.open("w", encoding="utf-8") as stream:
            process = subprocess.run(command, cwd=repo, env=environment, stdout=stream,
                                     stderr=subprocess.STDOUT, text=True, check=False)
        stage = {"stage": label, "argv": command, "cwd": str(repo), "return_code": process.returncode,
                 "seconds": time.perf_counter() - stage_start, "log": str(log.relative_to(output_root)),
                 "log_sha256": sha(log)}
        receipt["commands"].append(stage)
        dump(manifest, receipt)
        if process.returncode:
            receipt.update(status="FAILED", failed_stage=label, elapsed_seconds=time.perf_counter() - started)
            dump(manifest, receipt)
            raise RuntimeError(f"TN10 stage {label} failed; see {log}")
        print(f"PASS {label}: {stage['seconds']:.3f} seconds", flush=True)
    for name, source in inputs.items():
        assert sha(inputs_root / name) == source["sha256"], f"Raw input changed: {name}"
    for name, expected in receipt["source_sha256"].items():
        assert sha(repo / "src" / name) == expected, f"Source changed during run: {name}"
    primary = json.loads((card / "results" / "primary_results.json").read_text())
    boot = json.loads((card / "results" / "bootstrap_results.json").read_text())
    assert boot["repetitions"] == boot["valid"] == 10000 and boot["unusable"] == 0
    receipt.update(status="PASS", completed_utc=datetime.now(timezone.utc).isoformat(),
                   elapsed_seconds=time.perf_counter() - started, raw_inputs_unchanged=True,
                   source_unchanged=True, primary_hypothesis_status=primary["hypothesis_status"],
                   theta_weeks=primary["theta"]["estimate_weeks"],
                   theta_fixed_group_ci95=primary["theta"]["ci95_weeks"],
                   theta_full_score_bootstrap_ci95=boot["percentile_ci95_weeks"]["theta"])
    # Manifests inside individual stages bind their detailed outputs. This
    # receipt binds every generated scientific table/figure, excluding caches.
    receipt["output_sha256"] = {
        str(path.relative_to(output_root)): sha(path)
        for folder in (card / "results", card / "verification", revision, diagnostics)
        for path in sorted(folder.rglob("*")) if path.is_file()
        and "cache" not in path.relative_to(folder).parts
        and "matplotlib_cache" not in path.relative_to(folder).parts
    }
    dump(manifest, receipt)
    print(json.dumps({"status": "PASS", "manifest": str(manifest),
                      "seconds": receipt["elapsed_seconds"], "theta_weeks": receipt["theta_weeks"]}), flush=True)
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs-root", type=Path, required=True, help="Read-only directory of original source files")
    parser.add_argument("--output-root", type=Path, required=True, help="Fresh shared study output root")
    args = parser.parse_args()
    run(args.inputs_root, args.output_root)
