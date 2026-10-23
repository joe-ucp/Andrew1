"""Card003 entry point. Reports missing inputs honestly and never substitutes data.

Default mode rebuilds the completed report, or executes the complete pipeline
when original CLVer and independently frozen DCLP3 models are available.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import io
import json
import platform
import os
import re
import subprocess
import time
import zipfile
from datetime import datetime, timezone
from pathlib import Path
from card003_core import sha256_file, EXPECTED_CLVER_SHA256

REPO = Path(__file__).resolve().parents[1]
MEMBERS = {
    "DeviceCGM.txt": "8759338f647b063d57f4e2ca6a3f14d1fa3c087a46fb3495c301d8aee04ca39c",
    "DevicePumpCGMDATAGXB.txt": "9c9ae8a2cb1c313f7a16cc381fc8674c2501a132ad9c849f468272a3b7846b6c",
    "PtRoster.txt": "bd809e20652b5de7754e57559c516fe2b4db26b1c1b687038742c97934f1d372",
    "SampleResults.txt": "8cd1f5fdebe0910d7109439da2ed38014ebf8b0b4c2fbc3bc1d15a6baae3fddb",
    "InsulinUse.txt": "a29dc971c2cc4227900f8f5e4d0232393518d39df7f088d5bbb34a8e1d36ad16",
    "DiabPhysExam.txt": "ecdd29cc83408692162ef1bbcab1dbf16cd19c9f24695575832c2578a4d32cb2",
    "VisitInfo.txt": "5764abe5510308f91d18165fb5c114248ed35cf30b92d4c5d4e6eae0ceedd25c",
}


def dump(path, obj):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False, allow_nan=False), encoding="utf-8")


def archive_audit(path, data_root):
    """Hash originals and inspect headers; no eligibility inferred from headers."""
    import hashlib
    result = {"archive": str(path), "bytes": path.stat().st_size,
              "sha256": sha256_file(path), "expected_sha256": EXPECTED_CLVER_SHA256,
              "members": [], "audited_utc": datetime.now(timezone.utc).isoformat()}
    result["source_version_matches"] = result["sha256"] == EXPECTED_CLVER_SHA256
    with zipfile.ZipFile(path) as z:
        result["archive_inventory"] = [{"name": e.filename, "bytes": e.file_size} for e in z.infolist()]
        for basename, expected in MEMBERS.items():
            found = [n for n in z.namelist() if n.split("/")[-1] == basename]
            if len(found) != 1:
                result["members"].append({"name": basename, "status": "missing_or_ambiguous"})
                continue
            name = found[0]
            h = hashlib.sha256()
            with z.open(name) as f:
                for b in iter(lambda: f.read(1048576), b""):
                    h.update(b)
            with z.open(name) as f:
                bom = f.read(4)
            encoding = "utf-16" if bom.startswith((b"\xff\xfe",b"\xfe\xff")) else "utf-8-sig"
            with z.open(name) as f:
                header = io.TextIOWrapper(f,encoding=encoding).readline().strip()
            delimiters = {d: header.count(d) for d in ("\t", "|", ",")}
            delim = max(delimiters, key=delimiters.get)
            result["members"].append({"name": name, "bytes": z.getinfo(name).file_size,
                "sha256": h.hexdigest(), "expected_sha256": expected,
                "hash_matches": h.hexdigest() == expected,
                "encoding": encoding, "delimiter": repr(delim),
                "columns": next(csv.reader([header], delimiter=delim))})
        for name in z.namelist():
            if "read" in name.lower() and name.lower().endswith(".txt"):
                (data_root / "raw" / Path(name).name).write_bytes(z.read(name))
    dump(data_root / "results" / "archive_audit.json", result)
    return result


def verify_clock_closure(root):
    from verify_card003_clocks import verify
    evidence=verify(root)
    if not evidence.get("strict_primary_date_mask_pass") or not evidence.get("computational_checks_pass_under_midnight_interpretation"):
        raise ValueError("Independent source clock audit did not pass")


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--data-root", type=Path, required=True)
    ap.add_argument("--archive", type=Path)
    ap.add_argument("--audit-only", action="store_true")
    ap.add_argument("--rebuild",action="store_true",help="Reconstruct the source intervals and rerun all empirical analyses")
    ap.add_argument("--freeze-policy-from", type=Path,
                    help="Validate and freeze Card003 policy from a completed Card001 DCLP3 full validation grid")
    args = ap.parse_args()
    t0 = time.perf_counter()
    root = args.data_root.resolve()
    for name in ("raw", "processed", "results", "literature", "history"):
        (root / name).mkdir(parents=True, exist_ok=True)
    env = os.environ.copy()
    env["PYTHONDONTWRITEBYTECODE"] = "1"
    checks = subprocess.run([sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_card003*.py", "-v"],
                            cwd=REPO, capture_output=True, text=True, env=env)
    test_log = checks.stdout + checks.stderr
    (root / "results" / "software_tests.log").write_text(test_log, encoding="utf-8")
    count_match = re.search(r"Ran (\d+) tests", test_log)
    archive = args.archive or root / "raw" / "CLVerPublicDataset.zip"
    status = {
        "card": "003", "complete": False, "hypothesis_status": "not_tested",
        "statistical_result": None, "clinical_result": None,
        "data_root": str(root), "repo": str(REPO),
        "expected_source_sha256": EXPECTED_CLVER_SHA256,
        "python": platform.python_version(), "executed_utc": datetime.now(timezone.utc).isoformat(),
        "software_tests": {"passed": checks.returncode == 0,
                           "count": int(count_match[1]) if count_match else None,
                           "log": str(root / "results" / "software_tests.log")},
        "implemented": ["normalized MMTT AUC and assay gates", "clinical IDAA1c and date/unit gates",
            "prospective interval and overlap gates", "independent episode cross-check",
            "boundary-aware frozen-policy replay kernel", "prespecified PPML and child bootstrap",
            "qualified-H1 conjunction", "source archive/member hash and schema inventory"],
        "remaining": ["obtain authorized CLVer archive", "validate raw-to-canonical source adapter against glossary",
            "inspect fasting/status and date encoding", "audit chronological and source/clock exclusions",
            "obtain hash-locked DCLP3 policy selected by the complete Card003 procedure",
            "construct all eligible real clinical/CGM intervals and event/alarm ledgers",
            "evaluate transfer and residual-event gates", "estimate primary and companion associations",
            "execute all prespecified damage/sensitivity tests", "complete independent reproduction"],
        "no_substitution": "No fabricated or simulated patients, clinical estimates, policy-transfer results or efficacy claims.",
    }
    if args.freeze_policy_from:
        from card003_policy import freeze_from_card001
        frozen = freeze_from_card001(args.freeze_policy_from, root / "models")
        status["frozen_policy"] = str(frozen)
        status["frozen_policy_sha256"] = sha256_file(frozen)
        status["remaining"].remove("obtain hash-locked DCLP3 policy selected by the complete Card003 procedure")
    if archive.exists():
        audit = archive_audit(archive, root)
        status["archive_audit"] = str(root / "results" / "archive_audit.json")
        status["source_version_matches"] = audit["source_version_matches"]
        status["remaining"].remove("obtain authorized CLVer archive")
        status["status"] = "source_audited_analysis_not_executed" if audit["source_version_matches"] else "source_version_requires_audit"
    else:
        status["status"] = "awaiting_authorized_data"
        status["missing_input"] = str(archive)
    summary_path=root/"results"/"analysis_summary.json"
    if summary_path.exists() and json.loads(summary_path.read_text()).get("complete") and not args.rebuild:
        from card003_boundary import extend_families,finalize as finalize_boundaries
        if not (root/"results"/"boundary_extension"/"corrected_policy_results.json").exists():extend_families(root)
        finalize_boundaries(root)
        from card003_diagnostics import run as run_diagnostics
        run_diagnostics(root)
        subprocess.run([sys.executable,"-B",str(REPO/"src"/"verify_card003.py"),"--data-root",str(root),"--expected-result",str(root/"results"/"model_primary.json"),"--bootstrap-draws",str(root/"results"/"bootstrap_primary.csv.gz")],check=True,cwd=REPO,env=env)
        verify_clock_closure(root)
        from card003_report import build_report
        complete=build_report(root,REPO)
        print(json.dumps({"status":complete["status"],"complete":True,"hypothesis_status":complete["hypothesis_status"]}))
        return
    if archive.exists() and (root/"models"/"frozen_policy.json").exists() and not args.audit_only and not args.freeze_policy_from:
        if checks.returncode:raise RuntimeError("Scientific kernel tests failed")
        if args.rebuild or not (root/"processed"/"intervals_all_coverage.csv").exists():
            subprocess.run([sys.executable,"-B",str(REPO/"src"/"card003_data.py"),"--archive",str(archive),"--output",str(root/"processed"),"--prepare"],check=True,cwd=REPO,env=env)
        subprocess.run([sys.executable,"-B",str(REPO/"src"/"card003_analyze.py"),"--data-root",str(root),"--stage","complete"],check=True,cwd=REPO,env=env)
        from card003_boundary import extend_families,finalize as finalize_boundaries
        extend_families(root)
        finalize_boundaries(root)
        from card003_diagnostics import run as run_diagnostics
        run_diagnostics(root)
        subprocess.run([sys.executable,"-B",str(REPO/"src"/"verify_card003.py"),"--data-root",str(root),"--expected-result",str(root/"results"/"model_primary.json"),"--bootstrap-draws",str(root/"results"/"bootstrap_primary.csv.gz")],check=True,cwd=REPO,env=env)
        verify_clock_closure(root)
        from card003_report import build_report
        complete=build_report(root,REPO)
        print(json.dumps({"status":complete["status"],"complete":True,"hypothesis_status":complete["hypothesis_status"]}))
        return
    status["elapsed_seconds"] = round(time.perf_counter() - t0, 3)
    dump(root / "results" / "card003_status.json", status)
    report = f"""# Card 003 execution status

Status: **{status['status']}**. H1 is **not tested**. The clinical research run is not complete.

The source card supplies a research design and earlier feasibility counts; those counts are not newly reproduced results. No rate ratio, remission association, warning sensitivity, or clinical efficacy estimate is available from this execution.

Research data/results: `{root}`. Source code: `{REPO / 'src'}`.

The local assay, clinical-remission, temporal-window, independent episode, replay-boundary, PPML and uncertainty kernels are implemented. Software-fixture tests run: {status['software_tests']['count']}; passed: {status['software_tests']['passed']}. These fixtures exercise the implemented kernels; raw-source adapter inspection remains outstanding at this preparation stage.

## Outstanding work

""" + "\n".join("- " + x for x in status["remaining"]) + f"""

## Interpretation

The requested endpoint is the later observed rate of significant low-CGM episodes without a qualifying replay warning, conditional on a useful policy frozen in DCLP3. It cannot be estimated by adding or multiplying a TN10 disease-delay estimate and a CLVer hypoglycemia rate. IDAA1c at or below 9 and laboratory secretion remain separate measurements. The specified outcomes are observed sensor-low burden, measured secretion and the visit-level IDAA1c remission definition.

The source-data access requirement is recorded at [Jaeb CLVer](https://public.jaeb.org/dataset/591): a registrant supplies a name, email, institution/planned use and agrees to the terms before download. A documented publicly visible archive URL does not replace completing the site's required access step. Existing authorized local data can be provided at the missing-input path.

Machine-readable execution status and the source audit, if available, are stored under `{root / 'results'}`. Rerun: `python src/card003_run.py --data-root "{root}"`.
"""
    (root / "reports").mkdir(parents=True, exist_ok=True)
    (root / "reports" / "card003_execution.md").write_text(report, encoding="utf-8")
    print(json.dumps({"status": status["status"], "complete": False,
                      "hypothesis_status": "not_tested", "status_file": str(root / "results" / "card003_status.json")}))


if __name__ == "__main__":
    main()
