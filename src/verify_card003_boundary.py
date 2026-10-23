"""Independent finite-history/refractory proof audit under the supplied study root."""
from pathlib import Path
import argparse
import hashlib
import json
import numpy as np
import pandas as pd
from card003_analyze import policy_arrays
from card003_boundary import synchronized
from card001_core import ready


def inspect(path, choice, start, known):
    with np.load(path) as cache:
        p, cut, _ = policy_arrays(cache, choice)
        t, g = cache["time"], cache["glucose"]
        candidate = np.isfinite(g) & (g >= 54) & np.isfinite(p) & (p >= cut)
        end = np.searchsorted(t, start - 1800, side="right") - 1
        witnesses = [i for i in range(known, end - 4) if not candidate[i:i+6].any()]
        observed = [i for i in witnesses if np.isfinite(g[i:i+6]).all()]
        assert np.all(np.diff(t) == 300), "Irregular grid"
        assert bool(witnesses) == synchronized(t, g, p, cut, start, known)
        if choice["kind"] == "M15":
            assert np.array_equal(np.isfinite(p), ready(g, 4) & (g >= 54))
        return dict(fully_determined_six_bin_witness=bool(witnesses),
                    fully_observed_six_bin_witness=bool(observed),
                    earliest_witness_index=witnesses[0] if witnesses else None,
                    known_from=int(known))


def replay_from(candidates, next_legal):
    alarms = []
    for index in np.flatnonzero(candidates):
        if index >= next_legal:
            alarms.append(int(index))
            next_legal = int(index) + 6
    return alarms


def main(root):
    root = Path(root).resolve()
    frame = pd.read_csv(root / "processed" / "intervals_all_coverage.csv")
    frozen = json.loads((root / "models" / "frozen_policy.json").read_text())
    choice = frozen["choices"]["primary"]
    variants = ["primary", "phase1", "phase2", "phase3", "phase4", "delay5", "run3",
                "blank_units_excluded", "date_only_midnight_assumption"]
    rows = []
    for variant in variants:
        for row in frame.to_dict("records"):
            path = root / "processed" / "predictions" / variant / (str(row["pair_id"]) + ".npz")
            known = (24 if choice["kind"] in ("M120", "missing", "switch") else 3) + (variant == "delay5")
            result = inspect(path, choice, row["start_seconds"], known)
            rows.append(dict(variant=variant, pair_id=row["pair_id"], **result))
    output = pd.DataFrame(rows)
    output.to_csv(root / "results" / "independent_refractory_audit.csv", index=False)
    secondary = pd.read_csv(root / "results" / "boundary_extension" / "synchronization_audit.csv")
    assert len(secondary) == len(frame) * 3 and not secondary.duplicated(["pair_id", "policy"]).any()
    starts = frame.set_index("pair_id").start_seconds.to_dict()
    secondary_rows = []
    for row in secondary.to_dict("records"):
        choice2 = frozen["choices"]["frozen_family"][row["policy"]]
        path = root / "processed" / "predictions" / "extended_boundary" / f"{row['pair_id']}_{row['lookback_minutes']}.npz"
        result = inspect(path, choice2, starts[row["pair_id"]], 24)
        secondary_rows.append({**row, **result})
    output2 = pd.DataFrame(secondary_rows)
    output2.to_csv(root / "results" / "independent_secondary_refractory_audit.csv", index=False)
    t, g, p = np.arange(60) * 300., np.full(60, 100.), np.full(60, .5)
    p[:24] = np.nan
    prefix_rejected = not synchronized(t, g, p, .1, 50 * 300, 24)
    continuous = np.ones(60, bool)
    phase_distinct = replay_from(continuous, 0) != replay_from(continuous, 1)
    quiet = continuous.copy()
    quiet[24:30] = False
    suffixes = [[v for v in replay_from(quiet, i) if v >= 30] for i in range(6)]
    quiet_synchronizes = all(v == suffixes[0] for v in suffixes)
    summary = dict(n_windows=len(frame), n_interval_variants=len(output),
        expected_interval_variants=len(frame) * len(variants),
        all_fully_determined_witnesses=bool(output.fully_determined_six_bin_witness.all()),
        all_fully_observed_witnesses=bool(output.fully_observed_six_bin_witness.all()),
        fully_observed_missing_by_variant=output.groupby("variant").fully_observed_six_bin_witness.apply(lambda s: int((~s).sum())).to_dict(),
        secondary_policy_windows=len(output2),
        all_secondary_fully_determined_witnesses=bool(output2.fully_determined_six_bin_witness.all()),
        unknown_prefix_rejected=prefix_rejected,
        continuous_candidates_preserve_initial_phase=phase_distinct,
        six_known_quiet_bins_synchronize_all_possible_initial_timers=quiet_synchronizes,
        errors=[], source_sha256=hashlib.sha256(Path(__file__).with_name("card003_boundary.py").read_bytes()).hexdigest(),
        verification_source_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        claim="Independent computational audit of the primary policy's nine variants and extended secondary policies.")
    summary["pass"] = bool(summary["all_fully_determined_witnesses"] and summary["all_secondary_fully_determined_witnesses"] and prefix_rejected and phase_distinct and quiet_synchronizes)
    (root / "results" / "independent_refractory_audit.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    assert summary["pass"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True)
    main(parser.parse_args().data_root)
