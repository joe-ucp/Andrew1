"""Separate DCLP3 multi-source robustness adapter; never the primary stream.

Copies the already committed participant allocation and time bounds. No model,
threshold or event-dependent decision is made here. Aggregate source clocks lack
the transmitter counter needed to reproduce raw-stream clock barriers.
"""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

import numpy as np
import pandas as pd
from card001_core import causal_grid


def prepare_aggregate(archive, primary_root, output):
    archive, primary_root, output = map(Path, (archive, primary_root, output))
    original = json.loads((primary_root / 'participant_manifest.json').read_text())
    with archive.open('rb') as stream:
        archive_hash = hashlib.file_digest(stream, 'sha256').hexdigest()
    if archive_hash != original['archive_sha256']:
        raise ValueError('Aggregate source must be the exact committed DCLP3 archive')
    processed = output / 'processed'
    processed.mkdir(parents=True, exist_ok=True)
    rows = [dict(row, cohort='DCLP3_aggregate', file=str(processed / f"DCLP3_{row['id']}.npz"))
            for row in original['participants']]
    bounds = {r['id']: (r['start'], r['end']) for r in rows}
    parts = {pid: [] for pid in bounds}
    periods = {}; raw_rows = 0; invalid_times = 0; in_window_periods = {}
    with zipfile.ZipFile(archive) as z:
        with z.open('Data Files/cgm.txt') as source:
            for chunk in pd.read_csv(source, sep='|', encoding='utf-8', chunksize=250000):
                raw_rows += len(chunk)
                for key, value in chunk.Period.value_counts(dropna=False).items():
                    periods[str(key)] = periods.get(str(key), 0) + int(value)
                dt = pd.to_datetime(chunk.DataDtTm, format='%d%b%y:%H:%M:%S', errors='coerce')
                valid = dt.notna(); invalid_times += int((~valid).sum())
                chunk = chunk[valid].copy()
                chunk['time'] = dt[valid].astype('datetime64[ns]').astype('int64').to_numpy(dtype=float) / 1e9
                for pid, part in chunk.groupby('PtID', sort=False):
                    if int(pid) not in bounds:
                        continue
                    start, end = bounds[int(pid)]
                    part = part[(part.time >= start) & (part.time < end)]
                    for key, value in part.Period.value_counts(dropna=False).items():
                        in_window_periods[str(key)] = in_window_periods.get(str(key), 0) + int(value)
                    if len(part):
                        parts[int(pid)].append(part[['time', 'CGM']].to_numpy(dtype=float))
                print(f'aggregate_rows_parsed={raw_rows}', flush=True)
    if invalid_times:
        raise ValueError(f'{invalid_times} invalid aggregate timestamps; resolve before replay')
    # Both bounds and period identity must agree. Do not silently use mislabeled
    # baseline data if source dates place it inside a post-randomization window.
    if any('post' not in key.lower() or 'random' not in key.lower() for key in in_window_periods):
        raise ValueError(f'Unexpected within-window aggregate periods: {in_window_periods}')
    audit = []
    for row in rows:
        arrays = parts.pop(row['id'])
        a = np.concatenate(arrays) if arrays else np.empty((0, 2))
        if len(a):
            grid = causal_grid(a[:, 0], a[:, 1], row['start'], row['end'] - 1e-6)
            duplicates = len(a) - len(np.unique(a[:, 0]))
        else:
            t = np.arange(np.ceil(row['start'] / 300) * 300, row['end'], 300)
            grid = dict(time=t, glucose=np.full(len(t), np.nan), source_time=np.full(len(t), np.nan), conflict=np.zeros(len(t), bool))
            duplicates = 0
        np.savez_compressed(row['file'], **grid)
        audit.append(dict(id=row['id'], split=row['split'], arm=row['arm'],
                          raw_in_window=len(a), duplicate_timestamps=duplicates,
                          conflict_grid_bins=int(grid['conflict'].sum()),
                          observed_grid_bins=int(np.isfinite(grid['glucose']).sum()),
                          calendar_grid_bins=len(grid['time'])))
    manifest = dict(source_semantics_audited=True, archive=str(archive), archive_sha256=archive_hash,
                    parent_manifest_sha256=hashlib.sha256((primary_root / 'participant_manifest.json').read_bytes()).hexdigest(),
                    analysis_role='multi_source_aggregate_sensitivity', participants=rows,
                    exclusions=original['exclusions'],
                    rules=dict(clock='aggregate DataDtTm; same epoch-anchored latest-in-bin causal grid',
                               window='exact copy of primary manifest boundaries and allocations',
                               no_imputation=True, raw_primary_only=False,
                               clock_barrier='not identifiable from aggregate fields; no transmitter counter',
                               interpretation='retrospective multi-source union; does not establish application receipt availability'),
                    audit=dict(total_raw_rows=raw_rows, invalid_timestamps=invalid_times,
                               periods=periods, in_window_periods=in_window_periods,
                               participants=len(rows)))
    pd.DataFrame(audit).to_csv(output / 'source_time_audit.csv', index=False)
    (output / 'participant_manifest.json').write_text(json.dumps(manifest, indent=2), encoding='utf-8')
    print(json.dumps(manifest['audit'], indent=2))
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--archive', type=Path, required=True)
    parser.add_argument('--primary-root', type=Path, required=True)
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    prepare_aggregate(args.archive, args.primary_root, args.output)
