"""Reconstruct the warning study from the two original public source archives.

No result cache is an input. Logs, derived data and checkpoints stay below
--work-root. Run with the release's pinned environment.
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

SOURCE = Path(__file__).resolve().parent
EXPECTED_SOURCE_SHA256 = {
    'DCLP3': '6dbe8580fe38612622f0ebf5da58fed04b4ed556d1e06bf3a67813791222f5dc',
    'REPLACE_BG': '6d4bb4e47e71fff72cd3eb94c5a5e60c9f32e066606bb3883821debf528dc004',
}


def disjoint_output(output, protected):
    """Reject an output inside, equal to, or containing protected input trees."""
    output = Path(output).resolve()
    for path in protected:
        path = Path(path).resolve()
        if output == path or output.is_relative_to(path) or path.is_relative_to(output):
            raise ValueError(f'Output directory must be disjoint from source/input tree: {path}')
    return output


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--dclp3-archive', type=Path, required=True)
    parser.add_argument('--replacebg-archive', type=Path, required=True)
    parser.add_argument('--work-root', type=Path, required=True,
                        help='New warning-study directory, usually SOURCE_RUN/card001')
    parser.add_argument('--diagnostics-root', type=Path)
    parser.add_argument('--clinical-root', type=Path,
                        help='Optional independently reconstructed Card003 root for the CLVer containment check')
    parser.add_argument('--threads', type=int, default=8)
    parser.add_argument('--workers', type=int, default=2)
    parser.add_argument('--stage', choices=['prepare', 'train', 'core', 'primary',
                        'robustness', 'diagnostics', 'all'], default='all')
    args = parser.parse_args()
    if not 1 <= args.threads <= 64 or not 1 <= args.workers <= 8:
        parser.error('threads must be 1..64 and workers 1..8')
    dclp = args.dclp3_archive.resolve()
    replace = args.replacebg_archive.resolve()
    if not dclp.is_file() or not replace.is_file():
        parser.error('Both original source ZIP archives are required')
    protected = [SOURCE.parent, dclp.parent, replace.parent]
    root = disjoint_output(args.work_root, protected)
    if root.exists() and any(root.iterdir()) and not (root / 'original_inputs.json').is_file():
        raise ValueError('Nonempty work directory is not owned by this reconstruction; use a new directory')
    if args.diagnostics_root:
        diagnostics_path = disjoint_output(args.diagnostics_root, protected)
        if diagnostics_path.exists() and any(diagnostics_path.iterdir()) and not (root / 'execution_logs' / 'revision_diagnostics.json').is_file():
            raise ValueError('Nonempty diagnostics directory lacks a matching execution receipt')
    input_record = {'DCLP3': {'file': str(dclp), 'sha256': digest(dclp)},
                    'REPLACE_BG': {'file': str(replace), 'sha256': digest(replace)}}
    for name, item in input_record.items():
        if item['sha256'] != EXPECTED_SOURCE_SHA256[name]:
            raise ValueError(f'{name} archive differs from the exact published source release')
    root.mkdir(parents=True, exist_ok=True)
    logs = root / 'execution_logs'
    logs.mkdir(exist_ok=True)
    temp = root / 'temporary'
    temp.mkdir(exist_ok=True)
    commitment = root / 'original_inputs.json'
    if commitment.exists() and json.loads(commitment.read_text()) != input_record:
        raise ValueError('Original inputs changed; use a new work directory')
    if not commitment.exists():
        commitment.write_text(json.dumps(input_record, indent=2), encoding='utf-8')
    versions = {name: importlib.metadata.version(name) for name in
                ['numpy', 'pandas', 'scipy', 'scikit-learn', 'numba', 'matplotlib']}
    env = dict(os.environ, PYTHONPATH=str(SOURCE), PYTHONDONTWRITEBYTECODE='1',
               PYTHONNOUSERSITE='1',
               OMP_NUM_THREADS=str(args.threads), OPENBLAS_NUM_THREADS='1',
               MKL_NUM_THREADS='1', NUMBA_NUM_THREADS=str(args.threads),
               MPLCONFIGDIR=str(root / 'matplotlib_cache'),
               TMP=str(temp), TEMP=str(temp))
    env.pop('PYTHONHOME', None)
    environment_record = {
        'python': sys.version, 'executable': sys.executable, 'packages': versions,
        'threads': args.threads, 'workers': args.workers,
        'source_directory': str(SOURCE), 'original_inputs': input_record,
        'only_release_source_on_pythonpath': True}
    environment_path = logs / 'environment.json'
    if environment_path.exists() and json.loads(environment_path.read_text()) != environment_record:
        raise ValueError('Execution environment changed; use a new work directory')
    if not environment_path.exists():
        environment_path.write_text(json.dumps(environment_record, indent=2), encoding='utf-8')

    def stage_sources(script):
        names = {'card001_core.py', 'card001_models.py', 'reproduce_warning.py', script}
        if script == 'card001_dataflow.py':
            names.add('card001_data.py')
        if script == 'card001_development_audit.py':
            names.add('card001_run.py')
        if script in ['card001_run.py', 'revision_warning.py', 'card001_sensitivity.py',
                      'card001_conservative_clock.py']:
            names |= {'card001_run.py', 'card001_evaluate.py', 'card001_data.py'}
        if script == 'card001_run.py':
            names.add('card001_sensitivity.py')
        if script == 'revision_warning.py':
            names.add('card003_policy.py')
        return {name: digest(SOURCE / name) for name in sorted(names)}

    def file_state():
        return {str(p.resolve()): (p.stat().st_size, p.stat().st_mtime_ns)
                for folder in {root, args.diagnostics_root.resolve() if args.diagnostics_root else root / 'diagnostics'}
                if folder.exists() for p in folder.rglob('*') if p.is_file()
                and not any(x in p.parts for x in ['execution_logs', 'temporary', 'matplotlib_cache', '__pycache__'])}

    def run(name, script, *arguments, threads=None):
        command = [sys.executable, '-B', str(SOURCE / script), *map(str, arguments)]
        receipt = logs / f'{name}.json'
        sources = stage_sources(script)
        if receipt.exists():
            old = json.loads(receipt.read_text())
            if old.get('returncode') != 0 or old.get('source_sha256') != sources or old.get('command') != command:
                raise ValueError(f'{name} checkpoint does not match this execution; use a new work directory')
            for path, expected in old.get('output_sha256', {}).items():
                if not Path(path).is_file() or digest(path) != expected:
                    raise ValueError(f'{name} output changed: {path}')
            if not old.get('output_sha256'):
                raise ValueError(f'{name} checkpoint lacks verified outputs')
            print(f'REUSE VERIFIED {name}', flush=True)
            return
        if (logs / f'{name}.log').exists():
            raise ValueError(f'{name} has an incomplete execution log; use a new work directory')
        before = file_state()
        record = {'stage': name, 'command': command, 'started_utc':
                  datetime.now(timezone.utc).isoformat(),
                  'script_sha256': digest(SOURCE / script), 'source_sha256': sources}
        active_env = env if threads is None else dict(env, OMP_NUM_THREADS=str(threads),
                                                     NUMBA_NUM_THREADS=str(threads))
        started = time.perf_counter()
        print(f'START {name}', flush=True)
        with (logs / f'{name}.log').open('w', encoding='utf-8') as log:
            completed = subprocess.run(command, cwd=SOURCE.parent, env=active_env,
                                       stdout=log, stderr=subprocess.STDOUT, check=False)
        record.update(returncode=completed.returncode, seconds=time.perf_counter()-started,
                      finished_utc=datetime.now(timezone.utc).isoformat(),
                      log_sha256=digest(logs / f'{name}.log'))
        after = file_state()
        record['output_sha256'] = {path: digest(path) for path, state in after.items()
                                   if before.get(path) != state}
        receipt.write_text(json.dumps(record, indent=2), encoding='utf-8')
        print(f'END {name}: exit={completed.returncode}, seconds={record["seconds"]:.1f}', flush=True)
        if completed.returncode:
            raise RuntimeError(f'{name} failed; see {logs / (name + ".log")}')

    def prepared(path, archive):
        manifest = path / 'participant_manifest.json'
        if not manifest.exists():
            return False
        obj = json.loads(manifest.read_text())
        if obj['archive_sha256'] != digest(archive):
            raise ValueError('Prepared data came from a different source archive')
        if not all(Path(p['file']).is_file() for p in obj['participants']):
            raise ValueError('Incomplete source preparation; use a new work directory')
        stage_name, script = ('prepare_dclp3', 'card001_data.py') if path == root else (
            ('prepare_replacebg', 'card001_replacebg.py') if path == root / 'external' / 'replacebg'
            else ('prepare_aggregate', 'card001_aggregate.py'))
        receipt = logs / f'{stage_name}.json'
        if not receipt.exists():
            raise ValueError('Prepared source has no execution receipt; use a new work directory')
        prior = json.loads(receipt.read_text())
        if prior.get('returncode') != 0 or prior.get('source_sha256') != stage_sources(script) or not prior.get('output_sha256'):
            raise ValueError('Prepared source has no matching code/output commitment')
        for file, expected in prior['output_sha256'].items():
            if not Path(file).is_file() or digest(file) != expected:
                raise ValueError(f'Prepared source output changed: {file}')
        return True

    external = root / 'external' / 'replacebg'
    if args.stage in ['prepare', 'train', 'core', 'all']:
        if not prepared(root, dclp):
            run('prepare_dclp3', 'card001_data.py', '--archive', dclp, '--output', root)
    if args.stage in ['train', 'core', 'all']:
        run('train', 'card001_run.py', '--data-root', root, '--stage', 'train')
        run('calibrate', 'card001_run.py', '--data-root', root, '--stage', 'calibrate')
    if args.stage in ['prepare', 'core', 'primary', 'all']:
        if not prepared(external, replace):
            run('prepare_replacebg', 'card001_replacebg.py', '--archive', replace,
                '--output', external)
    if args.stage in ['core', 'primary', 'all']:
        run('internal_test', 'card001_run.py', '--data-root', root, '--stage', 'test')
        run('external_test', 'card001_run.py', '--data-root', root, '--stage', 'test',
            '--external-root', external)
    if args.stage in ['robustness', 'all']:
        run('internal_sensitivity', 'card001_run.py', '--data-root', root,
            '--stage', 'sensitivity', '--workers', args.workers,
            threads=max(1, args.threads // args.workers))
        run('external_sensitivity', 'card001_run.py', '--data-root', root,
            '--stage', 'sensitivity', '--external-root', external, '--workers', args.workers,
            threads=max(1, args.threads // args.workers))
        if not prepared(root / 'aggregate', dclp):
            run('prepare_aggregate', 'card001_aggregate.py', '--archive', dclp,
                '--primary-root', root, '--output', root / 'aggregate')
        run('aggregate_test', 'card001_run.py', '--data-root', root, '--stage', 'test',
            '--external-root', root / 'aggregate', '--cohort-split', 'test',
            '--result-prefix', 'aggregate')
        run('inclusive_endpoint', 'card001_run.py', '--data-root', root, '--stage', 'test',
            '--external-root', external / 'inclusive_day182',
            '--result-prefix', 'external_inclusive_day182')
        run('conservative_clock_fit', 'card001_conservative_clock.py', '--data-root', root,
            '--stage', 'fit')
        run('conservative_clock_internal', 'card001_conservative_clock.py', '--data-root', root,
            '--stage', 'evaluate')
        run('conservative_clock_external', 'card001_conservative_clock.py', '--data-root', root,
            '--stage', 'evaluate', '--external-root', external)
        run('mask_stress', 'card001_stress.py', '--data-root', root, '--threads', args.threads)
        run('source_dataflow', 'card001_dataflow.py', '--data-root', root)
        run('source_clock_audit', 'card001_clock_audit.py', '--data-root', root)
        run('development_opportunities', 'card001_development_audit.py', '--data-root', root)
    if args.stage in ['diagnostics', 'all']:
        run('result_accounting', 'card001_result_audit.py', '--data-root', root)
        diagnostics = args.diagnostics_root.resolve() if args.diagnostics_root else root / 'diagnostics'
        more = ['--clinical-root', args.clinical_root.resolve()] if args.clinical_root else []
        run('revision_diagnostics', 'revision_warning.py', '--warning-root', root,
            '--output-root', diagnostics, *more, '--stage', 'all')
    print(json.dumps({'status': 'REQUESTED_STAGES_COMPLETE', 'stage': args.stage,
                      'work_root': str(root), 'logs': str(logs)}), flush=True)


if __name__ == '__main__':
    main()
