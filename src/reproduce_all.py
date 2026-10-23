"""Execute the three frozen investigations from six original provider files."""
import argparse
from datetime import datetime, timezone
import json
from pathlib import Path
import subprocess
import sys
import time

sys.path.insert(0, str(Path(__file__).resolve().parent))
from release import PACKAGE, environment, check_inputs


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--inputs-root', type=Path, required=True)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    inputs, out = args.inputs_root.resolve(), args.output_root.resolve()
    for protected in (inputs, PACKAGE):
        if out == protected or out.is_relative_to(protected) or protected.is_relative_to(out):
            raise ValueError('Output, input and package trees must be disjoint')
    if out.exists() and any(p.name not in {'environment.json', 'runtime'} for p in out.iterdir()):
        raise FileExistsError('Use a new complete-run directory')
    check_inputs(inputs)
    env = environment(out)
    logs = out/'logs'; logs.mkdir(exist_ok=False)
    commands = []
    def run(name, script, *arguments):
        command = [sys.executable, '-B', str(PACKAGE/'src'/script), *map(str, arguments)]
        started = time.monotonic()
        print('START '+name, flush=True)
        with (logs/(name+'.log')).open('w', encoding='utf-8') as stream:
            result = subprocess.run(command, cwd=PACKAGE, env=env, stdout=stream, stderr=subprocess.STDOUT)
        commands.append(dict(name=name, command=command, returncode=result.returncode,
                             seconds=time.monotonic()-started, log='logs/'+name+'.log'))
        (out/'commands.json').write_text(json.dumps(commands, indent=2)+'\n', encoding='utf-8')
        if result.returncode:
            raise RuntimeError(name+' failed; inspect '+str(logs/(name+'.log')))
        print('PASS '+name, flush=True)
    run('scientific_tests', 'run_tests.py', '--output-root', out/'test_results')
    run('tn10', 'reproduce_tn10.py', '--inputs-root', inputs, '--output-root', out)
    warning_args = ['--dclp3-archive', inputs/'DCLP3 Public Dataset - Release 3 - 2022-08-04.zip',
                    '--replacebg-archive', inputs/'Replace-BG Dataset.zip', '--work-root', out/'card001',
                    '--diagnostics-root', out/'revision_20261022_r2/card001', '--threads', '8', '--workers', '2']
    run('warning_primary', 'reproduce_warning.py', *warning_args, '--stage', 'core')
    run('warning_robustness', 'reproduce_warning.py', *warning_args, '--stage', 'robustness')
    run('clinical', 'reproduce_clinical.py', '--run-root', out,
        '--archive', inputs/'CLVerPublicDataset.zip', '--warning-root', out/'card001', '--stage', 'complete')
    run('warning_diagnostics', 'reproduce_warning.py', *warning_args,
        '--clinical-root', out/'card003', '--stage', 'diagnostics')
    run('endpoint_decomposition', 'paper_endpoint_decomposition.py', '--data-root', out,
        '--output-root', out/'revision_20261022_r3/paper/endpoints')
    run('compare_authoritative', 'compare_results.py', '--run-root', out, '--output', out/'result_comparison.json')
    check_inputs(inputs)
    receipt = dict(status='PASS', completed_utc=datetime.now(timezone.utc).isoformat(),
                   scope='Original-source reproduction; six provider inputs and release source only',
                   original_inputs_unchanged=True, commands=commands)
    (out/'reproduction_result.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
    print(json.dumps(receipt, indent=2))


if __name__ == '__main__':
    main()
