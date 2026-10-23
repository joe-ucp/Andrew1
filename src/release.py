"""Portable release verification and scientific reproduction entry point."""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parents[1]
PRIVATE_ROOT_DIRECTORIES = {'_UCP', '_private'}
PRIVATE_PATH_COMPONENTS = {name.casefold() for name in PRIVATE_ROOT_DIRECTORIES}


def sha(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def load(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def contained(root, relative):
    root = Path(root).resolve()
    value = Path(relative)
    if value.is_absolute() or value.drive or '..' in value.parts:
        raise ValueError(f'Unsafe manifest path: {relative}')
    target = (root / value).resolve()
    if not target.is_relative_to(root):
        raise ValueError(f'Manifest path escapes root: {relative}')
    return target


def separate(output, *inputs):
    output = Path(output).resolve()
    for other in (PACKAGE, *inputs):
        other = Path(other).resolve()
        if output == other or output.is_relative_to(other) or other.is_relative_to(output):
            raise ValueError('Output, package and input trees must be disjoint')
    if output.exists():
        raise FileExistsError(f'Use a new output directory: {output}')
    return output


def environment(output):
    env = dict(os.environ)
    env.pop('PYTHONHOME', None)
    env.update(PYTHONPATH=str(PACKAGE / 'src'), PYTHONNOUSERSITE='1',
               PYTHONDONTWRITEBYTECODE='1', PYTHONUTF8='1', PYTHONIOENCODING='utf-8',
               OPENBLAS_NUM_THREADS='1', MKL_NUM_THREADS='1', OMP_NUM_THREADS='8',
               NUMBA_NUM_THREADS='8')
    for key, name in [('TMP', 'tmp'), ('TEMP', 'tmp'), ('MPLCONFIGDIR', 'matplotlib'),
                      ('NUMBA_CACHE_DIR', 'numba')]:
        folder = Path(output) / 'runtime' / name
        folder.mkdir(parents=True, exist_ok=True)
        env[key] = str(folder)
    return env


def verify_git_private_paths(private_directories, allowed):
    """Audit the index and all locally reachable commit trees without changing Git."""
    env = {key: value for key, value in os.environ.items() if not key.startswith('GIT_')}
    env.update(GIT_OPTIONAL_LOCKS='0', GIT_NO_REPLACE_OBJECTS='1',
               GIT_NO_LAZY_FETCH='1', GIT_CONFIG_NOSYSTEM='1', GIT_CONFIG_GLOBAL=os.devnull)

    def git(*arguments, allowed_returncodes=(0,)):
        command = ['git', '-c', 'core.fsmonitor=false', '-C', str(PACKAGE), *arguments]
        try:
            result = subprocess.run(command, env=env, capture_output=True, check=False)
        except OSError as error:
            raise ValueError('Git is required to verify checkout privacy') from error
        if result.returncode not in allowed_returncodes:
            detail = result.stderr.decode('utf-8', errors='replace').strip()
            raise ValueError('Git privacy audit failed: ' + detail)
        return result

    top = git('rev-parse', '--show-toplevel').stdout.decode('utf-8').strip()
    if Path(top).resolve() != PACKAGE.resolve():
        raise ValueError('Root Git metadata refers to a different working tree')
    if git('rev-parse', '--is-shallow-repository').stdout.strip() != b'false':
        raise ValueError('Complete local Git history is required for the privacy audit')
    for name in private_directories:
        if git('check-ignore', '--no-index', '-q', '--', name + '/',
               allowed_returncodes=(0, 1)).returncode:
            raise ValueError('Root private directory must be ignored by Git: ' + name)

    def inspect_paths(raw, location):
        paths = [item.decode('utf-8', errors='surrogateescape')
                 for item in raw.split(b'\0') if item]
        forbidden = [path for path in paths
                     if any(part.casefold() in PRIVATE_PATH_COMPONENTS for part in path.split('/'))]
        if forbidden:
            raise ValueError(f'Private paths in Git {location}: {forbidden}')
        unexpected = sorted(set(paths) - allowed)
        if unexpected:
            raise ValueError(f'Non-allowlisted paths in Git {location}: {unexpected}')
        return len(paths)

    indexed = inspect_paths(git('ls-files', '-z').stdout, 'index')
    commits = git('rev-list', '--all').stdout.decode('ascii').splitlines()
    history_paths = 0
    for commit in commits:
        history_paths += inspect_paths(
            git('ls-tree', '-r', '-z', '--name-only', commit).stdout, 'history ' + commit)
    return dict(status='PASS', index_paths_checked=indexed, commits_checked=len(commits),
                commit_tree_paths_checked=history_paths,
                scope='Index and every commit tree reachable from all local refs; replacement objects disabled')


def verify():
    # Inspect members before opening manifests. Root Git metadata and the two
    # private directories are opaque; neither they nor filesystem links are walked.
    actual = set()
    git_metadata = False
    private_directories = []
    directories = [PACKAGE]
    while directories:
        directory = directories.pop()
        for path in directory.iterdir():
            relative = path.relative_to(PACKAGE)
            if path.is_symlink() or path.is_junction():
                raise ValueError('Forbidden release member: ' + str(path))
            if relative.as_posix() == '.git':
                if path.is_dir():
                    git_metadata = True
                    continue
                if path.is_file():
                    lines = path.read_text(encoding='utf-8').splitlines()
                    if len(lines) == 1 and lines[0].startswith('gitdir: ') and lines[0][8:].strip():
                        git_metadata = True
                        continue
                raise ValueError('Invalid root Git metadata: ' + str(path))
            if len(relative.parts) == 1 and relative.name in PRIVATE_ROOT_DIRECTORIES:
                if not path.is_dir():
                    raise ValueError('Root private exclusion must be a directory: ' + str(path))
                private_directories.append(relative.name)
                continue
            if any(part.casefold() in {'.git', '__pycache__'} | PRIVATE_PATH_COMPONENTS
                   for part in relative.parts):
                raise ValueError('Forbidden release member: ' + str(path))
            if path.is_dir():
                directories.append(path)
            elif path.is_file():
                actual.add(relative.as_posix())
            else:
                raise ValueError('Unsupported release member: ' + str(path))
    expected = load(PACKAGE / 'PACKAGE_SHA256.json')
    allowed = set(load(PACKAGE / 'RELEASE_ALLOWLIST.json')['files'])
    if allowed != actual:
        raise ValueError(f'Allowlist mismatch: missing={sorted(allowed-actual)}, extra={sorted(actual-allowed)}')
    if set(expected) != allowed - {'PACKAGE_SHA256.json'}:
        raise ValueError('Hash manifest does not cover the exact release allowlist')
    for relative, digest in expected.items():
        if sha(contained(PACKAGE, relative)) != digest:
            raise ValueError('Changed package file: ' + relative)
    git_privacy = (verify_git_private_paths(private_directories, allowed) if git_metadata
                   else dict(status='NOT_APPLICABLE', scope='No root Git metadata exists'))
    return dict(status='PASS', scope='Public-file integrity and local Git privacy; no model or data reproduction',
                files=len(actual), hashes_checked=len(expected), exact_allowlist=True,
                root_git_metadata_excluded=True,
                root_private_directories_excluded=sorted(private_directories),
                git_privacy=git_privacy)


def check_inputs(root):
    root = Path(root).resolve()
    records = load(PACKAGE / 'provenance/source_data.json')
    for record in records:
        path = contained(root, record['filename'])
        if not path.is_file() or path.stat().st_size != record['bytes'] or sha(path) != record['sha256']:
            raise ValueError('Original input missing or changed: ' + record['filename'])
    return dict(status='PASS', original_inputs_verified=len(records), source_bytes_unchanged=True)


def check_environment():
    required = {}
    for line in (PACKAGE / 'requirements-lock.txt').read_text().splitlines():
        if line.strip() and not line.lstrip().startswith('#'):
            name, version = line.strip().split('==', 1)
            required[name] = version
    actual = {name: importlib.metadata.version(name) for name in required}
    if actual != required:
        raise ValueError('Pinned environment differs: ' + str({n: (required[n], actual[n]) for n in required if required[n] != actual[n]}))
    return dict(python=sys.version, packages=actual)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    sub = parser.add_subparsers(dest='command', required=True)
    sub.add_parser('verify', help='Verify the allowlist and all packaged bytes; no external inputs')
    inputs = sub.add_parser('inputs', help='Verify the six provider files before analysis')
    inputs.add_argument('--inputs-root', type=Path, required=True)
    tests = sub.add_parser('tests', help='Run the included scientific and release-safety tests')
    tests.add_argument('--output-root', type=Path, required=True)
    reproduce = sub.add_parser('reproduce', help='Reproduce analyses from the six original provider files')
    reproduce.add_argument('--inputs-root', type=Path, required=True)
    reproduce.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    if args.command == 'verify':
        result = verify()
    elif args.command == 'inputs':
        result = check_inputs(args.inputs_root)
    else:
        verify()
        versions = check_environment()
        out = separate(args.output_root, *([args.inputs_root] if args.command == 'reproduce' else []))
        if args.command == 'reproduce':
            check_inputs(args.inputs_root)
        out.mkdir(parents=True)
        env = environment(out)
        script = 'run_tests.py' if args.command == 'tests' else 'reproduce_all.py'
        command = [sys.executable, '-I', '-B', str(PACKAGE / 'src' / script), '--output-root', str(out)]
        if args.command == 'reproduce':
            command += ['--inputs-root', str(args.inputs_root.resolve())]
        receipt = dict(command=command, environment=versions, scope=args.command)
        (out / 'environment.json').write_text(json.dumps(receipt, indent=2)+'\n', encoding='utf-8')
        subprocess.run(command, cwd=PACKAGE, env=env, check=True)
        result = dict(status='PASS', scope=args.command, output_root=str(out))
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
