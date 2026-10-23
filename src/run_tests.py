"""Run the public tests without importing another checkout or user site."""
import argparse
import importlib
import io
import json
from pathlib import Path
import sys
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / 'src'))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-root', type=Path, required=True)
    args = parser.parse_args()
    out = args.output_root.resolve()
    if out == ROOT or out.is_relative_to(ROOT) or ROOT.is_relative_to(out):
        raise ValueError('Test output must be separate from the release')
    if out.exists() and any(p.name not in {'environment.json', 'runtime'} for p in out.iterdir()):
        raise FileExistsError('Use a new test output directory')
    out.mkdir(parents=True, exist_ok=True)
    suite = unittest.defaultTestLoader.discover(str(ROOT / 'tests'))
    stream = io.StringIO()
    result = unittest.TextTestRunner(stream=stream, verbosity=2).run(suite)
    (out / 'tests.log').write_text(stream.getvalue(), encoding='utf-8')
    modules = {}
    for name, module in sorted(sys.modules.items()):
        if name.startswith(('card001', 'card002', 'card003', 'revision_', 'test_')):
            path = Path(module.__file__).resolve()
            if not path.is_relative_to(ROOT):
                raise RuntimeError('Import escaped release: ' + name)
            modules[name] = path.relative_to(ROOT).as_posix()
    report = dict(status='PASS' if result.wasSuccessful() else 'FAIL', tests=result.testsRun,
                  failures=len(result.failures), errors=len(result.errors), skipped=len(result.skipped),
                  python=sys.version, imported_scientific_modules=modules,
                  scope='Implementation and numerical unit checks; source-data reproduction is separate')
    (out / 'tests.json').write_text(json.dumps(report, indent=2)+'\n', encoding='utf-8')
    print(stream.getvalue(), end='')
    print(json.dumps(report, indent=2))
    if not result.wasSuccessful():
        raise SystemExit(1)


if __name__ == '__main__':
    main()
