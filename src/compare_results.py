"""Compare reproduced scientific evidence and figures with the authoritative release."""
import argparse
import csv
import json
import math
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))
from release import PACKAGE, contained, load, sha

ATOL = 1e-10
RTOL = 1e-10
# These fields describe execution time or a manifest containing output locations.
# They are reported as exclusions; no scientific estimate or decision is excluded.
EXCLUSIONS = {
    'card001/results/internal_comparison.json': {'cohort_manifest_sha256'},
    'card001/results/external_comparison.json': {'cohort_manifest_sha256'},
    **{'card003/results/model_'+n+'.json': {'elapsed_seconds'} for n in
       ['primary', 'total_events', 'remission', 'quadratic_AUC']},
    # The portable verifier corrects the description of the old endpoint assertion.
    # Ascertainment status/counts and every numerical field are still compared.
    'revision_20261022_r2/card002/endpoint_provenance_revision.json': {'original_assertion_narrowed'},
}
# One aggregate reuses a full fitted-result object, including its runtime.
# Exclude only this exact nested field; neighboring scientific fields stay strict.
NESTED_EXCLUSIONS = {
    'revision_20261022_r2/card003/results/never_warning_identity.json':
        [('original_total_result_reused', 'elapsed_seconds')],
}
FIGURES = {
    'warning_diagnostics.png': 'revision_20261022_r2/card001/principal_warning_diagnostics.png',
    'secretion_rates.png': 'figures/secretion_rates.png',
    'tn10_curves.png': 'revision_20261022_r2/card002/tn10_outcome_curves.png',
    'tn10_uncertainty_profile.png': 'revision_20261022_r3/paper/tn10/tn10_bootstrap_followup_diagnostics.png',
}


def compare(expected, actual, trail=''):
    if isinstance(expected, bool) or expected is None:
        if expected != actual: raise AssertionError(trail)
    elif isinstance(expected, (int,float)):
        if not isinstance(actual,(int,float)): raise AssertionError(trail+' type')
        if float(expected).is_integer() and float(actual).is_integer():
            if expected != actual: raise AssertionError(trail+' integer')
        elif math.isnan(expected):
            if not math.isnan(actual): raise AssertionError(trail+' NaN')
        elif not math.isclose(expected, actual, rel_tol=RTOL, abs_tol=ATOL):
            raise AssertionError(f'{trail}: {expected} != {actual}')
    elif isinstance(expected, dict):
        if set(expected) != set(actual): raise AssertionError(trail+' keys')
        for k in expected: compare(expected[k],actual[k],trail+'/'+str(k))
    elif isinstance(expected,list):
        if len(expected)!=len(actual): raise AssertionError(trail+' length')
        for i,(x,y) in enumerate(zip(expected,actual)):compare(x,y,trail+'/'+str(i))
    elif expected != actual:
        raise AssertionError(f'{trail}: {expected!r} != {actual!r}')


def csv_value(value):
    if value=='':return value
    try:return float(value)
    except ValueError:return value


def compare_file(reference, actual, relative, fields=None):
    if sha(reference)==sha(actual):return dict(method='exact SHA-256',sha256=sha(actual))
    if reference.suffix=='.json':
        left,right=load(reference),load(actual)
        # Only declared release projections omit participant-specific parameters.
        if fields is not None:
            if not isinstance(left,dict) or set(left) != set(fields):
                raise AssertionError(relative+' declared projection fields')
            right={key:right[key] for key in left}
        omitted=EXCLUSIONS.get(relative,set())
        changes={key:dict(authoritative=left.pop(key),reproduced=right.pop(key)) for key in omitted}
        for path in NESTED_EXCLUSIONS.get(relative, []):
            left_parent, right_parent = left, right
            for key in path[:-1]:
                left_parent, right_parent = left_parent[key], right_parent[key]
            changes['.'.join(path)] = dict(authoritative=left_parent.pop(path[-1]),
                                            reproduced=right_parent.pop(path[-1]))
        compare(left,right,relative)
        return dict(method='structured comparison', excluded_execution_metadata=changes,
                    authoritative_sha256=sha(reference),reproduced_sha256=sha(actual))
    if reference.suffix in {'.csv','.tsv'}:
        separator='\t' if reference.suffix=='.tsv' else ','
        with reference.open(newline='',encoding='utf-8') as f: left=list(csv.reader(f,delimiter=separator))
        with actual.open(newline='',encoding='utf-8') as f: right=list(csv.reader(f,delimiter=separator))
        if left[0]!=right[0]:raise AssertionError(relative+' column headers')
        compare([[csv_value(v) for v in row] for row in left[1:]],
                [[csv_value(v) for v in row] for row in right[1:]],relative)
        return dict(method='ordered CSV comparison', rows=len(left)-1,
                    authoritative_sha256=sha(reference),reproduced_sha256=sha(actual))
    raise AssertionError('Byte difference: '+relative)


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--run-root',type=Path,required=True)
    parser.add_argument('--output',type=Path,required=True)
    args=parser.parse_args()
    out=args.output.resolve()
    if out.exists() or out.is_relative_to(PACKAGE):
        raise ValueError('Comparison output must be new and outside the package')
    records={}; failures={}
    for entry in load(PACKAGE/'provenance/evidence_sources.json'):
        relative=entry.get('generated_path',entry.get('external_cache_path'))
        try: records[relative]=compare_file(contained(PACKAGE,entry['package_path']),
                                           contained(args.run_root,relative),relative,entry.get('fields'))
        except Exception as error:failures[relative]=str(error)
    from PIL import Image, ImageChops
    images={}
    for name,relative in FIGURES.items():
        try:
            with Image.open(PACKAGE/'manuscript/figures'/name) as left, Image.open(contained(args.run_root,relative)) as right:
                if left.size!=right.size or ImageChops.difference(left.convert('RGB'),right.convert('RGB')).getbbox():
                    raise AssertionError('Figure pixels differ')
            images[name]=dict(method='exact RGB pixels',source=relative)
        except Exception as error:failures[name]=str(error)
    report=dict(status='PASS' if not failures else 'FAIL',evidence=records,figures=images,failures=failures,
                numeric_tolerance=dict(atol=ATOL,rtol=RTOL),
                rationale='Integer JSON counts and strings are exact; floats permit only binary/serialization roundoff. Row order and all scientific fields are preserved. Execution-only exclusions are explicit. Figures require exact pixels in the declared environment.')
    out.parent.mkdir(parents=True,exist_ok=True)
    out.write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps(dict(status=report['status'],evidence_checked=len(records),figures_checked=len(images),failures=failures),indent=2))
    if failures:raise SystemExit(1)


if __name__=='__main__':main()
