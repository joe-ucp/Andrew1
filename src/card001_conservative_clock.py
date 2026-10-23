"""Predeclared conservative-clock sensitivity, including training refits.

Mask stalled transmitter counters, ambiguous same-local-time counter groups,
and their next distinct-time boundary within 30 minutes. Retain original primary
hyperparameters and operating points. Never select this result over primary.
"""
from pathlib import Path
import argparse
import hashlib
import json
import pickle
import time

import numpy as np
from card001_core import SEED
from card001_models import design, participant_weights


def sha(path):
    with Path(path).open('rb') as stream: return hashlib.file_digest(stream,'sha256').hexdigest()


def prepare(root):
    out=root/'conservative_clock'; out.mkdir(exist_ok=True)
    target=out/'participant_manifest.json'
    if target.exists():
        previous=json.loads(target.read_text())
        if previous['parent_manifest_sha256']!=sha(root/'participant_manifest.json'):
            raise ValueError('Original manifest changed')
        return out
    manifest=json.loads((root/'participant_manifest.json').read_text())
    manifest['parent_manifest_sha256']=sha(root/'participant_manifest.json')
    manifest['analysis_role']='conservative_clock_sensitivity'
    manifest['sensitivity_code_sha256']=sha(__file__)
    manifest['predeclaration']='Rule specified before primary test metrics were opened; no outcome-based thresholds, tuning or participant selection'
    manifest['rules']['additional_clock_rule']=__doc__
    changed=[]; processed=out/'processed'; processed.mkdir(exist_ok=True)
    for row in manifest['participants']:
        path=root/'clocks'/f"DCLP3_{row['id']}.npz"
        if not path.exists(): continue
        with np.load(path) as raw:
            t=raw['time']; tx=raw['transmitter']; dt=np.diff(t)
            stalled=t[1:][(dt>0)&(dt<=1800)&(np.diff(tx)==0)]
            unique,first,count=np.unique(t,return_index=True,return_counts=True)
            boundary=[]; ambiguous=0
            for j in np.flatnonzero(count>1):
                values=tx[first[j]:first[j]+count[j]]
                if not np.all((values==values[0]) | (np.isnan(values)&np.isnan(values[0]))):
                    ambiguous+=1; boundary.append(unique[j])
                    if j+1<len(unique) and unique[j+1]-unique[j]<=1800:
                        boundary.append(unique[j+1])
            times=np.r_[stalled,boundary]
        if not len(times): continue
        original=Path(row['file'])
        with np.load(original) as z: a={k:z[k].copy() for k in z.files}
        mask=np.isin(a['time'],np.ceil(times/300)*300)
        removed=int((mask&np.isfinite(a['glucose'])).sum())
        a['glucose'][mask]=np.nan; a['clock_barrier']|=mask
        targetfile=processed/original.name
        np.savez_compressed(targetfile,**a)
        row['file']=str(targetfile)
        changed.append(dict(id=row['id'],split=row['split'],stalled_intervals=len(stalled),
                            ambiguous_counter_groups=ambiguous,newly_unavailable_bins=removed,
                            original_grid_sha256=sha(original),raw_clock_sha256=sha(path),
                            modified_grid_sha256=sha(targetfile)))
    manifest['clock_changes']=changed
    target.write_text(json.dumps(manifest,indent=2))
    print(json.dumps({split:sum(r['newly_unavailable_bins'] for r in changed if r['split']==split)
                      for split in ['train','validation','test']},indent=2),flush=True)
    return out


def fit(root):
    from sklearn.ensemble import HistGradientBoostingClassifier
    from card001_run import load_participants
    import sklearn
    out=prepare(root); folder=out/'models'; folder.mkdir(exist_ok=True)
    participants=load_participants(out,['train'])
    for kind in ['M120','M15','missing']:
        selected_path=root/'models'/f'{kind}_fit.json'
        selected=json.loads(selected_path.read_text())['selected']
        modelpath=folder/f'{kind}.pkl'; commitpath=folder/f'{kind}_commitment.json'
        commitment=dict(code_sha256=sha(__file__),manifest_sha256=sha(out/'participant_manifest.json'),
            primary_selection_sha256=sha(selected_path),sklearn=sklearn.__version__,
            source_sha256={name:sha(Path(__file__).parent/name) for name in ['card001_core.py','card001_models.py']},
            max_leaf_nodes=selected['max_leaf_nodes'],l2_regularization=selected['l2_regularization'],
            max_iter=200,min_samples_leaf=50,learning_rate=.05,early_stopping=False,seed=SEED)
        if modelpath.exists():
            if not commitpath.exists() or json.loads(commitpath.read_text())!=commitment:
                raise ValueError('Alternative model lacks matching commitment')
            continue
        commitpath.write_text(json.dumps(commitment,indent=2))
        started=time.perf_counter(); X,y,groups,omissions=design(participants,kind)
        model=HistGradientBoostingClassifier(max_leaf_nodes=selected['max_leaf_nodes'],l2_regularization=selected['l2_regularization'],
            max_iter=200,min_samples_leaf=50,learning_rate=.05,early_stopping=False,random_state=SEED)
        model.fit(X,y,sample_weight=participant_weights(groups))
        with modelpath.open('wb') as stream: pickle.dump(model,stream)
        result=dict(training_rows=len(y),omissions=omissions,model_sha256=sha(modelpath),seconds=time.perf_counter()-started)
        (folder/f'{kind}_fit.json').write_text(json.dumps(result,indent=2))
        print(f"conservative_clock_model={kind}; seconds={result['seconds']:.1f}",flush=True)
        del X,y,groups,model


def evaluate(root,external_root=None):
    from card001_evaluate import test_run
    out=prepare(root); cohort='external' if external_root else 'internal'
    if not (root/'results'/f'{cohort}_comparison.json').exists(): raise ValueError('Primary replay must precede sensitivity')
    source=root/'calibration.json'; frozen=json.loads(source.read_text())
    clean=lambda p:{k:p[k] for k in ['kind','threshold','short_threshold'] if k in p}
    alternate=dict(status='FROZEN_SENSITIVITY_ORIGINAL_OPERATING_POINTS',
        original_calibration_sha256=sha(source),model_sha256={k:sha(out/'models'/f'{k}.pkl') for k in ['M120','M15','missing']},
        manifest_sha256=sha(out/'participant_manifest.json'),comparator=clean(frozen['comparator']),
        comparators={k:clean(v) for k,v in frozen['comparators'].items()},switch=clean(frozen['switch']),
        development_precision=frozen['development_precision'],
        interpretation='Original operating points and original development precision audit; no alternate validation selection')
    target=out/'calibration.json'
    if target.exists() and json.loads(target.read_text())!=alternate: raise ValueError('Alternative replay commitment changed')
    target.write_text(json.dumps(alternate,indent=2))
    test_run(out,external_root,prefix_override=f'conservative_clock_{cohort}')


if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__); p.add_argument('--data-root',type=Path,required=True)
    p.add_argument('--stage',choices=['prepare','fit','evaluate'],required=True); p.add_argument('--external-root',type=Path)
    a=p.parse_args()
    if a.stage=='prepare': prepare(a.data_root)
    elif a.stage=='fit': fit(a.data_root)
    else: evaluate(a.data_root,a.external_root)
