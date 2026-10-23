"""Independent CLVer source arithmetic and Poisson likelihood reconstruction.

Does not import the production scientific kernels or use statsmodels' optimizer.
Clinical verification reads original archive tables, not the transformed ledger.
"""
from pathlib import Path
import argparse
import hashlib
import json
import zipfile

import numpy as np
import pandas as pd


def digest(path):
    with Path(path).open('rb') as stream:
        return hashlib.file_digest(stream, 'sha256').hexdigest()


def original_tables(root):
    with zipfile.ZipFile(root / 'raw' / 'CLVerPublicDataset.zip') as z:
        return {name: pd.read_csv(z.open('Data Tables/' + name + '.txt'), sep='|',
                                 dtype=str, keep_default_na=False, encoding='utf-8-sig')
                for name in ['SampleResults', 'InsulinUse', 'DiabPhysExam', 'VisitInfo']}


def clinical(root):
    tables = original_tables(root)
    samples = tables['SampleResults']
    derived = pd.read_csv(root / 'processed' / 'clinical_derivations.csv')
    checks = []; issues = []
    for row in derived.itertuples():
        pid = str(int(row.PtID)); visit = f'{int(row.index_visit)} Week'
        s = samples[(samples.PtID == pid) & (samples.Visit == visit)]
        if row.eligible:
            values = []
            for minute in [0,15,30,60,90,120]:
                slot = s[s.ResultName.eq(f'C-PEP-{minute}')]
                v = pd.to_numeric(slot.Value, errors='coerce')
                assert len(v) and v.notna().all() and v.ge(0).all() and v.nunique() == 1
                assert slot.Units.str.strip().eq('nmol/L').all()
                assert slot.STAResultStatus.str.strip().str.lower().eq('complete').all()
                assert pd.to_datetime(slot.CollectionDt).dt.normalize().eq(pd.Timestamp(row.collection_date)).all()
                if minute == 0:
                    assert slot.STASampleFastStatus.str.strip().eq('Fasting').all()
                values.append(float(v.iloc[0]))
            # Direct quadrature weights differ from production's adjacent-pair integral.
            auc = float(np.dot([7.5,15,22.5,30,30,15], values) / 120)
            np.testing.assert_allclose(auc, row.auc_nmol_l, rtol=0, atol=1e-12)
            np.testing.assert_allclose(max(values), row.peak_nmol_l, rtol=0, atol=1e-12)
            checks.append(('auc', int(row.PtID), int(row.index_visit)))
        if np.isfinite(row.remission):
            key = str(int(row.visit_key))
            def unique(name, fields):
                frame = tables[name]
                part = frame[frame.PtID.eq(pid) & frame.ParentLoginVisitID.eq(key)][fields].drop_duplicates()
                assert len(part) == 1, (name, pid, key)
                return part.iloc[0]
            weight = unique('DiabPhysExam', ['Weight','WeightUnits','WeightUnk'])
            dose = unique('InsulinUse', ['TDIUnits','TDIUnitsUnk','InsDoseInfoDt'])
            visitrow = unique('VisitInfo', ['Visit','VisitDt'])
            assert visitrow.Visit == visit
            assert weight.WeightUnk != '1' and dose.TDIUnitsUnk != '1'
            kg = float(weight.Weight)
            if weight.WeightUnits.strip() == 'lbs': kg *= 0.45359237
            else: assert weight.WeightUnits.strip() == 'kg'
            hb = s[s.ResultName.eq('GLYHB')][['Value','CollectionDt']].drop_duplicates()
            assert len(hb) == 1
            h = float(hb.iloc[0].Value); daily = float(dose.TDIUnits)
            anchor = pd.Timestamp(row.anchor_date).normalize()
            for date in (visitrow.VisitDt, dose.InsDoseInfoDt, hb.iloc[0].CollectionDt):
                assert abs((pd.Timestamp(date).normalize() - anchor).days) <= 7
            result = h + 4 * daily / kg
            np.testing.assert_allclose([kg, daily/kg, result],
                [row.weight_kg, row.insulin_u_kg_day, row.idaa1c], rtol=0, atol=1e-11)
            assert int(result <= 9) == int(row.remission)
            checks.append(('idaa1c', int(row.PtID), int(row.index_visit)))
    result = dict(status='PASS', implementation='Independent original-table joins and fixed quadrature weights; no production import',
                  archive_sha256=digest(root/'raw'/'CLVerPublicDataset.zip'),
                  canonical_clinical_sha256=digest(root/'processed'/'clinical_derivations.csv'),
                  source_visit_rows=len(derived),
                  auc_visits_verified=sum(x[0]=='auc' for x in checks),
                  clinical_remission_visits_verified=sum(x[0]=='idaa1c' for x in checks),
                  scope='Reconstruct every included laboratory value and known clinical state; does not certify all exclusion decisions or clinical truth')
    (root/'results').mkdir(exist_ok=True)
    (root/'results'/'independent_clinical_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


def poisson_newton(X, y, offset):
    if np.linalg.matrix_rank(X) != X.shape[1] or y.sum() == 0:
        raise ValueError('Unidentified design or no events')
    beta = np.zeros(X.shape[1]); beta[0] = np.log(y.sum() / np.exp(offset).sum())
    def objective(b):
        eta = offset + X @ b
        if np.max(eta) > 700: return np.inf
        return np.exp(eta).sum() - y @ eta
    for iteration in range(100):
        mu = np.exp(offset + X @ beta)
        gradient = X.T @ (mu-y)
        hessian = X.T @ (mu[:,None]*X)
        step = np.linalg.solve(hessian, gradient)
        if np.max(np.abs(step)) < 1e-10:
            return beta, mu, iteration
        alpha = 1.; before = objective(beta)
        for _ in range(45):
            trial = beta-alpha*step
            if objective(trial) <= before - 1e-4*alpha*(gradient @ step) + 1e-10:
                beta = trial; break
            alpha *= .5
        else: raise ValueError('Independent Newton line search did not converge')
        if np.linalg.norm(beta) > 1000:
            raise ValueError('Potential nonfinite optimum')
    raise ValueError('Independent Newton iteration limit')


def independent_design(f):
    assert f.target_visit.isin([26,39,52]).all()
    assert f.management.isin(['usual care','intensive/Tandem','intensive/Medtronic']).all()
    assert f.cohort.isin(['A','B']).all()
    return np.column_stack([np.ones(len(f)), f.auc_nmol_l/.20,
        f.target_visit.eq(39),f.target_visit.eq(52),(f.lag_days-60)/30,
        f.age-12,f.management.eq('intensive/Tandem'),f.management.eq('intensive/Medtronic'),f.cohort.eq('B')]).astype(float)


def model(root, expected_path, draws_path=None):
    path = root/'results'/'primary_intervals.csv'
    raw = pd.read_csv(path)
    required = ['unwarned','observed_days','auc_nmol_l','target_visit','lag_days','age','management','cohort','PtID']
    f = raw.dropna(subset=required).reset_index(drop=True)
    X = independent_design(f); y = f.unwarned.to_numpy(dtype=float); offset = np.log(f.observed_days.to_numpy(dtype=float))
    assert np.isfinite(X).all() and np.isfinite(offset).all() and np.all(y>=0) and np.all(y==np.floor(y))
    beta, mu, iterations = poisson_newton(X,y,offset)
    inverse = np.linalg.inv(X.T@(mu[:,None]*X))
    groups = f.PtID.to_numpy(); ids = np.unique(groups)
    scores = np.stack([X[groups==pid].T@(y[groups==pid]-mu[groups==pid]) for pid in ids])
    n,k = X.shape; m=len(ids)
    covariance = inverse @ (scores.T@scores) @ inverse * m/(m-1)*(n-1)/(n-k)
    se=float(np.sqrt(covariance[1,1])); rr=float(np.exp(beta[1]))
    expected=json.loads(expected_path.read_text())
    if 'rr' not in expected and 'primary' in expected: expected=expected['primary']
    np.testing.assert_allclose([beta[1],rr,se],[expected['beta'],expected['rr'],expected['robust_se']],rtol=1e-6,atol=1e-8)
    result=dict(status='PASS',implementation='Independent analytic Newton optimizer and explicitly assembled cluster sandwich; no production kernel or statsmodels optimizer',
        n_intervals=n,n_children=m,dropped_essential_covariates=len(raw)-n,iterations=iterations,
        beta=float(beta[1]),rr=rr,robust_se=se,
        robust_95_ci=np.exp([beta[1]-1.959963984540054*se,beta[1]+1.959963984540054*se]).tolist(),
        data_sha256=digest(path),expected_result_sha256=digest(expected_path))
    if draws_path:
        draws=pd.read_csv(draws_path,sep='\t' if draws_path.suffix=='.tsv' else ',')
        valid=draws[draws.rr.notna()].head(25)
        errors=[]
        for draw in valid.itertuples():
            picks=[int(float(p)) for p in json.loads(draw.sampled_child_ids)]
            ix=np.concatenate([np.flatnonzero(groups==pid) for pid in picks])
            b,_,_=poisson_newton(X[ix],y[ix],offset[ix])
            delta=abs(np.exp(b[1])-draw.rr); errors.append(float(delta))
            np.testing.assert_allclose(np.exp(b[1]),draw.rr,rtol=1e-5,atol=1e-7)
        result['first_valid_bootstrap_draws_reconstructed']=len(errors)
        result['max_bootstrap_RR_absolute_difference']=max(errors) if errors else None
    (root/'results'/'independent_primary_verification.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--data-root',type=Path,required=True)
    parser.add_argument('--clinical-only',action='store_true')
    parser.add_argument('--expected-result',type=Path)
    parser.add_argument('--bootstrap-draws',type=Path)
    args=parser.parse_args()
    clinical(args.data_root)
    if not args.clinical_only:
        if args.expected_result is None: parser.error('--expected-result required for primary fit verification')
        model(args.data_root,args.expected_result,args.bootstrap_draws)
