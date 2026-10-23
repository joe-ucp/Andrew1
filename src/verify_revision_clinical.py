"""Independent Newton and explicit cluster-covariance checks of revision fits.

Does not import the primary or revision clinical analysis implementation.
This verifies same-data arithmetic, not causal identification or external validity.
"""
from pathlib import Path
import argparse
import json
import numpy as np
import pandas as pd
from scipy.special import expit


def design(f):
    return np.column_stack([np.ones(len(f)), f.auc_nmol_l / .20,
        f.target_visit.eq(39), f.target_visit.eq(52), (f.lag_days-60)/30,
        f.age-12, f.management.eq('intensive/Tandem'),
        f.management.eq('intensive/Medtronic'), f.cohort.eq('B')]).astype(float)


def fit_newton(x, y, offset, binary=False):
    beta = np.zeros(x.shape[1])
    beta[0] = np.log(y.sum()/np.exp(offset).sum()) if not binary else 0
    def terms(b):
        eta = x @ b + offset
        if binary:
            mu = expit(eta)
            loss = np.logaddexp(0, eta).sum() - y @ eta
            w = mu * (1-mu)
        else:
            mu = np.exp(eta)
            loss = mu.sum() - y @ eta
            w = mu
        return loss, mu, w
    for _ in range(100):
        old, mu, w = terms(beta)
        score = x.T @ (y-mu)
        info = x.T @ (w[:, None]*x)
        step = np.linalg.solve(info, score)
        if np.max(np.abs(step)) < 1e-11:
            break
        rate = 1.
        while terms(beta+rate*step)[0] > old+1e-10:
            rate /= 2
            if rate < 1e-10:
                raise AssertionError('Newton line search failed')
        beta += rate*step
    else:
        raise AssertionError('Newton did not converge')
    _, mu, w = terms(beta)
    return beta, mu, x.T @ (w[:, None]*x)


def sandwich(x, y, mu, info, ids):
    unique = pd.unique(ids)
    scores = np.stack([x[ids == p].T @ (y[ids == p]-mu[ids == p]) for p in unique])
    bread = np.linalg.inv(info)
    n, p = x.shape
    return bread @ (scores.T @ scores) @ bread * len(unique)/(len(unique)-1)*(n-1)/(n-p)


def main(root):
    rev = root/'revision_20261022_r2'
    r = rev/'card003/results'
    f = pd.read_csv(root/'card003/results/primary_intervals.csv', dtype={'PtID':str})
    subset = f.loc[f.groupby('PtID').unwarned.transform('sum') > 0].copy()
    assert (len(subset), subset.PtID.nunique()) == (164, 58)
    x = design(subset); y = subset.unwarned.to_numpy(float)
    b, mu, info = fit_newton(x, y, np.log(subset.observed_days.to_numpy()))
    cov = sandwich(x,y,mu,info,subset.PtID.to_numpy())
    saved = json.loads((r/'same_sample_pooled.json').read_text())
    np.testing.assert_allclose(b, list(saved['coefficients'].values()), atol=1e-9, rtol=0)
    np.testing.assert_allclose(cov, saved['cluster_covariance'], atol=1e-9, rtol=0)

    # Expand interval binomial counts independently of the event-table adapter.
    event_rows = np.repeat(np.arange(len(f)), f.events.astype(int))
    ey = np.concatenate([np.r_[np.ones(int(row.unwarned)),np.zeros(int(row.matched))] for row in f.itertuples()])
    ef = f.iloc[event_rows]; ex = design(ef)
    assert len(ef) == 408 and ef.PtID.nunique() == 72
    eb, em, ei = fit_newton(ex,ey,np.zeros(len(ey)),True)
    ec = sandwich(ex,ey,em,ei,ef.PtID.to_numpy())
    logistic = json.loads((r/'conditional_warning_logistic.json').read_text())
    np.testing.assert_allclose(eb,list(logistic['coefficients'].values()),atol=1e-9,rtol=0)
    np.testing.assert_allclose(ec,logistic['cluster_covariance'],atol=1e-9,rtol=0)

    draws = pd.read_csv(r/'paired_child_bootstrap.csv.gz')
    assert len(draws) == 5000 and draws[['unwarned_beta','events_beta']].notna().all().all()
    np.testing.assert_allclose(draws.difference_log_rr,draws.unwarned_beta-draws.events_beta,atol=1e-14)
    contrast = json.loads((r/'paired_comparison.json').read_text())
    np.testing.assert_allclose(np.percentile(draws.difference_log_rr,[2.5,97.5]),contrast['difference_log_rr_95_ci'],atol=1e-14)
    rows = {p: np.flatnonzero(f.PtID.to_numpy()==p) for p in f.PtID.unique()}
    bx = design(f); offset = np.log(f.observed_days.to_numpy())
    max_error = 0.
    for draw in draws.iloc[np.linspace(0,4999,25,dtype=int)].itertuples():
        ix = np.concatenate([rows[str(p)] for p in json.loads(draw.sampled_child_ids)])
        for outcome in ('unwarned','events'):
            bb,_,_ = fit_newton(bx[ix],f[outcome].to_numpy(float)[ix],offset[ix])
            error = abs(bb[1]-getattr(draw,outcome+'_beta'))
            max_error=max(max_error,error)
            assert error < 1e-8
    bands = pd.read_csv(r/'joint_auc_bands.csv')
    assert bands.events.sum()==408 and bands.unwarned.sum()==180 and bands.intervals.sum()==276
    result = dict(status='PASS',same_sample_children=58,same_sample_intervals=164,
        independent_same_sample_rr=float(np.exp(b[1])),independent_conditional_or=float(np.exp(eb[1])),
        all_paired_draw_quantiles_verified=5000,independently_refit_draws=25,
        largest_refit_log_coefficient_error=max_error,method='Independent Newton solver and explicit cluster covariance; same-data arithmetic only')
    (rev/'qa').mkdir(parents=True,exist_ok=True)
    (rev/'qa/independent_clinical_revision.json').write_text(json.dumps(result,indent=2))
    print(json.dumps(result,indent=2))


if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--data-root',type=Path,required=True)
    main(p.parse_args().data_root)
