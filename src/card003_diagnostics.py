"""Actual-data invariance and quadratic interpretation checks; no primary retuning."""
from pathlib import Path
import argparse,json
import numpy as np
import pandas as pd
from card003_analyze import load_frozen,predict_grid,predicted_curve,dump
from card003_core import fit_association,fixed_design


def run(root):
    root=Path(root)
    f=pd.read_csv(root/"results"/"primary_intervals.csv",dtype={"PtID":str})
    frozen,models,_=load_frozen(root)
    primary=json.loads((root/"results"/"model_primary.json").read_text())
    ids=sorted(f.PtID.unique());rename={p:f"renamed_{len(ids)-j:03d}" for j,p in enumerate(ids)}
    renamed=f.copy();renamed["PtID"]=renamed.PtID.map(rename)
    check,_=fit_association(renamed,n_boot=0)
    diffs={k:abs(check[k]-primary[k]) for k in ("beta","rr","robust_se")}
    if max(diffs.values())>1e-10:raise AssertionError("Patient-ID relabeling changed the numerical estimate")
    prefixes=[]
    for source,group in f.groupby("source"):
        row=group.iloc[0]
        with np.load(row.grid_path) as data:g=data["glucose"]
        cut=len(g)//2
        changed=g.copy();changed[cut:]=np.where(np.arange(len(g)-cut)%2,350.,45.)
        original=predict_grid(models,g);perturbed=predict_grid(models,changed)
        for kind in original:
            a,b=original[kind][:cut],perturbed[kind][:cut]
            if not np.array_equal(a,b,equal_nan=True):raise AssertionError(f"Future glucose changed prior {kind} prediction")
            prefixes.append(dict(source=source,kind=kind,prior_bins=cut,identical=True))
    result=dict(status="PASS",patient_relabeling_absolute_differences=diffs,
        actual_model_prefix_invariance=prefixes,scope="Actual canonical cohort and frozen model arithmetic; no clinical efficacy inference")
    dump(root/"results"/"invariance_diagnostics.json",result)
    # Full uncertainty for the nonconstant quadratic exposure contrast.
    import statsmodels.api as sm
    from statsmodels.stats.sandwich_covariance import cov_cluster
    X=fixed_design(f,quadratic=True)
    fit=sm.GLM(f.unwarned.to_numpy(),X.to_numpy(),family=sm.families.Poisson(),offset=np.log(f.observed_days)).fit(maxiter=100)
    cov=cov_cluster(fit,f.PtID.to_numpy());k=list(X.columns).index("exposure");q=list(X.columns).index("exposure_squared")
    rows=[]
    for quantile in (.05,.25,.5,.75,.95):
        value=float(f.auc_nmol_l.quantile(quantile));contrast=np.zeros(len(X.columns));contrast[k]=1;contrast[q]=2*value/.20+1
        logrr=float(contrast@fit.params);se=float(np.sqrt(contrast@cov@contrast));ci=np.exp([logrr-1.95996398454*se,logrr+1.95996398454*se])
        rows.append(dict(quantile=quantile,auc_nmol_l=value,contrast_to=value+.20,
                         inside_observed_support=value+.20<=f.auc_nmol_l.max(),rr=float(np.exp(logrr)),robust_95_ci=ci.tolist()))
    diagnostic=dict(model="quadratic AUC",coefficients=dict(zip(X.columns,fit.params.tolist())),
        covariance_columns=list(X.columns),child_cluster_covariance=cov.tolist(),
        squared_term_beta=float(fit.params[q]),squared_term_robust_se=float(np.sqrt(cov[q,q])),
        local_per_0_20_contrasts=rows,
        meaning="exp(beta_linear) alone is not a finite +0.20 contrast. At AUC x, RR(x to x+0.20)=exp(beta_linear+beta_squared*(2*x/0.20+1)).")
    dump(root/"results"/"quadratic_diagnostics.json",diagnostic)
    predicted_curve(root,f)
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data-root",type=Path,required=True);a=p.parse_args()
    print(json.dumps(run(a.data_root)))
