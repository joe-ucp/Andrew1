"""Frozen TN10 analysis: public randomized delay contrast.

Generated artifacts remain beneath --data-root. Supply the three original
inputs with --raw-root; they are read-only and verified by SHA-256 and size.
"""
from __future__ import annotations
import argparse
import hashlib
import importlib.metadata
import json
import os
from pathlib import Path
import re
import sys
from datetime import datetime, timezone

SOURCES = {
    "TN10_metadata.csv": ("https://zenodo.org/records/17063348/files/TN10_metadata.csv?download=1", "dabd7beff35f82aef231f9c9e0fd5b91974cdbc9f175edd8831c475824517586", 44999),
    "TN10_raw_counts.csv": ("https://zenodo.org/records/17063348/files/TN10_raw_counts.csv?download=1", "74085ba79ef79c44d5e25679e36f3615562c0a644700f2332426d71b942e1ee8", 8508942),
    "JCI176403_Table6.pdf": ("https://www.jci.org/articles/view/176403/sd/pdf/render/6", "93b00f6ccb5688eea28f936cfc8e72117e65f70b612c656531cea11c5b51f90c", 210945),
}

# Preserve the historical analysis identifier used in frozen result provenance.
CONFIG = {
    "card": "STRIKE_CARD_002", "outcome_category": 1,
    "population": "published_TN10_baseline_RNA_subset", "expected_n": 73,
    "expected_treatment_n": 42, "expected_placebo_n": 31,
    "source_gene_n": 359, "primary_gene_n": 312,
    "gene_mapping": "exact_symbol_intersection",
    "expression": "log2(1+CPM); library sum over all 13021 supplied rows",
    "standardization": "pooled baseline gene mean and sample SD",
    "score": "unweighted mean standardized gene expression",
    "split": "high iff score strictly exceeds pooled median",
    "horizon_weeks": 156, "hypothesis_margin_weeks": 26,
    "estimand": "treatment effect low minus treatment effect high",
    "interval": "four-cell Welch-Satterthwaite 95% t interval",
    "bootstrap_repetitions": 10000, "bootstrap_seed": 20261022,
    "bootstrap_stratification": "randomized arm", "clinical_deployment": False,
    "cure_claim": False,
}

def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, allow_nan=False) + "\n", encoding="utf-8")

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def retrieve(root, download=False, raw_root=None):
    from urllib.request import urlopen, Request
    def fetch(url):
        with urlopen(Request(url, headers={"User-Agent": "UCP-Research/1.0 public-source-reproduction"}), timeout=120) as response:
            return response.read()
    raw = Path(raw_root).resolve() if raw_root is not None else root / "raw"
    if raw_root is not None and download:
        raise ValueError("--raw-root is read-only; acquire original inputs separately")
    if raw_root is None:
        raw.mkdir(parents=True, exist_ok=True)
    manifest = []
    for name, (url, expected, size) in SOURCES.items():
        target = raw / name
        if not target.exists() and download:
            target.write_bytes(fetch(url))
        if not target.exists():
            raise FileNotFoundError(f"{target}; rerun with --download")
        actual = sha256(target)
        if actual != expected or target.stat().st_size != size:
            raise ValueError(f"Source identity failure for {name}: {actual}, {target.stat().st_size}")
        manifest.append({"name": name, "url": url, "sha256": actual, "bytes": size, "verified": True})
    api = raw / "zenodo_record.json"
    if not api.exists() and download:
        api.write_bytes(fetch("https://zenodo.org/api/records/17063348"))
    if api.exists():
        record = json.loads(api.read_text(encoding="utf-8"))
        license_info = record.get("metadata", {}).get("license")
        manifest.append({"name": api.name, "url": "https://zenodo.org/api/records/17063348", "sha256": sha256(api), "license": license_info})
    write_json(root / "retrieval_manifest.json", {"verified_at_utc": datetime.now(timezone.utc).isoformat(), "sources": manifest})
    return raw

def parse_table(path):
    import pdfplumber
    import pandas as pd
    with pdfplumber.open(path) as pdf:
        pages = [page.extract_text(layout=False) or "" for page in pdf.pages]
    text = "\n\n".join(pages)
    # Source table has symbol / Pearson correlation / P-value rows.
    pattern = re.compile(r"^\s*(\S+)\s+(-?\d+\.\d+)\s+([\d.]+(?:[Ee][+-]?\d+)?)\s*$", re.MULTILINE)
    rows = [(m[0], float(m[1]), float(m[2])) for m in pattern.findall(text)]
    frame = pd.DataFrame(rows, columns=["symbol", "pearson_correlation", "p_value"])
    return frame[frame.pearson_correlation < 0].copy(), text

def score_expression(expression):
    import numpy as np
    sd = expression.std(axis=0, ddof=1)
    if not np.all(np.isfinite(sd) & (sd > 0)):
        raise ValueError("Nonfinite or zero-variance signature gene")
    score = ((expression - expression.mean(axis=0)) / sd).mean(axis=1)
    return score, score > np.median(score)

def mean_summary(values):
    import numpy as np
    from scipy.stats import t
    n = len(values)
    if n < 2:
        raise ValueError("Insufficient comparison cell")
    mean = float(np.mean(values))
    var = float(np.var(values, ddof=1))
    se = float(np.sqrt(var / n))
    width = float(t.ppf(.975, n-1) * se)
    return {"n": n, "mean_weeks": mean, "sd_weeks": float(np.sqrt(var)),
            "variance": var, "se_weeks": se, "ci95_weeks": [mean-width, mean+width]}

def welch_contrast(cells, weights):
    import numpy as np
    from scipy.stats import t
    estimate = sum(weight * cell["mean_weeks"] for cell, weight in zip(cells, weights))
    parts = [weight**2 * cell["variance"] / cell["n"] for cell, weight in zip(cells, weights)]
    variance = sum(parts)
    se = float(np.sqrt(variance))
    denom = sum(part**2 / (cell["n"] - 1) for part, cell in zip(parts, cells))
    df = variance**2/denom if denom > 0 else 1e12
    width = float(t.ppf(.975, df) * se)
    return {"estimate_weeks": estimate, "estimate_months": estimate*7/30.4375,
            "se_weeks": se, "df": df, "ci95_weeks": [estimate-width, estimate+width],
            "p_two_sided_zero": float(2*t.sf(abs(estimate/se), df)) if se > 0 else None}

def estimate(y, arm, high):
    cells = {f"{g}_{a}": mean_summary(y[(high == (g == 'high')) & (arm == (a == 'treated'))])
             for g in ["low", "high"] for a in ["treated", "placebo"]}
    low = welch_contrast([cells["low_treated"], cells["low_placebo"]], [1, -1])
    high_effect = welch_contrast([cells["high_treated"], cells["high_placebo"]], [1, -1])
    theta = welch_contrast([cells[k] for k in ["low_treated", "low_placebo", "high_treated", "high_placebo"]], [1,-1,-1,1])
    prognosis = welch_contrast([cells["low_placebo"], cells["high_placebo"]], [1,-1])
    lo, hi = theta["ci95_weeks"]
    status = "supported" if lo >= 26 else "refuted" if hi < 26 else "inconclusive"
    return {"cells": cells, "delta_low": low, "delta_high": high_effect, "theta": theta,
            "placebo_low_minus_high": prognosis, "hypothesis_status": status,
            "card_verdict": {"supported":"supported", "refuted":"damaged", "inconclusive":"unresolved"}[status],
            "interval_crosses_zero": bool(lo <= 0 <= hi), "margin_weeks": 26}

def continuous_models(y, arm, score, baseline):
    import numpy as np
    import pandas as pd
    import statsmodels.api as sm
    age = pd.to_numeric(baseline["AGE_AT_RAND"], errors="coerce").to_numpy(float)
    hba1c = pd.to_numeric(baseline["HBA1C"], errors="coerce").to_numpy(float)
    complete = np.isfinite(age) & np.isfinite(hba1c)
    z = (score-score.mean())/score.std(ddof=1)
    def model(mask, adjusted):
        x = {"intercept": np.ones(mask.sum()), "treatment": arm[mask].astype(float),
             "score_z": z[mask], "treatment_x_score_z": arm[mask]*z[mask]}
        if adjusted:
            x.update(age=age[mask], hba1c=hba1c[mask])
        fit = sm.OLS(y[mask], pd.DataFrame(x)).fit(cov_type="HC3", use_t=True)
        ci = fit.conf_int()
        return {"n": int(mask.sum()), "residual_df": int(fit.df_resid), "covariance": "HC3", "interval_reference": "t residual df",
                "coefficients": {term: {"estimate":float(fit.params[term]), "se":float(fit.bse[term]),
                    "ci95":ci.loc[term].tolist(), "p":float(fit.pvalues[term])} for term in x}}
    return {"score_unit": "one pooled baseline composite-score sample SD", "omitted_covariates_n": int((~complete).sum()),
            "adjusted_complete": model(complete, True), "unadjusted_complete":model(complete, False),
            "unadjusted_full": model(np.ones(len(y), dtype=bool), False),
            "dichotomized_complete_fixed_primary_split":estimate(y[complete], arm[complete], (score>np.median(score))[complete])}

def bootstrap(expression, y, arm, n=10000, seed=20261022):
    import numpy as np
    rng = np.random.default_rng(seed)
    arm_indices = [np.flatnonzero(arm == a) for a in [False, True]]
    values = []
    unusable = 0
    for _ in range(n):
        idx = np.concatenate([rng.choice(i, len(i), replace=True) for i in arm_indices])
        _, high = score_expression(expression[idx])
        yy, aa = y[idx], arm[idx]
        parts = [yy[(aa == a) & (high == h)] for h in [False,True] for a in [True,False]]
        if min(map(len, parts)) < 2:
            unusable += 1
            continue
        means = [p.mean() for p in parts]
        low, high_effect = means[0]-means[1], means[2]-means[3]
        values.append([low, high_effect, low-high_effect])
    arr = np.asarray(values)
    return arr, {"repetitions":n, "seed":seed, "valid":len(values), "unusable":unusable,
        "recomputed_in_each_resample":["gene means", "gene sample SDs", "composite scores", "pooled median", "group labels"],
        "percentile_ci95_weeks":{k:np.percentile(arr[:,j],[2.5,97.5]).tolist() for j,k in enumerate(["delta_low","delta_high","theta"])},
        "theta_bootstrap_sd":float(arr[:,2].std(ddof=1)), "theta_negative_fraction":float((arr[:,2]<0).mean())}

def make_figures(root, primary):
    import numpy as np
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":11, "axes.spines.top":False, "axes.spines.right":False, "savefig.dpi":180})
    fig, ax = plt.subplots(figsize=(9,5.3), layout="constrained")
    names = ["Low / teplizumab", "Low / placebo", "High / teplizumab", "High / placebo"]
    cells = [primary["cells"][k] for k in ["low_treated","low_placebo","high_treated","high_placebo"]]
    means = np.array([c["mean_weeks"] for c in cells])
    cis = np.array([c["ci95_weeks"] for c in cells])
    ax.bar(range(4),means,color=["#166b8b","#96b4bf","#166b8b","#96b4bf"],width=.65)
    ax.errorbar(range(4),means,yerr=np.array([means-cis[:,0],cis[:,1]-means]),fmt="none",color="#172b36",capsize=5)
    ax.set_xticks(range(4),[f"{name}\nn = {c['n']}" for name,c in zip(names,cells)])
    ax.set_ylabel("Mean stage-3-free time within 156 weeks")
    ax.set_title("TN10 baseline RNA subset: observed restricted means")
    ax.set_ylim(0,max(170,float(cis[:,1].max()+12)))
    ax.text(.01,.98,"Error bars: 95% t intervals for cell means", transform=ax.transAxes, va="top",fontsize=9)
    fig.savefig(root/"four_cell_outcomes.png")
    plt.close(fig)
    fig, ax = plt.subplots(figsize=(9,5), layout="constrained")
    rows = [primary[k] for k in ["delta_low","delta_high","theta"]]
    labels = ["Low score: treatment - placebo (22 / 15)", "High score: treatment - placebo (20 / 16)", "Interaction: low benefit - high benefit (n = 73)"]
    for j,(row,label) in enumerate(zip(rows,labels)):
        est=row["estimate_weeks"]; lo,hi=row["ci95_weeks"]
        ax.errorbar(est,2-j,xerr=[[est-lo],[hi-est]],fmt="o",color="#166b8b",capsize=5)
        ax.text(hi+3,2-j,f"{est:.1f} [{lo:.1f}, {hi:.1f}]",va="center",fontsize=9)
    ax.axvline(0,color="#64717a",linewidth=1)
    ax.axvline(26,color="#bc6624",linestyle="--",label="Interaction margin: 26 weeks")
    ax.set_yticks([2,1,0],labels)
    ax.set_xlim(min(r["ci95_weeks"][0] for r in rows)-10,max(r["ci95_weeks"][1] for r in rows)+50)
    ax.set_ylim(-.6,2.7)
    ax.set_xlabel("Additional stage-3-free weeks over a 156-week horizon")
    ax.set_title("Randomized arm contrasts with Welch 95% intervals")
    ax.legend(loc="upper left",fontsize=9)
    fig.savefig(root/"forest_delay.png")
    plt.close(fig)

def run_analysis(root, metadata, counts, genes):
    import numpy as np
    import pandas as pd
    from scipy.stats import norm
    results = root / "results"
    results.mkdir(exist_ok=True)
    baseline = metadata.loc[metadata["Time.Point"] == "Baseline"].reset_index(drop=True)
    assert len(metadata)==195 and counts.shape==(13021,196)
    assert len(baseline)==73 and baseline.libraryID.is_unique
    assert counts.index.is_unique and counts.columns.is_unique
    assert baseline.libraryID.isin(counts.columns).all()
    assert set(baseline.TREATMENTDESC)=={"Teplizumab","Placebo"}
    arm = (baseline.TREATMENTDESC == "Teplizumab").to_numpy()
    assert int(arm.sum()) == 42
    assert len(genes)==359 and genes.symbol.is_unique
    matched = sorted(set(genes.symbol) & set(counts.index))
    gene_text = "\n".join(matched)+"\n"
    gene_hash = hashlib.sha256(gene_text.encode()).hexdigest()
    assert len(matched)==312 and gene_hash=="09c7b17296f6ec5b1b90357b555ec62748821592f042070a6f0ec525cd106f11"
    (results/"primary_genes.txt").write_bytes(gene_text.encode())
    genes["in_primary_projection"] = genes.symbol.isin(matched)
    genes.to_csv(results/"gene_coverage.tsv",sep="\t",index=False)
    genes.loc[~genes.in_primary_projection].to_csv(results/"missing_genes.tsv",sep="\t",index=False)
    raw = counts.loc[:,baseline.libraryID].to_numpy(float).T
    assert np.isfinite(raw).all() and (raw>=0).all()
    libsize=raw.sum(axis=1)
    gene_positions=counts.index.get_indexer(matched)
    expression=np.log2(1+1e6*raw[:,gene_positions]/libsize[:,None])
    score, high=score_expression(expression)
    assert len(np.unique(score))==73
    assert [int(((~high)&arm).sum()),int(((~high)&(~arm)).sum()),int((high&arm).sum()),int((high&(~arm)).sum())]==[22,15,20,16]
    survival=pd.to_numeric(baseline.Survival, errors="raise").to_numpy(float)
    event=pd.to_numeric(baseline.DIAG_T1D, errors="raise").to_numpy()
    assert set(event)=={0,1} and int(event.sum())==54
    assert np.isfinite(survival).all() and (survival>=0).all()
    assert (survival[event==0]>=156).all()
    treatment_date=pd.to_datetime(baseline.TREATMENTDATE,format="%m/%d/%Y")
    diagnosis_date=pd.to_datetime(baseline.IDDMDATE.where(event==1),format="%m/%d/%Y")
    elapsed=(diagnosis_date-treatment_date).dt.days.to_numpy()/7
    assert np.allclose(elapsed[event==1],survival[event==1],atol=5e-5,rtol=0)
    assert (baseline.IDDMDATE[event==0]==".").all()
    cohort={"metadata_records":len(metadata),"count_gene_rows":counts.shape[0],"count_libraries":counts.shape[1],
        "baseline_n":73,"treated_n":42,"placebo_n":31,"original_randomized_n":76,
        "absent_RNA_n":3,"absent_RNA_status":"unclassified; two treatment and one placebo versus original allocation",
        "metadata_libraries_absent_from_counts":sorted(set(metadata.libraryID)-set(counts.columns)),
        "count_libraries_absent_from_metadata":sorted(set(counts.columns)-set(metadata.libraryID)),
        "baseline_library_unique":True,"participant_identity_basis":"one baseline library per participant per published cohort; no separate participant identifier in public metadata",
        "source_genes":359,"matched_genes":312,"missing_genes":47,"primary_gene_sha256":gene_hash}
    checks={"diagnosis_n":54,"non_diagnosis_n":19,"encoded_early_censor_n_156":0,
        "non_event_min_encoded_duration_weeks":float(survival[event==0].min()),"event_dates_match_elapsed_weeks":True,
        "max_date_duration_rounding_error_weeks":float(abs(elapsed[event==1]-survival[event==1]).max()),
        "all_non_events_have_dot_diagnosis_date":True,"all_73_restricted_outcomes_numerically_defined":True,
        "clinical_complete_followup_independently_verified":False,"unverified_non_event_last_contact_n":19,
        "scope":"Released numerical encoding only; no individual last-contact dates for 19 non-events",
        "all_invariants_pass":True}
    write_json(results/"cohort_flow.json",cohort)
    write_json(results/"endpoint_checks.json",checks)
    # Only after all input and score invariants pass do outcome-by-score contrasts begin.
    y=np.minimum(survival,156)
    primary=estimate(y,arm,high)
    primary.update(CONFIG)
    primary["analysis_status"]="executed"
    primary["cohort_n"]=73
    primary["source_checks_pass"]=True
    primary["score_median"]=float(np.median(score))
    write_json(results/"primary_results.json",primary)
    continuous=continuous_models(y,arm,score,baseline)
    write_json(results/"continuous_models.json",continuous)
    boot_values,boot=bootstrap(expression,y,arm)
    pd.DataFrame(boot_values,columns=["delta_low","delta_high","theta"]).to_csv(results/"bootstrap_draws.tsv",sep="\t",index=False)
    write_json(results/"bootstrap_results.json",boot)
    sensitivity={"horizon_104":estimate(np.minimum(survival,104),arm,high)}
    sf=pd.to_numeric(baseline.sizeFactor).to_numpy(float)
    assert np.isfinite(sf).all() and (sf>0).all()
    sf_score,sf_high=score_expression(np.log2(1+raw[:,gene_positions]/sf[:,None]))
    sensitivity["supplied_sizeFactor"]=estimate(y,arm,sf_high)
    sensitivity["supplied_sizeFactor"]["processing"]="log2(1+count/sizeFactor); upstream size factors may include nonbaseline samples"
    sensitivity["supplied_sizeFactor"]["group_changes_vs_primary"]=int((sf_high!=high).sum())
    positive=(raw>0).all(axis=0)
    gmeans=np.exp(np.log(raw[:,positive]).mean(axis=0))
    ratios=np.median(raw[:,positive]/gmeans,axis=1)
    ratios/=np.exp(np.log(ratios).mean())
    mr_score,mr_high=score_expression(np.log2(1+raw[:,gene_positions]/ratios[:,None]))
    sensitivity["baseline_median_ratio"]=estimate(y,arm,mr_high)
    sensitivity["baseline_median_ratio"]["processing"]="baseline-only median count/geometric-mean ratios, all-positive genes; factors geometric-mean scaled to 1; log2(1+normalized counts)"
    sensitivity["baseline_median_ratio"]["size_factor_gene_n"]=int(positive.sum())
    sensitivity["baseline_median_ratio"]["group_changes_vs_primary"]=int((mr_high!=high).sum())
    loo=[]
    for drop in range(len(y)):
        keep=np.arange(len(y))!=drop
        _,loo_high=score_expression(expression[keep])
        result=estimate(y[keep],arm[keep],loo_high)
        loo.append({"dropped_row_zero_based":drop,"theta_weeks":result["theta"]["estimate_weeks"],"delta_low_weeks":result["delta_low"]["estimate_weeks"],"delta_high_weeks":result["delta_high"]["estimate_weeks"]})
    pd.DataFrame(loo).to_csv(results/"leave_one_out.tsv",sep="\t",index=False)
    loo_values=np.array([v["theta_weeks"] for v in loo])
    sensitivity["leave_one_out"]={"repetitions":73,"theta_min_weeks":float(loo_values.min()),"theta_max_weeks":float(loo_values.max()),"sign_changes_vs_primary":int((np.sign(loo_values)!=np.sign(primary["theta"]["estimate_weeks"])).sum()),"recomputed_score_and_split":True}
    write_json(results/"sensitivities.json",sensitivity)
    # De-identified local derivation file excludes public source dates; never committed.
    pd.DataFrame({"libraryID":baseline.libraryID,"treatment":arm.astype(int),"score":score,"high":high,"restricted_weeks":y,"survival_weeks":survival,"event":event}).to_csv(results/"analysis_dataset.tsv",sep="\t",index=False)
    planning={"design":"four equally sized independent cells, normal approximation, common within-cell SD",
        "formula":"N = 16 * sigma^2 * (z_(.975)+z_(.80))^2 / 26^2, rounded up to multiple of 4",
        "target":"80% power for 26-week interaction versus zero; equivalently true interaction 52 versus margin 26; not post hoc power",
        "does_not_model":["estimated expression threshold","assay measurement error","unequal cell sizes","multiple testing","independent assay transport"],
        "scenarios":[{"assumed_sd_weeks":s,"total_complete_n":int(4*np.ceil((16*s*s*(norm.ppf(.975)+norm.ppf(.8))**2/26**2)/4)),"bounded_range_max_sd":s==78} for s in [40,50,60,78]]}
    write_json(results/"validation_design.json",planning)
    make_figures(results,primary)
    packages=["numpy","pandas","scipy","statsmodels","matplotlib","pdfplumber","pypdf"]
    write_json(root/"environment.json",{"python":sys.version,"executable":sys.executable,"packages":{p:importlib.metadata.version(p) for p in packages},"analysis_script_sha256":sha256(Path(__file__))})
    write_report(root,primary,continuous,boot,sensitivity,planning)
    print(json.dumps({"status":primary["hypothesis_status"],"theta":primary["theta"],"delta_low":primary["delta_low"],"delta_high":primary["delta_high"],"bootstrap":boot,"output":str(results)},indent=2))

def write_report(root,primary,continuous,boot,sensitivity,planning):
    def fmt(row):
        e=row['estimate_weeks']; lo,hi=row['ci95_weeks']
        return f"{e:.2f} weeks (95% CI {lo:.2f} to {hi:.2f})"
    ci=continuous['adjusted_complete']['coefficients']['treatment_x_score_z']
    rows="\n".join(f"| {k.replace('_',' ')} | {v['n']} | {v['mean_weeks']:.2f} | {v['ci95_weeks'][0]:.2f} to {v['ci95_weeks'][1]:.2f} |" for k,v in primary['cells'].items())
    sensrows="\n".join(f"| {name.replace('_',' ')} | {fmt(sensitivity[name]['theta'])} | {sensitivity[name]['hypothesis_status']} |" for name in ['horizon_104','supplied_sizeFactor','baseline_median_ratio'])
    report=f"""# TN10 randomized disease-delay analysis

Measured outcome: category (1), delay of clinical stage 3 onset over 156 weeks.

**Execution complete; magnitude hypothesis {primary['hypothesis_status']} (analysis status: {primary['card_verdict']}).** This is an exploratory secondary analysis of the 73 participants with published baseline RNA, with randomized treatment contrasts and a fixed molecular score.

## Frozen question and reconstruction

The prespecified analysis plan specified a fixed AbATE-derived gene list, exact-symbol projection, pooled baseline log-CPM gene standardization, median split, complete 156-week restricted event-free outcome, and the low-minus-high difference in randomized treatment benefit. Its threshold was 26 weeks. The three retrieved inputs match the prespecified byte lengths and SHA-256 hashes exactly. All 73 baseline libraries join uniquely; allocation is 42 teplizumab and 31 placebo; 312 of 359 signature symbols match. The sorted primary list hash matches the prespecified value. Fifty-four participants eventually have recorded diagnosis; 19 are non-events. No non-event is encoded as censored before week 156; individual last-contact dates are unavailable for the 19 non-events, so the estimates remain conditional on the released encoding. Date-derived durations match recorded event durations up to rounding. The 196th count column has no corresponding metadata; it is not used. Gene selection and score construction use no TN10 outcomes.

The original trial randomized 76 people (44/32); the three absent RNA participants remain unclassified. The metadata contain library identity rather than an independent participant identifier; the one-baseline-library-per-person interpretation rests additionally on the published cohort description. Randomization blocks are unavailable. Causal interpretation in the RNA subset assumes baseline-sample inclusion has not induced material treatment-outcome selection.

## Primary results

| Cell | n | Mean restricted weeks | 95% t interval |
|---|---:|---:|---|
{rows}

Low-score treatment-added delay: **{fmt(primary['delta_low'])}**. High-score treatment-added delay: **{fmt(primary['delta_high'])}**.

The interaction, low benefit minus high benefit, is **{fmt(primary['theta'])}**; Welch-Satterthwaite df {primary['theta']['df']:.2f}, SE {primary['theta']['se_weeks']:.2f} weeks, two-sided test against zero p={primary['theta']['p_two_sided_zero']:.4f}. The interval {'crosses' if primary['interval_crosses_zero'] else 'does not cross'} zero. It spans the prespecified 26-week margin, so neither the proposed magnitude nor equal subgroup benefit is established. This verdict follows the original interval rule rather than post-result hypothesis selection.

Placebo prognosis, low-score minus high-score restricted time, is {fmt(primary['placebo_low_minus_high'])}. Prognostic separation and benefit modification are distinct estimands.

![Four cell means](four_cell_outcomes.png)

![Treatment-added delay and interaction](forest_delay.png)

## Required sensitivities

| Analysis | Interaction estimate and 95% interval | Magnitude status |
|---|---|---|
{sensrows}

The arm-stratified participant bootstrap performed all {boot['repetitions']:,} prespecified draws, with {boot['unusable']} unusable draws. Each draw recomputed gene means/SDs, scores, median, and groups. Its theta percentile interval is {boot['percentile_ci95_weeks']['theta'][0]:.2f} to {boot['percentile_ci95_weeks']['theta'][1]:.2f} weeks. The negative-theta fraction is {boot['theta_negative_fraction']:.3f}; this is a resampling diagnostic, not a posterior probability. Complete participant-level draws are saved in the configured analysis output directory.

The adjusted continuous-score model omitted {continuous['omitted_covariates_n']} participant for a missing baseline covariate (n={continuous['adjusted_complete']['n']}). Its treatment-by-score coefficient is {ci['estimate']:.2f} weeks per pooled composite-score SD, HC3/t 95% CI {ci['ci95'][0]:.2f} to {ci['ci95'][1]:.2f}, p={ci['p']:.4f}. A negative coefficient corresponds to less treatment benefit with higher expression; it is not numerically the dichotomized theta. Full and complete-case unadjusted continuous models and the complete-case four-cell contrast are retained in continuous_models.json.

Supplied size factors may have been learned using later samples. The baseline-only median-ratio analysis instead used {sensitivity['baseline_median_ratio']['size_factor_gene_n']} genes positive in every baseline library. Primary results remain the log-CPM result regardless of sensitivity agreement. Normalization changed {sensitivity['supplied_sizeFactor']['group_changes_vs_primary']} and {sensitivity['baseline_median_ratio']['group_changes_vs_primary']} group labels respectively.

Leave-one-out resampling recomputed the score/split all 73 times: theta range {sensitivity['leave_one_out']['theta_min_weeks']:.2f} to {sensitivity['leave_one_out']['theta_max_weeks']:.2f} weeks, with {sensitivity['leave_one_out']['sign_changes_vs_primary']} sign changes versus the primary estimate. Valid influential observations were retained.

## Precision and validation planning

The interval is the measured precision limit; lack of significance is not evidence of equal effects. For a hypothetical four-cell balanced validation study, detecting a true 26-week interaction versus zero with 80% power at two-sided 5% significance gives the following normal-approximation complete-sample totals. These same totals apply to a true 52-week interaction versus a 26-week null margin. They do not mean that a true effect exactly at 26 can be supported above 26 with 80% power. Assumed within-cell SD is deliberately varied, not selected from the observed effect.

| Assumed within-cell SD, weeks | Complete participants |
|---|---:|
"""
    report += "\n".join(f"| {r['assumed_sd_weeks']} | {r['total_complete_n']} |" for r in planning['scenarios'])
    report += """

These planning calculations omit threshold-estimation error, assay transport, imbalance, missing data and multiplicity. An independent validation protocol requires a transferable assay and fixed threshold; the pooled-cohort median used here is not such an assay.

## Interpretation and boundaries

This computation recovers the clinical-benefit comparison omitted by a treated-outcome classifier. The primary dichotomized interval supports a positive direction within this analysis, but it does not resolve the six-month magnitude; the adjusted continuous interaction remains imprecise across zero. The published molecular observations remain compatible with these results. Small cells, selection into the RNA subset, prior published use of TN10, bulk-cell mixture, the 47 missing genes, cohort-relative scoring, and lack of individualized safety data limit conclusions. The outcome scale is additional stage-3-free time within the 156-week horizon.

## Reproduction and provenance

From the repository, run `python -B src/reproduce_tn10.py --inputs-root INPUT_DIRECTORY --output-root NEW_OUTPUT_DIRECTORY` in the pinned Python environment. The three original files are read-only and verified by exact SHA-256 and byte count; this complete runner performs no downloads. It generates primary and sensitivity analyses, all 10,000 score-reconstruction bootstrap draws, independent arithmetic checks, outcome curves, ascertainment bounds, and label-instability diagnostics. Commands, timings, source hashes and environment versions are recorded under the chosen output directory. See `docs/tn10_reproduction.md`. Raw participant files are not committed.

Dataset: Lemaitre, Mathieu, Gysemans and Type 1 Diabetes TrialNet, [Zenodo v1 / DOI 10.5281/zenodo.17063348](https://zenodo.org/records/17063348), CC BY-SA 4.0 as verified from the archive API. Preserve attribution/share-alike for redistributed dataset adaptations. Gene source: Sassi et al., [JCI 2025;135:e176403](https://www.jci.org/articles/view/176403), [Supplemental Table 6](https://www.jci.org/articles/view/176403/sd/6), CC BY 4.0 article license. Trial interpretation: Herold et al., [NEJM 2019, DOI 10.1056/NEJMoa1902226](https://doi.org/10.1056/NEJMoa1902226). Prior methodological boundary: Stene, [Diabetologia 2025, DOI 10.1007/s00125-025-06512-5](https://doi.org/10.1007/s00125-025-06512-5). Related empirical prior art: Galderisi et al., [Diabetes 2026, DOI 10.2337/db25-0801](https://doi.org/10.2337/db25-0801), uses baseline secretion, subgroup restricted survival and interaction analysis for abatacept in stage 1 T1D. These studies establish the preceding methods and biological context. See literature_audit.md for the source-specific comparison.
"""
    (root/"results"/"scientific_report.md").write_text(report,encoding="utf-8")
    (root/"results"/"literature_audit.md").write_text(LITERATURE_AUDIT,encoding="utf-8")


LITERATURE_AUDIT = """# TN10 source and novelty audit

Date: 2026-10-22, America/New_York. This source update records the inspected literature and executed queries.

## Sources inspected

- [Zenodo record 17063348 v1](https://zenodo.org/records/17063348) and its API supplied open download locations, authorship, publication/version identity, and CC BY-SA 4.0 metadata. Both CSVs were downloaded and matched the protocol SHA-256 and exact byte count. No controlled-access boundary was encountered for these inputs. The public record's description says 1:1 randomization, while actual published/randomized allocation is 44/32 and RNA allocation 42/31; allocation counts come from the trial and actual metadata, not an inference from the record's prose.
- [Sassi et al. JCI article](https://www.jci.org/articles/view/176403) was fetched in full as HTML. Methods explicitly identify 73 baseline individuals (42 teplizumab/31 placebo), a 359-gene AbATE baseline negative-correlation list, and placebo-associated gene filtering for the separate TN10 list. The article/supplement is CC BY 4.0. Table 6 was parsed across all 11 pages; the first page was rendered and visually checked against the symbol/correlation/P-value columns. The 359 parsed symbols and sorted exact 312-gene projection match the frozen protocol hash. The 26-gene combined classifier was not used.
- [Herold et al., NEJM 2019](https://doi.org/10.1056/NEJMoa1902226) is the randomized-trial source; its overall clinical delay result is established prior knowledge. Later follow-up in this computation is sourced to the public molecular release, not represented as identical to the original report's cutoff.
- [Stene, Diabetologia 2025](https://doi.org/10.1007/s00125-025-06512-5) discusses risk heterogeneity versus treatment-effect heterogeneity, small-trial precision and multiplicity. The conceptual distinction and interaction methodology are prior art.
- [Galderisi et al., Diabetes 2026](https://doi.org/10.2337/db25-0801), [primary institutional PDF](https://scholarworks.indianapolis.iu.edu/server/api/core/bitstreams/48340b3b-1983-4731-ac97-e3e2c914d471/content), studies baseline insulin secretion, subgroup restricted mean progression-free survival, and treatment interaction in stage 1 abatacept. It directly precludes claiming novelty for biomarker-stratified time gained or randomized interaction testing in T1D generally. It differs in intervention, stage, assay and marker from this exact TN10 fixed-list contrast.
- [Suomi et al., EBioMedicine 2026](https://doi.org/10.1016/j.ebiom.2026.106374), [institutional source](https://www.repository.cam.ac.uk/items/ba648822-382f-4936-8c6a-c80970bf2b34), was identified as adjacent prognostic transcriptomics in newly diagnosed T1D. The institutional full-page fetch timed out; only source metadata were verified in this update. It is not used to assert a particular numerical or causal finding here.

## Search log and inference

Queries executed included `teplizumab neutrophil signature restricted mean treatment interaction Sassi TN10 2026`, `teplizumab restricted mean gene`, `Sassi teplizumab interaction`, `teplizumab AbATE treatment effect signature`, `TN10 restricted signature`, the exact Stene title, the exact Suomi title, and the exact Galderisi title with `restricted`.

The inspected primary-source studies did not report this fixed AbATE-only 312-gene TN10 156-week absolute treatment-benefit contrast. Its contribution is the specified randomized four-cell estimate and its score-reconstruction uncertainty, within the combined research program. The executed queries and inspected sources define this search record.

"""

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--data-root", type=Path, required=True, help="Generated analysis directory")
    parser.add_argument("--raw-root", type=Path, help="Read-only directory containing the three original inputs")
    parser.add_argument("--download", action="store_true")
    parser.add_argument("--verify-only", action="store_true")
    args = parser.parse_args()
    root = args.data_root.resolve()
    root.mkdir(parents=True, exist_ok=True)
    os.environ["MPLCONFIGDIR"] = str(root / "cache" / "matplotlib")
    write_json(root / "frozen_config.json", CONFIG)
    raw = retrieve(root, args.download, args.raw_root)
    import pandas as pd
    metadata = pd.read_csv(raw / "TN10_metadata.csv")
    counts = pd.read_csv(raw / "TN10_raw_counts.csv", index_col=0)
    genes, table_text = parse_table(raw / "JCI176403_Table6.pdf")
    (root / "table6_extracted.txt").write_text(table_text, encoding="utf-8")
    print(json.dumps({"metadata_shape": metadata.shape, "columns": metadata.columns.tolist(), "counts_shape": counts.shape, "table_rows": len(genes), "table_head": genes.head().to_dict("records")}))
    if args.verify_only:
        return
    run_analysis(root, metadata, counts, genes)

if __name__ == "__main__":
    main()
