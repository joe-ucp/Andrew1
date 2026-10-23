"""Reproduce TN10 companion tables/PNG figures without modifying original outputs."""
from __future__ import annotations
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import sys
sys.dont_write_bytecode = True
import numpy as np
import pandas as pd
from card002_run import SOURCES, score_expression, estimate

CELLS = [(False, 1, "Low / teplizumab", "low_treated"),
         (False, 0, "Low / placebo", "low_placebo"),
         (True, 1, "High / teplizumab", "high_treated"),
         (True, 0, "High / placebo", "high_placebo")]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path, value):
    Path(path).write_text(json.dumps(value, indent=2, allow_nan=False) + "\n", encoding="utf-8")


def curve(cell, horizon=156.):
    """Exact empirical survival conditional on the supplied no-early-censor encoding."""
    event_times = cell.loc[(cell.event == 1) & (cell.survival_weeks <= horizon), "survival_weeks"].to_numpy()
    times = np.unique(np.r_[0., event_times, horizon])
    survival = np.array([1 - np.sum(event_times <= t) / len(cell) for t in times])
    area = float(np.dot(np.diff(times), survival[:-1]))
    return times, survival, area


def point_interval(ax, estimate_value, bounds, y, color, label):
    ax.errorbar(estimate_value, y, xerr=[[estimate_value-bounds[0]], [bounds[1]-estimate_value]],
                fmt="o", color=color, capsize=4, markersize=5, linewidth=1.6)
    ax.text(.99, y, label, transform=ax.get_yaxis_transform(), ha="right", va="center", fontsize=9,
            bbox={"facecolor":"white", "edgecolor":"none", "pad":2, "alpha":.9})


def figures(output, data, primary, boot, continuous, identification_bounds):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":10, "axes.spines.top":False, "axes.spines.right":False,
                         "savefig.dpi":200, "font.family":"DejaVu Sans"})
    colors = ["#145f86", "#145f86", "#ad4a23", "#ad4a23"]
    styles = ["-", "--", "-", "--"]
    fig, (ax, tableax) = plt.subplots(2, 1, figsize=(9.5, 6.4), gridspec_kw={"height_ratios":[4.2,1.7]})
    ticktimes = [0,26,52,78,104,130,156]
    riskrows, curve_rows = [], []
    for i, (high, arm, label, key) in enumerate(CELLS):
        cell = data[(data.high == high) & (data.treatment == arm)]
        t, s, area = curve(cell)
        assert np.isclose(area, primary["cells"][key]["mean_weeks"], rtol=0, atol=1e-10)
        ax.step(t, s, where="post", color=colors[i], linestyle=styles[i], linewidth=2, label=f"{label} (n={len(cell)})")
        for tt, ss in zip(t, s):
            curve_rows.append(dict(group=key, weeks=float(tt), stage3_free=float(ss), restricted_mean_weeks=area))
        for tt in ticktimes:
            riskrows.append(dict(group=key, label=label, week=tt,
                                 encoded_at_risk=int((cell.survival_weeks >= tt).sum())))
    ax.set(xlim=(0,156), ylim=(0,1.03), xticks=ticktimes,
           ylabel="Proportion without recorded stage 3 diagnosis", xlabel="Weeks from treatment initiation")
    ax.legend(loc="lower left", fontsize=9, frameon=False)
    ax.set_title("Four randomized-arm / baseline-score outcome curves", loc="left", fontsize=12, weight="bold")
    ax.grid(axis="y", alpha=.2)
    tableax.set_axis_off()
    celltext = [[label]+[str(v["encoded_at_risk"]) for v in riskrows if v["group"] == key] for _,_,label,key in CELLS]
    table = tableax.table(cellText=celltext, colLabels=["Group / week"]+[str(t) for t in ticktimes],
                          colWidths=[.27]+[.73/7]*7, loc="center", cellLoc="center", bbox=[0,.02,1,.90])
    table.auto_set_font_size(False); table.set_fontsize(9)
    for (row,col), cell in table.get_celld().items():
        cell.set_edgecolor("#dddddd"); cell.set_linewidth(.35)
        if row == 0: cell.set_facecolor("#f0f3f5")
    tableax.text(0,1.08,"Encoded number at risk immediately before each week", transform=tableax.transAxes, fontsize=10)
    fig.subplots_adjust(left=.11, right=.98, top=.92, bottom=.13, hspace=.36)
    fig.text(.11,.035,"Conditional on released durations: no early censoring is encoded. Individual last-contact dates are unavailable.\nAll 19 non-events share a common 551.4286-week duration; these are not independently verified clinical risk sets.",fontsize=9)
    fig.savefig(output / "tn10_outcome_curves.png"); plt.close(fig)
    pd.DataFrame(riskrows).to_csv(output / "tn10_curve_risk_counts.csv", index=False)
    pd.DataFrame(curve_rows).to_csv(output / "tn10_curve_steps.csv", index=False)

    fig, axes = plt.subplots(3, 1, figsize=(10.5, 9.5), gridspec_kw={"height_ratios":[1,1.5,1.25]})
    a,b,c = axes
    for y,key in [(1,"delta_low"),(0,"delta_high")]:
        row=primary[key]; bounds=row["ci95_weeks"]; value=row["estimate_weeks"]
        point_interval(a,value,bounds,y,"#145f86",f"{value:.2f} ({bounds[0]:.2f}, {bounds[1]:.2f})")
    a.set(yticks=[1,0],yticklabels=["Low-score benefit", "High-score benefit"],ylim=(-.6,1.7),xlim=(-55,190),
          xlabel="Treatment minus placebo: restricted weeks")
    a.set_title("A  Subgroup treatment effects — fixed-group Welch intervals",loc="left",fontsize=11,weight="bold")
    theta=primary["theta"]["estimate_weeks"]
    for y,interval,label in [(2,primary["theta"]["ci95_weeks"],"Fixed-group Welch"),
                            (1,boot["percentile_ci95_weeks"]["theta"],"Full-score bootstrap")]:
        point_interval(b,theta,interval,y,"#145f86",f"{theta:.2f} ({interval[0]:.2f}, {interval[1]:.2f})")
    b.plot(identification_bounds["theta_range_weeks"],[0,0],color="#555555",linewidth=2,marker="|",markersize=10)
    b.text(.99,0,f"Range {identification_bounds['theta_range_weeks'][0]:.2f} to {identification_bounds['theta_range_weeks'][1]:.2f}",
           transform=b.get_yaxis_transform(),ha="right",va="center",fontsize=9,bbox={"facecolor":"white","edgecolor":"none","pad":2})
    b.set(yticks=[2,1,0],yticklabels=["Fixed-group Welch", "Full-score bootstrap", "Unknown follow-up range"],ylim=(-.6,2.7),xlim=(-55,190),
          xlabel="Low-score benefit minus high-score benefit: restricted weeks")
    b.axvline(26,color="#ad4a23",linestyle="--",linewidth=1.3)
    b.text(28,2.42,"26-week interaction margin",fontsize=9,color="#ad4a23")
    b.set_title("B  Interaction — sampling intervals and follow-up identification range",loc="left",fontsize=11,weight="bold")
    for y,key in [(2,"unadjusted_full"),(1,"unadjusted_complete"),(0,"adjusted_complete")]:
        row=continuous[key]["coefficients"]["treatment_x_score_z"]
        point_interval(c,row["estimate"],row["ci95"],y,"#ad4a23",f"{row['estimate']:.2f} ({row['ci95'][0]:.2f}, {row['ci95'][1]:.2f})")
    c.set(yticks=[2,1,0],yticklabels=["Unadjusted, n=73","Unadjusted, n=72","Age/HbA1c adjusted, n=72"],
          ylim=(-.6,2.7),xlim=(-65,85),xlabel="Treatment × score interaction: weeks per pooled score SD")
    c.set_title("C  Continuous-score interactions — HC3 / t intervals",loc="left",fontsize=11,weight="bold")
    for ax in axes:
        ax.axvline(0,color="#666666",linewidth=.8)
        ax.grid(axis="x",alpha=.15)
    fig.subplots_adjust(left=.26,right=.985,top=.95,bottom=.13,hspace=.66)
    fig.text(.26,.025,"Blue/orange points: existing estimates with 95% sampling intervals.\nGray range: endpoint-identification sensitivity, not a confidence interval.\nThe 26-week margin applies only to panel B. Panel C has a different estimand and scale.",fontsize=9)
    fig.savefig(output / "tn10_effect_uncertainty.png");plt.close(fig)


def main(original, output, raw_root=None):
    original,output=Path(original).resolve(),Path(output).resolve()
    raw_root = Path(raw_root).resolve() if raw_root is not None else original / "raw"
    if original==output or original in output.parents or raw_root==output or raw_root in output.parents:
        raise ValueError("Revision output must be separate from its read-only inputs")
    output.mkdir(parents=True,exist_ok=True)
    filenames=[raw_root/name for name in SOURCES]
    filenames += [original/"results"/name for name in ["analysis_dataset.tsv","primary_results.json", "bootstrap_results.json",
        "bootstrap_draws.tsv","continuous_models.json","sensitivities.json","endpoint_checks.json","primary_genes.txt"]]
    before={str(p):digest(p) for p in filenames}
    for name,(_,expected,size) in SOURCES.items():
        p=raw_root/name
        assert digest(p)==expected and p.stat().st_size==size
    load=lambda name:json.loads((original/"results"/name).read_text())
    primary,boot,continuous,sens=[load(n) for n in ["primary_results.json","bootstrap_results.json","continuous_models.json","sensitivities.json"]]
    data=pd.read_csv(original/"results"/"analysis_dataset.tsv",sep="\t")
    metadata=pd.read_csv(raw_root/"TN10_metadata.csv")
    baseline=metadata[metadata["Time.Point"]=="Baseline"].reset_index(drop=True)
    assert len(data)==73 and data.libraryID.is_unique and baseline.libraryID.tolist()==data.libraryID.tolist()
    assert np.array_equal(data.treatment.to_numpy(),baseline.TREATMENTDESC.eq("Teplizumab").astype(int))
    assert np.array_equal(data.event.to_numpy(),pd.to_numeric(baseline.DIAG_T1D))
    assert np.allclose(data.survival_weeks,pd.to_numeric(baseline.Survival),rtol=0,atol=1e-10)
    counts=pd.read_csv(raw_root/"TN10_raw_counts.csv",index_col=0)
    genes=(original/"results"/"primary_genes.txt").read_text().splitlines()
    raw=counts.loc[:,baseline.libraryID].to_numpy(float).T
    expr=np.log2(1+1e6*raw[:,counts.index.get_indexer(genes)]/raw.sum(axis=1)[:,None])
    scores,high=score_expression(expr)
    assert np.allclose(scores,data.score,rtol=0,atol=1e-12) and np.array_equal(high,data.high)
    y=np.minimum(data.survival_weeks.to_numpy(),156)
    assert np.allclose(y,data.restricted_weeks,rtol=0,atol=1e-10)
    reproduced=estimate(y,data.treatment.to_numpy(bool),high)
    for key in ["theta","delta_low","delta_high","placebo_low_minus_high"]:
        assert np.isclose(reproduced[key]["estimate_weeks"],primary[key]["estimate_weeks"],atol=1e-10,rtol=0)
        assert np.allclose(reproduced[key]["ci95_weeks"],primary[key]["ci95_weeks"],atol=1e-10,rtol=0)
    draws=pd.read_csv(original/"results"/"bootstrap_draws.tsv",sep="\t")
    assert len(draws)==10000 and np.isfinite(draws.to_numpy()).all()
    assert np.allclose(np.percentile(draws.theta,[2.5,97.5]),boot["percentile_ci95_weeks"]["theta"],atol=1e-10,rtol=0)
    non_event=data[data.event==0]
    endpoint=dict(n=73,events_by_156=int(((data.event==1)&(data.survival_weeks<=156)).sum()),
        recorded_events_after_156=int(((data.event==1)&(data.survival_weeks>156)).sum()),
        non_event_n=len(non_event),encoded_early_censor_n=int((non_event.survival_weeks<156).sum()),
        unique_non_event_encoded_durations=sorted(non_event.survival_weeks.unique().tolist()),
        maximum_event_duration=float(data.loc[data.event==1,"survival_weeks"].max()),
        clinical_complete_followup_independently_verified=False,
        limitation="All19 non-events share551.4286 weeks, equal to maximum event duration rounded; no individual last-contact field. Encoded completeness does not verify actual clinical ascertainment.",
        numerical_result_status="unchanged conditional reproduction",original_assertion_narrowed="All 73 restricted outcomes are numerically defined by the encoding; this does not independently verify clinical ascertainment")
    assert endpoint["events_by_156"]==32 and endpoint["recorded_events_after_156"]==22 and endpoint["non_event_n"]==19
    dump(output/"endpoint_provenance_revision.json",endpoint)
    # Finite-sample linear contrast extrema: retain all recorded diagnosis times,
    # vary only the 19 non-event records with unverifiable follow-up over [0,156].
    # No censoring model or population sampling uncertainty is imposed here.
    extreme_y_low=y.copy();extreme_y_high=y.copy();bound_rows=[]
    for h,a,label,key in CELLS:
        mask=(data.high.to_numpy()==h)&(data.treatment.to_numpy()==a)
        uncertain=mask&(data.event.to_numpy()==0)
        sign=1 if key in ("low_treated","high_placebo") else -1
        extreme_y_low[uncertain]=0 if sign>0 else 156
        extreme_y_high[uncertain]=156 if sign>0 else 0
        fixed_sum=float(y[mask&~uncertain].sum());n=int(mask.sum());u=int(uncertain.sum())
        bound_rows.append(dict(group=key,contrast_sign=sign,n=n,unverified_non_events=u,
            cell_mean_min=fixed_sum/n,cell_mean_max=(fixed_sum+156*u)/n))
    contrast=lambda v:sum(r["contrast_sign"]*v[(data.high.to_numpy()==h)&(data.treatment.to_numpy()==a)].mean()
                          for r,(h,a,_,_) in zip(bound_rows,CELLS))
    bound_values=[float(contrast(extreme_y_low)),float(contrast(extreme_y_high))]
    assert bound_values[0] <= primary["theta"]["estimate_weeks"] <= bound_values[1]
    bound_result=dict(theta_range_weeks=bound_values,uncertain_endpoint_n=19,
        uncertain_endpoint_range_weeks=[0,156],recorded_event_durations_held_fixed=54,
        interpretation="Conservative finite-sample contrast identification range, not confidence interval or population treatment-effect bound",
        assumption="No additional follow-up information imposed on19 non-event records; later RNA sampling not used to tighten this deliberately conservative box",
        clinical_hypothesis_new=False,post_results_provenance_sensitivity=True,cells=bound_rows)
    dump(output/"followup_identification_range.json",bound_result)
    pd.DataFrame(bound_rows).to_csv(output/"followup_identification_cells.csv",index=False)
    cellrows=[]
    for h,a,label,key in CELLS:
        x=data[(data.high==h)&(data.treatment==a)];r=primary["cells"][key]
        cellrows.append(dict(group=key,label=label,n=len(x),events_by_156=int(((x.event==1)&(x.survival_weeks<=156)).sum()),
            at_ceiling_156=int((x.restricted_weeks==156).sum()),non_event_duration_unverified=int((x.event==0).sum()),
            restricted_mean_weeks=r["mean_weeks"],lower95=r["ci95_weeks"][0],upper95=r["ci95_weeks"][1]))
    pd.DataFrame(cellrows).to_csv(output/"four_group_outcomes.csv",index=False)
    rows=[]
    for key,label in [("unadjusted_full","Unadjusted full"),("unadjusted_complete","Unadjusted complete case"),("adjusted_complete","Age and HbA1c adjusted complete case")]:
        m=continuous[key];r=m["coefficients"]["treatment_x_score_z"]
        rows.append(dict(model=label,n=m["n"],estimate_weeks_per_score_sd=r["estimate"],lower95=r["ci95"][0],upper95=r["ci95"][1],p_zero=r["p"],interval="HC3 t"))
    pd.DataFrame(rows).to_csv(output/"continuous_interactions.csv",index=False)
    cc=continuous["dichotomized_complete_fixed_primary_split"]
    dump(output/"complete_case_four_cell.json",cc)
    rows=[]
    for name,n,horizon,r,changes in [("Primary CPM",73,156,primary,0),("Complete case, primary split",72,156,cc,0),
        ("CPM,104-week horizon",73,104,sens["horizon_104"],0),
        ("Supplied size factors",73,156,sens["supplied_sizeFactor"],sens["supplied_sizeFactor"]["group_changes_vs_primary"]),
        ("Baseline median ratios",73,156,sens["baseline_median_ratio"],sens["baseline_median_ratio"]["group_changes_vs_primary"])]:
        theta=r["theta"]
        rows.append(dict(analysis=name,n=n,horizon_weeks=horizon,theta_weeks=theta["estimate_weeks"],lower95=theta["ci95_weeks"][0],upper95=theta["ci95_weeks"][1],changed_group_labels=changes))
    pd.DataFrame(rows).to_csv(output/"normalization_horizon_completecase.csv",index=False)
    figures(output,data,primary,boot,continuous,bound_result)
    after={str(p):digest(p) for p in filenames}
    assert before==after,"Original input changed during revision"
    result=dict(completed_utc=datetime.now(timezone.utc).isoformat(),status="complete_conditional_reproduction",
        original_files_unchanged=True,original_input_sha256=before,source_sha256=digest(__file__),
        exact_score_and_group_reconstruction=True,primary_effects_and_intervals_reproduced=True,
        bootstrap_percentiles_reconciled=True,curve_areas_reconcile_to_original_means=True,
        clinical_followup_provenance_status="unverified19non_events",primary_magnitude_verdict_unchanged="inconclusive",
        output_sha256={p.name:digest(p) for p in output.iterdir() if p.is_file() and p.name!="revision_manifest.json"})
    dump(output/"revision_manifest.json",result)
    print(json.dumps({k:v for k,v in result.items() if not k.endswith("sha256")},indent=2))


if __name__=="__main__":
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--original-root",type=Path,required=True)
    ap.add_argument("--output-root",type=Path,required=True)
    ap.add_argument("--raw-root",type=Path,help="Read-only original input directory")
    args=ap.parse_args();main(args.original_root,args.output_root,args.raw_root)
