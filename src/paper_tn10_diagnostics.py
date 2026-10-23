"""Portable follow-up and saved-bootstrap diagnostics; no new clinical fits.

Provide the research data directory with --data-root. Original and revision-2
inputs can be relocated independently with --original-root and --revision2-root.
All new artifacts are written beneath --output-root; inputs are read-only.
"""
from __future__ import annotations
import argparse,hashlib,json,os,shutil
from datetime import datetime,timezone
from pathlib import Path
import numpy as np
import pandas as pd
from card002_run import score_expression

SEED=20261022
N_DRAWS=10000
HORIZON=156.
CELLS=[(False,True,"low_treated",1),(False,False,"low_placebo",-1),(True,True,"high_treated",-1),(True,False,"high_placebo",1)]


def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def dump(path,data):
    def safe(v):
        if isinstance(v,dict):return {str(k):safe(x) for k,x in v.items()}
        if isinstance(v,(tuple,list,np.ndarray)):return [safe(x) for x in v]
        if isinstance(v,np.generic):return safe(v.item())
        if isinstance(v,float) and not np.isfinite(v):return None
        return v
    Path(path).write_text(json.dumps(safe(data),indent=2,allow_nan=False)+"\n",encoding="utf-8")


def protocol():
    return dict(status="POST_RESULTS_DESCRIPTIVE_DIAGNOSTICS",scope="Exactly two descriptive diagnostics: follow-up identification and saved-bootstrap distribution/assignment stability; no new clinical models, no replacement of original estimates or inference.",
        followup=dict(horizon_weeks=156,uncertain_records=19,fixed_recorded_diagnoses=54,
            assumption="All 19 non-event records have stage-3-free restricted outcomes in [L,156], for a hypothetical common documented minimum L between 0 and 156. No value of L is established by this calculation.",
            formulas="lower=theta0-(156-L)*(9/22+2/16); upper=theta0+(156-L)*(2/15+6/20)",
            evaluations="Integer weeks 0..156 plus exact lower-bound crossings of 0 and 26 weeks; report algebraic crossings and customary 26-week grid values.",
            interpretation="Finite-sample identification profile with fixed groups and recorded diagnosis endpoints. Not a confidence interval, population effect bound, or satisfaction of the original sampling-interval decision rule. Later RNA samples supply no follow-up information here."),
        bootstrap=dict(source="Original saved 10000 full-score bootstrap draws",seed=SEED,repetitions=N_DRAWS,
            distribution="Mean, median, sample SD, adjusted sample skewness, minimum, maximum and percentiles 2.5,5,25,50,75,95,97.5. Compare centers and interval widths with original point and fixed-group interval.",
            replay="Saved draws contain effects only. Reproduce original placebo-then-treated stratified sampling and complete gene standardization/score/pooled-median construction solely to recover group labels. Require all three original effect columns to agree draw by draw within 1e-10 weeks.",
            draw_weighted_flips="Each of the 73 sampled slots per draw contributes once; repeated source participants count repeatedly. Compare reconstructed label with that participant's original full-sample label.",
            unique_inbag_flips="Each distinct source participant appearing in a draw contributes once; absent participants are not classified. Pooled rate divides summed participant-draw flips by summed distinct participant-draw appearances.",
            participant_flips="For each participant, flips divided by the number of draws in which that participant appears at least once. Repeated copies must have identical labels.",
            interpretation="Label changes and distribution displacement are descriptive properties of the existing bootstrap. They do not isolate a causal contribution of relabeling from sample composition or establish superiority of one interval method."),
        outputs="Aggregate summaries, full retained label/index arrays in the output directory, one supplementary PNG, concise report, original-input and output hashes.")


def freeze(out,original,revision2):
    out=Path(out).resolve()
    for source_root in (original,revision2):
        source_root=Path(source_root).resolve()
        if out==source_root or source_root in out.parents:raise ValueError("Output must be separate from the preserved input directories")
    out.mkdir(parents=True,exist_ok=True);path=out/"diagnostic_protocol.json";spec=protocol()
    if path.exists():
        previous=json.loads(path.read_text());previous.pop("frozen_utc",None)
        if previous!=spec:raise ValueError("Protocol differs from prior frozen version")
    else:dump(path,dict(spec,frozen_utc=datetime.now(timezone.utc).isoformat()))
    snapshot=out/"source_snapshot"/sha(__file__);snapshot.mkdir(parents=True,exist_ok=True)
    copy=snapshot/Path(__file__).name
    if copy.exists():assert sha(copy)==sha(__file__)
    else:shutil.copy2(__file__,copy)


def followup_profile(data,primary,out,revision2):
    rows=[];y=data.restricted_weeks.to_numpy();high=data.high.to_numpy(bool);arm=data.treatment.to_numpy(bool);event=data.event.to_numpy()
    theta=float(primary["theta"]["estimate_weeks"])
    for h,a,name,sign in CELLS:
        mask=(high==h)&(arm==a);uncertain=mask&(event==0)
        rows.append(dict(group=name,sign=sign,n=int(mask.sum()),uncertain=int(uncertain.sum()),fixed_sum=float(y[mask&~uncertain].sum())))
    assert [r["n"] for r in rows]==[22,15,20,16] and [r["uncertain"] for r in rows]==[9,2,6,2]
    positive=sum(r["uncertain"]/r["n"] for r in rows if r["sign"]>0)
    negative=sum(r["uncertain"]/r["n"] for r in rows if r["sign"]<0)
    cross0=HORIZON-theta/positive;cross26=HORIZON-(theta-26)/positive
    levels=np.unique(np.r_[np.arange(157,dtype=float),cross0,cross26]);records=[]
    for level in levels:
        lower=theta-(HORIZON-level)*positive;upper=theta+(HORIZON-level)*negative
        # Independent cell-wise extrema, avoiding algebraic reuse.
        direct_lower=sum(r["sign"]*(r["fixed_sum"]+r["uncertain"]*(level if r["sign"]>0 else HORIZON))/r["n"] for r in rows)
        direct_upper=sum(r["sign"]*(r["fixed_sum"]+r["uncertain"]*(HORIZON if r["sign"]>0 else level))/r["n"] for r in rows)
        np.testing.assert_allclose([lower,upper],[direct_lower,direct_upper],rtol=0,atol=1e-10)
        records.append(dict(assumed_minimum_stage3_free_weeks=level,theta_lower_weeks=lower,theta_upper_weeks=upper,range_width_weeks=upper-lower))
    frame=pd.DataFrame(records);frame.to_csv(out/"followup_profile.csv",index=False)
    old=json.loads((revision2/"followup_identification_range.json").read_text())
    np.testing.assert_allclose(frame.iloc[0][["theta_lower_weeks","theta_upper_weeks"]],old["theta_range_weeks"],rtol=0,atol=1e-10)
    np.testing.assert_allclose(frame.iloc[-1][["theta_lower_weeks","theta_upper_weeks"]],[theta,theta],rtol=0,atol=1e-10)
    result=dict(theta_original_weeks=theta,positive_uncertain_weight=positive,negative_uncertain_weight=negative,
        lower_crosses_zero_at_L_weeks=cross0,lower_crosses_26_at_L_weeks=cross26,
        inequality="Strict lower>0 requires L strictly greater than first crossing; strict lower>26 requires L strictly greater than second crossing. Equality at the crossing is not strict support.",
        lower_formula=f"{theta:.14g} - (156 - L) * {positive:.14g}",upper_formula=f"{theta:.14g} + (156 - L) * {negative:.14g}",
        clinical_minimum_followup_verified=False,sampling_uncertainty_included=False,original_decision_rule_satisfied_by_profile=False,
        cells=rows,interpretation=protocol()["followup"]["interpretation"])
    dump(out/"followup_profile_summary.json",result)
    return frame,result


def distribution(draws,primary,out):
    from scipy.stats import skew
    values=draws.theta.to_numpy();original=float(primary["theta"]["estimate_weeks"])
    levels=[0,.025,.05,.25,.5,.75,.95,.975,1]
    quantiles={f"{q*100:g}":float(np.quantile(values,q)) for q in levels}
    fixed=primary["theta"]["ci95_weeks"]
    result=dict(n=len(values),original_theta_weeks=original,mean_weeks=float(values.mean()),median_weeks=float(np.median(values)),
        sample_sd_weeks=float(values.std(ddof=1)),adjusted_sample_skewness=float(skew(values,bias=False)),percentiles_weeks=quantiles,
        mean_minus_original_weeks=float(values.mean()-original),median_minus_original_weeks=float(np.median(values)-original),
        fixed_group_interval_weeks=fixed,fixed_group_interval_width_weeks=float(fixed[1]-fixed[0]),
        bootstrap_percentile_interval_weeks=[quantiles["2.5"],quantiles["97.5"]],bootstrap_interval_width_weeks=quantiles["97.5"]-quantiles["2.5"],
        bootstrap_interval_midpoint_weeks=(quantiles["97.5"]+quantiles["2.5"])/2,
        interpretation="Descriptive distribution of the existing arm-stratified full-score bootstrap, conditional on original released outcomes. It is not a posterior distribution or a newly fitted clinical model.")
    dump(out/"bootstrap_distribution_summary.json",result);return result


def replay_labels(data,draws,out,original,raw_root=None):
    raw_root=Path(raw_root).resolve() if raw_root is not None else original/"raw"
    meta=pd.read_csv(raw_root/"TN10_metadata.csv");baseline=meta.loc[meta["Time.Point"]=="Baseline"].reset_index(drop=True)
    assert baseline.libraryID.tolist()==data.libraryID.tolist()
    counts=pd.read_csv(raw_root/"TN10_raw_counts.csv",index_col=0)
    genes=(original/"results"/"primary_genes.txt").read_text().splitlines()
    raw=counts.loc[:,baseline.libraryID].to_numpy(float).T;expression=np.log2(1+1e6*raw[:,counts.index.get_indexer(genes)]/raw.sum(axis=1)[:,None])
    score,original_high=score_expression(expression)
    np.testing.assert_allclose(score,data.score,rtol=0,atol=1e-12);np.testing.assert_array_equal(original_high,data.high)
    arm=data.treatment.to_numpy(bool);y=data.restricted_weeks.to_numpy();rng=np.random.default_rng(SEED)
    arm_indices=[np.flatnonzero(arm==a) for a in [False,True]]
    selected=np.empty((N_DRAWS,len(data)),np.uint8);slot_labels=np.empty(selected.shape,bool);unique_labels=np.full(selected.shape,-1,np.int8)
    values=[];rows=[];valid=0;unusable=0
    for draw in range(N_DRAWS):
        idx=np.concatenate([rng.choice(i,len(i),replace=True) for i in arm_indices]);_,high=score_expression(expression[idx]);yy,aa=y[idx],arm[idx]
        parts=[yy[(aa==a)&(high==h)] for h in [False,True] for a in [True,False]]
        if min(map(len,parts))<2:unusable+=1;continue
        means=[p.mean() for p in parts];delta_low,delta_high=means[0]-means[1],means[2]-means[3];values.append([delta_low,delta_high,delta_low-delta_high])
        selected[valid]=idx;slot_labels[valid]=high;unique,first=np.unique(idx,return_index=True)
        for person in unique:
            labels=high[idx==person]
            if np.any(labels!=labels[0]):raise AssertionError("Repeated participant received inconsistent reconstructed labels")
        unique_labels[valid,unique]=high[first]
        slot_flips=high!=original_high[idx];unique_flips=high[first]!=original_high[unique]
        rows.append(dict(draw=draw,valid_row=valid,sampled_slots=len(idx),draw_weighted_flips=int(slot_flips.sum()),
            unique_participants=len(unique),unique_inbag_flips=int(unique_flips.sum()),
            draw_weighted_flip_fraction=float(slot_flips.mean()),unique_inbag_flip_fraction=float(unique_flips.mean()),
            **{f"{name}_n":len(part) for (_,_,name,_),part in zip(CELLS,parts)}))
        valid+=1
    assert valid==len(draws)==N_DRAWS and unusable==0
    reproduced=np.asarray(values);saved=draws[["delta_low","delta_high","theta"]].to_numpy();maxdiff=float(np.max(np.abs(reproduced-saved)))
    np.testing.assert_allclose(reproduced,saved,rtol=0,atol=1e-10)
    np.savez_compressed(out/"bootstrap_assignments.npz",source_participant_indices=selected,reconstructed_high_by_slot=slot_labels,
        unique_inbag_high=unique_labels,original_high=original_high,library_ids=data.libraryID.to_numpy(str))
    bydraw=pd.DataFrame(rows);bydraw.to_csv(out/"bootstrap_assignment_by_draw.csv",index=False)
    perperson=[]
    for i,pid in enumerate(data.libraryID):
        included=unique_labels[:,i]>=0;flips=included&(unique_labels[:,i]!=original_high[i])
        perperson.append(dict(library_id=pid,original_high=bool(original_high[i]),draws_included=int(included.sum()),draws_flipped=int(flips.sum()),
            fraction_flipped_given_included=float(flips.sum()/included.sum())))
    persons=pd.DataFrame(perperson);persons.to_csv(out/"bootstrap_assignment_by_participant.csv",index=False)
    result=dict(seed=SEED,draws_requested=N_DRAWS,valid_draws=valid,unusable=unusable,max_absolute_original_effect_discrepancy_weeks=maxdiff,
        all_original_effect_draws_reproduced=True,
        draw_weighted_flipped_slots=int(bydraw.draw_weighted_flips.sum()),draw_weighted_denominator=int(bydraw.sampled_slots.sum()),
        draw_weighted_flip_fraction=float(bydraw.draw_weighted_flips.sum()/bydraw.sampled_slots.sum()),
        unique_inbag_flipped_participant_draws=int(bydraw.unique_inbag_flips.sum()),unique_inbag_denominator=int(bydraw.unique_participants.sum()),
        pooled_unique_inbag_flip_fraction=float(bydraw.unique_inbag_flips.sum()/bydraw.unique_participants.sum()),
        mean_unique_inbag_flip_fraction_per_draw=float(bydraw.unique_inbag_flip_fraction.mean()),
        median_unique_inbag_flip_fraction_per_draw=float(bydraw.unique_inbag_flip_fraction.median()),
        unique_inbag_flip_fraction_per_draw_q25_q75=np.quantile(bydraw.unique_inbag_flip_fraction,[.25,.75]),
        draws_with_any_flip=int((bydraw.unique_inbag_flips>0).sum()),
        median_unique_participants_per_draw=float(bydraw.unique_participants.median()),
        participant_conditional_flip_fraction_median=float(persons.fraction_flipped_given_included.median()),
        participant_conditional_flip_fraction_max=float(persons.fraction_flipped_given_included.max()),
        participants_ever_flipped=int((persons.draws_flipped>0).sum()),
        denominator_definitions={k:protocol()["bootstrap"][k] for k in ["draw_weighted_flips","unique_inbag_flips","participant_flips"]},
        interpretation=protocol()["bootstrap"]["interpretation"])
    assert result["unique_inbag_flipped_participant_draws"]==int(persons.draws_flipped.sum())
    assert result["unique_inbag_denominator"]==int(persons.draws_included.sum())
    dump(out/"bootstrap_assignment_summary.json",result);return result


def figure(draws,profile,followup,summary,out):
    os.environ.setdefault("MPLCONFIGDIR",str(out/"matplotlib_cache"))
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size":12,"axes.titlesize":13,"axes.spines.top":False,"axes.spines.right":False})
    fig,axes=plt.subplots(2,1,figsize=(9.5,8.1),layout="constrained")
    ax=axes[0];ax.hist(draws.theta,bins=45,color="#24758b",alpha=.85,edgecolor="white",linewidth=.35)
    ax.axvline(summary["original_theta_weeks"],color="#a14a29",ls="--",lw=2,label=f"Original {summary['original_theta_weeks']:.2f} weeks")
    ax.axvline(summary["mean_weeks"],color="#172f45",lw=1.7,label=f"Bootstrap mean {summary['mean_weeks']:.2f}")
    ax.set(title="A  Existing full-score bootstrap distribution",xlabel="Low-score benefit minus high-score benefit (weeks)",ylabel="Bootstrap draws")
    ax.legend(frameon=False,fontsize=11,loc="upper left")
    ax.text(.98,.96,f"n = 10,000\nMedian {summary['median_weeks']:.2f}\nSD {summary['sample_sd_weeks']:.2f}",transform=ax.transAxes,ha="right",va="top",fontsize=11)
    ax=axes[1];x=profile.assumed_minimum_stage3_free_weeks
    ax.fill_between(x,profile.theta_lower_weeks,profile.theta_upper_weeks,color="#c4dfe7",label="Finite-sample identification range")
    ax.plot(x,profile.theta_lower_weeks,color="#24758b",lw=2);ax.plot(x,profile.theta_upper_weeks,color="#24758b",lw=2)
    ax.axhline(0,color="#666",lw=.8);ax.axhline(26,color="#a14a29",ls="--",lw=1)
    for key,target,label in [("lower_crosses_zero_at_L_weeks",0,"Lower bound = 0"),("lower_crosses_26_at_L_weeks",26,"Lower bound = 26")]:
        cross=followup[key];ax.plot(cross,target,"o",color="#a14a29",ms=5)
        ax.annotate(f"{label}\nL = {cross:.2f}",(cross,target),xytext=(cross+3,target+38),fontsize=10.5,arrowprops={"arrowstyle":"-","color":"#666"})
    ax.set(xlim=(0,156),xticks=[0,26,52,78,104,130,156],ylabel="Interaction range (weeks)",
        xlabel="Assumed common minimum stage-3-free follow-up L (weeks)",title="B  Follow-up identification profile — not a confidence interval")
    ax.legend(frameon=False,loc="upper right",fontsize=10.5)
    fig.savefig(out/"tn10_bootstrap_followup_diagnostics.png",dpi=220);plt.close(fig)


def report(distribution,assignment,followup,profile,out):
    s=distribution;a=assignment;f=followup
    lines=["# TN10: bounded follow-up and saved-bootstrap diagnostics","",
        "The original clinical models, thresholds, participant groups, outcomes and bootstrap draws are unchanged. The protocol was saved before these descriptive computations. No minimum follow-up duration is inferred from the common encoded endpoint or later RNA sampling.","",
        "## Follow-up identification profile","",f"With the 54 recorded diagnosis outcomes and original four-cell group membership fixed, the 19 uncertain outcomes are allowed to range from L to 156 weeks. The profile is:","",
        f"- Lower: {f['lower_formula']}",f"- Upper: {f['upper_formula']}","",
        f"The lower endpoint equals zero at L={f['lower_crosses_zero_at_L_weeks']:.6f} weeks and equals 26 at L={f['lower_crosses_26_at_L_weeks']:.6f} weeks. Strict inequalities require L strictly above those crossings.","",
        "| Assumed minimum L (weeks) | Lower (weeks) | Upper (weeks) |","|---:|---:|---:|"]
    for level in [0,26,52,78,104,130,156]:
        r=profile.loc[profile.assumed_minimum_stage3_free_weeks==level].iloc[0];lines.append(f"| {level} | {r.theta_lower_weeks:.4f} | {r.theta_upper_weeks:.4f} |")
    lines += ["",f["interpretation"],"","No actual common L has been verified. The calculation shows what additional documented stage-3-free ascertainment would imply for this finite dataset; sampling uncertainty remains a separate question.","",
        "## Existing bootstrap distribution","",f"The 10,000 retained interaction draws have mean {s['mean_weeks']:.4f}, median {s['median_weeks']:.4f} and sample SD {s['sample_sd_weeks']:.4f} weeks. Relative to the original {s['original_theta_weeks']:.4f}-week estimate, their mean is displaced by {s['mean_minus_original_weeks']:.4f} weeks and median by {s['median_minus_original_weeks']:.4f} weeks. Adjusted sample skewness is {s['adjusted_sample_skewness']:.4f}.","",
        f"The fixed-group interval width is {s['fixed_group_interval_width_weeks']:.4f} weeks; the full-score percentile interval width is {s['bootstrap_interval_width_weeks']:.4f} weeks. Both displacement and width are therefore reported directly rather than describing the difference solely as increased uncertainty.","",
        "| Percentile | Interaction (weeks) |","|---:|---:|"]
    lines += [f"| {q} | {v:.4f} |" for q,v in s["percentiles_weeks"].items()]
    lines += ["","## Reconstructed assignment instability","",
        f"The original bootstrap saved effects, not labels. Exact replay of the original seed, arm-order sampling and score construction reproduced all three effect columns in all 10,000 draws; maximum absolute discrepancy was {a['max_absolute_original_effect_discrepancy_weeks']:.3g} weeks. No new outcome model was fitted.","",
        f"Counting repeated sampled copies separately, {a['draw_weighted_flipped_slots']:,} of {a['draw_weighted_denominator']:,} sampled slots changed group versus the source participant's original label ({a['draw_weighted_flip_fraction']:.2%}). Counting each distinct participant only once per draw, {a['unique_inbag_flipped_participant_draws']:,} of {a['unique_inbag_denominator']:,} participant-draw appearances changed ({a['pooled_unique_inbag_flip_fraction']:.2%}). Absent participants are excluded from this latter denominator.","",
        f"At least one label changed in {a['draws_with_any_flip']:,} draws. The median per-draw unique-participant flip fraction was {a['median_unique_inbag_flip_fraction_per_draw']:.2%}; its interquartile range was {a['unique_inbag_flip_fraction_per_draw_q25_q75'][0]:.2%}–{a['unique_inbag_flip_fraction_per_draw_q25_q75'][1]:.2%}. Participant-specific conditional flip frequencies, complete draw indices and labels are retained in the output directory.","",a["interpretation"],"",
        "The figure displays the existing bootstrap distribution and the assumption-indexed follow-up profile. It does not merge identification uncertainty with sampling uncertainty or change the original decision rule.",""]
    (out/"tn10_diagnostics_report.md").write_text("\n".join(lines),encoding="utf-8")


def run(out,data_root,original,revision2,raw_root=None):
    # Logical identifiers bind source bytes without publishing machine paths.
    raw_root=Path(raw_root).resolve() if raw_root is not None else original/"raw"
    inputs={f"original/raw/{n}":raw_root/n for n in ["TN10_metadata.csv","TN10_raw_counts.csv"]}
    inputs.update({f"original/results/{n}":original/"results"/n for n in ["analysis_dataset.tsv","primary_results.json","bootstrap_results.json","bootstrap_draws.tsv","primary_genes.txt"]})
    inputs.update({f"revision2/{n}":revision2/n for n in ["followup_identification_range.json","followup_identification_cells.csv"]})
    inputs["source/card002_run.py"]=Path(__file__).parent/"card002_run.py"
    before={name:sha(p) for name,p in inputs.items()};data=pd.read_csv(original/"results"/"analysis_dataset.tsv",sep="\t")
    primary=json.loads((original/"results"/"primary_results.json").read_text());bootstrap=json.loads((original/"results"/"bootstrap_results.json").read_text())
    draws=pd.read_csv(original/"results"/"bootstrap_draws.tsv",sep="\t")
    assert len(data)==73 and len(draws)==N_DRAWS and list(draws)==["delta_low","delta_high","theta"]
    np.testing.assert_allclose(np.quantile(draws.theta,[.025,.975]),bootstrap["percentile_ci95_weeks"]["theta"],rtol=0,atol=1e-10)
    profile,f=followup_profile(data,primary,out,revision2);s=distribution(draws,primary,out);a=replay_labels(data,draws,out,original,raw_root)
    figure(draws,profile,f,s,out);report(s,a,f,profile,out)
    after={name:sha(p) for name,p in inputs.items()};assert before==after,"Original source changed during diagnostics"
    def logical_root(root,flag):
        try:return Path(root).relative_to(data_root).as_posix()
        except ValueError:return f"<{flag}>"
    dump(out/"diagnostic_manifest.json",dict(status="COMPLETE_BOUNDED_DIAGNOSTICS",completed_utc=datetime.now(timezone.utc).isoformat(),
        source_sha256=sha(__file__),protocol_sha256=sha(out/"diagnostic_protocol.json"),original_inputs_unchanged=True,
        original_input_sha256=before,input_roots={"original":logical_root(original,"original-root"),"revision2":logical_root(revision2,"revision2-root"),"source":"src"},
        all_bootstrap_effect_draws_reproduced=True,clinical_models_refitted=False,actual_L_verified=False,
        artifacts={p.name:sha(p) for p in out.iterdir() if p.is_file() and p.name!="diagnostic_manifest.json"}))
    print(json.dumps(dict(status="complete",followup=f,distribution=s,assignments=a),default=lambda x:x.tolist() if isinstance(x,np.ndarray) else x),flush=True)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root",type=Path,required=True,help="Research data root containing card002 and revision_20261022_r2")
    parser.add_argument("--original-root",type=Path,help="Override original card002 inputs (default: DATA_ROOT/card002)")
    parser.add_argument("--revision2-root",type=Path,help="Override prior bounds (default: DATA_ROOT/revision_20261022_r2/card002)")
    parser.add_argument("--output-root",type=Path,help="New output directory (default: DATA_ROOT/revision_20261022_r3/paper/reconstruction_tn10)")
    parser.add_argument("--raw-root",type=Path,help="Read-only original input directory")
    parser.add_argument("--stage",choices=["freeze","run"],default="run")
    args=parser.parse_args();data_root=args.data_root.resolve()
    original=(args.original_root or data_root/"card002").resolve()
    revision2=(args.revision2_root or data_root/"revision_20261022_r2"/"card002").resolve()
    out=(args.output_root or data_root/"revision_20261022_r3"/"paper"/"reconstruction_tn10").resolve()
    freeze(out,original,revision2)
    if args.stage=="run":run(out,data_root,original,revision2,args.raw_root)
