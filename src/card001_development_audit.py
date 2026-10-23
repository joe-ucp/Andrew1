"""Development-only opportunity decomposition; never opens held-out outcomes."""
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from card001_core import episodes,episodes_independent,ready
from card001_run import load_participants


def run(root):
    rows=[]
    for p in load_participants(root,["train","validation"]):
        g=p["glucose"]; onsets,censored=episodes(g)
        other,unknown=episodes_independent(g)
        np.testing.assert_array_equal(onsets,other); np.testing.assert_array_equal(censored,unknown)
        opportunities={}
        for n in [4,25]:
            available=ready(g,n)&(g>=54)
            opportunities[n]=np.array([available[max(0,e-6):max(0,e-2)].any() for e in onsets],bool)
        if np.any(opportunities[25]&~opportunities[4]): raise AssertionError("Long history not subset of short history")
        rows.append(dict(id=p["id"],split=p["split"],arm=p["arm"],events=len(onsets),censored_low_runs=len(censored),M120_opportunities=int(opportunities[25].sum()),M15_opportunities=int(opportunities[4].sum()),additional_short_opportunities=int((opportunities[4]&~opportunities[25]).sum()),observed_days=np.isfinite(g).sum()/288,calendar_days=len(g)/288))
    frame=pd.DataFrame(rows); frame.to_csv(root/"development_opportunity_audit.csv",index=False)
    report={}
    for split,part in frame.groupby("split"):
        total=part[["events","censored_low_runs","M120_opportunities","M15_opportunities","additional_short_opportunities","observed_days","calendar_days"]].sum().to_dict()
        total["participants"]=len(part)
        total["M120_coverage"]=total["M120_opportunities"]/total["events"] if total["events"] else None
        total["M15_coverage"]=total["M15_opportunities"]/total["events"] if total["events"] else None
        report[split]=total
    report["interpretation"]="Opportunity coverage bounds sensitivity; this does not estimate predictive sensitivity or clinical benefit. Unknown low onsets remain outside event denominator."
    report["timing_note"]="Executed during prespecified model fitting but before validation threshold selection or any held-out outcome replay; no model/endpoint choice was changed in response."
    (root/"development_opportunity_audit.json").write_text(json.dumps(report,indent=2)); print(json.dumps(report,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data-root",type=Path,required=True)
    run(p.parse_args().data_root)
