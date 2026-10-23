"""Read-only duplicate-internal-clock audit; computes no low-event outcomes."""
import argparse
import json
from pathlib import Path
import numpy as np


def audit(root):
    root=Path(root); manifest=json.loads((root/"participant_manifest.json").read_text())
    pmap={r["id"]:r for r in manifest["participants"]}
    rows=[];total=duplicates=ambiguities=0
    for path in sorted((root/"clocks").glob("*.npz")):
        with np.load(path) as z:
            t=z["time"];x=z["transmitter"]
        if not len(t): continue
        if (np.diff(t)<0).any(): raise AssertionError("Source clocks must be sorted")
        u,j,c=np.unique(t,return_index=True,return_counts=True)
        low=np.minimum.reduceat(x,j);high=np.maximum.reduceat(x,j)
        bad=np.isfinite(low)&np.isfinite(high)&(low!=high)
        n=int(bad.sum());pid=int(path.stem.split("_")[-1])
        total+=len(t);duplicates+=int((c-1).sum());ambiguities+=n
        if n:
            with np.load(pmap[pid]["file"]) as g:
                bins=np.ceil(u[bad]/300)*300
                idx=np.searchsorted(g["time"],bins)
                idx=idx[idx<len(g["time"])]
                group=np.flatnonzero(bad)
                both=np.unique(np.r_[group,group[group+1<len(u)]+1])
                conservative_bins=np.unique(np.ceil(u[both]/300)*300)
                both_idx=np.searchsorted(g["time"],conservative_bins)
                both_idx=both_idx[both_idx<len(g["time"])]
                delta=np.diff(u); first=x[j];last=x[j+c-1]
                def flags(increment):
                    return (delta>0)&(delta<=1800)&np.isfinite(increment)&((increment<0)|((increment>0)&(np.abs(delta-increment)>60)))
                original=flags(first[1:]-last[:-1]);reversed_order=flags(last[1:]-first[:-1])
                changed=np.unique(np.ceil(u[1:][original!=reversed_order]/300)*300)
                changed_idx=np.searchsorted(g["time"],changed)
                changed_idx=changed_idx[changed_idx<len(g["time"])]
                rows.append(dict(id=pid,split=pmap[pid]["split"],ambiguous_timestamps=n,
                    already_unavailable_bins=int((~np.isfinite(g["glucose"][idx])).sum()),
                    still_observed_bins=int(np.isfinite(g["glucose"][idx]).sum()),
                    group_plus_following_distinct_time_bins=len(both_idx),
                    group_plus_following_still_observed_bins=int(np.isfinite(g["glucose"][both_idx]).sum()),
                    within_duplicate_order_reversal_changed_boundary_bins=len(changed_idx),
                    changed_boundary_still_observed_bins=int(np.isfinite(g["glucose"][changed_idx]).sum())))
    result=dict(raw_rows=total,duplicate_rows=duplicates,ambiguous_internal_times=ambiguities,
        participants=rows,outcomes_evaluated=False,raw_arrays_modified=False)
    output=root/"results";output.mkdir(exist_ok=True)
    (output/"duplicate_internal_clock_audit.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":
    parser=argparse.ArgumentParser();parser.add_argument("--data-root",required=True,type=Path)
    args=parser.parse_args()
    print(json.dumps(audit(args.data_root),indent=2))
