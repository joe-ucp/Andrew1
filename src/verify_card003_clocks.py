"""Independent CLVer literal-format clock audit; no clinical outcome fitting."""
import argparse
from collections import Counter
import json
from pathlib import Path
import zipfile
import numpy as np
import pandas as pd


def verify(root):
    root=Path(root);raw={};date_only_times={};counts=Counter();issues=[]
    roster=pd.read_csv(root/"processed"/"randomized_roster.csv",dtype={"PtID":str})
    sources=dict(zip(roster.PtID,roster.source))
    for pid in sources:
        with np.load(root/"processed"/"raw_streams"/f"{pid}.npz") as z:
            order=np.argsort(z["recid"])
            explicit=z["explicit_time"] if "explicit_time" in z else np.ones(len(order),bool)
            raw[pid]=(z["recid"][order],z["time"][order],explicit[order])
            date_only_times[pid]=np.setdiff1d(z["time"][~explicit],z["time"][explicit])
            counts["primary_explicit_time_rows"]+=int(explicit.sum())
            counts["date_only_physical_candidates"]+=int((~explicit).sum())
            t,x=z["time"],z["transmitter"]
            counts["retained_rows"]+=len(t)
            if len(t):
                unique,j,c=np.unique(t,return_index=True,return_counts=True)
                lo=np.minimum.reduceat(x,j);hi=np.maximum.reduceat(x,j)
                bad=np.isfinite(lo)&np.isfinite(hi)&(lo!=hi)
                counts["same_time_internal_conflicts"]+=int(bad.sum())
                counts["duplicate_timestamp_rows"]+=len(t)-len(unique)
                if bad.any():issues.append(dict(kind="same_time_internal_conflict",id=pid,count=int(bad.sum())))
            counts["raw_barrier_times"]+=len(z["barrier_times"])
    original={}
    with zipfile.ZipFile(root/"raw"/"CLVerPublicDataset.zip") as archive:
        for table,source in [("DeviceCGM.txt","DeviceCGM"),("DevicePumpCGMDATAGXB.txt","GXB")]:
            member=next(n for n in archive.namelist() if n.endswith("/"+table))
            subtotal=Counter();failed_strings=Counter()
            for frame in pd.read_csv(archive.open(member),sep="|",encoding="utf-8-sig",dtype=str,
                                     keep_default_na=False,usecols=["PtID","RecID","DeviceDtTm"],chunksize=250000):
                # Original table dates are explicit US month/day/year with AM/PM.
                parsed=pd.to_datetime(frame.DeviceDtTm,format="%m/%d/%Y %I:%M:%S %p",errors="coerce")
                bad=parsed.isna();failed_strings.update(frame.loc[bad,"DeviceDtTm"])
                date_only=frame.DeviceDtTm.str.fullmatch(r"\d{1,2}/\d{1,2}/\d{4}",na=False)
                subtotal["source_rows"]+=len(frame);subtotal["strict_fulltimestamp_format_unparsed"]+=int(bad.sum())
                subtotal["date_only_requires_midnight_interpretation"]+=int(date_only.sum())
                parsed.loc[date_only]=pd.to_datetime(frame.loc[date_only,"DeviceDtTm"],format="%m/%d/%Y",errors="coerce")
                subtotal["unrecognized_literal_date_formats"]+=int(parsed.isna().sum())
                frame["seconds"]=parsed.astype("datetime64[ns]").astype("int64").to_numpy()/1e9
                frame["explicit_time"]=(~date_only)&parsed.notna()
                for pid,part in frame.groupby("PtID",sort=False):
                    if sources.get(pid)!=source:continue
                    ids,times,explicit=raw[pid]
                    if not len(ids):continue
                    rid=pd.to_numeric(part.RecID).to_numpy()
                    loc=np.searchsorted(ids,rid);valid=loc<len(ids)
                    at=np.flatnonzero(valid);at=at[ids[loc[at]]==rid[at]]
                    if not len(at):continue
                    subset=part.iloc[at]
                    matched=np.isclose(times[loc[at]],subset.seconds.to_numpy(),rtol=0,atol=1e-6)
                    subtotal["retained_rows_crosschecked"]+=len(at)
                    subtotal["epoch_mismatches_under_explicit_midnight_interpretation"]+=int((~matched).sum())
                    subtotal["explicit_time_flag_mismatches"]+=int((explicit[loc[at]]!=subset.explicit_time.to_numpy()).sum())
            original[table]={**subtotal,"unparsed_literal_values":dict(failed_strings)}
    for path in (root/"processed"/"interval_grids").glob("*.npz"):
        with np.load(path) as z:
            t=z["time"];s=z["source_time"];valid=np.isfinite(s)
            counts["interval_grids"]+=1
            counts["noncausal_source_bins"]+=int((valid&((s>t)|(s<=t-300))).sum())
            counts["observed_conflict_bins"]+=int((z["conflict"]&np.isfinite(z["glucose"])).sum())
            counts["observed_barrier_bins"]+=int((z["clock_barrier"]&np.isfinite(z["glucose"])).sum())
            pid=path.stem.split("_")[0]
            counts["observed_bins_from_date_only_source"]+=int((np.isfinite(z["glucose"])&np.isin(s,date_only_times[pid])).sum())
    result=dict(source_counts=original,canonical_counts=dict(counts),issues=issues,
        scope="Independent literal timestamp format, retained epoch, identity and causal grid audit")
    result["computational_checks_pass_under_midnight_interpretation"]=all(v.get("epoch_mismatches_under_explicit_midnight_interpretation",0)==0 and v["unrecognized_literal_date_formats"]==0 for v in original.values()) and all(counts[k]==0 for k in ["noncausal_source_bins","observed_conflict_bins","observed_barrier_bins"])
    result["strict_primary_date_mask_pass"]=all(v.get("explicit_time_flag_mismatches",0)==0 for v in original.values()) and counts["observed_bins_from_date_only_source"]==0
    result["date_only_source_resolution_requires_explicit_adjudication"]=any(v["date_only_requires_midnight_interpretation"] for v in original.values()) and not result["strict_primary_date_mask_pass"]
    result["primary_resolution"]="Date-only candidates excluded; original physical records preserved" if result["strict_primary_date_mask_pass"] else "Primary time-resolution review required"
    result["interpretation_boundary"]="Date-only literals have no printed clock time. Midnight matches their stored candidate epochs but remains an interpretation confined to a separately labeled sensitivity. The strict primary excludes those records; date_only_source_records.csv retains adjacent-clock evidence."
    (root/"results"/"independent_clock_verification.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data-root",type=Path,required=True)
    a=p.parse_args();print(json.dumps(verify(a.data_root),indent=2))
