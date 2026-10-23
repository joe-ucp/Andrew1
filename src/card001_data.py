"""DCLP3 original archive adapter. Study semantics are explicitly recorded."""
from __future__ import annotations
import argparse,hashlib,json,zipfile,time
from pathlib import Path
import numpy as np
import pandas as pd
from card001_core import causal_grid,participant_split


def seconds(values):
    dt=pd.to_datetime(values,format="mixed",errors="coerce")
    return dt.astype("datetime64[ns]").astype("int64").to_numpy(dtype=float)/1e9


def prepare_dclp3(archive,out):
    out=Path(out); out.mkdir(parents=True,exist_ok=True)
    processed=out/"processed"; processed.mkdir(exist_ok=True)
    clockdir=out/"clocks"; clockdir.mkdir(exist_ok=True)
    z=zipfile.ZipFile(archive)
    def table(name,**kw): return pd.read_csv(z.open("Data Files/"+name),sep="|",encoding="utf-16",**kw)
    roster=table("PtRoster_a.txt")
    visits=table("VisitInfo_a.txt")
    randomized=roster[roster.RandDt.notna() & roster.trtGroup.isin(["CLC","SAP"])].copy()
    final=visits[visits.Visit.eq("26-Week")].set_index("PtID").VisitDt.to_dict()
    rows=[]; bounds={}
    for r in randomized.itertuples():
        start=pd.Timestamp(r.RandDt).timestamp()
        end=start+182*86400
        # Visit dates have day resolution. Retain the whole final observed day.
        if r.PtID in final: end=min(end,pd.Timestamp(final[r.PtID]).timestamp()+86400)
        split,reason=participant_split(r.PtID)
        bounds[int(r.PtID)]=(start,end)
        rows.append(dict(id=int(r.PtID),cohort="DCLP3",split=split,split_reason=reason,arm=r.trtGroup,start=start,end=end,file=str(processed/f"DCLP3_{r.PtID}.npz")))
    manifest=dict(source_semantics_audited=True,archive=str(archive),archive_sha256=hashlib.file_digest(open(archive,"rb"),"sha256").hexdigest(),participants=rows,exclusions=[dict(id=int(r.PtID),reason="not_randomized") for r in roster[~roster.PtID.isin(randomized.PtID)].itertuples()],rules={"clock":"DataDtTm_adj when nonempty else DataDtTm; epoch anchored local five-minute grid", "window":"randomization through 182 days; whole final 26-week visit day when earlier", "clock_barrier":"transmitter reset or disagreement >60sec between adjacent positive local/transmitter increments with local gap <=1800sec; first grid bin at barrier unavailable", "raw_primary_only":True,"no_outlier_deletion":True})
    # Freeze participant allocation before reading CGM outcomes.
    (out/"participant_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    parts={int(p):[] for p in randomized.PtID}; total=0; invalid_time=0
    for chunk in table("DexcomClarityCGM_a.txt",chunksize=250000,usecols=["PtID","RecID","DataDtTm","DataDtTm_adj","CGM","TransmitterTime"],dtype={"DataDtTm":"str","DataDtTm_adj":"str"}):
        total+=len(chunk)
        effective=chunk.DataDtTm_adj.fillna(chunk.DataDtTm)
        dt=pd.to_datetime(effective,format="mixed",errors="coerce")
        valid=dt.notna(); invalid_time+=int((~valid).sum())
        chunk=chunk[valid].copy()
        chunk["time"]=seconds(effective[valid])
        for pid,part in chunk.groupby("PtID",sort=False):
            if int(pid) not in bounds: continue
            start,end=bounds[int(pid)]
            part=part[(part.time>=start)&(part.time<end)]
            if len(part): parts[int(pid)].append(part[["time","CGM","TransmitterTime","RecID"]].to_numpy(dtype=float))
        print(f"parsed_raw_rows={total}",flush=True)
    audit=[]
    for row in rows:
        pid=row["id"]
        a=np.concatenate(parts.pop(pid)) if parts[pid] else np.empty((0,4))
        if len(a):
            a=a[np.argsort(a[:,0],kind="stable")]
            dt=np.diff(a[:,0]); dtx=np.diff(a[:,2])
            reset=np.isfinite(dtx)&(dtx<0)&(dt>0)&(dt<=1800)
            inconsistent=np.isfinite(dtx)&(dtx>0)&(dt>0)&(dt<=1800)&(np.abs(dt-dtx)>60)
            barrier_times=a[1:,0][reset|inconsistent]
            grid=causal_grid(a[:,0],a[:,1],row["start"],row["end"]-1e-6)
            barriers=np.isin(grid["time"],np.ceil(barrier_times/300)*300)
            original_grid=grid["glucose"].copy()
            grid["glucose"][barriers]=np.nan
            # Preserve source arrays so every sensitivity can be rebuilt honestly.
            np.savez_compressed(clockdir/f"DCLP3_{pid}.npz",time=a[:,0],glucose=a[:,1],transmitter=a[:,2],recid=a[:,3],barrier_times=barrier_times)
            duplicate_count=len(a)-len(np.unique(a[:,0])); conflict_count=int(grid["conflict"].sum())
        else:
            t=np.arange(np.ceil(row["start"]/300)*300,row["end"],300)
            grid=dict(time=t,glucose=np.full(len(t),np.nan),source_time=np.full(len(t),np.nan),conflict=np.zeros(len(t),bool))
            original_grid=grid["glucose"].copy(); barriers=np.zeros(len(t),bool); duplicate_count=conflict_count=0
        np.savez_compressed(row["file"],**grid,clock_barrier=barriers,glucose_without_clock_barriers=original_grid)
        audit.append(dict(id=pid,split=row["split"],arm=row["arm"],raw_in_window=len(a),duplicate_timestamps=duplicate_count,conflict_grid_bins=conflict_count,clock_barrier_bins=int(barriers.sum()),observed_grid_bins=int(np.isfinite(grid["glucose"]).sum()),calendar_grid_bins=len(grid["time"])))
    pd.DataFrame(audit).to_csv(out/"source_time_audit.csv",index=False)
    manifest["audit"]=dict(total_raw_rows=total,invalid_timestamps=invalid_time,participants=len(rows),schema="original raw DexcomClarityCGM_a.txt; no merged pump stream")
    (out/"participant_manifest.json").write_text(json.dumps(manifest,indent=2),encoding="utf-8")
    return manifest


if __name__=="__main__":
    p=argparse.ArgumentParser(); p.add_argument("--archive",required=True); p.add_argument("--output",required=True)
    a=p.parse_args(); prepare_dclp3(Path(a.archive),Path(a.output))
