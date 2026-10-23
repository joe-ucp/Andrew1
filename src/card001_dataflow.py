"""Reconstruct source-only reason counts without examining held-out outcomes."""
import argparse,json,zipfile
from pathlib import Path
import numpy as np
import pandas as pd
from card001_data import seconds


def run(root):
    manifest=json.loads((root/"participant_manifest.json").read_text())
    rows={p["id"]:p for p in manifest["participants"]}
    counts=dict(source_rows=0,not_randomized=0,invalid_timestamp=0,before_randomization=0,at_or_after_followup_end=0,in_randomized_window=0)
    with zipfile.ZipFile(manifest["archive"]) as z:
        with z.open("Data Files/DexcomClarityCGM_a.txt") as f:
            for chunk in pd.read_csv(f,sep="|",encoding="utf-16",usecols=["PtID","DataDtTm","DataDtTm_adj"],dtype={"DataDtTm":"str","DataDtTm_adj":"str"},chunksize=250000):
                counts["source_rows"]+=len(chunk)
                for pid,p in chunk.groupby("PtID"):
                    if int(pid) not in rows: counts["not_randomized"]+=len(p); continue
                    effective=p.DataDtTm_adj.fillna(p.DataDtTm)
                    valid=pd.to_datetime(effective,format="mixed",errors="coerce").notna()
                    counts["invalid_timestamp"]+=int((~valid).sum())
                    ts=seconds(effective[valid]); row=rows[int(pid)]
                    before=ts<row["start"]; after=ts>=row["end"]
                    counts["before_randomization"]+=int(before.sum());counts["at_or_after_followup_end"]+=int(after.sum())
                    counts["in_randomized_window"]+=int((~before&~after).sum())
    if counts["source_rows"]!=sum(counts[k] for k in counts if k!="source_rows"): raise AssertionError("Source-flow reasons do not partition original rows")
    (root/"source_data_flow.json").write_text(json.dumps(counts,indent=2))
    pd.DataFrame([dict(stage=k,records=v) for k,v in counts.items()]).to_csv(root/"source_data_flow.csv",index=False)
    clock=[]
    for path in (root/"clocks").glob("DCLP3_*.npz"):
        a=np.load(path); dt=np.diff(a["time"]); dx=np.diff(a["transmitter"])
        short=(dt>0)&(dt<=1800)
        for reason,which in [("transmitter_reset",short&(dx<0)),("clock_increment_disagreement",short&(dx>0)&(np.abs(dt-dx)>60)),("stalled_counter_sensitivity_only",short&(dx==0))]:
            for i in np.flatnonzero(which)+1:
                clock.append(dict(id=int(path.stem.split("_")[1]),time=float(a["time"][i]),reason=reason,primary_barrier=reason!="stalled_counter_sensitivity_only"))
    pd.DataFrame(clock).to_csv(root/"source_clock_exclusion_ledger.csv",index=False)
    print(json.dumps(counts,indent=2))

if __name__=="__main__":
    p=argparse.ArgumentParser();p.add_argument("--data-root",type=Path,required=True)
    run(p.parse_args().data_root)
