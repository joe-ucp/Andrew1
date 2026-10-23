"""REPLACE-BG original-archive adapter and recorded-clock audit.

Only an authorized local source archive is accepted. Inspect the archive's
original README and roster semantics before deriving any external test cohort.
The source card describes candidate field names; it does not verify their data.
"""
from __future__ import annotations
import argparse
import csv
import hashlib
import io
import json
import re
import zipfile
import time
from collections import Counter
from pathlib import Path

import numpy as np
import pandas as pd
from card001_core import causal_grid

EXPECTED_BYTES = 174942355
EXPECTED_SHA256 = "6d4bb4e47e71fff72cd3eb94c5a5e60c9f32e066606bb3883821debf528dc004"
REQUIRED_CGM_COLUMNS = {
    "PtID", "DeviceDtTmDaysFromEnroll", "DeviceTm",
    "DexInternalDtTmDaysFromEnroll", "DexInternalTm",
    "RecordType", "GlucoseValue", "RecID",
}


def file_hash(path):
    with Path(path).open("rb") as f:
        return hashlib.file_digest(f, "sha256").hexdigest()


def relative_clock_seconds(days, clocks):
    """Enrollment-relative integer date plus strictly valid local clock time.

    No artificial date or host timezone enters this pseudo-time coordinate.
    Invalid day/time values remain NaN. Numeric clocks without documented units
    are deliberately unsupported. Accepted textual forms are HH:MM[:SS[.fff]]
    and h:MM[:SS] AM/PM, interpreted only as time of day.
    """
    day = pd.to_numeric(pd.Series(days), errors="coerce").to_numpy(float)
    raw = pd.Series(clocks).astype("string").str.strip()
    pattern = r"^(\d{1,2}):(\d{2})(?::(\d{2}(?:\.\d+)?))?\s*([AaPp][Mm])?$"
    fields = raw.str.extract(pattern)
    hour = pd.to_numeric(fields[0], errors="coerce").to_numpy(float)
    minute = pd.to_numeric(fields[1], errors="coerce").to_numpy(float)
    second = pd.to_numeric(fields[2].fillna("0"), errors="coerce").to_numpy(float)
    ampm = fields[3].str.upper()
    is12 = ampm.notna().to_numpy()
    valid = np.isfinite(day) & (day == np.floor(day)) & np.isfinite(hour)
    valid &= (minute >= 0) & (minute < 60) & (second >= 0) & (second < 60)
    valid &= np.where(is12, (hour >= 1) & (hour <= 12), (hour >= 0) & (hour < 24))
    hour[is12] = np.mod(hour[is12],12) + np.where(ampm[is12].eq("PM"),12,0)
    out = day*86400+hour*3600+minute*60+second
    out[~valid] = np.nan
    return out


def clock_barriers(device_seconds, internal_seconds, *, max_local_gap=1800, discrepancy=60):
    """Conservative local/internal-clock mismatch boundaries, not corrections.

    Caller supplies records in stable device-time order after duplicate audit.
    A gap marker breaks event and feature continuity; clocks are never shifted to
    improve coverage. Missing internal clocks cannot establish clock agreement.
    """
    device=np.asarray(device_seconds,float)
    internal=np.asarray(internal_seconds,float)
    if device.ndim!=1 or device.shape!=internal.shape:
        raise ValueError("Device and internal clocks must be aligned vectors")
    if len(device)<2:
        return np.array([],float),np.array([],dtype="U32")
    local_delta=np.diff(device)
    internal_delta=np.diff(internal)
    adjacent=np.isfinite(local_delta)&(local_delta>0)&(local_delta<=max_local_gap)&np.isfinite(internal_delta)
    reset=adjacent&(internal_delta<0)
    frozen=adjacent&(internal_delta==0)
    disagreement=adjacent&(internal_delta>0)&(np.abs(local_delta-internal_delta)>discrepancy)
    flagged=reset|frozen|disagreement
    reasons=np.where(reset,"internal_clock_reset",np.where(frozen,"internal_clock_stalled","clock_increment_disagreement"))
    return device[1:][flagged],reasons[flagged]


def canonical_arrays(device, glucose, internal, recid, start, end):
    """Shared canonical grid contract; end is exclusive, all source rows retained.

    Source-specific record types/units/cohort/date windows must be validated by
    the archive adapter before this pure recorded-time transformation.
    """
    device=np.asarray(device,float); glucose=np.asarray(glucose,float)
    internal=np.asarray(internal,float); recid=np.asarray(recid)
    if not (device.ndim==1 and device.shape==glucose.shape==internal.shape==recid.shape):
        raise ValueError("Raw arrays must have identical one-dimensional shape")
    keep=np.isfinite(device)&(device>=start)&(device<end)
    order=np.argsort(device[keep],kind="stable")
    device,glucose,internal,recid=(x[keep][order] for x in [device,glucose,internal,recid])
    if len(device):
        # Group clock evidence before comparing increments. Source export order
        # among same-time rows must never decide where continuity is asserted.
        unique,first=np.unique(device,return_index=True)
        smallest=np.minimum.reduceat(internal,first)
        largest=np.maximum.reduceat(internal,first)
        ambiguous=np.isfinite(smallest)&np.isfinite(largest)&(smallest!=largest)
        consistent=smallest.copy(); consistent[ambiguous]=np.nan
        barrier_times,reasons=clock_barriers(unique,consistent)
        barrier_times=np.r_[barrier_times,unique[ambiguous]]
        reasons=np.r_[reasons,np.full(int(ambiguous.sum()),"same_time_internal_disagreement",dtype="U32")]
        sort=np.argsort(barrier_times,kind="stable")
        barrier_times,reasons=barrier_times[sort],reasons[sort]
        grid=causal_grid(device,glucose,start,end-1e-6)
        barriers=np.isin(grid["time"],np.ceil(barrier_times/300)*300)
    else:
        gridtime=np.arange(np.ceil(start/300)*300,end,300)
        grid={"time":gridtime,"glucose":np.full(len(gridtime),np.nan),"source_time":np.full(len(gridtime),np.nan),"conflict":np.zeros(len(gridtime),bool)}
        barriers=np.zeros(len(gridtime),bool); barrier_times=np.array([],float); reasons=np.array([],dtype="U32")
    grid["glucose_without_clock_barriers"]=grid["glucose"].copy()
    grid["glucose"][barriers]=np.nan
    grid["clock_barrier"]=barriers
    raw={"time":device,"glucose":glucose,"transmitter":internal,"recid":recid,"barrier_times":barrier_times,"barrier_reasons":reasons}
    return grid,raw


def inspect_archive(archive,out):
    """Inventory/source-schema evidence only; never compute low-event outcomes."""
    archive=Path(archive).resolve(); out=Path(out).resolve()
    out.mkdir(parents=True,exist_ok=True)
    docs=out/"documentation"; docs.mkdir(exist_ok=True)
    result={"archive":str(archive),"bytes":archive.stat().st_size,"sha256":file_hash(archive),
            "source_card_expected_bytes":EXPECTED_BYTES,"same_size_as_source_card":archive.stat().st_size==EXPECTED_BYTES,
            "semantics_verified":False,"external_outcomes_evaluated":False,"tables":{},"documentation":[]}
    with zipfile.ZipFile(archive) as z:
        result["inventory"]=[{"name":m.filename,"bytes":m.file_size} for m in z.infolist()]
        for member in z.infolist():
            name=Path(member.filename).name
            if re.search(r"readme|glossary|dictionary",name,re.I):
                content=z.read(member)
                target=docs/name
                if target.exists() and target.read_bytes()!=content:
                    raise ValueError("Documentation basename collision")
                target.write_bytes(content)
                result["documentation"].append({"member":member.filename,"file":str(target),"sha256":hashlib.sha256(content).hexdigest()})
        for base in ["HDeviceCGM.txt","HPtRoster.txt"]:
            matches=[n for n in z.namelist() if Path(n).name==base]
            if len(matches)!=1:
                raise ValueError(f"Expected exactly one {base}; found {matches}")
            with z.open(matches[0]) as f:
                prefix=f.read(65536)
            encoding="utf-16" if prefix[:2] in [b"\xff\xfe",b"\xfe\xff"] else "utf-8-sig"
            first=prefix.decode(encoding,errors="strict").splitlines()[0]
            delimiter=max(["|","\t",","],key=first.count)
            columns=next(csv.reader([first],delimiter=delimiter))
            result["tables"][base]={"member":matches[0],"encoding":encoding,"delimiter":delimiter,"columns":columns}
        if not REQUIRED_CGM_COLUMNS.issubset(result["tables"]["HDeviceCGM.txt"]["columns"]):
            raise ValueError("CGM schema differs from source card; new source audit required")
    (out/"archive_schema_audit.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    return result


def cohort_rows(roster, visits, final, out, *, days=182):
    """Literal enrollment-time card cohort, never selected by CGM outcomes.

    No withdrawal date is published except a possible death date. The final
    documented non-missed visit is therefore the follow-up truncation, including
    its entire day. This is not asserted to be the actual withdrawal instant.
    """
    if roster.PtID.duplicated().any():
        raise ValueError("Participant roster has duplicate identities")
    randomized=roster.RandDtDaysAfterEnroll.notna() & roster.TrtGroup.isin(["CGM Only","CGM+BGM"])
    dated=visits[visits.VisitDtDaysFromEnroll.notna() & ~visits.VisitMiss.eq(1)]
    final_days=dated.groupby("PtID").VisitDtDaysFromEnroll.max().to_dict()
    deaths=final[final.DeathDtDaysAfterEnroll.notna()].set_index("PtID").DeathDtDaysAfterEnroll.to_dict()
    rows=[]; exclusions=[]
    for row in roster.sort_values("PtID").itertuples():
        pid=int(row.PtID)
        if not bool(randomized.loc[row.Index]):
            exclusions.append(dict(id=pid,reason="not_randomized")); continue
        if pid==183:
            exclusions.append(dict(id=pid,reason="generation_schema_exposure")); continue
        if pid not in final_days:
            raise ValueError(f"No documented final follow-up for {pid}")
        final_day=float(final_days[pid])
        if pid in deaths: final_day=min(final_day,float(deaths[pid]))
        if not np.isfinite(final_day) or final_day<0:
            raise ValueError(f"Invalid final follow-up for {pid}")
        end=min(days*86400,(final_day+1)*86400)
        rows.append(dict(id=pid,cohort="REPLACE-BG",split="external_test",
            split_reason="independent_external_cohort",arm=row.TrtGroup,
            status=row.PtStatus,start=0.0,end=float(end),
            randomization_day=float(row.RandDtDaysAfterEnroll),
            final_documented_visit_day=final_day,
            endpoint_reason="final_documented_followup_whole_day" if end<days*86400 else "fixed_enrollment_horizon",
            file=str(Path(out)/"processed"/f"REPLACE-BG_{pid}.npz")))
    return rows,exclusions


def prepare_replacebg(archive,out):
    """Reconstruct original recorded clocks/masks without computing outcomes."""
    begun=time.monotonic(); archive=Path(archive).resolve(); out=Path(out).resolve()
    inspection=inspect_archive(archive,out)
    if inspection["sha256"]!=EXPECTED_SHA256:
        raise ValueError("Archive differs from audited original; inspect new source before analysis")
    with zipfile.ZipFile(archive) as z:
        def table(name,**kw):
            return pd.read_csv(z.open("Data Tables/"+name),sep="|",encoding="utf-8-sig",**kw)
        roster=table("HPtRoster.txt"); visits=table("HVisitInfo.txt"); final=table("HPostRandPtFinalStatus.txt")
        roots=[out,out/"inclusive_day182"]
        manifests=[]
        for target,days in zip(roots,[182,183]):
            for directory in [target,target/"processed",target/"clocks"]: directory.mkdir(parents=True,exist_ok=True)
            rows,exclusions=cohort_rows(roster,visits,final,target,days=days)
            manifest=dict(source_semantics_audited=True,archive=str(archive),archive_sha256=inspection["sha256"],
                source_url="https://public.jaeb.org/dataset/546",source_provided_by_user=True,
                participants=rows,exclusions=exclusions,external_outcomes_evaluated=False,
                rules={"clock":"DeviceDtTmDaysFromEnroll*86400 + DeviceTm in seconds; internal clock cross-check only; no clock correction or timezone inference",
                    "window":f"enrollment [0,{days}*86400), truncated at whole final dated non-missed HVisitInfo visit day (or earlier death day)",
                    "endpoint_convention":"primary elapsed 182 days" if days==182 else "prespecified inclusive integer day 182 sensitivity",
                    "includes_prerandomization_runin":True,
                    "withdrawal_date_limitation":"HPostRandPtFinalStatus publishes reason but no withdrawal date; final documented visit is not the exact withdrawal instant",
                    "record_type":"CGM only; calibration and other record types excluded",
                    "units":"mg/dL, original README HDeviceCGM.GlucoseValue",
                    "clock_barrier":"internal reset, stall, or positive increment discrepancy >60 sec when adjacent positive local gap <=1800 sec; first barrier grid bin unavailable",
                    "same_time_clock_audit":"Internal clocks grouped by device time; differing finite internal values mark that time unavailable; raw export row order cannot choose continuity",
                    "no_outlier_deletion":True,"raw_primary_only":True,
                    "grid":"5 min enrollment-relative midnight phase; latest source in (t-300,t]; exact repeats collapsed, conflicting same-time glucose unavailable"})
            # Cohort and endpoint are committed before reading CGM values.
            (target/"participant_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
            manifests.append(manifest)
        rows=manifests[1]["participants"]; bounds={r["id"]:(r["start"],r["end"]) for r in rows}
        parts={pid:[] for pid in bounds}; counts=Counter(); record_types=Counter(); per_pid={pid:Counter() for pid in bounds}
        columns=["PtID","RecID","DeviceDtTmDaysFromEnroll","DeviceTm","DexInternalDtTmDaysFromEnroll","DexInternalTm","RecordType","GlucoseValue"]
        for chunk in table("HDeviceCGM.txt",chunksize=250000,usecols=columns,dtype={"DeviceTm":"str","DexInternalTm":"str","RecordType":"str"}):
            counts["all_source_rows"]+=len(chunk);record_types.update(chunk.RecordType.fillna("MISSING").value_counts().to_dict())
            counts["non_CGM_rows"]+=int((~chunk.RecordType.eq("CGM")).sum())
            chunk=chunk[chunk.RecordType.eq("CGM")].copy()
            counts["outside_cohort_rows"]+=int((~chunk.PtID.isin(bounds)).sum())
            chunk=chunk[chunk.PtID.isin(bounds)].copy()
            chunk["time"]=relative_clock_seconds(chunk.DeviceDtTmDaysFromEnroll,chunk.DeviceTm)
            chunk["internal"]=relative_clock_seconds(chunk.DexInternalDtTmDaysFromEnroll,chunk.DexInternalTm)
            chunk["glucose"]=pd.to_numeric(chunk.GlucoseValue,errors="coerce")
            for pid,part in chunk.groupby("PtID",sort=False):
                pid=int(pid);start,end=bounds[pid]
                per_pid[pid]["source_CGM_rows"]+=len(part)
                valid=np.isfinite(part.time)
                per_pid[pid]["invalid_device_time"]+=int((~valid).sum())
                in_window=valid&(part.time>=start)&(part.time<end)
                per_pid[pid]["outside_inclusive_window"]+=int((valid&~in_window).sum())
                part=part[in_window]
                per_pid[pid]["invalid_internal_time_in_inclusive_window"]+=int((~np.isfinite(part.internal)).sum())
                per_pid[pid]["missing_or_nonnumeric_glucose_in_inclusive_window"]+=int((~np.isfinite(part.glucose)).sum())
                if len(part): parts[pid].append(part[["time","glucose","internal","RecID"]].to_numpy(dtype=float))
            print(f"replacebg_parsed_rows={counts['all_source_rows']}",flush=True)
        for target,manifest in zip(roots,manifests):
            audit=[]; hashes=[]
            for row in manifest["participants"]:
                pid=row["id"]
                a=np.concatenate(parts[pid]) if parts[pid] else np.empty((0,4))
                grid,raw=canonical_arrays(a[:,0],a[:,1],a[:,2],a[:,3],row["start"],row["end"])
                np.savez_compressed(row["file"],**grid)
                rawfile=target/"clocks"/f"REPLACE-BG_{pid}.npz"
                np.savez_compressed(rawfile,**raw)
                audit.append(dict(id=pid,split=row["split"],arm=row["arm"],**per_pid[pid],
                    raw_in_window=len(raw["time"]),duplicate_timestamps=len(raw["time"])-len(np.unique(raw["time"])),
                    conflict_grid_bins=int(grid["conflict"].sum()),clock_barrier_bins=int(grid["clock_barrier"].sum()),
                    ambiguous_internal_clock_times=int((raw["barrier_reasons"]=="same_time_internal_disagreement").sum()),
                    observed_grid_bins=int(np.isfinite(grid["glucose"]).sum()),calendar_grid_bins=len(grid["time"])))
                hashes.extend([dict(file=str(path),sha256=file_hash(path)) for path in [Path(row["file"]),rawfile]])
            pd.DataFrame(audit).to_csv(target/"source_time_audit.csv",index=False)
            manifest["audit"]=dict(**counts,record_types=dict(record_types),participants=len(manifest["participants"]),
                original_roster_participants=len(roster),dropped_participants=int(roster.PtStatus.eq("Dropped").sum()),
                schema="original HDeviceCGM, HPtRoster, HVisitInfo, HPostRandPtFinalStatus and README",elapsed_seconds=time.monotonic()-begun)
            (target/"participant_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
            (target/"canonical_file_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
        inspection["semantics_verified"]=True
        inspection["adapter_sha256"]=file_hash(Path(__file__))
        inspection["cohort_endpoint_frozen_before_CGM_scan"]=True
        (out/"archive_schema_audit.json").write_text(json.dumps(inspection,indent=2)+"\n",encoding="utf-8")
    return manifests[0]


def rebuild_canonical_from_preserved_raw(out):
    """Apply a pre-outcome clock audit repair without rereading the large ZIP."""
    out=Path(out).resolve()
    for target in [out,out/"inclusive_day182"]:
        manifest=json.loads((target/"participant_manifest.json").read_text())
        if manifest.get("external_outcomes_evaluated"):
            raise ValueError("Pre-outcome repair cannot relabel an evaluated cohort")
        audit=pd.read_csv(target/"source_time_audit.csv").set_index("id");hashes=[]
        for row in manifest["participants"]:
            rawfile=target/"clocks"/f"REPLACE-BG_{row['id']}.npz"
            with np.load(rawfile) as z:
                grid,raw=canonical_arrays(z["time"],z["glucose"],z["transmitter"],z["recid"],row["start"],row["end"])
            np.savez_compressed(rawfile,**raw);np.savez_compressed(row["file"],**grid)
            audit.loc[row["id"],"clock_barrier_bins"]=int(grid["clock_barrier"].sum())
            audit.loc[row["id"],"observed_grid_bins"]=int(np.isfinite(grid["glucose"]).sum())
            audit.loc[row["id"],"ambiguous_internal_clock_times"]=int((raw["barrier_reasons"]=="same_time_internal_disagreement").sum())
            hashes.extend([dict(file=str(p),sha256=file_hash(p)) for p in [Path(row["file"]),rawfile]])
        audit.reset_index().to_csv(target/"source_time_audit.csv",index=False)
        manifest["rules"]["same_time_clock_audit"]="Internal clocks grouped by device time; differing finite internal values mark that time unavailable; raw export row order cannot choose continuity"
        (target/"participant_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n",encoding="utf-8")
        (target/"canonical_file_hashes.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
    inspection=json.loads((out/"archive_schema_audit.json").read_text())
    inspection["adapter_sha256"]=file_hash(Path(__file__))
    inspection["pre_outcome_clock_grouping_repair"]=True
    (out/"archive_schema_audit.json").write_text(json.dumps(inspection,indent=2)+"\n",encoding="utf-8")
    return {"status":"canonical_clock_audit_complete","external_outcomes_evaluated":False}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--archive",type=Path,required=True)
    parser.add_argument("--output",type=Path,required=True)
    parser.add_argument("--inspect-only",action="store_true")
    parser.add_argument("--rebuild-clocks",action="store_true")
    args=parser.parse_args()
    result=rebuild_canonical_from_preserved_raw(args.output) if args.rebuild_clocks else (inspect_archive(args.archive,args.output) if args.inspect_only else prepare_replacebg(args.archive,args.output))
    print(json.dumps({k:v for k,v in result.items() if k not in ["participants","inventory"]},indent=2))
