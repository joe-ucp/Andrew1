"""Audited original-archive CLVer adapter for the pinned public release.

Literal encodings and excluded measurements are preserved. Protocol-based date
and implicit-unit assumptions are explicit in the source and clinical ledgers.
"""
from __future__ import annotations
import sys
sys.dont_write_bytecode = True
import argparse
import csv
import io
import json
import re
import zipfile
from decimal import Decimal, InvalidOperation
from pathlib import Path

import numpy as np
import pandas as pd
from card003_core import EXPECTED_CLVER_SHA256, sha256_file
from card003_core import derive_mmt_visit,clinical_remission,later_interval,overlapping_pairs,independent_episode_labeler

REQUIRED_COLUMNS = {
    "PtRoster.txt": {"PtID","TrtGroup","Cohort","DeviceTrtGroup","AgeAsofEnrollDt","RandDt"},
    "SampleResults.txt": {"PtID","Visit","CollectionDt","CollectionTime","ResultName","Value","Units","STAResultStatus","STASampleFastStatus","TempIndicator"},
    "InsulinUse.txt": {"PtID","Visit","ParentLoginVisitID","InsDoseInfoDt","TDIUnits","TDIUnitsUnk"},
    "DiabPhysExam.txt": {"PtID","Visit","ParentLoginVisitID","Weight","WeightUnits","WeightUnk"},
    "VisitInfo.txt": {"PtID","Visit","ParentLoginVisitID","VisitDt","VisitMiss"},
    "DeviceCGM.txt": {"PtID","DeviceDtTm","RecordType","Value","Units","TransTime","Unusable","UnusableReason","RecID"},
    "DevicePumpCGMDATAGXB.txt": {"PtID","DeviceDtTm","CGMValue","HighLowIndicator","Unusable","UnusableReason","RecID"},
}
DOMAIN_FIELDS = {"TrtGroup","Cohort","DeviceTrtGroup","Visit","Units","STAResultStatus",
    "STASampleFastStatus","TempIndicator","WeightUnits","WeightUnk","TDIUnitsUnk","VisitMiss",
    "RecordType","Unusable","UnusableReason","HighLowIndicator","ResultName"}


class ArchiveTables:
    """No archive extraction: members are read by unique basename as streams."""
    def __init__(self, path):
        self.path = Path(path)
        self.archive = zipfile.ZipFile(self.path)

    def __enter__(self):
        return self

    def __exit__(self, *_):
        self.archive.close()

    def resolve(self, basename):
        members = [n for n in self.archive.namelist() if n.replace("\\","/").split("/")[-1] == basename]
        if len(members) != 1:
            raise ValueError(f"Expected unique archive member {basename}; found {len(members)}")
        return members[0]

    def format(self, basename):
        member = self.resolve(basename)
        with self.archive.open(member) as f:
            bom = f.read(4)
        if bom.startswith((b"\xff\xfe\x00\x00",b"\x00\x00\xfe\xff")):
            raise ValueError("Unexpected UTF-32 source; inspect and explicitly support before parsing")
        encoding = "utf-16" if bom.startswith((b"\xff\xfe",b"\xfe\xff")) else "utf-8-sig"
        with self.archive.open(member) as f:
            header = io.TextIOWrapper(f,encoding=encoding).readline().rstrip("\r\n")
        counts = {d:header.count(d) for d in ("|","\t",",")}
        delimiter = max(counts,key=counts.get)
        if counts[delimiter] == 0:
            raise ValueError(f"No recognized multi-column delimiter in {member}")
        columns = next(csv.reader([header],delimiter=delimiter))
        if len(set(columns)) != len(columns):
            raise ValueError(f"Duplicate column name in {member}")
        return dict(member=member,encoding=encoding,delimiter=delimiter,columns=columns)

    def read(self, basename, *, columns=None, nrows=None):
        spec = self.format(basename)
        if columns is not None and not set(columns).issubset(spec["columns"]):
            raise ValueError(f"Missing requested columns for {basename}")
        with self.archive.open(spec["member"]) as f:
            return pd.read_csv(f,sep=spec["delimiter"],encoding=spec["encoding"],
                dtype=str,keep_default_na=False,na_filter=False,usecols=columns,nrows=nrows)

    def chunks(self, basename, *, columns=None, chunksize=250000):
        spec = self.format(basename)
        if columns is not None and not set(columns).issubset(spec["columns"]):
            raise ValueError(f"Missing requested columns for {basename}")
        with self.archive.open(spec["member"]) as f:
            for chunk in pd.read_csv(f,sep=spec["delimiter"],encoding=spec["encoding"],
                dtype=str,keep_default_na=False,na_filter=False,usecols=columns,chunksize=chunksize):
                yield chunk

    def inventory(self):
        return [{"member":x.filename,"uncompressed_bytes":x.file_size,"compressed_bytes":x.compress_size}
                for x in self.archive.infolist()]


def canonical_id(value):
    text = str(value).strip()
    try:
        d = Decimal(text)
    except InvalidOperation as exc:
        raise ValueError(f"Invalid participant ID encoding: {text!r}") from exc
    if not d.is_finite() or d != d.to_integral_value() or d < 0:
        raise ValueError(f"Nonintegral/negative/unknown participant ID: {text!r}")
    return str(int(d))


def numeric(values):
    """Only literal finite numbers; inequality assay values are not imputed."""
    out = pd.to_numeric(values,errors="coerce").astype(float)
    return out.where(np.isfinite(out))


def source_dates(values, *, date_format=None):
    """No numeric-to-epoch coercion or guessed study-day origin."""
    s = pd.Series(values,dtype="string").str.strip()
    if date_format is None and s.str.fullmatch(r"[+-]?\d+(?:\.\d+)?",na=False).any():
        raise ValueError("Numeric date encoding needs an explicitly documented format/origin")
    return pd.to_datetime(s,format=date_format or "mixed",errors="coerce")


def measurement_times(values):
    """A CGM measurement needs an explicit time; date-only is not midnight.

    The pinned archive uses month/day/year and a 12-hour clock with seconds.
    Date-only export strings have plausible midnight neighbors but no explicit
    time-of-day field, so their primary measurement time remains unknown.
    """
    return source_dates(values,date_format="%m/%d/%Y %I:%M:%S %p")


def collapse_and_flag(frame, key, semantic_columns):
    """Analytical duplicate collapse, with all conflicting keys excluded.

    Identity/provenance columns such as RecID must not make otherwise identical
    analytical rows appear to be distinct measurements. The caller explicitly
    supplies the measurement, date, unit and quality columns that define identity.
    """
    key, semantic_columns = list(key),list(semantic_columns)
    missing = (set(key)|set(semantic_columns)) - set(frame.columns)
    if missing:
        raise ValueError(f"Missing duplicate-audit fields: {sorted(missing)}")
    cols = list(dict.fromkeys(key+semantic_columns))
    unique = frame.drop_duplicates(cols).copy()
    conflict_mask = unique.duplicated(key,keep=False)
    conflicts = unique.loc[conflict_mask,key].drop_duplicates()
    valid = unique.loc[~conflict_mask].copy()
    return valid,conflicts,dict(raw_rows=len(frame),analytical_duplicate_rows=len(frame)-len(unique),
                               conflicting_keys=len(conflicts),retained_rows=len(valid))


def inspect_archive(path, output):
    """Schema/domain inspection; never marks the source adapter validated."""
    path, output = Path(path),Path(output)
    output.mkdir(parents=True,exist_ok=True)
    digest = sha256_file(path)
    result = dict(archive=str(path),bytes=path.stat().st_size,sha256=digest,
        expected_sha256=EXPECTED_CLVER_SHA256,source_version_matches=digest==EXPECTED_CLVER_SHA256,
        schema_validated=False,clinical_derivations_executed=False,tables={})
    with ArchiveTables(path) as tables:
        result["inventory"] = tables.inventory()
        for name, expected in REQUIRED_COLUMNS.items():
            try:
                spec = tables.format(name)
                info = dict(spec,missing_expected_columns=sorted(expected-set(spec["columns"])))
                # The full small clinical tables can be audited without sampling.
                # CGM inspection is deliberately bounded; later complete stream
                # audits must count all source flags, conflicts and clock problems.
                limit = 10000 if name.startswith("Device") else None
                frame = tables.read(name,nrows=limit)
                info["inspected_rows"] = len(frame)
                info["sample_only"] = limit is not None
                info["domains"] = {c:{str(k):int(v) for k,v in frame[c].value_counts(dropna=False).items()}
                    for c in frame.columns if c in DOMAIN_FIELDS}
                info["blank_counts"] = {c:int(frame[c].str.strip().eq("").sum()) for c in frame.columns}
                if name == "PtRoster.txt" and "PtID" in frame:
                    info["unique_participant_ids"] = int(frame.PtID.nunique())
                    info["duplicate_participant_id_rows"] = int(frame.PtID.duplicated(keep=False).sum())
                result["tables"][name] = info
            except Exception as exc:
                result["tables"][name] = {"inspection_error":type(exc).__name__+": "+str(exc)}
    target = output / "schema_inventory.json"
    target.write_text(json.dumps(result,indent=2,ensure_ascii=False),encoding="utf-8")
    return result


def _visit_number(value):
    match = re.fullmatch(r"(13|26|39|52) Week",str(value).strip())
    return int(match[1]) if match else None


def _unique_row(frame, columns):
    dedup = frame.drop_duplicates(columns)
    return dedup.iloc[0] if len(dedup) == 1 else None


def reconstruct_clinical(path, root):
    """Validated CLVer release adapter; derived rows use the specified output root."""
    root=Path(root); out=root/"processed";out.mkdir(exist_ok=True)
    if sha256_file(Path(path)) != EXPECTED_CLVER_SHA256:
        raise ValueError("CLVer release hash differs from audited version")
    with ArchiveTables(path) as z:
        raw={name:z.read(name) for name in ["PtRoster.txt","SampleResults.txt","InsulinUse.txt","DiabPhysExam.txt","VisitInfo.txt"]}
    for frame in raw.values():
        frame["PtID"]=frame.PtID.map(canonical_id)
    roster=raw["PtRoster.txt"]
    randomized=roster[roster.RandDt.ne("") & roster.TrtGroup.isin(["A.HCL","A.Usual Care","B.HCL","B.Usual Care"])].copy()
    if randomized.PtID.duplicated().any():raise ValueError("Duplicate randomized roster identity")
    randomized["age"]=numeric(randomized.AgeAsofEnrollDt)
    randomized["delivery_arm"]=np.where(randomized.TrtGroup.str.endswith("Usual Care"),"usual care","intensive")
    randomized["management"]=np.where(randomized.delivery_arm.eq("usual care"),"usual care",
        np.where(randomized.DeviceTrtGroup.eq("Tandem"),"intensive/Tandem",np.where(randomized.DeviceTrtGroup.eq("MedTronic"),"intensive/Medtronic","unknown")))
    randomized["source"]=np.where(randomized.DeviceTrtGroup.eq("Tandem"),"GXB","DeviceCGM")
    randomized["cohort"]=randomized.Cohort
    randomized.to_csv(out/"randomized_roster.csv",index=False)
    lab=raw["SampleResults.txt"];ins=raw["InsulinUse.txt"];weight=raw["DiabPhysExam.txt"];vis=raw["VisitInfo.txt"]
    for f in [lab,ins,weight,vis]:f["week"]=f.Visit.map(_visit_number)
    lab["date"]=source_dates(lab.CollectionDt)
    lab["value"]=numeric(lab.Value)
    vis["date"]=source_dates(vis.VisitDt)
    ins["date"]=source_dates(ins.InsDoseInfoDt)
    clinical=[]
    for pid in randomized.PtID:
        for week in [13,26,39,52]:
            sl=lab[(lab.PtID==pid)&(lab.week==week)]
            cp=sl[sl.ResultName.str.startswith("C-PEP-")].copy()
            cp["minute"]=cp.ResultName.str.extract(r"C-PEP-(\d+)$",expand=False).astype(float)
            samples=cp.rename(columns={"Units":"units","STAResultStatus":"status","STASampleFastStatus":"fasting","TempIndicator":"temperature"})
            assay=derive_mmt_visit(samples[["minute","value","units","date","status","fasting","temperature"]])
            # Source audit shows visit-level fasting status repeated on samples;
            # explicitly discordant nonbaseline status cannot establish fasting.
            if len(cp) and not cp.STASampleFastStatus.eq("Fasting").all():
                assay.update(eligible=False,reason="fasting_not_established_on_all_assay_rows",auc_nmol_l=np.nan,peak_nmol_l=np.nan)
            row=dict(PtID=pid,index_visit=week,**assay)
            base=cp[cp.minute.eq(0)]
            anchors=base.date.dropna().dt.normalize().unique()
            row["anchor_date"]=pd.Timestamp(anchors[0]).isoformat() if len(anchors)==1 and base.date.notna().all() else None
            vv=vis[(vis.PtID==pid)&(vis.week==week)&vis.VisitMiss.ne("1")]
            vr=_unique_row(vv,["ParentLoginVisitID","date"])
            row["visit_date"]=vr.date.isoformat() if vr is not None and pd.notna(vr.date) else None
            row["visit_key"]=vr.ParentLoginVisitID if vr is not None else None
            remission=dict(idaa1c=np.nan,remission=np.nan,insulin_u_kg_day=np.nan,weight_kg=np.nan,reason="missing_or_conflicting_clinical_join")
            hba=sl[sl.ResultName.eq("GLYHB") & sl.STAResultStatus.eq("Complete") & sl.Units.eq("%")]
            hr=_unique_row(hba,["value","date","Units","STAResultStatus"])
            row["hba1c"]=float(hr.value) if hr is not None and np.isfinite(hr.value) else np.nan
            if vr is not None:
                ii=ins[(ins.PtID==pid)&(ins.week==week)&(ins.ParentLoginVisitID==vr.ParentLoginVisitID)]
                ww=weight[(weight.PtID==pid)&(weight.week==week)&(weight.ParentLoginVisitID==vr.ParentLoginVisitID)]
                ir=_unique_row(ii,["TDIUnits","TDIUnitsUnk","date"])
                wr=_unique_row(ww,["Weight","WeightUnits","WeightUnk"])
                if ir is not None and wr is not None and hr is not None:
                    remission=clinical_remission(hr.value,ir.TDIUnits,wr.Weight,wr.WeightUnits,
                        row["anchor_date"],hr.date,ir.date,vr.date,
                        dose_known=ir.TDIUnitsUnk!="1",weight_known=wr.WeightUnk!="1")
                    row.update(insulin_date=ir.date.isoformat() if pd.notna(ir.date) else None,
                               hba1c_date=hr.date.isoformat() if pd.notna(hr.date) else None,
                               exam_date_proxy=vr.date.isoformat() if pd.notna(vr.date) else None,
                               daily_insulin_units=ir.TDIUnits,source_weight=wr.Weight,source_weight_units=wr.WeightUnits)
            row.update({("remission_reason" if k=="reason" else k):v for k,v in remission.items()})
            row["exam_date_provenance"]="VisitInfo.VisitDt through participant+ParentLoginVisitID; visit-level proxy per CRF/protocol"
            clinical.append(row)
    cf=pd.DataFrame(clinical)
    cf.to_csv(out/"clinical_derivations.csv",index=False)
    # Preserve source assay values, including excluded/non-numeric ones.
    lab[lab.PtID.isin(randomized.PtID)&lab.week.notna()].to_csv(out/"source_assay_ledger.csv",index=False)
    pairs=[]
    for r in randomized.to_dict("records"):
        pid=r["PtID"];person=cf[cf.PtID.eq(pid)].set_index("index_visit")
        for index,target in [(13,26),(26,39),(39,52)]:
            a=person.loc[index];b=person.loc[target]
            row={**r,**a.to_dict(),"PtID":pid,"index_visit":index,"target_visit":target,"pair_id":f"{pid}_{target}","target_anchor":b.anchor_date}
            interval=later_interval(a.anchor_date,[b.anchor_date] if b.anchor_date else [],b.visit_date)
            row["interval_eligible"]=bool(a.eligible and interval["eligible"])
            row["interval_reason"]=a.reason if not a.eligible else interval["reason"]
            for k in ["start","end","lag_days"]:
                row[k]=interval.get(k)
            row["contemporaneous_auc_nmol_l"]=b.auc_nmol_l
            row["contemporaneous_assay_eligible"]=b.eligible
            row["contemporaneous_remission"]=b.remission
            pairs.append(row)
    pf=pd.DataFrame(pairs)
    for pid,group in pf[pf.interval_eligible].groupby("PtID"):
        bad=overlapping_pairs(list(zip(group.start,group.end)))
        for i in bad:
            ix=group.index[i];pf.loc[ix,"interval_eligible"]=False;pf.loc[ix,"interval_reason"]="overlapping_outcome_intervals"
    pf.to_csv(out/"interval_derivations.csv",index=False)
    summary=dict(roster_total=len(roster),randomized=len(randomized),assay_visits=len(cf),eligible_assay_visits=int(cf.eligible.sum()),
        attempted_pairs=len(pf),temporally_eligible_pairs=int(pf.interval_eligible.sum()),
        eligible_children=int(pf.loc[pf.interval_eligible,"PtID"].nunique()),
        assay_exclusion_reasons=cf.reason.value_counts().to_dict(),interval_exclusion_reasons=pf.interval_reason.value_counts().to_dict(),
        known_remission_visits=int(cf.remission.notna().sum()),remission_visits=int(cf.remission.eq(1).sum()),
        exam_date_limitation="Same-visit VisitInfo date is a documented protocol-level proxy; no independent weight timestamp is available.")
    (root/"results"/"clinical_audit.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return randomized,cf,pf,summary


def reconstruct_sources(path,root,roster):
    """Audit chosen sources and preserve ungridded arrays for every sensitivity."""
    root=Path(root);folder=root/"processed"/"raw_streams";folder.mkdir(exist_ok=True)
    choices=dict(zip(roster.PtID,roster.source));parts={pid:[] for pid in choices}
    counters={};exclusions=[]
    with ArchiveTables(path) as z:
        for table,source in [("DeviceCGM.txt","DeviceCGM"),("DevicePumpCGMDATAGXB.txt","GXB")]:
            counts={"raw_rows":0,"chosen_source_rows":0,"unusable":0,"flagged_placeholders":0,"nonnumeric":0,"invalid_time":0,"unsupported_units":0,"unit_domain":{}}
            for chunk in z.chunks(table):
                counts["raw_rows"]+=len(chunk)
                chunk=chunk[chunk.PtID.isin(choices)].copy()
                chunk=chunk[chunk.PtID.map(choices).eq(source)].copy()
                counts["chosen_source_rows"]+=len(chunk)
                if source=="DeviceCGM":
                    chunk=chunk[chunk.RecordType.eq("CGM")].copy()
                    for k,v in chunk.Units.value_counts().items():counts["unit_domain"][k]=counts["unit_domain"].get(k,0)+int(v)
                    unitok=chunk.Units.isin(["","mg/dL"])
                    counts["unsupported_units"]+=int((~unitok).sum());chunk=chunk[unitok].copy()
                    glucose=numeric(chunk.Value);trans=numeric(chunk.TransTime);explicit=chunk.Units.eq("mg/dL")
                else:
                    flags=chunk.HighLowIndicator.ne("0");counts["flagged_placeholders"]+=int(flags.sum())
                    chunk=chunk[~flags].copy();glucose=numeric(chunk.CGMValue)
                    trans=pd.Series(np.nan,index=chunk.index);explicit=pd.Series(True,index=chunk.index)
                usable=chunk.Unusable.eq("False");counts["unusable"]+=int((~usable).sum())
                dates=measurement_times(chunk.DeviceDtTm)
                counts["nonnumeric"]+=int(glucose.isna().sum());counts["invalid_time"]+=int(dates.isna().sum())
                for ix in chunk.index[dates.isna()]:
                    reason="missing_time_of_day" if ":" not in chunk.loc[ix,"DeviceDtTm"] else "unparsed_timestamp"
                    exclusions.append(dict(source=source,PtID=chunk.loc[ix,"PtID"],RecID=chunk.loc[ix,"RecID"],reason=reason,value=chunk.loc[ix,"DeviceDtTm"]))
                # Retain explicitly labeled date-only midnight candidates for a
                # separate sensitivity, never as primary timed measurements.
                explicit_time=dates.notna()
                date_only=chunk.DeviceDtTm.str.fullmatch(r"\d{1,2}/\d{1,2}/\d{4}")
                nominal_dates=dates.copy()
                nominal_dates.loc[date_only]=source_dates(chunk.loc[date_only,"DeviceDtTm"])
                good=usable&glucose.notna()&nominal_dates.notna()
                selected=chunk[good].copy()
                selected["time"]=nominal_dates[good].astype("datetime64[ns]").astype("int64")/1e9
                selected["glucose"]=glucose[good];selected["transmitter"]=trans[good];selected["explicit_unit"]=explicit[good].astype(int)
                selected["explicit_time"]=explicit_time[good].astype(int)
                for pid,group in selected.groupby("PtID",sort=False):
                    parts[pid].append(group[["time","glucose","transmitter","explicit_unit","RecID","explicit_time"]].to_numpy(dtype=float))
            counters[source]=counts
    audits=[]
    for pid in choices:
        a=np.concatenate(parts[pid]) if parts[pid] else np.empty((0,6))
        if len(a):a=a[np.argsort(a[:,0],kind="stable")]
        def barriers(values):
            dt=np.diff(values[:,0]);dtx=np.diff(values[:,2])
            barrier=np.isfinite(dtx)&(dt>0)&(dt<=1800)&((dtx<0)|((dtx>0)&(np.abs(dt-dtx)>60)))
            return values[1:,0][barrier]
        primary=a[a[:,5].astype(bool)]
        barrier_times=barriers(primary);midnight_barriers=barriers(a)
        unique,first,count=np.unique(a[:,0],return_index=True,return_counts=True)
        conflicts=sum(len(set(a[i:i+n,1]))>1 for i,n in zip(first[count>1],count[count>1]))
        counter_conflicts=[a[i,0] for i,n in zip(first[count>1],count[count>1]) if len(set(a[i:i+n,2][np.isfinite(a[i:i+n,2])]))>1]
        np.savez_compressed(folder/f"{pid}.npz",time=a[:,0],glucose=a[:,1],transmitter=a[:,2],explicit_unit=a[:,3].astype(bool),recid=a[:,4],explicit_time=a[:,5].astype(bool),barrier_times=barrier_times,midnight_barrier_times=midnight_barriers,counter_conflict_times=np.asarray(counter_conflicts))
        audits.append(dict(PtID=pid,source=choices[pid],valid_source_rows=len(primary),date_only_sensitivity_rows=len(a)-len(primary),duplicate_timestamp_rows=len(a)-len(unique),counter_conflict_groups=len(counter_conflicts),
            conflicting_same_time_groups=conflicts,clock_barriers=len(barrier_times),transmitter_observed_rows=int(np.isfinite(a[:,2]).sum()),
            glucose_min=float(a[:,1].min()) if len(a) else None,glucose_max=float(a[:,1].max()) if len(a) else None))
    pd.DataFrame(audits).to_csv(root/"results"/"source_time_audit.csv",index=False)
    pd.DataFrame(exclusions).to_csv(root/"results"/"source_exclusion_ledger.csv",index=False)
    payload=dict(counts=counters,clock_rule="No invented offset; insert unknown bin at transmitter resets/incompatible positive increments >60 seconds within 30-minute adjacent gaps, retaining subsequent observed source segment. No transmitter clock exists in GXB.",
        implicit_units="Blank DeviceCGM units interpreted mg/dL based on source value domain/protocol; repeat excluding blanks. GXB unit is implicit mg/dL by source/protocol, not an explicit Units field.",
        source_rule="Tandem roster uses GXB only; all others use DeviceCGM only",
        timestamp_rule="Primary requires explicit time-of-day; date-only records retained only in labeled midnight-assumption sensitivity. Duplicate device times with conflicting finite transmitter counters are unavailable.")
    (root/"results"/"source_audit.json").write_text(json.dumps(payload,indent=2),encoding="utf-8")
    return audits,payload


def window_grid(root,pid,start,end,phase_minutes=0,delay_minutes=0,exclude_blank=False,clock_barriers=True,include_date_only=False,lookback_minutes=160):
    from card001_core import causal_grid
    raw=np.load(Path(root)/"processed"/"raw_streams"/f"{pid}.npz")
    t=raw["time"];g=raw["glucose"]
    keep=(t>=start-(lookback_minutes+60)*60)&(t<=end+2*3600)
    if not include_date_only:keep &= raw["explicit_time"]
    if exclude_blank:keep &= raw["explicit_unit"]
    if keep.any():
        grid=causal_grid(t[keep],g[keep],start-lookback_minutes*60,end+50*60,phase_seconds=phase_minutes*60,delay_seconds=delay_minutes*60)
    else:
        phase=phase_minutes*60;tt=np.arange(np.ceil((start-lookback_minutes*60-phase)/300)*300+phase,np.floor((end+3000-phase)/300)*300+phase+1,300)
        grid=dict(time=tt,glucose=np.full(len(tt),np.nan),source_time=np.full(len(tt),np.nan),conflict=np.zeros(len(tt),bool))
    barrier_key="midnight_barrier_times" if include_date_only else "barrier_times"
    barrier=np.isin(grid["time"],np.ceil((raw[barrier_key]+delay_minutes*60-phase_minutes*60)/300)*300+phase_minutes*60)
    grid["clock_barrier"]=barrier
    grid["counter_conflict"]=np.isin(grid["source_time"],raw["counter_conflict_times"])
    grid["glucose"][grid["counter_conflict"]]=np.nan
    if clock_barriers:grid["glucose"][barrier]=np.nan
    return grid


def construct_intervals(root,pairs,phase_minutes=0,delay_minutes=0,exclude_blank=False,run_length=4,save_grids=True):
    from card001_core import episodes
    root=Path(root);folder=root/"processed"/"interval_grids";folder.mkdir(exist_ok=True)
    rows=[];events=[]
    for row in pairs[pairs.interval_eligible].to_dict("records"):
        start=pd.Timestamp(row["start"]).timestamp();end=pd.Timestamp(row["end"]).timestamp()
        grid=window_grid(root,row["PtID"],start,end,phase_minutes,delay_minutes,exclude_blank)
        g=grid["glucose"];t=grid["time"];inside=(t>=start)&(t<end)
        onset,censored=episodes(g,run_length)
        aa,bb=independent_episode_labeler(g,run_length)
        if not np.array_equal(aa,onset) or not np.array_equal(bb,censored):raise AssertionError("Independent event label mismatch")
        selected=onset[(t[onset]>=start)&(t[onset]<end)]
        row.update(start_seconds=start,end_seconds=end,observed_bins=int(np.isfinite(g[inside]).sum()),grid_bins=int(inside.sum()),
            observed_days=float(np.isfinite(g[inside]).sum()/288),coverage=float(np.isfinite(g[inside]).mean()),events=len(selected),
            censored_runs=int(np.sum((t[censored]>=start)&(t[censored]<end))),conflict_bins=int(grid["conflict"][inside].sum()),
            clock_barrier_bins=int(grid["clock_barrier"][inside].sum()))
        row["counter_conflict_bins"]=int(grid["counter_conflict"][inside].sum())
        if save_grids:
            path=folder/f"{row['pair_id']}.npz";np.savez_compressed(path,**grid);row["grid_path"]=str(path)
        rows.append(row)
        events += [dict(PtID=row["PtID"],pair_id=row["pair_id"],onset_seconds=float(t[e])) for e in selected]
    frame=pd.DataFrame(rows)
    if save_grids:
        frame.to_csv(root/"processed"/"intervals_all_coverage.csv",index=False)
        primary=frame[frame.coverage>=.7]
        primary.to_csv(root/"processed"/"intervals_primary.csv",index=False)
        pd.DataFrame(events).to_csv(root/"results"/"independent_event_ledger.csv",index=False)
        summary=dict(intervals=len(frame),children=int(frame.PtID.nunique()),primary_intervals=len(primary),primary_children=int(primary.PtID.nunique()),
            primary_events=int(primary.events.sum()),primary_event_children=int(primary.loc[primary.events>0,"PtID"].nunique()),
            primary_days=float(primary.observed_days.sum()),primary_known_remission=int(primary.remission.notna().sum()),primary_remission=int(primary.remission.eq(1).sum()))
        (root/"results"/"interval_audit.json").write_text(json.dumps(summary,indent=2),encoding="utf-8")
    return frame


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--archive",type=Path,required=True)
    ap.add_argument("--output",type=Path,required=True)
    ap.add_argument("--prepare",action="store_true",help="Reconstruct actual clinical and CGM intervals using the audited CLVer release")
    args = ap.parse_args()
    result = inspect_archive(args.archive,args.output)
    if args.prepare:
        root=args.output.parent if args.output.name=="processed" else args.output
        roster,clinical,pairs,audit=reconstruct_clinical(args.archive,root)
        print(json.dumps(audit),flush=True)
        reconstruct_sources(args.archive,root,roster)
        frame=construct_intervals(root,pairs)
        print((root/"results"/"interval_audit.json").read_text(),flush=True)
    print(json.dumps({"sha256":result["sha256"],"source_version_matches":result["source_version_matches"],
        "schema_inventory_only":not args.prepare,"clinical_derivations_executed":args.prepare,
        "output":str(args.output/"schema_inventory.json")}))


if __name__=="__main__": main()
