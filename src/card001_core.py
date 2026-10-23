"""Prespecified recorded-time CGM core, shared by Cards 001 and 003.

All indices are FIVE MINUTE bins. No outcome interpolation is performed.
Synthetic fixtures validate implementation only; they are never efficacy evidence.
"""
from __future__ import annotations
import hashlib
import warnings
import numpy as np

SEED = 20261022
EXPOSED = {1,3,4,5,6,7,8,10,12,13,14,15,19,20,22,32,36,39,50,51,53,100,101,102,103,104,105,106,107,108,109}


def participant_split(ptid):
    ptid = str(int(ptid))
    if int(ptid) in EXPOSED:
        return "train", "generation_exposed"
    bucket = int(hashlib.sha256(("t1d-card-001:" + ptid).encode()).hexdigest()[:8],16) % 10000
    return ("train" if bucket < 6000 else "validation" if bucket < 8000 else "test"), str(bucket)


def causal_grid(times, glucose, start=None, end=None, phase_seconds=0, delay_seconds=0):
    """Latest reading in (t-300,t], never later; conflicts stay unavailable.

    times/start/end are numeric seconds in a consistent origin. Input is one
    participant and one reconciled clock segment. Boundaries are inclusive.
    Delay shifts availability, leaving original source_time intact. A fixed
    epoch/day anchor (rather than first retained sample) makes prefix invariance
    meaningful. NaN glucose is retained as unavailable, not silently dropped.
    """
    ts = np.asarray(times, dtype=float)
    g = np.asarray(glucose, dtype=float)
    if ts.shape != g.shape or ts.ndim != 1:
        raise ValueError("times and glucose must be equal one-dimensional arrays")
    valid = np.isfinite(ts)
    ts, g = ts[valid], g[valid]
    if not len(ts):
        return {k: np.array([],dtype=bool if k=="conflict" else float) for k in ["time","glucose","source_time","conflict"]}
    order = np.argsort(ts,kind="stable")
    ts, g = ts[order], g[order]
    unique, first, counts = np.unique(ts, return_index=True, return_counts=True)
    values = g[first].copy()
    conflict = np.zeros(len(unique),dtype=bool)
    for j in np.flatnonzero(counts > 1):
        v = g[first[j]:first[j]+counts[j]]
        # Numeric vs unavailable at same clock position is also conflict.
        if not np.all((v==v[0]) | (np.isnan(v)&np.isnan(v[0]))):
            values[j] = np.nan
            conflict[j] = True
    available = unique + delay_seconds
    if start is None: start = available[0]
    if end is None: end = available[-1]
    first_bin = np.ceil((start-phase_seconds)/300)*300+phase_seconds
    last_bin = np.floor((end-phase_seconds)/300)*300+phase_seconds
    grid = np.arange(first_bin,last_bin+1,300,dtype=float)
    idx = np.searchsorted(available,grid,side="right")-1
    legal = (idx>=0)
    clipped = idx.clip(0)
    legal &= (available[clipped] > grid-300)
    out = np.full(len(grid),np.nan)
    source = np.full(len(grid),np.nan)
    conflicts = np.zeros(len(grid),bool)
    out[legal] = values[idx[legal]]
    source[legal] = unique[idx[legal]]
    conflicts[legal] = conflict[idx[legal]]
    return dict(time=grid,glucose=out,source_time=source,conflict=conflicts)


def ready(glucose,n_bins):
    observed=np.isfinite(glucose)
    if n_bins<1: raise ValueError("n_bins must be positive")
    if len(observed)==0: return np.array([],dtype=bool)
    count=np.convolve(observed.astype(np.int16),np.ones(n_bins,dtype=np.int16),mode="full")[:len(observed)]
    return count==n_bins


def episodes(glucose,run_length=4):
    """Stateful endpoint: recovery arms; missingness resets to unknown.

    censored contains sustained low-run starts whose onset cannot be identified.
    Confirmed events remain recorded if a later missing bin censors duration.
    """
    g=np.asarray(glucose,float)
    onsets=[]; censored=[]
    state="unknown"; high=0; low=0
    for i,value in enumerate(g):
        if not np.isfinite(value):
            state="unknown"; high=low=0
        elif value>=54:
            high+=1; low=0
            if high>=run_length: state="armed"
        else:
            low+=1; high=0
            if low==run_length:
                if state=="armed": onsets.append(i-run_length+1)
                elif state=="unknown": censored.append(i-run_length+1)
                state="event"
    return np.array(onsets,dtype=np.int64),np.array(censored,dtype=np.int64)


def episodes_independent(glucose,run_length=4):
    """Independent run-encoding implementation with no shared detector state."""
    x=np.asarray(glucose,float)
    if len(x)==0: return np.array([],dtype=np.int64),np.array([],dtype=np.int64)
    category=np.where(~np.isfinite(x),-1,np.where(x>=54,1,0))
    edges=np.r_[0,np.flatnonzero(np.diff(category))+1,len(x)]
    latest_recovery=-1; latest_barrier=-1; barrier_kind="gap"; events=[]; unknown=[]
    for start,stop in zip(edges[:-1],edges[1:]):
        kind=category[start]
        if kind==-1: latest_barrier=stop-1; barrier_kind="gap"
        elif kind==1 and stop-start>=run_length: latest_recovery=start+run_length-1
        elif kind==0 and stop-start>=run_length:
            if latest_recovery>latest_barrier: events.append(start)
            elif barrier_kind=="gap": unknown.append(start)
            latest_barrier=start+run_length-1; barrier_kind="event"
    return np.asarray(events,dtype=np.int64),np.asarray(unknown,dtype=np.int64)


def labels(glucose,onsets=None):
    """Future is used only for target construction. Unknown negatives are NaN."""
    x=np.asarray(glucose,float); n=len(x)
    if onsets is None: onsets=episodes(x)[0]
    y=np.full(n,np.nan)
    # A fully observed current-through-45-minute window suffices to establish
    # absence of qualifying onsets 15-30 minutes ahead, including confirmation.
    future_observed=ready(x[::-1],10)[::-1]
    y[future_observed]=0
    for e in onsets:
        y[max(0,e-6):max(0,e-2)]=1
    return y


def features(glucose,kind="M120"):
    """Strict causal glucose lags and OLS summaries; no global normalization."""
    x=np.asarray(glucose,float)
    if kind not in {"M120","M15","missing"}: raise ValueError("Unknown model kind")
    nlag=4 if kind=="M15" else 25
    lag=np.full((len(x),nlag),np.nan,dtype=np.float32)
    for j in range(nlag):
        if j<len(x): lag[j:,j]=x[:len(x)-j]
    windows=[4] if kind=="M15" else [4,7,13,25]
    columns=[lag]
    names=[f"lag_{j*5}min" for j in range(nlag)]
    if kind=="missing":
        columns.append(np.isfinite(lag).astype(np.float32))
        names += [f"observed_{j*5}min" for j in range(nlag)]
        last=np.maximum.accumulate(np.where(np.isfinite(x),np.arange(len(x)),-1))
        since=np.where(last>=0,(np.arange(len(x))-last)*5,np.nan)
        columns.append(since[:,None]); names.append("minutes_since_observation")
    with warnings.catch_warnings():
        warnings.simplefilter("ignore",RuntimeWarning)
        for window in windows:
            a=lag[:,:window]; mask=np.isfinite(a)
            t=np.broadcast_to(-5*np.arange(window,dtype=float),a.shape)
            av=np.where(mask,a,0); tv=np.where(mask,t,0)
            count=mask.sum(axis=1); denom=np.maximum(count,1)
            avg=av.sum(axis=1)/denom; tavg=tv.sum(axis=1)/denom
            centered=np.where(mask,t-tavg[:,None],0)
            ss=(centered**2).sum(axis=1)
            slope=np.divide((centered*av).sum(axis=1),ss,out=np.full(len(x),np.nan),where=ss>0)
            summaries=np.column_stack([avg,np.nanstd(a,axis=1),np.nanmin(a,axis=1),np.nanmax(a,axis=1),slope])
            summaries[count==0]=np.nan
            columns.append(summaries)
            names += [f"{s}_{(window-1)*5}min" for s in ["mean","sd","minimum","maximum","slope"]]
            if kind=="missing": columns.append(count[:,None]); names.append(f"count_{(window-1)*5}min")
    return np.hstack(columns).astype(np.float32),names


def linear_forecast(glucose):
    f,names=features(glucose,"M15")
    out=np.asarray(glucose)+30*f[:,names.index("slope_15min")]
    out[~ready(glucose,4)]=np.nan
    return out


def replay(glucose,probability,threshold,refractory_bins=6):
    """Threshold may be per-bin (switch). NaN is explicit unavailability."""
    x=np.asarray(glucose,float); p=np.asarray(probability,float)
    threshold=np.broadcast_to(threshold,p.shape)
    candidates=np.flatnonzero(np.isfinite(x)&(x>=54)&np.isfinite(p)&(p>=threshold))
    alarms=[]; next_legal=-1
    for i in candidates:
        if i>=next_legal:
            alarms.append(i); next_legal=i+refractory_bins
    return np.asarray(alarms,dtype=np.int64)


def match_alarms(alarms,onsets,glucose):
    """Chronological one-to-one matching. Unobservable outcomes remain cost."""
    alarms=np.asarray(alarms,dtype=int); onsets=np.asarray(onsets,dtype=int)
    used=set(); ledger=[]
    for a in np.sort(alarms):
        eligible=np.flatnonzero((onsets>=a+3)&(onsets<=a+6))
        available=[int(j) for j in eligible if int(j) not in used]
        event=available[0] if available else None
        if event is not None:
            used.add(event); category="matched"
        elif len(eligible): category="duplicate"
        elif np.any((onsets>=a)&(onsets<a+3)): category="late"
        elif a+10>len(glucose) or not np.isfinite(glucose[a:a+10]).all(): category="unobservable_outcome"
        else: category="unmatched_observed"
        ledger.append(dict(alarm_index=int(a),event_index=event,category=category))
    return ledger, np.asarray(sorted(used),dtype=int)


def independent_match_count(alarms,onsets):
    """Separate pair-list scorer checks maximum chronological matching."""
    taken=set(); count=0
    for a in sorted(alarms):
        for j,e in enumerate(onsets):
            if j not in taken and 15 <= (e-a)*5 <= 30:
                taken.add(j); count+=1; break
    return count


def participant_score(glucose,onsets,alarms,availability):
    ledger,matched=match_alarms(alarms,onsets,glucose)
    if len(matched)!=independent_match_count(alarms,onsets):
        raise AssertionError("Independent alarm matching disagreement")
    opportunity=np.array([np.any(availability[max(0,e-6):max(0,e-2)]) for e in onsets],bool)
    if any(not opportunity[j] for j in matched):
        raise AssertionError("A warning was matched without a legal prediction opportunity")
    return dict(events=len(onsets),opportunities=int(opportunity.sum()),matched=len(matched),alarms=len(alarms),unconfirmed=len(alarms)-len(matched),observed_days=float(np.isfinite(glucose).sum()/288),ledger=ledger,matched_event_indices=matched.tolist())


def paired_bootstrap(rows,n_boot=10000,seed=SEED):
    """rows columns: events, matched_switch, matched_comparator,
    unconfirmed_switch, unconfirmed_comparator, observed_days. Clusters are people.
    """
    rows=np.asarray(rows,float)
    if rows.ndim!=2 or rows.shape[1]!=6: raise ValueError("Expected N x 6 participant rows")
    if len(rows)==0: raise ValueError("No participants available for bootstrap")
    rng=np.random.default_rng(seed); samples=[]; undefined=0
    for _ in range(n_boot):
        s=rows[rng.integers(0,len(rows),len(rows))].sum(axis=0)
        if s[0]==0 or s[5]==0: undefined+=1; continue
        samples.append([(s[1]-s[2])/s[0],(s[3]-s[4])/s[5]])
    a=np.asarray(samples)
    return dict(replicates=n_boot,undefined=undefined,sensitivity_difference_ci=np.quantile(a[:,0],[.025,.975]).tolist() if len(a) else None,unconfirmed_difference_ci=np.quantile(a[:,1],[.025,.975]).tolist() if len(a) else None)
