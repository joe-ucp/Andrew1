"""Frozen learner and calibration machinery. No patient fit has been run yet.

Functions consume audited canonical five-minute participant timelines. This is
not a Jaeb source adapter. Hyperparameters/thresholds follow Card 001 section12.4.
"""
from __future__ import annotations
import hashlib,time
import numpy as np
from card001_core import SEED,features,ready,labels,episodes,replay,match_alarms


def participant_weights(groups):
    _,inverse,counts=np.unique(groups,return_inverse=True,return_counts=True)
    weights=1/counts[inverse]
    return weights/weights.mean()


def design(participants,kind):
    xs=[]; ys=[]; groups=[]; omission=[]
    for p in participants:
        g=np.asarray(p["glucose"],float)
        x,_=features(g,kind); y=labels(g)
        legal=ready(g,25 if kind=="M120" else 4)
        selected=legal&np.isfinite(y)
        xs.append(x[selected]); ys.append(y[selected]); groups.extend([p["id"]]*selected.sum())
        omission.append(dict(id=p["id"],ready_bins=int(legal.sum()),unknown_labels=int((legal&~np.isfinite(y)).sum()),training_bins=int(selected.sum())))
    if not xs or not sum(len(x) for x in xs): raise ValueError("No eligible labelled training data")
    return np.concatenate(xs),np.concatenate(ys),np.asarray(groups),omission


def fit_model(participants,kind):
    """Five participant-grouped folds, equal person weight, no early stopping."""
    from sklearn.ensemble import HistGradientBoostingClassifier
    from sklearn.metrics import log_loss
    X,y,groups,omission=design(participants,kind)
    if len(np.unique(groups))<5: raise ValueError("Five grouped folds require at least five training participants")
    if len(np.unique(y))<2: raise ValueError("Training labels do not contain both classes")
    # Card001 allows grouped five-fold CV without prescribing assignment;
    # use Card003's stronger exact commitment to share the fitted learners.
    ids=sorted({p["id"] for p in participants},key=lambda x:hashlib.sha256(("t1d-card-003-fold:"+str(int(x))).encode()).hexdigest())
    folds={pid:j%5 for j,pid in enumerate(ids)}
    fold_index=np.array([folds[g] for g in groups])
    options=[]
    for leaves in [7,15]:
        for l2 in [0,1]:
            fold_losses=[]
            for fold in range(5):
                train=np.flatnonzero(fold_index!=fold); validation=np.flatnonzero(fold_index==fold)
                if len(np.unique(y[train]))<2: raise ValueError("A CV training fold has no positive or negative outcomes")
                model=HistGradientBoostingClassifier(max_leaf_nodes=leaves,l2_regularization=l2,learning_rate=.05,max_iter=200,min_samples_leaf=50,early_stopping=False,random_state=SEED)
                t0=time.perf_counter()
                model.fit(X[train],y[train],sample_weight=participant_weights(groups[train]))
                pred=model.predict_proba(X[validation])[:,1]
                loss=log_loss(y[validation],pred,sample_weight=participant_weights(groups[validation]),labels=[0,1])
                fold_losses.append((float(loss),len(np.unique(groups[validation]))))
                print(f"fit {kind} leaves={leaves} l2={l2} fold={fold} rows={len(train)} sec={time.perf_counter()-t0:.1f} loss={loss:.8f}",flush=True)
            weighted_loss=sum(a*b for a,b in fold_losses)/sum(b for a,b in fold_losses)
            options.append(dict(max_leaf_nodes=leaves,l2_regularization=l2,weighted_log_loss=weighted_loss,fold_losses=fold_losses))
    selected=min(options,key=lambda o:(o["weighted_log_loss"],o["max_leaf_nodes"],-o["l2_regularization"]))
    model=HistGradientBoostingClassifier(max_leaf_nodes=selected["max_leaf_nodes"],l2_regularization=selected["l2_regularization"],learning_rate=.05,max_iter=200,min_samples_leaf=50,early_stopping=False,random_state=SEED)
    model.fit(X,y,sample_weight=participant_weights(groups))
    return model,dict(kind=kind,selected=selected,candidates=options,omissions=omission,training_participants=len(np.unique(groups)),training_rows=len(y),folds={str(k):int(v) for k,v in folds.items()})


def predict(model,glucose,kind):
    X,_=features(glucose,kind)
    can=ready(glucose,25 if kind=="M120" else 4)&(np.asarray(glucose)>=54)
    out=np.full(len(glucose),np.nan)
    if can.any(): out[can]=model.predict_proba(X[can])[:,1]
    return out


def calibration_score(participants,kind,threshold,short_threshold=None):
    totals=np.zeros(4) # events, matched, issued, observed days
    for p in participants:
        g=p["glucose"]; probs=p["predictions"]
        if kind=="switch":
            long=np.isfinite(probs["M120"])
            score=np.where(long,probs["M120"],probs["M15"])
            cut=np.where(long,threshold,short_threshold)
        else: score=probs[kind]; cut=threshold
        a=replay(g,score,cut)
        onsets=p.get("onsets")
        if onsets is None: onsets=episodes(g)[0]
        _,matched=match_alarms(a,onsets,g)
        totals+=np.array([len(onsets),len(matched),len(a),np.isfinite(g).sum()/288])
    events,matched,alarms,days=totals
    return dict(events=int(events),matched=int(matched),alarms=int(alarms),observed_days=float(days),sensitivity=float(matched/events) if events else None,unconfirmed_per_day=float((alarms-matched)/days) if days else None,total_per_day=float(alarms/days) if days else None)


def admissible(score,margin_cap=None):
    return score["sensitivity"] is not None and score["unconfirmed_per_day"] is not None and score["unconfirmed_per_day"]<=1 and score["total_per_day"]<=2 and (margin_cap is None or score["unconfirmed_per_day"]<=margin_cap+1e-12)


def ranking(row):
    s=row["score"]
    threshold=row["threshold"]
    physical=(-threshold if row["kind"]=="linear" and np.isfinite(threshold) else threshold)
    return (-s["sensitivity"],s["unconfirmed_per_day"],s["total_per_day"],-physical,-row.get("short_threshold",0))


def calibrate(participants):
    """Exact full-grid replay; never-alarm is encoded by +infinity internally.

    No test participants are accepted. Caller must preserve frozen model identity
    and persist this result before loading internal or external test outcomes.
    """
    if any(p.get("split")!="validation" for p in participants):
        raise ValueError("Calibration may only read DCLP3 validation participants")
    probabilities=np.r_[np.arange(101)/100,np.inf]
    table=[]; best={}
    for kind in ["M120","M15","linear","missing"]:
        thresholds=np.r_[-np.arange(40,101),np.inf] if kind=="linear" else probabilities
        candidates=[]
        for cut in thresholds:
            row=dict(kind=kind,threshold=float(cut),score=calibration_score(participants,kind,cut))
            table.append(row)
            if admissible(row["score"]): candidates.append(row)
        if not candidates: raise ValueError(f"No scorable validation operating point for {kind}")
        best[kind]=min(candidates,key=ranking)
    comparator=min(best.values(),key=ranking)
    margin_cap=comparator["score"]["unconfirmed_per_day"]+.10
    switch=[]
    for long_cut in probabilities:
        for short_cut in probabilities:
            row=dict(kind="switch",threshold=float(long_cut),short_threshold=float(short_cut),score=calibration_score(participants,"switch",long_cut,short_cut))
            table.append(row)
            if admissible(row["score"],margin_cap): switch.append(row)
    if not switch: raise ValueError("No scorable switch operating point")
    return dict(comparator=comparator,comparators=best,switch=min(switch,key=ranking),validation_grid=table,threshold_tie_rule="higher long threshold then higher short threshold; exact residual policy ties retain listed comparator order")
