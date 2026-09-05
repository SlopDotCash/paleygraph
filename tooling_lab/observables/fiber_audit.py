#!/usr/bin/env python3
"""Exact target-separation audit on actual prime-field character translates.

Equal features / unequal M6 refutes determination, not a uniform upper bound.
The engine retains finite fiber maxima and explicitly disallows target leakage.
"""
from __future__ import annotations

import argparse
from collections import Counter
from hashlib import sha256
from itertools import combinations, islice
import json
from math import comb
from pathlib import Path
import random
import time

import numpy as np

HERE = Path(__file__).resolve().parent
FEATURE_GROUPS = {
    "M2": ("m2",),
    "M2_M4": ("m2", "m4"),
    "plus_boundary": ("m2", "m4", "b2", "b4"),
    "plus_energy": ("m2", "m4", "b2", "b4", "energy"),
    "plus_quartic_deck": ("m2", "m4", "b2", "b4", "quartic_deck"),
    "all_admissible": ("m2", "m4", "b2", "b4", "energy", "quartic_deck"),
    "leakage_control_M6": ("m2", "m4", "m6"),
}


def quadratic_character(p):
    assert p >= 3 and all(p % d for d in range(2, int(p**0.5) + 1))
    return np.array([0] + [1 if pow(x, (p-1)//2, p) == 1 else -1
                          for x in range(1, p)], dtype=np.int64)


def enumerate_sets(p, n, samples, seed):
    """All normalized sets, or a reproducible uniform sample of that slice.

    All n>=2 affine orbits meet {C: 0,1 in C}; the slice is NOT uniform on
    affine orbits. No frequency statement for all sets is inferred.
    """
    if not (isinstance(n,int) and isinstance(p,int) and 4<=n<=p):
        raise ValueError("Supported domain: integer 4 <= n <= prime p")
    if samples is not None and not (isinstance(samples,int) and 0<=samples<=comb(p-2,n-2)):
        raise ValueError("Sample count must fit the normalized slice")
    if samples is None:
        for tail in combinations(range(2, p), n-2):
            yield (0, 1) + tail
    else:
        rng = random.Random(seed)
        seen = set()
        while len(seen) < samples:
            c = (0, 1) + tuple(sorted(rng.sample(range(2, p), n-2)))
            if c not in seen:
                seen.add(c)
                yield c


def batch_observables(sets, p, chi):
    c = np.array(sets, dtype=np.int64)
    b, n = c.shape
    if n<4 or 1000*p*n**6>np.iinfo(np.int64).max:
        raise ValueError("Outside conservative int64-exact arithmetic range")
    rows = chi[(np.arange(p)[None, None, :] - c[:, :, None]) % p]
    f = rows.sum(axis=1)
    f2 = f*f
    m2 = f2.sum(axis=1)
    m4 = (f2*f2).sum(axis=1)
    m6 = (f2*f2*f2).sum(axis=1)
    assert np.all(m2 == n*(p-n))
    at_boundary = np.take_along_axis(f, c, axis=1)
    b2 = (at_boundary**2).sum(axis=1)
    b4 = (at_boundary**4).sum(axis=1)
    pairs = ((c[:, :, None] + c[:, None, :]) % p).reshape(b, -1)
    pair_counts = np.zeros((b,p), dtype=np.int64)
    np.add.at(pair_counts, (np.arange(b)[:,None], pairs), 1)
    energy = (pair_counts*pair_counts).sum(axis=1)
    deck = []
    for q in combinations(range(n), 4):
        deck.append(np.prod(rows[:, q, :], axis=1).sum(axis=1))
    decks = np.sort(np.array(deck, dtype=np.int64).T, axis=1)
    # Directly check the pass-21 sixth-moment identity. T6 is diagnostic
    # output only, forbidden as a proposed independent predictor of M6.
    t6 = np.zeros(b, dtype=np.int64)
    for q in combinations(range(n), 6):
        t6 += np.prod(rows[:, q, :], axis=1).sum(axis=1)
    a0 = 15*n**3 - 30*n**2 + 16*n
    a2 = 90*n*n - 300*n + 272
    a4 = 360*n - 960
    rhs = p*a0 - a2*n*(n-1)//2 + a4*decks.sum(axis=1) + 720*t6 - 15*b4 - 15*b2 - n
    assert np.array_equal(m6, rhs)
    return [{"C": tuple(map(int,c[i])), "m2":int(m2[i]), "m4":int(m4[i]),
             "m6":int(m6[i]), "b2":int(b2[i]), "b4":int(b4[i]),
             "energy":int(energy[i]), "quartic_deck":tuple(map(int,decks[i])),
             "t6_excluded":int(t6[i])} for i in range(b)]


def summarize_fibers(records, feature_names, threshold):
    fibers = {}
    for i,r in enumerate(records):
        key = tuple(r[f] for f in feature_names)
        y = r["m6"]
        if key not in fibers:
            fibers[key] = [y,y,i,i,1]
        else:
            f=fibers[key]; f[4]+=1
            if y<f[0]: f[0],f[2]=y,i
            if y>f[1]: f[1],f[3]=y,i
    ambiguous=[f for f in fibers.values() if f[0]!=f[1]]
    worst=max(fibers.values(),key=lambda f:f[1]-f[0])
    return {"feature_names": list(feature_names), "fibers":len(fibers),
            "ambiguous_fibers":len(ambiguous),
            "records_in_ambiguous_fibers":sum(f[4] for f in ambiguous),
            "total_record_weighted_range":sum(f[4]*(f[1]-f[0]) for f in fibers.values()),
            "max_range":worst[1]-worst[0],
            "threshold":threshold,
            "fibers_straddling_threshold":sum(f[0]<=threshold<f[1] for f in fibers.values()),
            "witness": [records[worst[2]],records[worst[3]]] if worst[0]!=worst[1] else None}


def direct_record(c, p):
    """Independent scalar implementation for exact witness replay."""
    def ch(x):
        x%=p
        return 0 if x==0 else 1 if pow(x,(p-1)//2,p)==1 else -1
    f=[sum(ch(x-y) for y in c) for x in range(p)]
    pc=Counter((a+b)%p for a in c for b in c)
    def corr(q):
        from math import prod
        return sum(prod(ch(x-y) for y in q) for x in range(p))
    return {"C":tuple(c), "m2":sum(x*x for x in f),
            "m4":sum(x**4 for x in f), "m6":sum(x**6 for x in f),
            "b2":sum(f[x]**2 for x in c), "b4":sum(f[x]**4 for x in c),
            "energy":sum(v*v for v in pc.values()),
            "quartic_deck":tuple(sorted(corr(q) for q in combinations(c,4))),
            "t6_excluded":sum(corr(q) for q in combinations(c,6))}


def run_case(p,n,samples=None,seed=20260905):
    started=time.monotonic()
    chi=quadratic_character(p)
    source=iter(enumerate_sets(p,n,samples,seed))
    records=[]
    while batch:=list(islice(source, 128)):
        records.extend(batch_observables(batch,p,chi))
    # Gaussian sixth moment is a baseline only, NOT an asserted inequality.
    threshold=15*p*n**3
    groups={name:summarize_fibers(records,fs,threshold) for name,fs in FEATURE_GROUPS.items()}
    replayed=0
    for group in groups.values():
        if group["witness"]:
            left,right=group["witness"]
            assert left["m6"] != right["m6"]
            assert all(left[f]==right[f] for f in group["feature_names"])
            for w in group["witness"]:
                assert direct_record(w["C"],p)==w
                replayed+=1
    assert groups["leakage_control_M6"]["ambiguous_fibers"]==0
    # Feature refinement cannot increase aggregate finite ambiguity.
    for a,b in [("M2","M2_M4"),("M2_M4","plus_boundary"),
                ("plus_boundary","plus_energy"),("plus_boundary","plus_quartic_deck"),
                ("plus_energy","all_admissible"),("plus_quartic_deck","all_admissible")]:
        assert groups[b]["total_record_weighted_range"] <= groups[a]["total_record_weighted_range"]
    # Affine invariance of all chosen observables is independently tested;
    # use even moments and sign-invariant quartic correlations.
    affine_checks=0
    for r in records[::max(1,len(records)//7)][:7]:
        for a,b in [(2,3),(p-1,1)]:
            transformed=direct_record(sorted((a*x+b)%p for x in r["C"]),p)
            assert all(r[f]==transformed[f] for f in r if f!="C")
            affine_checks+=1
    best=min(records,key=lambda r:r["m6"])
    worst=max(records,key=lambda r:r["m6"])
    return {"p":p,"n":n,"sampling":"all normalized sets" if samples is None else "seeded normalized sample",
            "seed":seed,"records":len(records),"affine_orbits_are_not_uniformly_weighted":True,
            "m6_range":[best["m6"],worst["m6"]],"groups":groups,
            "checks":{"sixth_identity":len(records),"second_identity":len(records),
                      "scalar_witness_replays":replayed,"affine_replays":affine_checks,
                      "leakage_control":True,"refinement_monotonicity":True},
            "elapsed_seconds":round(time.monotonic()-started,3)}


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--suite", choices=["toy","scale","holdout"], default="toy")
    ap.add_argument("--output",type=Path)
    args=ap.parse_args()
    configurations={"toy":[(13,6,None),(17,6,None),(29,6,None)],
                    "scale":[(41,6,None),(61,6,None)],
                    "holdout":[(257,6,8192),(1297,6,4096),(4099,8,1024)]}[args.suite]
    results=[]
    for p,n,samples in configurations:
        r=run_case(p,n,samples)
        results.append(r)
        print(json.dumps({"p":p,"n":n,"records":r["records"],"seconds":r["elapsed_seconds"],
                          "ambiguity":{k:v["ambiguous_fibers"] for k,v in r["groups"].items()}}),flush=True)
    report={"schema_version":1,"suite":args.suite,
            "status":"finite exact diagnostic; no asymptotic bound or global novelty claim",
            "script_sha256":sha256(Path(__file__).read_bytes()).hexdigest(),
            "excluded_predictors":{"m6":"target itself; positive leakage control only",
                                   "t6_excluded":"equivalent target after known lower and boundary terms"},
            "cases":results}
    output=args.output or HERE/f"results_{args.suite}.json"
    output.write_text(json.dumps(report,indent=2)+"\n")


if __name__ == "__main__":
    main()
