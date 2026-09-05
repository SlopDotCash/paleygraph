#!/usr/bin/env python3
"""Iteration 2: exact integer moment trades and character-realizability gaps.

For six columns, fixed M2/M4 and boundary B2/B4 leave one primitive
histogram move (-10,15,-6,1). We replay its occurrence in actual kernels
and compare the entire nonnegative histogram fiber to actual-set fibers.
The move is elementary integer-kernel / Markov-basis mathematics, not a
new general mathematical invention.
"""
from collections import Counter
from hashlib import sha256
from itertools import islice
import json
from pathlib import Path

import numpy as np
import sympy as sp

from fiber_audit import batch_observables, direct_record, enumerate_sets, quadratic_character

HERE=Path(__file__).resolve().parent
MOVE=(-10,15,-6,1)
LEVELS=(0,2,4,6)
STEP=23040


def histogram(c,p):
    chi=quadratic_character(p)
    f=[int(sum(chi[(x-y)%p] for y in c)) for x in range(p)]
    interior=Counter(abs(f[x]) for x in range(p) if x not in c)
    boundary=Counter(abs(f[x]) for x in c)
    return [interior[a] for a in LEVELS], [boundary[a] for a in (1,3,5)]


def interval(h):
    # h + t*MOVE >= 0; coefficients of either sign are nonzero.
    lower=max(-(h[i]//v) for i,v in enumerate(MOVE) if v>0)
    upper=min(h[i]//(-v) for i,v in enumerate(MOVE) if v<0)
    assert lower<=0<=upper
    return lower,upper


def circuit_certificate():
    a=sp.Matrix([[x**(2*k) for x in LEVELS] for k in range(3)])
    v=sp.Matrix(MOVE)
    assert a.rank()==3 and a*v==sp.zeros(3,1)
    assert sum(vv*x**6 for x,vv in zip(LEVELS,MOVE))==STEP
    # Since v[-1]=1, all integral kernel vectors are integer multiples of v.
    b=sp.Matrix([[x**(2*k) for x in (1,3,5)] for k in range(3)])
    assert b.det()!=0
    return {"constraint_matrix": [list(map(int,a.row(i))) for i in range(3)],
            "integer_kernel_generator":list(MOVE),"rank":3,
            "sixth_moment_step":STEP,"boundary_matrix_determinant":int(b.det()),
            "integer_kernel_argument":"rank one over Q and final coordinate of generator is 1"}


def replay_collision(c,p,features):
    l,r=c
    left,right=direct_record(l['C'],p),direct_record(r['C'],p)
    assert all(left[f]==right[f] for f in features)
    lh,lb=histogram(l['C'],p); rh,rb=histogram(r['C'],p)
    t=rh[-1]-lh[-1]
    assert lb==rb and [b-a for a,b in zip(lh,rh)]==[t*v for v in MOVE]
    assert right['m6']-left['m6']==STEP*t
    lo,hi=interval(lh)
    return {"p":p,"features_preserved":features,"left":left,"right":right,
            "interior_histograms":[lh,rh],"boundary_histogram":lb,
            "primitive_steps":t,"nonnegative_histogram_t_interval":[lo,hi],
            "relaxed_m6_interval":[left['m6']+STEP*lo,left['m6']+STEP*hi]}


def realizability_gap(p):
    """Exhaustive over normalized sets; all affine orbits are covered."""
    source=iter(enumerate_sets(p,6,None,0)); chi=quadratic_character(p)
    fibers={}; count=0
    while batch:=list(islice(source,128)):
        for r in batch_observables(batch,p,chi):
            key=tuple(r[f] for f in ('m2','m4','b2','b4'))
            if key not in fibers:
                fibers[key]={'min':r['m6'],'max':r['m6'],'example':r,'count':0}
            f=fibers[key]
            f['min']=min(f['min'],r['m6']);f['max']=max(f['max'],r['m6']);f['count']+=1
            count+=1
    gaps=[]; tight=0
    for key,f in fibers.items():
        h,b=histogram(f['example']['C'],p)
        lo,hi=interval(h)
        relaxed_min=f['example']['m6']+STEP*lo
        relaxed_max=f['example']['m6']+STEP*hi
        assert relaxed_min<=f['min']<=f['max']<=relaxed_max
        if relaxed_max==f['max']:tight+=1
        else:gaps.append({'feature_key':key,'example':f['example'],'actual_m6_interval':[f['min'],f['max']],
                         'relaxed_m6_interval':[relaxed_min,relaxed_max],
                         'upper_gap':relaxed_max-f['max'],'interior_histogram':h,
                         'boundary_histogram':b,'normalized_sets_in_fiber':f['count']})
    gaps.sort(key=lambda f:f['upper_gap'],reverse=True)
    return {'p':p,'n':6,'normalized_sets':count,'actual_fibers':len(fibers),
            'upper_relaxation_tight_fibers':tight,'upper_realizability_gap_fibers':len(gaps),
            'max_upper_gap':max([f['upper_gap'] for f in gaps],default=0),
            'largest_gap_witnesses':gaps[:4]}


def main():
    collisions=[]
    for filename in ['results_toy.json','results_scale.json','results_holdout.json']:
        path=HERE/filename
        if not path.exists():continue
        for case in json.loads(path.read_text())['cases']:
            if case['n']!=6:continue
            for group in ('plus_boundary','plus_quartic_deck','all_admissible'):
                g=case['groups'][group]
                if g['witness']:
                    collisions.append(replay_collision(g['witness'],case['p'],g['feature_names']))
    report={'status':'exact finite mechanism and exhaustive toy realizability audit; not a uniform bound',
            'circuit':circuit_certificate(),'replayed_collisions':collisions,
            'realizability':[realizability_gap(p) for p in (29,41,61)],
            'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'results_trades.json').write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'collision_replays':len(collisions),
                      'realizability':[{k:v for k,v in r.items() if k!='largest_gap_witnesses'} for r in report['realizability']]},indent=2))


if __name__=='__main__':main()
