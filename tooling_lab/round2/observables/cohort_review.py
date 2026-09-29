#!/usr/bin/env python3
"""Iteration 2: use affine-orbit-disjoint cohorts; audit exact features."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import prod
from pathlib import Path
import numpy as np
from critical_audit import HERE,design,metrics,contractions


def orbit_key(c,p):
    return min(tuple(sorted((x-a)*pow(b-a,-1,p)%p for x in c)) for a in c for b in c if a!=b)


def stable_split(key):
    return int.from_bytes(sha256(','.join(map(str,key)).encode()).digest()[:4],'big')%2


def exact_replay(r,p):
    c=r['C'];n=len(c);squares={x*x%p for x in range(1,p)}
    chi=lambda x:0 if x%p==0 else 1 if x%p in squares else -1
    f=[sum(chi(x-y) for y in c) for x in range(p)]
    k=[sum(prod(chi(x-y) for y in q) for x in range(p)) for q in combinations(c,4)]
    t6=sum(sum(prod(chi(x-y) for y in q) for x in range(p)) for q in combinations(c,6))
    assert k==r['quartic_deck'] and t6==r['t6'] and sum(x**6 for x in f)==r['m6']
    graph=[[0]*n for _ in range(n)];qs=list(combinations(range(n),4))
    for i,j in combinations(range(n),2):
        graph[i][j]=graph[j][i]=sum(v for q,v in zip(qs,k) if i not in q and j not in q)
    w=np.array(graph,dtype=object);q=n*np.eye(n,dtype=object)-np.ones((n,n),dtype=object)
    aq=q@w@q
    direct_center=int(np.trace(aq@aq@aq))
    assert direct_center==n**3*r['centered_triangle_numerator']
    assert contractions(graph)==(r['triangle'],r['centered_triangle_numerator'])
    return True


def main():
    results=[];toy_models=None;replays=0
    for p in [61,1297,2437,4129]:
        data=json.loads((HERE/f'data_{p}.json').read_text());n=data['n']
        keyed={};duplicates=0
        for r in data['random']:
            key=orbit_key(r['C'],p)
            if key in keyed:duplicates+=1
            else:keyed[key]=r
        train=[r for k,r in keyed.items() if stable_split(k)==0]
        test=[r for k,r in keyed.items() if stable_split(k)==1]
        train_keys={k for k in keyed if stable_split(k)==0}
        planted=[r for r in data['planted'] if orbit_key(r['C'],p) not in train_keys]
        assert len(train)>20 and len(test)>20
        local={};transfer={};models={}
        for mode in ['baseline','actual','shuffled']:
            x,y=design(train,p,n,mode);beta=np.linalg.lstsq(x,y,rcond=None)[0];models[mode]=beta
            local[mode]={}
            for name,rows in [('random_holdout',test),('three_anchor_holdout',planted)]:
                x,y=design(rows,p,n,mode);local[mode][name]=metrics(y,x@beta)
            if toy_models is not None:
                x,y=design(test,p,n,mode);transfer[mode]=metrics(y,x@toy_models[mode])
        if toy_models is None:toy_models=models
        # Replay extreme target inputs, including structured stress cases.
        witnesses=[min(test,key=lambda r:r['t6']),max(test,key=lambda r:r['t6']),max(planted,key=lambda r:r['t6'])]
        for r in witnesses:exact_replay(r,p);replays+=1
        results.append({'p':p,'n':n,'raw_random_count':len(data['random']),
                        'distinct_affine_orbits':len(keyed),'duplicate_orbit_records_removed':duplicates,
                        'training_orbits':len(train),'evaluation_orbits':len(test),'planted_evaluated':len(planted),
                        'model_comparison':local,'toy_model_transfer':transfer,
                        'exact_witnesses':witnesses})
        print(p,len(train),len(test),duplicates,{k:round(v['random_holdout']['r2_against_cohort_mean'],6) for k,v in local.items()},flush=True)
    (HERE/'cohort_review.json').write_text(json.dumps({'status':'orbit-disjoint finite statistical audit; no universal lower/upper or impossibility theorem',
                                                     'scalar_extreme_replays':replays,'cases':results,
                                                     'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest()},indent=2)+'\n')


if __name__=='__main__':main()
