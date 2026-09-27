#!/usr/bin/env python3
"""Rank-three pinning with zeros and parallel classes, via corner enumeration.

Pair-transfer concavity reduces occupancies to box-simplex corners. This is a
standard rounding-style argument, not an asserted new mathematical invention.
"""
from collections import defaultdict
from fractions import Fraction
from itertools import combinations
from math import comb,isqrt
import numpy as np


def determinant(a,b,c,p):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])+a[2]*(b[0]*c[1]-b[1]*c[0]))%p


def partition(p,columns):
    if not (type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))):raise ValueError('prime field required')
    groups=defaultdict(list);zeros=[]
    for i,v in enumerate(columns):
        if len(v)!=3 or any(type(x) is not int or not 0<=x<p for x in v):raise ValueError('canonical three-entry columns required')
        first=next((x for x in v if x),None)
        if first is None:zeros.append(i)
        else:
            inv=pow(first,-1,p);groups[tuple(x*inv%p for x in v)].append(i)
    records=[{'direction':list(v),'coordinates':indices} for v,indices in sorted(groups.items())]
    return zeros,records


def optimize(p,columns,s):
    n=len(columns)
    if not (3<=n<=4096 and type(s) is int and 0<=s<=n):raise ValueError('3<=n<=4096 and integer agreement threshold required')
    zeros,groups=partition(p,columns);m=len(groups)
    if not 3<=m<=22:raise ValueError('corner compiler requires between3 and22 nonzero projective classes')
    basis_triples=[t for t in combinations(range(m),3) if determinant(*(groups[i]['direction'] for i in t),p)]
    if not basis_triples:raise ValueError('evaluation rank must be three')
    capacities=[len(g['coordinates']) for g in groups];size=1<<m;t=s-min(s,len(zeros))
    # Every count is bounded by C(n,3), safely inside signed int64 at this cap.
    assert comb(n,3)<2**63
    weights=np.zeros(size,dtype=np.int64);scores=np.zeros(size,dtype=np.int64)
    for i,c in enumerate(capacities):weights[1<<i:1<<(i+1)]=weights[:1<<i]+c
    for i,j,k in basis_triples:scores[(1<<i)|(1<<j)|(1<<k)]=capacities[i]*capacities[j]*capacities[k]
    for i in range(m):
        step=1<<i;view=scores.reshape(-1,2*step);view[:,step:]+=view[:,:step]
    indices=np.arange(size,dtype=np.int64);best=comb(n,3)+1;best_mask=None;partial=None;candidate_count=0
    full=np.flatnonzero(weights==t);candidate_count+=len(full)
    if len(full):
        choice=int(full[np.argmin(scores[full])]);best=int(scores[choice]);best_mask=choice
    for i,c in enumerate(capacities):
        remaining=t-weights;valid=np.flatnonzero(((indices>>i)&1==0)&(remaining>0)&(remaining<c));candidate_count+=len(valid)
        if not len(valid):continue
        derivative=(scores[valid|(1<<i)]-scores[valid])//c
        values=scores[valid]+remaining[valid]*derivative;j=int(np.argmin(values));value=int(values[j])
        if value<best:best=value;best_mask=int(valid[j]);partial=(i,int(remaining[best_mask]))
    assert best_mask is not None
    occupancies=[c if best_mask>>i&1 else (partial[1] if partial is not None and partial[0]==i else 0) for i,c in enumerate(capacities)]
    witness=zeros[:min(s,len(zeros))]
    for g,a in zip(groups,occupancies):witness+=g['coordinates'][:a]
    assert len(witness)==s
    probability=Fraction(best,comb(n,3))
    return {'schema':'rank_three_capacitated_pinning_v1','p':p,'n':n,'s':s,'zero_coordinates':zeros,'projective_classes':groups,
            'independent_class_triples':[list(t) for t in basis_triples],
            'minimum_injective_triples':best,'uniform_triple_success_probability':[probability.numerator,probability.denominator],
            'worst_agreement_set':sorted(witness),'class_occupancies':occupancies,'partial_class':list(partial) if partial else None,
            'corner_candidates_evaluated':candidate_count,'subset_transform_states':size,
            'scope':'Exact declared rank-three evaluation-space query, at most22 nonzero projective classes. Pair-transfer concavity supplies a corner minimizer; no general cluster coverage or historical novelty claim.'}


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    data=json.loads(Path(sys.argv[1]).read_text())
    print(json.dumps(optimize(data['p'],data['columns'],data['s']),indent=2))
