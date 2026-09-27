#!/usr/bin/env python3
"""Initial exact worst-case pinning compiler for a declared rank-two cluster."""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations,product
import json
from math import comb,isqrt
from pathlib import Path

HERE=Path(__file__).resolve().parent


def compile_pinning(p,columns,s):
    n=len(columns)
    if not (type(p) is int and p>=2 and all(p%d for d in range(2,isqrt(p)+1)) and 0<=s<=n and n>=2):
        raise ValueError('prime field and valid agreement size required')
    groups=defaultdict(list);zeros=[]
    for i,(a,b) in enumerate(columns):
        if not all(type(x) is int and 0<=x<p for x in (a,b)):raise ValueError('noncanonical field entry')
        if not a and not b:zeros.append(i)
        else:groups[(1,b*pow(a,-1,p)%p) if a else (0,1)].append(i)
    if len(groups)<2:raise ValueError('evaluation space must have rank two')
    ordered=sorted(groups.items(),key=lambda x:(-len(x[1]),x[0]))
    chosen=zeros[:min(s,len(zeros))];remaining=s-len(chosen);counts=[]
    for direction,indices in ordered:
        take=min(remaining,len(indices));chosen.extend(indices[:take]);remaining-=take;counts.append(take)
    assert remaining==0
    t=sum(counts);minimum=(t*t-sum(x*x for x in counts))//2
    probability=F(minimum,comb(n,2))
    return {'p':p,'n':n,'agreement_size':s,'zero_coordinates':zeros,
            'projective_classes':[{'direction':list(k),'coordinates':v} for k,v in ordered],
            'class_counts_in_witness':counts,'worst_agreement_set':sorted(chosen),
            'minimum_injective_pairs':minimum,'uniform_pair_success_probability':[probability.numerator,probability.denominator]}


def brute(p,columns):
    n=len(columns);edges=[(i,j) for i,j in combinations(range(n),2) if (columns[i][0]*columns[j][1]-columns[i][1]*columns[j][0])%p]
    minima={s:comb(n,2)+1 for s in range(n+1)};counts=0
    for s in range(n+1):
        for A in combinations(range(n),s):
            a=set(A);v=sum(i in a and j in a for i,j in edges);minima[s]=min(minima[s],v);counts+=1
    return minima,edges,counts


def main():
    # All rank-two length-five evaluation patterns over F3, up to independent
    # nonzero scaling of columns. This includes non-RS abstract linear spaces.
    symbols=[(0,0),(0,1),(1,0),(1,1),(1,2)];patterns=checks=subsets=0
    for indices in product(range(5),repeat=5):
        if len(set(indices)-{0})<2:continue
        columns=[symbols[i] for i in indices];mins,edges,count=brute(3,columns);subsets+=count;patterns+=1
        for s in range(6):
            r=compile_pinning(3,columns,s);A=set(r['worst_agreement_set'])
            assert r['minimum_injective_pairs']==mins[s]==sum(i in A and j in A for i,j in edges);checks+=1
    out={'status':'toy_passed','scope':'Exact rank-two evaluation-pattern audit; no cluster-discovery algorithm or full proximity theorem.',
         'projective_patterns':patterns,'threshold_checks':checks,'agreement_sets_enumerated':subsets,
         'source_sha256':{'rank_two_pinning.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'preflight_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if not k.endswith('sha256')}),flush=True)


if __name__=='__main__':main()
