#!/usr/bin/env python3
"""Constructor-free certificate replay with determinant classes and integer DP."""
from fractions import Fraction as F
from itertools import combinations
from math import comb,isqrt


def validate_certificate(c):
    p,k,s=c['p'],c['k'],c['s'];domain=c['domain'];n=len(domain);basis=c['basis'];out=c['pinning']
    assert c['schema']=='rank_two_polynomial_pinning_v1'
    assert type(p) is int and 2<=p<=2**31-1 and all(p%d for d in range(2,isqrt(p)+1))
    assert c['n']==n and type(k) is int and 2<=k<=n<=min(p,100_000) and type(s) is int and 0<=s<=n
    assert len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==2 and all(len(b)<=k and all(type(x) is int and 0<=x<p for x in b) for b in [*basis,c['origin']])
    cols=[tuple(sum(v*pow(x,j,p) for j,v in enumerate(poly))%p for poly in basis) for x in domain]
    groups,zeros,determinants=determinant_partition(cols,p)
    assert len(groups)>=2
    exported=[r['coordinates'] for r in out['projective_classes']]
    assert sorted(map(tuple,map(sorted,exported)))==sorted(map(tuple,map(sorted,groups)))
    assert out['zero_coordinates']==zeros
    for record in out['projective_classes']:
        a,b=record['direction'];assert type(a) is int and type(b) is int and 0<=a<p and 0<=b<p
        assert a==1 or (a==0 and b==1)
        assert all((a*cols[i][1]-b*cols[i][0])%p==0 for i in record['coordinates'])
    optimum,dp_states=integer_optimum([len(g) for g in groups],len(zeros),s)
    assert optimum==out['minimum_injective_pairs']
    assert out['p']==p and out['n']==n and out['agreement_size']==s
    assert out['uniform_pair_success_probability']==[F(optimum,comb(n,2)).numerator,F(optimum,comb(n,2)).denominator]
    witness=out['worst_agreement_set'];assert witness==sorted(set(witness)) and len(witness)==s and all(0<=i<n for i in witness)
    assert out['class_counts_in_witness']==[len(set(witness)&set(g)) for g in exported]
    actual=sum((cols[i][0]*cols[j][1]-cols[i][1]*cols[j][0])%p!=0 for i,j in combinations(witness,2))
    assert actual==optimum
    return {'p':p,'n':n,'k':k,'s':s,'minimum':optimum,'full_coordinate_pairs_checked':determinants,
            'witness_pairs_checked':comb(s,2),'dynamic_program_transitions':dp_states,'projective_classes':len(groups),'zero_columns':len(zeros)}


def determinant_partition(columns,p):
    n=len(columns);zeros=[i for i,x in enumerate(columns) if x==(0,0)];nonzero=[i for i,x in enumerate(columns) if x!=(0,0)]
    parent=list(range(n))
    def find(i):
        while parent[i]!=i:parent[i]=parent[parent[i]];i=parent[i]
        return i
    comparisons=0
    for i,j in combinations(range(n),2):
        comparisons+=1
        if columns[i]!=(0,0) and columns[j]!=(0,0) and (columns[i][0]*columns[j][1]-columns[j][0]*columns[i][1])%p==0:
            parent[find(i)]=find(j)
    groups={}
    for i in nonzero:groups.setdefault(find(i),[]).append(i)
    return list(groups.values()),zeros,comparisons


def integer_optimum(capacities,zeros,s):
    # State t stores the fewest independent pairs among t chosen nonzero columns.
    # Taking a columns from the next class creates exactly a*t new pairs.
    dp={0:0};transitions=0
    for capacity in capacities:
        new={}
        for t,value in dp.items():
            for a in range(min(capacity,s-t)+1):
                candidate=value+a*t;key=t+a
                if key not in new or candidate<new[key]:new[key]=candidate
                transitions+=1
        dp=new
    answer=min(dp[s-z] for z in range(min(zeros,s)+1) if s-z in dp)
    return answer,transitions


if __name__=='__main__':
    import json,sys
    from pathlib import Path
    print(json.dumps(validate_certificate(json.loads(Path(sys.argv[1]).read_text())),indent=2))
