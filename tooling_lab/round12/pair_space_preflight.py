#!/usr/bin/env python3
"""Quadratic complete pair-space audit of the full saved simple rank-three input.

Two arithmetic constructions: cross-product normal forms and direct reduced row
spaces. This preflight is not yet a standalone general certificate interface.
"""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import heapq
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter

HERE=Path(__file__).resolve().parent


def evaluate(poly,x,p):
    result=0
    for a in poly[::-1]:result=(result*x+a)%p
    return result


def cross_key(a,b,p):
    v=((a[1]*b[2]-a[2]*b[1])%p,(a[2]*b[0]-a[0]*b[2])%p,(a[0]*b[1]-a[1]*b[0])%p)
    first=next((x for x in v if x),None);assert first is not None
    inv=pow(first,-1,p);return tuple(x*inv%p for x in v)


def row_space_key(a,b,p):
    rows=[list(a),list(b)];pivot=0
    for column in range(3):
        found=next((j for j in range(pivot,2) if rows[j][column]%p),None)
        if found is None:continue
        rows[pivot],rows[found]=rows[found],rows[pivot]
        inv=pow(rows[pivot][column]%p,-1,p);rows[pivot]=[x*inv%p for x in rows[pivot]]
        for j in range(2):
            if j==pivot:continue
            value=rows[j][column];rows[j]=[(x-value*y)%p for x,y in zip(rows[j],rows[pivot])]
        pivot+=1
        if pivot==2:break
    assert pivot==2
    return tuple(rows[0]+rows[1])


def catalogue(columns,p,key):
    groups={};n=len(columns)
    for i,j in combinations(range(n),2):
        code=key(columns[i],columns[j],p);mask=(1<<i)|(1<<j)
        if code in groups:
            old,count=groups[code];groups[code]=(old|mask,count+1)
        else:groups[code]=(mask,1)
    ordinary=0;blocks=[]
    for mask,count in groups.values():
        size=mask.bit_count();assert count==comb(size,2)
        if size==2:ordinary+=1
        else:blocks.append(tuple(i for i in range(n) if mask>>i&1))
    assert ordinary+sum(comb(len(L),2) for L in blocks)==comb(n,2)
    return ordinary,sorted(blocks)


def bound_and_witness(n,s,lines):
    incidence=[[] for _ in range(n)];occupancy=list(map(len,lines));active=set(range(n))
    degrees=[0]*n
    for j,L in enumerate(lines):
        for i in L:incidence[i].append(j);degrees[i]+=comb(len(L)-1,2)
    total=sum(comb(len(L),3) for L in lines)
    assert sum(degrees)==3*total
    upper=min(total,sum(sorted(degrees,reverse=True)[:s])//3,comb(s,3))
    loss=degrees[:];queue=[(d,i) for i,d in enumerate(loss)];heapq.heapify(queue)
    while len(active)>s:
        while True:
            value,i=heapq.heappop(queue)
            if i in active and value==loss[i]:break
        active.remove(i)
        for j in incidence[i]:
            before=occupancy[j];occupancy[j]-=1
            for v in lines[j]:
                if v in active:loss[v]-=max(0,before-2);heapq.heappush(queue,(loss[v],v))
    attained=sum(comb(a,3) for a in occupancy if a>=3)
    assert attained==sum(comb(len(active&set(L)),3) for L in lines if len(active&set(L))>=3)
    assert attained<=upper
    interval=[comb(s,3)-upper,comb(s,3)-attained]
    return {'total_dependent_triples':total,'dependency_vertex_degrees':degrees,
            'maximum_dependent_triples_interval':[attained,upper],'minimum_injective_triples_interval':interval,
            'uniform_triple_success_probability_interval':[[F(x,comb(n,3)).numerator,F(x,comb(n,3)).denominator] for x in interval],
            'worst_known_agreement_set':sorted(active),'status':'exact' if attained==upper else 'bounded'}


def main():
    path=HERE.parent/'round11/rank_three_actual.json';old=json.loads(path.read_text());source=old['large_declared_space']
    p,n,k,s=(source[x] for x in ('p','n','k','s'));domain=source['domain'];basis=source['basis'];start=perf_counter()
    columns=[tuple(evaluate(poly,x,p) for poly in basis) for x in domain]
    ordinary,lines=catalogue(columns,p,cross_key);producer_seconds=perf_counter()-start
    print(json.dumps({'phase':'cross_product_complete','all_pairs':comb(n,2),'nontrivial_lines':len(lines),'seconds':producer_seconds}),flush=True)
    start=perf_counter()
    other_columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(poly))%p for poly in basis) for x in domain]
    assert other_columns==columns
    other_ordinary,other_lines=catalogue(other_columns,p,row_space_key)
    assert (ordinary,lines)==(other_ordinary,other_lines);review_seconds=perf_counter()-start
    bound=bound_and_witness(n,s,lines)
    # Both catalogues prove all pair-span classes are complete. A dependent
    # triple belongs to a unique such space, so this covers all C(n,3) triples
    # without literally iterating through them.
    out={'status':'full_input_preflight_passed','scope':'Complete quadratic pair-space agreement for the saved simple rank-three polynomial space. Initial integer degree-bound interval with an attaining feasible set; exact optimum and standalone interface not claimed.',
         'p':p,'n':n,'k':k,'s':s,'domain':domain,'basis':basis,'ordinary_two_point_lines':ordinary,
         'nontrivial_maximal_lines':[list(L) for L in lines],'line_size_histogram':[[a,b] for a,b in sorted(Counter(map(len,lines)).items())],
         'all_pairs_per_implementation':comb(n,2),'triple_dependencies_covered_by_pair_partition':comb(n,3),
         'literal_triple_enumeration_performed':False,'bound':bound,'universal_degree_bound':comb(s-k+3,3),
         'producer_seconds':producer_seconds,'independent_arithmetic_seconds':review_seconds,
         'source_sha256':{'pair_space_preflight.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'../round11/rank_three_actual.json':sha256(path.read_bytes()).hexdigest()}}
    (HERE/'pair_space_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('domain','basis','nontrivial_maximal_lines','bound','source_sha256','input_sha256')}
                     | {'pinning_interval':bound['minimum_injective_triples_interval'],'dependent_triples':bound['total_dependent_triples']}),flush=True)


if __name__=='__main__':main()
