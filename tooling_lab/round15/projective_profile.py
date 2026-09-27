#!/usr/bin/env python3
"""Complete-block size profiles constrained by zero and projective capacities.

An arc block uses at most one coordinate of each nonzero projective class. With
t queried blocks, at least z+sum(max(c_i-t,0)) coordinates must be unqueried.
The compiler enforces this necessary condition before any higher-rank search.
"""
from math import comb


def optimize(capacities,zeros,s,d,max_block=64):
    if not capacities or any(type(c) is not int or c<1 for c in capacities):raise ValueError('positive integer projective capacities required')
    n=sum(capacities)+zeros
    if not all(type(x) is int for x in (zeros,s,d,max_block)) or zeros<0 or not 2<=d<=s<=n or not d<=max_block<=64:raise ValueError('invalid dimensions, zeros or threshold')
    rows=[]
    for t in range(1,min(n//d,(s-1)//(d-1))+1):
        u=min(s-1-t*(d-1),n-t*d);q,r=divmod(n-u,t)
        mandatory=zeros+sum(max(c-t,0) for c in capacities)
        admissible=mandatory<=u and q+(r>0)<=max_block
        rows.append({'queried_blocks':t,'unqueried_coordinates':u,'mandatory_unqueried':mandatory,
                     'queried_block_sizes':[q+1]*r+[q]*(t-r),
                     'queries':r*comb(q+1,d)+(t-r)*comb(q,d),
                     'maximum_query_avoiding_set':t*(d-1)+u,'admissible':admissible})
    feasible=[r for r in rows if r['admissible']]
    if not feasible:return {'status':'no_parallel_feasible_profile','n':n,'s':s,'dimension':d,'candidate_profiles':rows}
    best=min(feasible,key=lambda r:(r['queries'],r['queried_blocks']))
    return {'status':'parallel_profile_found','n':n,'s':s,'dimension':d,'zeros':zeros,'max_block':max_block,
            'best':best,'candidate_profiles':rows,
            'scope':'Optimal complete-block query count under the parallel-class capacity relaxation and block-size cap. Higher-rank dependencies still require separate arithmetic checks.'}


def pack(groups,zeros,profile):
    best=profile['best'];t=best['queried_blocks'];target=profile['n']-best['unqueried_coordinates'];blocks=[[] for _ in range(t)]
    unqueried=list(zeros);remaining=target;slot=0
    for group in groups:
        take=min(len(group),t,remaining)
        for i in group[:take]:blocks[slot%t].append(i);slot+=1
        unqueried.extend(group[take:]);remaining-=take
    assert remaining==0 and slot==target
    assert [len(B) for B in blocks]==best['queried_block_sizes']
    assert len(unqueried)==best['unqueried_coordinates']
    return [sorted(B) for B in blocks]+[[i] for i in sorted(unqueried)]
