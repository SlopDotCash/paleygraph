#!/usr/bin/env python3
"""Count actual scalar pairs realizing a specified centered difference box."""
from math import prod
import numpy as np


def count_pairs(p,g,delta,fixed,limit=100,chunk=262144):
    if not fixed:raise ValueError('at least one fixed coordinate is required')
    if type(limit)is not int or limit<0 or type(chunk)is not int or chunk<1:raise ValueError('invalid limits')
    m=p//2;intervals=[]
    for i,H in fixed.items():
        if type(i)is not int or i<0 or type(H)is not int:raise ValueError('invalid coordinate or difference')
        if H%p!=delta*pow(g,i,p)%p:raise ValueError('inconsistent scalar difference')
        lo=max(-m,H-m);hi=min(m,H+m)
        intervals.append((hi-lo+1,i,lo,hi))
    if any(size<=0 for size,*_ in intervals):
        return {'status':'complete','count':0,'scalars_a':[],'truncated':False,'seed_candidates':0,'survivors_after_coordinates':[]}
    intervals.sort();size,seed,lo,hi=intervals[0];inverse=pow(g,-seed,p)
    assert (p-1)**2+m<np.iinfo(np.int64).max and m*(p-1)<np.iinfo(np.int64).max
    survivors=[0]*len(intervals);found=[];total=0
    for start in range(lo,hi+1,chunk):
        y=np.arange(start,min(hi+1,start+chunk),dtype=np.int64)
        a=(y*inverse)%p;survivors[0]+=len(a)
        for j,(_,i,L,U) in enumerate(intervals[1:],1):
            values=(a*pow(g,i,p)+m)%p-m
            a=a[(values>=L)&(values<=U)];survivors[j]+=len(a)
            if not len(a):break
        total+=len(a)
        if len(found)<limit:found.extend(map(int,a[:limit-len(found)]))
    return {'status':'complete','count':total,'scalars_a':found,'truncated':total>len(found),
            'seed_coordinate':seed,'seed_interval':[lo,hi],'seed_candidates':size,
            'coordinate_order':[i for _,i,_,_ in intervals],'survivors_after_coordinates':survivors}
