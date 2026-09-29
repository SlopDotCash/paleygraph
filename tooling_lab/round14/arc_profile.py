#!/usr/bin/env python3
"""Optimal size profile within the complete-arc-block partition construction.

This is elementary separable convex optimization, not a global Turan optimum.
It assumes every d-subset of a queried block can be made independent; actual
arithmetic feasibility must be checked separately.
"""
from math import comb


def profile(n,s,d,max_block=64):
    if not all(type(v) is int for v in (n,s,d,max_block)) or not 2<=d<=s<=n or not d<=max_block<=64:raise ValueError('require2<=d<=s<=n and d<=max_block<=64')
    candidates=[]
    for t in range(1,min(n//d,(s-1)//(d-1))+1):
        unused=min(s-1-t*(d-1),n-t*d)
        q,r=divmod(n-unused,t)
        if q+(r>0)>max_block:continue
        sizes=[q]*(t-r)+[q+1]*r+[1]*unused
        queries=(t-r)*comb(q,d)+r*comb(q+1,d)
        capacity=t*(d-1)+unused
        candidates.append({'queried_blocks':t,'unqueried_coordinates':unused,'sizes':sizes,
                           'queries':queries,'no_query_capacity':capacity})
    if not candidates:return {'status':'no_profile_within_block_cap','n':n,'s':s,'dimension':d,'max_block':max_block}
    best=min(candidates,key=lambda row:(row['queries'],len(row['sizes']),row['sizes']))
    return {'status':'profile_found','n':n,'s':s,'dimension':d,'max_block':max_block,
            'queried_block_counts_considered':len(candidates),**best,
            'scope':'Minimum query count within partitions into complete d-uniform arc blocks and unqueried coordinates, subject to max_block. Feasibility in an actual evaluation matroid is not implied; global query-family optimality is not claimed.'}


if __name__=='__main__':
    import json,sys
    if len(sys.argv)!=4:raise SystemExit('Usage: arc_profile.py n s dimension')
    print(json.dumps(profile(*(int(v) for v in sys.argv[1:])),indent=2))
