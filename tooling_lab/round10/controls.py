#!/usr/bin/env python3
"""Affine reindexing, input rejection and explicit integer range audit."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path
from query import query

HERE=Path(__file__).resolve().parent


def main():
    rows=[];entries=0
    for q,c,degrees in [(17,[0,1,2,3,4,6,10],range(7)),(29,[0,1,3,4,7,11,18,20],[6])]:
        for d in degrees:
            original=query(q,c,d);lookup={tuple(p['deleted']):p for p in original['pairs']}
            for multiplier,shift in [(1,5),(2,3),(3,7)]:
                mapping={x:(multiplier*x+shift)%q for x in c};new=query(q,list(mapping.values()),d)
                inverse={y:x for x,y in mapping.items()};sign=1 if pow(multiplier,(q-1)//2,q)==1 else (-1)**d
                for pair in new['pairs']:
                    key=tuple(sorted(inverse[x] for x in pair['deleted']))
                    assert pair['target_sum']==sign*lookup[key]['target_sum']
                    assert F(*pair['mean'])==sign*F(*lookup[key]['mean']);entries+=1
                rows.append({'q':q,'n':len(c),'degree':d,'multiplier':multiplier,'shift':shift,'entries':len(new['pairs'])})
    invalid=[(9,[0,1],2),(13,[0,1,1],2),(13,list(range(12)),2),(13,[0,1,13],2),
             (257,list(range(65)),6),(17,list(range(7)),7),(13,[0,1],2.0),(True,[0,1],2)]
    for args in invalid:
        try:query(*args)
        except ValueError:pass
        else:raise AssertionError(('accepted invalid input',args))
    q=10_000_000;n=64
    bound=lambda j:sum((i+1)*comb(n,j-i) for i in range(j+1))
    phi=q*q*bound(6)+2*q*n*bound(5)+(n*n+q)*bound(4)
    accumulator=16*q*phi+24*n*phi
    assert accumulator<2**127 and q*comb(n,6)<2**63
    out={'status':'passed','affine_cases':rows,'all_pair_entries_reindexed':entries,'invalid_queries_rejected':len(invalid),
         'integer_bounds':{'q_max':q,'n_max':n,'degree_max':6,'hypothetical_two_division_e6_bound':bound(6),
                           'row_phi_bound':phi,'conservative_accumulator_bound':accumulator,'signed128_limit':2**127-1,
                           'inward_target_bound':q*comb(n,6),'signed64_limit':2**63-1},
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('controls.py','query.py','pair_means_backend.cpp','pair_means_backend')}}
    (HERE/'controls_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','affine_cases':len(rows),'reindexed_entries':entries,'invalid_inputs':len(invalid),'integer_bounds_passed':True}),flush=True)


if __name__=='__main__':main()
