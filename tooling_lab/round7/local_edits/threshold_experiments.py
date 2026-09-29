#!/usr/bin/env python3
"""Validate the usable API and finite one-sided bounds, then apply saved exact moments."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations_with_replacement
import json
from pathlib import Path
import local_edits
from query import query, summarize, threshold_bound

HERE=Path(__file__).resolve().parent


def main():
    checks=0
    for n in range(1,7):
        for values in combinations_with_replacement(range(-2,3),n):
            mean=F(sum(values),n);var=F(sum(x*x for x in values),n)-mean*mean
            summary={'neighbour_count':n,'neighbour_mean':[mean.numerator,mean.denominator],
                     'neighbour_variance':[var.numerator,var.denominator]}
            for threshold in (F(-3),F(-2),F(-1),F(-1,2),F(0),F(1,2),F(1),F(2),F(3)):
                for direction in ('upper','lower'):
                    r=threshold_bound(summary,threshold,direction)
                    count=sum(x>=threshold if direction=='upper' else x<=threshold for x in values)
                    assert count<=r['neighbour_count_upper']<=n
                    checks+=1
    invalid=[(9,[0,1],1),(21,[0,1],1),(13,[0,0],1),(13,[],0),(13,[-1,0],1),(13,[0,13],1),
             (13,[0,1],3),(13,[0,1],-1),(10000001,[0,1],1),(13,[0,1.5],1),(13,list(range(13)),6)]
    for args in invalid:
        try:query(*args)
        except ValueError:pass
        else:raise AssertionError(('invalid query accepted',args))
    boundaries=[]
    for d in (0,1,3,6):
        fast,_=query(101,list(range(64)),d)
        other=local_edits.from_prime(101,list(range(64)),d,chunk_size=17)
        for key in ('target','neighbour_count','delta_mean','delta_second_moment','neighbour_mean','neighbour_variance'):
            assert fast[key]==other[key],key
        boundaries.append({'q':101,'n':64,'degree':d,'exact_agreement':True})
    saved=json.loads((HERE/'results.json').read_text());gram=json.loads((HERE/'gram_results.json').read_text())
    rows=[]
    for r,g in zip(saved['scale_cases'],gram['scale_cases']):
        fast=summarize(g)
        for key in ('target','neighbour_count','delta_mean','delta_second_moment','neighbour_mean','neighbour_variance'):
            assert r[key]==fast[key],key
        direction='upper' if r['target']<0 else 'lower'
        half=threshold_bound(fast,F(r['target'],2),direction)
        zero=threshold_bound(fast,0,direction)
        rows.append({'q':r['q'],'n':r['n'],'family':r['family'],'target':r['target'],
                     'same_sign_failure_bound':zero,'half_magnitude_failure_bound':half,
                     'half_magnitude_persistence_lower_float':half['complement_count_lower']/r['neighbour_count'],
                     'scope':'All one-swap neighbours of this saved input; no multi-step or uniform-input assertion.'})
    out={'status':'passed','finite_distribution_event_bounds_checked':checks,'invalid_queries_rejected':len(invalid),
         'maximum_n_boundary_cases':boundaries,'scale_query_consumer_comparisons':8,'scale_threshold_bounds':rows,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('threshold_experiments.py','query.py','local_edits.py','gram_oracle.cpp','gram_oracle')},
         'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('results.json','gram_results.json')}}
    (HERE/'threshold_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('scale_threshold_bounds','source_sha256','input_sha256')}),flush=True)


if __name__=='__main__':main()
