#!/usr/bin/env python3
"""A completion query with a checked coupled uniqueness certificate."""
import argparse
from fractions import Fraction as F
from itertools import product
from math import floor
from pathlib import Path
import json
from erasure_query import query
from coupled_review import setup,norm_certificate,target_certificate,rational
from conditioned_review import certify


def verify_uniqueness(c, proof):
    certify(c); _,V,rows,Q=setup(c)
    bounds=[floor(2*min(rational(c['coordinate_radii'][j]),rational(c['conditioning']['cuts'][j]['radius']))) for j in V]
    assert proof['visible_directions']==V and proof['individual_difference_bounds']==bounds
    pending={t for t in product(*(range(-b,b+1) for b in bounds)) if any(t) and next(v for v in t if v)>0}
    initial=len(pending)
    for cut in proof['pair_cuts']:
        direction=cut['direction']; assert len(direction)==len(V) and all(type(x) is int for x in direction)
        q=[sum(direction[j]*Q[j][i] for j in range(len(V))) for i in range(c['N'])]
        norm=norm_certificate(q,rows,cut['certificate'])
        cap=floor(F(c['p']-1,c['p'])*norm)
        pending={t for t in pending if abs(sum(x*y for x,y in zip(direction,t)))<=cap}
    after_pairs=len(pending)
    for cut in proof['target_separators']:
        norm=target_certificate(c,cut,rows,Q)
        assert F(c['p']-1,c['p'])*norm<1
        target=tuple(cut['target'])
        pending.discard(target); pending.discard(tuple(-v for v in target))
    assert not pending, 'uncovered integer difference'
    return {'nonzero_sign_representatives':initial,'remaining_after_pair_cuts':after_pairs,
            'uniqueness_method':'coupled_visible_difference_separation'}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,required=True)
    parser.add_argument('--joint-certificate',type=Path,required=True)
    parser.add_argument('--input',type=Path,required=True)
    parser.add_argument('--cyclic-start',type=int)
    parser.add_argument('--candidate-budget',type=int,default=50000)
    args=parser.parse_args()
    try:
        data=json.loads(args.certificate.read_text()); c=data.get('certificate',data)
        proof=verify_uniqueness(c,json.loads(args.joint_certificate.read_text()))
        out=query(c,json.loads(args.input.read_text()),args.candidate_budget,args.cyclic_start)
        out.update(proof); out['universal_unique_completion']=True
        if out['status']=='complete': assert out['count']<=1
    except (ValueError,AssertionError,KeyError,IndexError,TypeError) as error:
        parser.exit(2,'invalid input or certificate: '+str(error)+'\n')
    print(json.dumps(out,indent=2))


if __name__=='__main__': main()
