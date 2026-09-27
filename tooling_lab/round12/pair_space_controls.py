#!/usr/bin/env python3
"""Literal triple controls for the pair-space grouping and degree interval."""
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'round11'))
from rank_three_verifier import rank
from pair_space_preflight import catalogue,cross_key,row_space_key,bound_and_witness

HERE=Path(__file__).resolve().parent


def audit(data,exact=None):
    p=data['p'];domain=data['domain'];n=len(domain);s=data['s']
    columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(b))%p for b in data['basis']) for x in domain]
    ordinary,lines=catalogue(columns,p,cross_key);assert (ordinary,lines)==catalogue(columns,p,row_space_key)
    literal={t for t in combinations(range(n),3) if rank([columns[i] for i in t],p)<3}
    from_lines={t for L in lines for t in combinations(L,3)};assert literal==from_lines
    bound=bound_and_witness(n,s,lines);A=set(bound['worst_known_agreement_set'])
    assert len(A)==s and bound['maximum_dependent_triples_interval'][0]==sum(set(t)<=A for t in literal)
    if exact is None:
        exact=comb(s,3)-max(sum(set(t)<=set(A) for t in literal) for A in combinations(range(n),s))
    low,high=bound['minimum_injective_triples_interval'];assert low<=exact<=high
    return {'p':p,'n':n,'s':s,'all_pairs_per_implementation':comb(n,2),'literal_triples':comb(n,3),'interval':[low,high],'literal_minimum':exact}


def main():
    path=HERE.parent/'round11/rank_three_actual.json';source=json.loads(path.read_text());rows=[]
    for case in source['small_cases']:rows.append(audit(case['certificate']['input'],case['literal_review']['minimum']))
    for f in ([0,0,1],[0,0,0,1]):
        for s in range(9):rows.append(audit({'p':11,'k':4,'s':s,'domain':list(range(8)),'basis':[[1],[0,1],f]}))
    for key in (cross_key,row_space_key):
        for a,b in [((0,0,0),(1,2,3)),((1,2,3),(2,4,6))]:
            try:key(a,b,17)
            except AssertionError:pass
            else:raise AssertionError('dependent pair was accepted')
    out={'status':'passed','cases':rows,'dependent_pairs_rejected':4,
         'source_sha256':{'pair_space_controls.py':sha256(Path(__file__).read_bytes()).hexdigest(),'pair_space_preflight.py':sha256((HERE/'pair_space_preflight.py').read_bytes()).hexdigest(),
                          '../round11/rank_three_verifier.py':sha256((HERE.parent/'round11/rank_three_verifier.py').read_bytes()).hexdigest()},
         'input_sha256':{'../round11/rank_three_actual.json':sha256(path.read_bytes()).hexdigest()}}
    (HERE/'pair_space_controls.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(rows),'literal_triples':sum(r['literal_triples'] for r in rows),'dependent_pairs_rejected':4}),flush=True)


if __name__=='__main__':main()
