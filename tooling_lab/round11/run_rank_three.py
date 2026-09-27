#!/usr/bin/env python3
"""Compile and separately replay actual polynomial examples at increasing size."""
from copy import deepcopy
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter
import numpy as np
from rank_three_pinning import certify
from rank_three_verifier import validate,rank

HERE=Path(__file__).resolve().parent


def literal(data):
    p=data['p'];domain=data['domain'];n=len(domain);s=data['s']
    columns=[tuple(sum(a*pow(x,j,p) for j,a in enumerate(b))%p for b in data['basis']) for x in domain]
    dependent=[sum(1<<i for i in t) for t in combinations(range(n),3) if rank([columns[i] for i in t],p)<3]
    masks=np.fromiter((sum(1<<i for i in A) for A in combinations(range(n),s)),dtype=np.uint64,count=comb(n,s))
    counts=np.zeros(len(masks),dtype=np.int64)
    for t in dependent:counts+=(masks&t)==t
    return {'agreement_sets':len(masks),'minimum':comb(s,3)-int(counts.max()),'dependent_triples':len(dependent)}


def main():
    rows=[];certs=[];names=[]
    polynomials={'quartic':[0,0,1,0,1],'sextic':[0,0,1,0,0,0,1]}
    for p,n,s in [(17,16,11),(29,24,16),(37,32,22),(67,48,33),(67,64,44)]:
        for label,f in polynomials.items():
            name=f'rank3_p{p}_n{n}_{label}'
            data={'p':p,'k':8,'s':s,'domain':list(range(1,n+1)),'basis':[[1],[0,1],f]}
            start=perf_counter();certificate=certify(data,node_limit=20_000);seconds=perf_counter()-start
            path=HERE/f'{name}.certificate.json';path.write_text(json.dumps(certificate,separators=(',',':'))+'\n')
            start=perf_counter();review=validate(certificate);review_seconds=perf_counter()-start
            exhaustive=literal(data) if n<=24 else None
            if exhaustive is not None:assert review['interval']==[exhaustive['minimum']]*2
            row={'name':name,**review,'compile_seconds':seconds,'review_seconds':review_seconds,'literal_review':exhaustive,
                 'universal_degree_bound':comb(s-8+3,3),'certificate':path.name}
            rows.append(row);certs.append(certificate);names.append(path.name)
            print(json.dumps(row),flush=True)
    # The bounded case must survive validation with its uncertainty explicit.
    limited=certify({'p':17,'k':8,'s':11,'domain':list(range(1,17)),'basis':[[1],[0,1],polynomials['quartic']]},node_limit=1)
    limited_review=validate(limited);assert limited_review['status']=='bounded'
    (HERE/'rank3_frontier_control.certificate.json').write_text(json.dumps(limited,separators=(',',':'))+'\n');names.append('rank3_frontier_control.certificate.json')
    invalid=[]
    for label in ('missing_line','added_coordinate','inflated_minimum','missing_tree_leaf','invalid_branch','bad_witness','false_exact'):
        bad=deepcopy(limited if label=='false_exact' else certs[0])
        if label=='missing_line':bad['nontrivial_maximal_lines'].pop()
        elif label=='added_coordinate':bad['nontrivial_maximal_lines'][0].append(16)
        elif label=='inflated_minimum':bad['pinning']['minimum_injective_triples_interval'][0]+=1
        elif label=='missing_tree_leaf':bad['pinning']['proof_tree'].pop()
        elif label=='invalid_branch':bad['pinning']['proof_tree'][0]=16
        elif label=='bad_witness':bad['pinning']['worst_known_agreement_set'].pop()
        else:bad['pinning']['status']='exact'
        try:validate(bad)
        except (AssertionError,ValueError,IndexError):invalid.append(label)
        else:raise AssertionError(('corruption accepted',label))
    out={'status':'passed','scope':'Exact or explicitly bounded worst agreement-set pinning for declared polynomial subspaces. Search trees and all evaluation triples separately checked.',
         'rows':rows,'frontier_control':limited_review,'corrupted_certificates_rejected':invalid,
         'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('run_rank_three.py','rank_three_pinning.py','rank_three_verifier.py')},
         'certificate_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in names}}
    (HERE/'rank_three_scale.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(rows),'exact_cases':sum(r['status']=='exact' for r in rows),'invalid_certificates':len(invalid)}),flush=True)


if __name__=='__main__':main()
