#!/usr/bin/env python3
"""Apply the capacitated query to every rank-three saved candidate quadruple."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import sys
from time import perf_counter
import numpy as np
sys.path.insert(0,str(Path(__file__).resolve().parent.parent/'round11'))
from rank_three_verifier import rank
from capacitated_pinning import optimize

HERE=Path(__file__).resolve().parent


def main():
    path=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json'
    source=json.loads(path.read_text());p,n,s,k=(source[x] for x in ('p','n','s','k'));nodes=source['nodes']
    for node in nodes:
        poly=node['coefficients'];assert len(poly)<=k
        assert node['codeword']==[sum(a*pow(x,j,p) for j,a in enumerate(poly))%p for x in source['domain']]
        support=[i for i,(v,u,w) in enumerate(zip(node['codeword'],source['u0'],source['u1'])) if v==(u+node['scalar']*w)%p]
        assert support==node['agreement_support'] and len(support)>=s
    agreement=np.fromiter((sum(1<<i for i in A) for A in combinations(range(n),s)),dtype=np.uint64,count=comb(n,s))
    rows=[];dependent=0;catalogues={};class_hist=Counter();started=perf_counter()
    for ids in combinations(range(len(nodes)),4):
        words=[nodes[i]['codeword'] for i in ids]
        columns=[tuple((word[j]-words[0][j])%p for word in words[1:]) for j in range(n)]
        if rank(columns,p)<3:dependent+=1;continue
        certificate=optimize(p,columns,s);class_hist[len(certificate['projective_classes'])]+=1
        # The verifier uses all raw coordinate triples and every raw agreement
        # set, without projective normalization, corners or subset transforms.
        dependent_triples=tuple(sum(1<<i for i in t) for t in combinations(range(n),3) if rank([columns[i] for i in t],p)<3)
        if dependent_triples not in catalogues:
            counts=np.zeros(len(agreement),dtype=np.int64)
            for triple in dependent_triples:counts+=(agreement&triple)==triple
            catalogues[dependent_triples]=comb(s,3)-int(counts.max())
        actual=catalogues[dependent_triples];assert certificate['minimum_injective_triples']==actual
        witness=set(certificate['worst_agreement_set'])
        assert actual==sum(rank([columns[i] for i in t],p)==3 for t in combinations(witness,3))
        assert len(witness)==s
        rows.append({'node_indices':list(ids),'certificate':certificate,'literal_minimum':actual})
        if len(rows)%500==0:print(json.dumps({'reviewed':len(rows),'seconds':perf_counter()-started}),flush=True)
    assert len(rows)==2894 and dependent==166
    old=json.loads((HERE.parent/'round11/rank_three_actual.json').read_text())
    byids={tuple(r['node_indices']):r for r in rows}
    for r in old['small_cases']:assert byids[tuple(r['node_indices'])]['literal_minimum']==r['independent_review']['interval'][0]
    out={'status':'actual_inputs_passed','scope':'All rank-three spaces spanned by quadruples of the saved18 actual nodes, including zeros and parallel columns. Nodes can belong to different received-word scalars; this is not a decoder or a cluster coverage theorem.',
         'rank_three_quadruples':len(rows),'dependent_quadruples':dependent,'unique_labeled_dependency_catalogues':len(catalogues),
         'agreement_sets_per_catalogue':len(agreement),'literal_agreement_sets_evaluated':len(catalogues)*len(agreement),
         'coordinate_triples_checked':len(rows)*comb(n,3),'previous_simple_cases_matched':len(old['small_cases']),
         'minimum_over_all_rank_three_spaces':min(r['literal_minimum'] for r in rows),'nonzero_class_count_histogram':[[a,b] for a,b in sorted(class_hist.items())],
         'rows':rows,'seconds':perf_counter()-started,
         'source_sha256':{'actual_preflight.py':sha256(Path(__file__).read_bytes()).hexdigest(),'capacitated_pinning.py':sha256((HERE/'capacitated_pinning.py').read_bytes()).hexdigest(),
                          '../round11/rank_three_verifier.py':sha256((HERE.parent/'round11/rank_three_verifier.py').read_bytes()).hexdigest()},
         'input_sha256':{'../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(path.read_bytes()).hexdigest(),
                         '../round11/rank_three_actual.json':sha256((HERE.parent/'round11/rank_three_actual.json').read_bytes()).hexdigest()}}
    (HERE/'actual_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('rows','source_sha256','input_sha256')}),flush=True)


if __name__=='__main__':main()
