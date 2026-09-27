#!/usr/bin/env python3
"""Apply the query-cover interface to every preceding actual rank-three space."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
from time import perf_counter
import numpy as np
from query_cover import find_partition
from cover_verifier import validate,evaluate

HERE=Path(__file__).resolve().parent


def main():
    sourcepath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';source=json.loads(sourcepath.read_text())
    previouspath=HERE.parent/'round12/actual_results.json';previous=json.loads(previouspath.read_text())
    p,n,k,s=(source[x] for x in ('p','n','k','s'));nodes=source['nodes']
    for node in nodes:
        assert [evaluate(node['coefficients'],x,p) for x in source['domain']]==node['codeword']
    agreement=np.fromiter((sum(1<<i for i in A) for A in combinations(range(n),s)),dtype=np.uint64,count=comb(n,s))
    rows=[];failed=[];histogram=Counter();started=perf_counter();queries=local_subsets=0
    for index,old in enumerate(previous['rows']):
        ids=old['node_indices'];polys=[nodes[i]['coefficients']+[0]*(k-len(nodes[i]['coefficients'])) for i in ids]
        data={'p':p,'n':n,'k':k,'s':s,'domain':source['domain'],'origin':polys[0],
              'basis':[[(a-b)%p for a,b in zip(poly,polys[0])] for poly in polys[1:]]}
        searches=[];found=None
        for sizes in ([4]*4,[5,5,6]):
            search=find_partition(data,sizes,seed=1306+index,max_attempts=64)
            certificate=search.pop('certificate',None);searches.append({'sizes':sizes,**search})
            if certificate is not None:found=certificate;break
        if found is None:failed.append({'node_indices':ids,'searches':searches});continue
        review=validate(found)
        # Separate global subset coverage, without using local capacities or
        # the partition proof. The verifier above checked every query's rank.
        hits=np.zeros(len(agreement),dtype=bool)
        for block in found['blocks']:
            for Q in block['queries']:
                mask=sum(1<<i for i in Q);hits|=(agreement&mask)==mask
        assert hits.all()
        queries+=review['queries_rank_checked'];local_subsets+=review['literal_local_subsets_checked']
        histogram[found['query_count']]+=1
        rows.append({'node_indices':ids,'searches':searches,'certificate':found,'review':review,'literal_agreement_sets_checked':len(agreement)})
        if len(rows)%500==0:print(json.dumps({'covered':len(rows),'not_found':len(failed),'seconds':perf_counter()-started}),flush=True)
    assert len(rows)+len(failed)==previous['rank_three_quadruples']==2894
    out={'status':'passed','scope':'Bounded partition searches on all preceding actual rank-three candidate spaces. Every found cover is independently validated and checked against all11-subsets. Failure to find a partition, if any, is not an impossibility result. Source nodes can have different received scalars; no cluster-cover theorem is implied.',
         'spaces_attempted':2894,'spaces_covered':len(rows),'spaces_not_found':len(failed),
         'query_count_histogram':sorted(histogram.items()),'queries_rank_checked':queries,
         'literal_local_subsets_checked':local_subsets,'literal_agreement_sets_checked':len(rows)*len(agreement),
         'rows':rows,'not_found':failed,'seconds':perf_counter()-started,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['family_preflight.py','query_cover.py','cover_verifier.py']},
         'input_sha256':{'../round12/actual_results.json':sha256(previouspath.read_bytes()).hexdigest(),
                         '../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json':sha256(sourcepath.read_bytes()).hexdigest()}}
    (HERE/'family_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('rows','not_found','source_sha256','input_sha256')}),flush=True)


if __name__=='__main__':main()
