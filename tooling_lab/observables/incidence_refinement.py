#!/usr/bin/env python3
"""Iteration 3: retain incidence between quartic correlations.

For |C|=6 each quartic subset is indexed by its omitted pair. Its
correlation defines a weighted K6 edge. The multiset of sorted incident
edge lists is an inexpensive permutation invariant. This is familiar
graph-invariant / color-refinement machinery specialized to the exact
arithmetic target, not an invented graph algorithm.
"""
from hashlib import sha256
from itertools import combinations,islice
import json
from pathlib import Path
import time

import numpy as np
from fiber_audit import batch_observables,direct_record,enumerate_sets,quadratic_character,summarize_fibers

HERE=Path(__file__).resolve().parent
BASE=('m2','m4','b2','b4','energy','quartic_deck')


def incidence_signature(c,p,chi):
    assert len(c)==6
    row=chi[(np.arange(p)[None,:]-np.array(c)[:,None])%p]
    edge=np.zeros((6,6),dtype=np.int64)
    for i,j in combinations(range(6),2):
        q=[k for k in range(6) if k not in (i,j)]
        edge[i,j]=edge[j,i]=np.prod(row[q,:],axis=0).sum()
    signature=tuple(sorted(tuple(sorted(int(edge[i,j]) for j in range(6) if j!=i)) for i in range(6)))
    return signature,edge.tolist()


def case(p,samples):
    start=time.monotonic();chi=quadratic_character(p)
    source=iter(enumerate_sets(p,6,samples,20260905)); fibers={};count=0
    while batch:=list(islice(source,128)):
        for r in batch_observables(batch,p,chi):
            fibers.setdefault(tuple(r[f] for f in BASE),[]).append(r);count+=1
    unresolved=[r for records in fibers.values() if len({r['m6'] for r in records})>1 for r in records]
    before=summarize_fibers(unresolved,BASE,15*p*6**3) if unresolved else None
    for r in unresolved:
        r['quartic_incidence'],graph=incidence_signature(r['C'],p,chi)
        # Degree-three graph contraction of already evaluated quartics.
        # Python ints avoid any hidden fixed-width overflow in the products.
        r['quartic_triangle']=sum(graph[i][j]*graph[j][k]*graph[k][i]
                                  for i in range(6) for j in range(6) for k in range(6))
    after=summarize_fibers(unresolved,BASE+('quartic_incidence',),15*p*6**3) if unresolved else None
    after_triangle=summarize_fibers(unresolved,BASE+('quartic_incidence','quartic_triangle'),15*p*6**3) if unresolved else None
    separated_witnesses=[]
    if before:
        assert after['total_record_weighted_range']<=before['total_record_weighted_range']
        # Store a deterministic original collision and its separating graphs.
        l,r=before['witness']
        lg=incidence_signature(l['C'],p,chi);rg=incidence_signature(r['C'],p,chi)
        assert all(direct_record(w['C'],p)['m6']==w['m6'] for w in (l,r))
        assert incidence_signature(tuple(reversed(l['C'])),p,chi)[0]==lg[0]
        separated_witnesses.append({'sets':[l['C'],r['C']],'m6':[l['m6'],r['m6']],
                                    'omitted_pair_weighted_graphs':[lg[1],rg[1]],
                                    'rooted_signatures':[lg[0],rg[0]],
                                    'separated':lg[0]!=rg[0]})
    return {'p':p,'n':6,'sampling':'all normalized sets' if samples is None else 'seeded normalized sample',
            'records':count,'original_fibers':len(fibers),
            'unresolved_records':len(unresolved),
            'before':before,'after':after,'after_triangle':after_triangle,'witnesses':separated_witnesses,
            'elapsed_seconds':round(time.monotonic()-start,3),
            'interpretation':'No remaining finite collision is not a theorem of determination; rich fingerprints may effectively identify inputs.'}


def main():
    cases=[]
    for p,samples in ((61,None),(257,8192)):
        r=case(p,samples);cases.append(r)
        print(json.dumps({'p':p,'records':r['records'],'unresolved_records':r['unresolved_records'],
                          'before_ambiguous':r['before']['ambiguous_fibers'],
                          'after_ambiguous':r['after']['ambiguous_fibers'],
                          'after_triangle_ambiguous':r['after_triangle']['ambiguous_fibers'],
                          'seconds':r['elapsed_seconds']}),flush=True)
    (HERE/'results_incidence.json').write_text(json.dumps({'status':'finite refinement only',
                                                         'script_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
                                                         'cases':cases},indent=2)+'\n')


if __name__=='__main__':main()
