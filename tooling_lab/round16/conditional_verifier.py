#!/usr/bin/env python3
"""Constructor-free verification of a single-class conditional query cover.

Verifies direct polynomial evaluations, the exact received-label histogram,
every residual arc query via integer determinants, and the complete branch
partition of parameter space. Does not certify optimality of the cutoff search.
"""
from collections import Counter
from hashlib import sha256
import json
from math import isqrt
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'round14'))
from arc_verifier import validate,determinant,evaluate


def validate_conditional(c):
    if not __debug__:raise RuntimeError('verifier must run without Python-O')
    assert c['status']=='found' and c['schema']=='single_projective_class_conditional_cover_v1'
    plan=c['plan'];data=plan['input'];p,n,k,s,d=(data[x] for x in ('p','n','k','s','dimension'))
    assert all(type(v) is int for v in (p,n,k,s,d)) and d==3 and 3<=k<=n<=min(p,4096) and 3<=s<=n
    assert 2<=p<=2**31-1 and all(p%j for j in range(2,isqrt(p)+1))
    domain=data['domain'];origin=data['origin'];basis=data['basis'];received=plan['received']
    assert len(domain)==len(set(domain))==n and all(type(x) is int and 0<=x<p for x in domain)
    assert len(basis)==3 and all(len(b)<=k and all(type(v) is int and 0<=v<p for v in b) for b in [origin,*basis])
    assert len(received)==n and all(type(v) is int and 0<=v<p for v in received)
    cols=[tuple(evaluate(b,x,p) for b in basis) for x in domain];orig=[evaluate(origin,x,p) for x in domain]
    a=next(v for v in cols if any(v));b=next(v for v in cols if any(a[i]*v[j]%p!=a[j]*v[i]%p for i in range(3) for j in range(i+1,3)))
    assert any(determinant([a,b,v])%p for v in cols)
    w=plan['projective_direction'];pivot=plan['pivot'];F=plan['class_coordinates'];outside=plan['outside_coordinates']
    assert type(pivot) is int and 0<=pivot<3 and len(w)==3 and all(type(v) is int and 0<=v<p for v in w)
    assert w[pivot]==1 and not any(w[:pivot])
    actual=[i for i,v in enumerate(cols) if v[pivot]!=0 and all(v[j]==v[pivot]*w[j]%p for j in range(3))]
    assert F==actual and F and outside==[i for i in range(n) if i not in set(F)]
    labels=[[i,(received[i]-orig[i])*pow(cols[i][pivot],-1,p)%p] for i in F];counts=Counter(b for _,b in labels)
    assert labels==plan['normalized_received_labels'] and sorted(counts.items())==[tuple(row) for row in plan['label_counts']]
    best=plan['best'];h=best['cutoff'];assert type(h) is int and 0<=h<=max(counts.values())
    heavy=sorted(b for b,m in counts.items() if m>h)
    assert heavy==best['heavy_values'] and best['light_agreement_threshold']==s-h
    free=[j for j in range(3) if j!=pivot];seen=set();queries=0
    for part in c['parts']:
        cert=part['certificate'];review=validate(cert);target=cert['input'];queries+=review['queries_rank_checked']
        assert target['p']==p and target['n']==len(outside) and target['k']==min(k,len(outside))
        assert target['domain']==[domain[i] for i in outside]
        if part['kind']=='light':
            key=('light',);assert s-h<=len(outside) and part['excluded_values']==heavy and target['s']==s-h and target['dimension']==3
            expected_origin=[orig[i] for i in outside];expected_basis=[[cols[i][j] for i in outside] for j in range(3)]
        else:
            assert part['kind']=='fibre';value=part['value'];key=('fibre',value)
            assert value in heavy and part['class_agreements']==counts[value] and part['free_parameter_indices']==free
            assert target['s']==s-counts[value] and target['dimension']==2
            expected_origin=[(orig[i]+value*cols[i][pivot])%p for i in outside]
            expected_basis=[[(cols[i][j]-w[j]*cols[i][pivot])%p for i in outside] for j in free]
        assert key not in seen;seen.add(key)
        assert [evaluate(target['origin'],x,p) for x in target['domain']]==expected_origin
        assert [[evaluate(poly,x,p) for x in target['domain']] for poly in target['basis']]==expected_basis
    expected={('fibre',b) for b in heavy if s-counts[b]<=len(outside)}
    if s-h<=len(outside):expected.add(('light',))
    assert seen==expected and queries==c['query_count'] and queries<=250000
    # For all b outside heavy, m_b<=h, including field values absent from the
    # histogram (m_b=0). Thus the light residual threshold is sufficient. For
    # each heavy b, m_b is exact, and the affine lift fixes theta.w=b. These
    # disjoint cases cover every parameter; impossible residual sizes are empty.
    return {'status':'verified','queries_rank_checked':queries,'branches_checked':len(seen),
            'projective_class_size':len(F),'outside_coordinates':len(outside),'routing_complete_in_supplied_space':True,
            'cutoff_optimality_checked':False}


def main():
    source=HERE/'conditional_results.json';results=json.loads(source.read_text());rows=[];bindings={source.name:sha256(source.read_bytes()).hexdigest()}
    for row in results['cases']:
        path=HERE/row['file'];cert=json.loads(path.read_text());checked=validate_conditional(cert)
        rows.append({'case':row['case'],'review':checked});bindings[path.name]=sha256(path.read_bytes()).hexdigest();print(json.dumps(rows[-1]),flush=True)
    out={'status':'passed','scope':'Every constituent query independently rank-checked, all affine lifts and received histograms checked, and the exhaustive heavy/light routing certificate validated. No candidate-decoding experiment or globally optimal conditional family is certified here.',
         'cases':rows,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['conditional_verifier.py','../round14/arc_verifier.py']},'input_sha256':bindings}
    (HERE/'conditional_review.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
