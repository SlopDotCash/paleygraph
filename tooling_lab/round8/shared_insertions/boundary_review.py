#!/usr/bin/env python3
"""Degree/boundary checks via literal deletion DP and full insertion transforms."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
import sys
from covariance import query

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'round7/local_edits'))
from review_results import literal_signs


def coefficient(values,d):
    if d<0:return 0
    result=[1]+[0]*d
    for x in values:
        for j in range(d,0,-1):result[j]+=x*result[j-1]
    return result[d]


def main():
    cases=[(q,list(range(n)),d) for q,n in ((5,2),(13,6),(29,8),(101,64)) for d in range(min(6,n)+1)]
    residues=sorted({x*x%257 for x in range(1,257)})[:64]
    cases.append((257,residues,6));rows=[]
    for q,c,d in cases:
        n=len(c);S=literal_signs(q);r,g=query(q,c,d);outside=sorted(set(range(q))-set(c));m=len(outside)
        R=[[coefficient([S[x][b] for b in c if b!=a],d-1) for x in range(q)] for a in c]
        T=sum(coefficient([row[a] for a in c],d) for row in S)
        assert T==r['target']
        if q==257:assert max(map(abs,R[0]))==7028847
        V=[[sum(R[a][x]*S[x][b] for x in range(q)) for b in c] for a in range(n)]
        assert V==r['internal_contraction']
        edits=[[sum(R[a][x]*S[x][b] for x in range(q))-V[a][a] for b in outside] for a in range(n)]
        means=[F(sum(row),m) for row in edits]
        for a in range(n):
            for b in range(n):
                assert sum(x*y for x,y in zip(R[a],R[b]))==r['derivative_gram'][a][b]
                value=F(sum(x*y for x,y in zip(edits[a],edits[b])),m)-means[a]*means[b]
                assert value==F(r['covariance_numerator'][a][b],r['covariance_denominator'])
        rows.append({'q':q,'n':n,'degree':d,'derivative_gram_entries':n*n,'covariance_entries':n*n,'outside_insertion_transforms':n*m})
    invalid=[(9,[0,1],1),(13,[0,0],1),(13,[],0),(13,[0,13],1),(13,[0,1],3),(13,[0,1.5],1)]
    for args in invalid:
        try:query(*args)
        except ValueError:pass
        else:raise AssertionError(('invalid accepted',args))
    out={'status':'passed','scope':'Separate literal deletion DP and all outside insertion transforms; includes all degrees0..6 and an actual all-positive64-point row. No production Gram formula reuse in the oracle.',
         'cases':rows,'invalid_queries_rejected':len(invalid),
         'source_sha256':{f'round8/shared_insertions/{name}':sha256((HERE/name).read_bytes()).hexdigest() for name in ('boundary_review.py','covariance.py','covariance_backend.cpp','covariance_backend')},
         'dependency_sha256':{'round7/local_edits/review_results.py':sha256((LAB/'round7/local_edits/review_results.py').read_bytes()).hexdigest()}}
    (HERE/'boundary_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'cases':len(rows),'exact_covariance_entries':sum(x['covariance_entries'] for x in rows),'invalid_queries_rejected':len(invalid)}),flush=True)


if __name__=='__main__':main()
