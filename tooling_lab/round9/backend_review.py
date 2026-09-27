#!/usr/bin/env python3
"""Separate count-formula/dense-matrix readback of C++ DP and NTT statistics."""
from hashlib import sha256
import json
from pathlib import Path
from preflight import literal
from two_insertions import from_matrix,query

HERE=Path(__file__).resolve().parent


def main():
    cases=[];entries=0
    inputs=[(5,[0,1]),(5,[0,1,2]),(13,list(range(6))),(17,[0,1,2,3,4,6,10]),
            (29,[0,1,3,4,7,11,18,20]),(101,list(range(12))),
            (257,sorted({x*x%257 for x in range(1,257)})[:64])]
    for q,c in inputs:
        S=literal.literal_signs(q)
        for d in range(min(6,len(c))+1):
            a=c[:2];r=query(q,c,a,d);v=from_matrix(S,c,a,d)
            for k,value in v['statistics'].items():assert r['statistics'][k]==value,(q,d,k)
            for k,value in v.items():
                if k!='statistics':assert r[k]==value,(q,d,k)
            entries+=len(c)**2+5*len(c)+6
            cases.append({'q':q,'n':len(c),'degree':d,'convolution_skipped':r['statistics']['convolution']['skipped_exact_zero_Q']})
    invalid=[(9,[0,1], [0,1],2),(13,[0,1,1],[0,1],2),(13,[0,1,2],[0,3],2),
             (13,[0,1,2],[0,0],2),(13,list(range(12)),[0,1],2),(13,[0,1],[0,1],3),
             (13,[0,1],[0,1],-1)]
    for args in invalid:
        try:query(*args)
        except ValueError:pass
        else:raise AssertionError(('accepted invalid query',args))
    invalid_matrices=[]
    S=literal.literal_signs(13)
    for i,j,value in [(0,0,1),(0,1,1.5),(0,1,0)]:
        bad=[row[:] for row in S];bad[i][j]=value;invalid_matrices.append(bad)
    for bad in invalid_matrices:
        try:from_matrix(bad,list(range(6)),[0,1],6)
        except ValueError:pass
        else:raise AssertionError('accepted invalid matrix')
    out={'status':'passed','scope':'Complete dense direct arithmetic comparison of selected statistics, Q and final moments, with coefficient counts distinct from C++ DP.',
         'cases':cases,'scalar_statistics_compared':entries,'invalid_queries_rejected':len(invalid),
         'invalid_matrices_rejected':len(invalid_matrices),
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('backend_review.py','two_insertions.py','two_insertions_backend.cpp','two_insertions_backend','preflight.py')},
         'input_sha256':{'../round7/local_edits/review_results.py':sha256((HERE.parent/'round7/local_edits/review_results.py').read_bytes()).hexdigest()}}
    (HERE/'backend_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('source_sha256','input_sha256','cases')}),flush=True)


if __name__=='__main__':main()
