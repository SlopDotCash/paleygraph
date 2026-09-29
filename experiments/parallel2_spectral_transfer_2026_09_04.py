#!/usr/bin/env python3
"""Exact quadratic-field matrix checks for the aggregate trace reduction."""
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def mul(X,Y,p):
    a,b=X; c,d=Y
    return a@c+p*(b@d), a@d+b@c


def sub(X,Y):
    return X[0]-Y[0],X[1]-Y[1]


def scale(X,n):
    return n*X[0],n*X[1]


def trace(X):
    return int(np.trace(X[0])),int(np.trace(X[1]))


def field(p):
    assert p%4==1 and all(p%d for d in range(2,isqrt(p)+1))
    chi=[0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
    S=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    eye=np.eye(p,dtype=object); J=np.ones((p,p),dtype=object)
    assert np.array_equal(S@S,p*eye-J)
    assert np.all(S@J==0)
    # 2pP=pI-J+sqrt(p)S; check P^2=P after clearing denominators.
    N=(p*eye-J,S)
    diff=sub(mul(N,N,p),scale(N,2*p))
    assert all(np.all(x==0) for x in diff)
    assert trace(N)==(p*(p-1),0)
    return chi,S,eye,J


def one_case(p,R,b,S,eye,J,anchors=None,chi=None):
    m=len(R);r=(p-1)//2
    assert m<=r and b>=2
    e=np.array([int(i in R) for i in range(p)],dtype=object)
    diag=b*e-1
    W=(-J*diag,S*diag)
    smallS=S[np.ix_(R,R)]
    smallJ=np.ones((m,m),dtype=object)
    smallI=np.eye(m,dtype=object)
    H0=(2*smallI,np.zeros((m,m),dtype=object))
    H1=(-b*smallJ,b*smallS)
    H=[H0,H1]
    Wpower=(eye,np.zeros((p,p),dtype=object))
    traces=[]
    for k in range(1,9):
        Wpower=mul(Wpower,W,p)
        if k>=2:
            H.append(sub(mul(H1,H[-1],p),scale(H[-2],p*p*(b-1))))
        left=trace(Wpower); right=trace(H[k])
        residual=p**k*((p-r-m)+(-1)**k*(r-m))
        assert left==(right[0]+residual,right[1])
        traces.append({'k':k,'trace_W_integer_part':left[0],
                       'trace_W_sqrt_part':left[1],
                       'residual_integer':residual})
    rayleigh_checks=0
    if anchors is not None:
        assert b==2**len(anchors)
        # Build the sum of all nonempty mask diagonals literally.
        F=np.zeros(p,dtype=object)
        for size in range(1,len(anchors)+1):
            for Z in combinations(anchors,size):
                for x in range(p):
                    value=1
                    for z in Z:value*=chi[(x-z)%p]
                    F[x]+=value
        for z in anchors:F[z]-=b//2
        assert np.array_equal(F,diag)
        # Every checked extension is an actual clique inside R.
        for size in range(1,min(4,m)+1):
            count=0
            for C in combinations(R,size):
                if not all(chi[(x-y)%p]==1 for x,y in combinations(C,2)):continue
                ell=len(C)
                rational=sum(p*int(x==y)-1 for x in C for y in C)
                radical=sum(int(S[x,y]) for x in C for y in C)
                assert (rational,radical)==(p*ell-ell*ell,ell*(ell-1))
                rayleigh_checks+=1;count+=1
                if count==3:break
    return {'p':p,'b':b,'anchors':anchors,'R':R,'m':m,
            'trace_identity_checks':8,'rayleigh_numerator_checks':rayleigh_checks,
            'traces':traces}


def main():
    cases=[]
    for p in (5,13,17,29,41):
        chi,S,eye,J=field(p)
        for a in range(1,4):
            count=0
            for I in combinations(range(p),a):
                if not all(chi[(x-y)%p]==1 for x,y in combinations(I,2)):continue
                R=[x for x in range(p) if all(chi[(x-z)%p]==1 for z in I)]
                cases.append(one_case(p,R,2**a,S,eye,J,list(I),chi))
                count+=1
                if count==2:break
        for R,b in [([],2),([0],4),(list(range((p-1)//2)),8)]:
            cases.append(one_case(p,R,b,S,eye,J))
    report={'status':'Exact algebra checks only; aggregate cancellation hypothesis unproved.',
            'arithmetic':'Integer pairs u+v sqrt(p), using NumPy object matrices.',
            'case_count':len(cases),
            'projection_identity_fields':5,
            'aggregate_trace_identity_checks':sum(x['trace_identity_checks'] for x in cases),
            'anchor_expansion_checks':sum(x['anchors'] is not None for x in cases),
            'clique_rayleigh_numerator_checks':sum(x['rayleigh_numerator_checks'] for x in cases),
            'cases':cases}
    files=['research/parallel2-spectral-transfer-2026-09-04.md',
           'experiments/parallel2_spectral_transfer_2026_09_04.py',
           'sources/kunisky-2303.16475v1.html']
    report['source_sha256']={x:sha256((ROOT/x).read_bytes()).hexdigest() for x in files}
    target=ROOT/'results/parallel2_spectral_transfer_2026_09_04.json'
    target.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({k:v for k,v in report.items() if k not in ('cases','source_sha256')},indent=2))


if __name__=='__main__':main()
