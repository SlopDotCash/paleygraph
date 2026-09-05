#!/usr/bin/env python3
"""Integer and cyclotomic checks of the two-anchor conic representation."""
from collections import Counter
from functools import lru_cache
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def trim(a):
    while len(a)>1 and a[-1]==0:a.pop()
    return a


def divide(a,b):
    a=a.copy(); q=[0]*max(1,len(a)-len(b)+1)
    assert b[-1]==1
    while len(a)>=len(b) and a!=[0]:
        k=len(a)-len(b);c=a[-1];q[k]=c
        for j,v in enumerate(b):a[k+j]-=c*v
        trim(a)
    return trim(q),trim(a)


@lru_cache(None)
def cyclotomic(n):
    a=[-1]+[0]*(n-1)+[1]
    for d in range(1,n):
        if n%d==0:
            a,r=divide(a,list(cyclotomic(d)));assert r==[0]
    return tuple(a)


def reduce_poly(a,n):
    return divide(trim(a.copy()),list(cyclotomic(n)))[1]


def cyclic_product(a,b,n):
    out=[0]*n
    for i,x in enumerate(a):
        if x:
            for j,y in enumerate(b):
                if y:out[(i+j)%n]+=x*y
    return out


def field(p):
    assert p%4==1 and all(p%d for d in range(2,isqrt(p)+1))
    chi=[0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
    n=p-1
    g=next(g for g in range(2,p) if len({pow(g,j,p) for j in range(n)})==n)
    ts=[pow(g,j,p) for j in range(n)]
    logs={t:j for j,t in enumerate(ts)}
    return chi,ts,logs,g


def one_field(p,do_fourier):
    chi,ts,logs,g=field(p);n=p-1
    inv=lambda t:pow(t,-1,p)
    xs=[pow((t+inv(t))*inv(2)%p,2,p) for t in ts]
    R=[x for x in range(p) if chi[x]==1 and chi[(x-1)%p]==1]
    assert len(R)==(p-5)//4
    fibers=Counter(xs)
    assert set(fibers)==set(R)|{0,1}
    assert fibers[0]==fibers[1]==2 and all(fibers[x]==4 for x in R)
    good=[j for j,x in enumerate(xs) if x not in (0,1)]
    assert len(good)==n-4
    f=[chi[(t*t-1)%p] for t in ts]
    M=np.array([[f[(i-j)%n]*f[(i+j)%n] for j in range(n)] for i in range(n)],dtype=object)
    literal=np.array([[chi[(x-y)%p] for y in xs] for x in xs],dtype=object)
    assert np.array_equal(M,literal)
    # Integer columns of norm 2: each R fiber and the +/- anchor combinations.
    cols=[]
    for x in R:cols.append([int(y==x) for y in xs])
    cols.append([int(y==0)+int(y==1) for y in xs])
    cols.append([int(y==0)-int(y==1) for y in xs])
    Q=np.array(cols,dtype=object).T
    m=len(R)
    S=np.array([[chi[(x-y)%p] for y in R] for x in R],dtype=object).reshape((m,m))
    B=np.zeros((m+2,m+2),dtype=object)
    B[:m,:m]=4*S;B[:m,m]=4;B[m,:m]=4;B[m,m]=2;B[m+1,m+1]=-2
    assert np.array_equal(Q.T@Q,4*np.eye(m+2,dtype=object))
    assert np.array_equal(M@Q,Q@B)
    assert np.array_equal(4*M,Q@B@Q.T)
    F=np.array([[int(xs[i]==x) for x in R] for i in good],dtype=object).reshape((len(good),m))
    Mg=M[np.ix_(good,good)]
    assert np.array_equal(F.T@F,4*np.eye(m,dtype=object))
    assert np.array_equal(F.T@Mg@F,16*S)
    assert np.array_equal(F.T@np.ones((len(good),len(good)),dtype=object)@F,
                          16*np.ones((m,m),dtype=object))
    # All six anharmonic transformations, tagged by their permutation signs.
    maps=[(1,lambda x:x),(-1,lambda x:(1-x)%p),(-1,lambda x:inv(x)),
          (-1,lambda x:x*inv(x-1)%p),(1,lambda x:inv(1-x)),
          (1,lambda x:(x-1)*inv(x)%p)]
    index={x:i for i,x in enumerate(R)}
    Ps=[];fixed=[]
    for sign,fun in maps:
        perm=[index[fun(x)%p] for x in R]
        P=np.zeros((m,m),dtype=object)
        for j,i in enumerate(perm):P[i,j]=1
        assert np.array_equal(P.T@S@P,S)
        Ps.append((sign,P));fixed.append(int(np.trace(P)))
    t=int(chi[2]==1);s=2*int(p%3==1)
    assert fixed==[m,t,t,t,s,s]
    T=sum((P for sign,P in Ps),np.zeros((m,m),dtype=object))
    A=sum((sign*P for sign,P in Ps),np.zeros((m,m),dtype=object))
    C=6*np.eye(m,dtype=object)-T-A
    for Z in (T,A,C):
        assert np.array_equal(Z@Z,6*Z)
        assert np.array_equal(S@Z,Z@S)
    assert np.all(T@A==0) and np.all(T@C==0) and np.all(A@C==0)
    multiplicities=[(m+3*t+2*s)//6,(m-3*t+2*s)//6,(m-s)//3]
    assert [int(np.trace(T)),int(np.trace(A)),int(np.trace(C))]==[6*multiplicities[0],6*multiplicities[1],12*multiplicities[2]]
    out={'p':p,'generator':g,'m':m,'conic_entry_checks':n*n,
         'exact_quotient_checks':6,'anharmonic_symmetries':6,
         'S3_multiplicities_trivial_sign_standard':multiplicities}
    if do_fourier:
        # F_j=sum_r f(g^r) zeta^(-jr), represented exactly in Z[zeta_n].
        Fs=[]
        for j in range(n):
            a=[0]*n
            for r,val in enumerate(f):a[(-j*r)%n]+=val
            Fs.append(a)
            if j%2:assert reduce_poly(a,n)==[0]
        # Independent Jacobi-sum identity F_(2j)=J(rho^-j,chi)+J(chi*rho^-j,chi).
        for j in range(n//2):
            right=[0]*n
            for v in ts:right[(-j*logs[v])%n]+=(1+chi[v])*chi[(1-v)%p]
            assert reduce_poly([a-b for a,b in zip(Fs[2*j],right)],n)==[0]
        assert reduce_poly(Fs[0],n)==[-2]
        offdiag=[];checks=0
        for a in range(n):
            for b in range(n):
                left=[0]*n
                for i in range(n):
                    for j in range(n):left[(-a*i+b*j)%n]+=int(M[i,j])
                right=[0]*n
                for j in range(n):
                    if (2*j-a-b)%n==0:
                        prod=cyclic_product(Fs[j],Fs[(j-b)%n],n)
                        right=[x+y for x,y in zip(right,prod)]
                assert reduce_poly([x-y for x,y in zip(left,right)],n)==[0]
                value=reduce_poly(left,n)
                if a%2 or b%2 or (a-b)%4:assert value==[0]
                if a!=b and value!=[0] and len(offdiag)<3:
                    offdiag.append({'a':a,'b':b,'N_times_Mellin_entry_coefficients':value})
                checks+=1
        out.update({'cyclotomic_DFT_entry_checks':checks,'jacobi_identity_checks':n//2,
                    'off_diagonal_witnesses':offdiag})
    return out


def main():
    cases=[one_field(p,p<=41) for p in (5,13,17,29,41,53,61,73,89,97,101,109,137,149)]
    sources=['research/parallel3-conic-2026-09-04.md',
             'research/parallel2-spectral-transfer-2026-09-04.md']
    result={'status':'Exact representations and symmetry checks; no improved uniform spectral bound.',
            'arithmetic':'Integer matrices and polynomial reduction modulo cyclotomic polynomials; no floating-point acceptance.',
            'cases':cases,'input_sha256':{s:sha256((ROOT/s).read_bytes()).hexdigest() for s in sources}}
    out=ROOT/'results/parallel3_conic_2026_09_04.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'fields':len(cases),'conic_entries':sum(c['conic_entry_checks'] for c in cases),
                      'DFT_entries':sum(c.get('cyclotomic_DFT_entry_checks',0) for c in cases),
                      'result':str(out)}))


if __name__=='__main__':main()
