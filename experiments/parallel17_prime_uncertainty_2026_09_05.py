#!/usr/bin/env python3
"""Exact prime-field determinant gaps and their quantitative limitation.

No floating-point values are used for acceptance. This is not a proof of
the Paley conjecture; the uniform arguments are in the accompanying note.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, permutations
from math import comb, isqrt
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


def qmul(x,y,p):
    a,b=x;c,d=y
    return a*c+p*b*d,a*d+b*c


def qsub(x,y):return x[0]-y[0],x[1]-y[1]


def qdiv(x,y,p):
    a,b=qmul(x,(y[0],-y[1]),p)
    norm=y[0]**2-p*y[1]**2
    assert norm
    aa,ra=divmod(a,norm);bb,rb=divmod(b,norm)
    assert ra==rb==0,('inexact quadratic division',x,y,p)
    return aa,bb


def both_positive(x,p):
    a,b=x
    return a>0 and a*a>p*b*b


def qdet(matrix,p,positive=False):
    a=[row[:] for row in matrix];n=len(a)
    if not n:return (1,0),0
    old=(1,0);sign=1;leading=0
    for k in range(n):
        if a[k][k]==(0,0):
            if positive:raise AssertionError('zero positive pivot')
            j=next((j for j in range(k+1,n) if a[j][k]!=(0,0)),None)
            if j is None:return (0,0),leading
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        if positive:
            assert both_positive(pivot,p),('nonpositive conjugate minor',p,k,pivot)
            leading+=1
        if k==n-1:return (sign*pivot[0],sign*pivot[1]),leading
        for i in range(k+1,n):
            for j in range(k+1,n):
                a[i][j]=qdiv(qsub(qmul(pivot,a[i][j],p),qmul(a[i][k],a[k][j],p)),old,p)
        old=pivot


def idet(matrix):
    a=[[int(x) for x in row] for row in matrix];n=len(a)
    if not n:return 1
    old=1;sign=1
    for k in range(n-1):
        if not a[k][k]:
            j=next((j for j in range(k+1,n) if a[j][k]),None)
            if j is None:return 0
            a[k],a[j]=a[j],a[k];sign=-sign
        pivot=a[k][k]
        for i in range(k+1,n):
            for j in range(k+1,n):
                value,rem=divmod(pivot*a[i][j]-a[i][k]*a[k][j],old)
                assert rem==0
                a[i][j]=value
        old=pivot
    return sign*a[-1][-1]


def character(p):
    assert p%4==1 and all(p%d for d in range(2,isqrt(p)+1))
    return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]


def paley(p):
    chi=character(p)
    S=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=np.int64)
    assert np.array_equal(S@S,p*np.eye(p,dtype=np.int64)-np.ones((p,p),dtype=np.int64))
    return S


COUNTS=Counter()


def compression(p,S,C,shift=False):
    m=len(C);assert 1<=m<=(p-1)//2 and len(set(C))==m
    T=S[np.ix_(C,C)];B=[[(p-1,0) if i==j else (-1,int(T[i,j]))
                        for j in range(m)] for i in range(m)]
    (a,b),leading=qdet(B,p,positive=True)
    COUNTS['positive_conjugate_leading_minors']+=2*leading
    norm=a*a-p*b*b;factor=4**m*p**(m-1)
    k,rem=divmod(norm,factor)
    assert rem==0 and k>0
    Delta=Fraction(p-m,p)
    theta=Fraction(4**m*k,p**(m+1))/Delta**2
    assert 0<theta<=1
    H=(p-m)*np.eye(m,dtype=np.int64)+np.ones((m,m),dtype=np.int64)
    TH=T@H
    independent=idet(p*(p-m)**2*np.eye(m,dtype=object)-TH.astype(object)@TH.astype(object))
    assert theta==Fraction(independent,(p*(p-m)**2)**m)
    COUNTS['independent_integer_determinants']+=1
    row_sums=T@np.ones(m,dtype=np.int64)
    exact_trace=Fraction(int(np.trace(TH@TH)),(p-m)**2)
    predicted=m*(m-1)+Fraction(2*int(row_sums@row_sums),p-m)+Fraction(int(row_sums.sum())**2,(p-m)**2)
    assert exact_trace==predicted
    t=exact_trace/p
    assert 0<=t<m
    amgm=(1-t/m)**m
    assert theta<=amgm
    COUNTS['second_moment_and_rational_decay_checks']+=1
    delta=Fraction(4**(m-1),p**m*(p-m));data_gap=k*delta
    assert data_gap==Delta*theta/4
    if shift:
        # 2p*(X-data_gap*I), with the rational denominator cleared.
        d=data_gap.denominator;n=data_gap.numerator
        shifted=[[(d*x[0]-(2*p*n if i==j else 0),d*x[1])
                  for j,x in enumerate(row)] for i,row in enumerate(B)]
        _,checked=qdet(shifted,p,positive=True)
        COUNTS['shifted_positive_conjugate_minors']+=2*checked
    if m<=3:
        direct=(0,0)
        for perm in permutations(range(m)):
            sign=(-1)**sum(perm[i]>perm[j] for i in range(m) for j in range(i+1,m))
            term=(sign,0)
            for i,j in enumerate(perm):term=qmul(term,B[i][j],p)
            direct=(direct[0]+term[0],direct[1]+term[1])
        assert direct==(a,b)
        COUNTS['permutation_determinants']+=1
    COUNTS['compressions']+=1
    return dict(p=p,m=m,determinant_B=[str(a),str(b)],integer_k=str(k),theta=str(theta),
                trace_A_squared=str(exact_trace),theta_amgm_upper=str(amgm),
                uniform_rational_gap=str(delta),data_rational_gap=str(data_gap),
                certificate_gap_upper=str(Delta*theta/2),shifted_gap_checked=shift)


def poly_reduce(a,p):
    b=[0]*p
    for j,x in enumerate(a):b[j%p]+=x
    return [x-b[-1] for x in b[:-1]]


def poly_mul(a,b,p):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):out[i+j]+=x*y
    return poly_reduce(out,p)


def interval_projection(p):
    h=(p-1)//4;freq=set(range(h+1,3*h+1));r=(p-1)//2
    assert len(freq)==r and {(-x)%p for x in freq}==freq and 0 not in freq
    kernels=[]
    for j in range(p):
        coeff=[0]*p
        for x in freq:coeff[j*x%p]+=1
        kernels.append(poly_reduce(coeff,p))
    assert kernels[0]==[r]+[0]*(p-2)
    assert [sum(k[j] for k in kernels) for j in range(p-1)]==[0]*(p-1)
    for x in range(p):
        total=[0]*(p-1)
        for j in range(p):
            term=poly_mul(kernels[j],kernels[(x-j)%p],p)
            total=[a+b for a,b in zip(total,term)]
        assert total==[p*a for a in kernels[x]]
        assert kernels[x]==kernels[-x%p]
    chi=character(p);gauss=poly_reduce(chi,p)
    assert poly_mul(gauss,gauss,p)==[p]+[0]*(p-2)
    off=[2*x for x in kernels[1]];off[0]+=1
    scaled_sign_entry=poly_mul(gauss,off,p)
    has_paley_entry=scaled_sign_entry in ([p]+[0]*(p-2),[-p]+[0]*(p-2))
    assert has_paley_entry==(p==5)
    COUNTS['cyclotomic_projection_convolution_rows']+=p
    return dict(p=p,rank=r,all_convolution_rows=p,real=True,annihilates_constants=True,
                zero_diagonal_sign_transform=True,gauss_square_checked=True,
                transformed_off_diagonal_entry_is_plus_or_minus_one=has_paley_entry)


def main():
    # Quadratic arithmetic and pivot edge cases, independent of Paley inputs.
    for p in (5,13):
        for x in [(0,0),(1,0),(0,1),(-7,3)]:
            for y in [(1,0),(2,1),(-3,2)]:assert qdiv(qmul(x,y,p),y,p)==x
        assert qdet([[(0,0),(1,0)],[(1,0),(0,0)]],p)[0]==(-1,0)
        assert qdet([[(1,0),(1,0)],[(1,0),(1,0)]],p)[0]==(0,0)
    assert idet([[0,1],[1,0]])==-1 and idet([[1,1],[1,1]])==0
    exhaustive=[]
    for p,maxm in [(5,2),(13,6),(17,5)]:
        S=paley(p);sizes=Counter();minimum_k={}
        for m in range(1,maxm+1):
            for rest in combinations(range(1,p),m-1):
                C=[0,*rest]
                result=compression(p,S,C,shift=C==list(range(m)))
                sizes[m]+=1
                minimum_k[m]=min(minimum_k.get(m,int(result['integer_k'])),int(result['integer_k']))
        exhaustive.append(dict(p=p,all_subsets_containing_zero_through_size=maxm,
                               count_by_size=dict(sizes),minimum_k_by_size={k:str(v) for k,v in minimum_k.items()}))
        print(json.dumps(dict(exhaustive_p=p,compressions=sum(sizes.values()))),flush=True)
    selected=[];boundaries=[]
    for p in [29,41,73,97,193,257]:
        S=paley(p);anchors=[0]
        for a in range(1,5):
            if len(anchors)<a:
                anchors.append(next(x for x in range(p) if all(S[x,y]==1 for y in anchors)))
            C=[x for x in range(p) if all(S[x,y]==1 for y in anchors)]
            if not C:continue
            row=compression(p,S,C,shift=len(C)<=32)
            row['anchors']=anchors[:];selected.append(row)
        if p==257:
            anchors=[0,1,62]
            assert all(S[x,y]==1 for x,y in combinations(anchors,2))
            C=[x for x in range(p) if all(S[x,y]==1 for y in anchors)]
            row=compression(p,S,C,shift=len(C)<=32)
            row['anchors']=anchors;selected.append(row)
        if p<=41:
            m=(p+1)//2;T=S[:m,:m]
            B=[[(p-1,0) if i==j else (-1,int(T[i,j])) for j in range(m)] for i in range(m)]
            det,_=qdet(B,p)
            assert det==(0,0)
            boundaries.append(dict(p=p,m=m,determinant=list(det)))
        print(json.dumps(dict(selected_p=p,last_m=selected[-1]['m'])),flush=True)
    intervals=[interval_projection(p) for p in (5,13,29)]
    leakage=[]
    for p in (257,1021,4093,12289):
        character(p);m=p//8;n=m-1
        # Integer binomial identities and exact rational certificate only;
        # the uniform trigonometric inequality is proved in the note.
        norm=comb(2*n,n)
        if p<=4093:assert sum(comb(n,j)**2 for j in range(n+1))==norm
        assert norm*(2*n+1)>=4**n
        bound=Fraction(2**n,norm);simple=Fraction(2*n+1,2**n)
        assert 0<bound<=simple<1
        leakage.append(dict(p=p,m=m,binomial_leakage_bound=str(bound),simple_leakage_bound=str(simple),
                            bound_less_than_one_millionth=bound<Fraction(1,10**6)))
    inputs=['research/parallel17-prime-uncertainty-2026-09-05.md',
            'experiments/parallel17_prime_uncertainty_2026_09_05.py',
            'sources/tao-uncertainty-math0308286v6.pdf','sources/tao-uncertainty-math0308286v6.html',
            'research/parallel16-square-field-barrier-2026-09-05.md']
    out=dict(status='Exact prime-field determinant gaps and quantitative limitation passed; full conjecture remains unproved.',
             arithmetic='Integer and quadratic-integer Bareiss determinants, both real embeddings, rational inequalities, and exact cyclotomic identities. No floating-point acceptance.',
             counts=dict(COUNTS),exhaustive=exhaustive,selected_localizations=selected,
             rank_boundary_singularities=boundaries,interval_projection_checks=intervals,
             interval_leakage_certificates=leakage,
             input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest=ROOT/'results/parallel17_prime_uncertainty_2026_09_05.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),counts=dict(COUNTS))),flush=True)


if __name__=='__main__':main()
