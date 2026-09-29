#!/usr/bin/env python3
"""Exact checks of corrected pair kernels and aggregates across ranks."""
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import json
import numpy as np

from parallel_necklace_2026_09_04 import matrices
from parallel2_necklace_2026_09_04 import positive_principal_minors
from parallel3_necklace_2026_09_04 import convolution_kernels

ROOT=Path(__file__).resolve().parents[1]


def rational_rank(rows):
    a=[[Fraction(x) for x in row] for row in rows]
    n=len(a);m=len(a[0]) if n else 0;r=0
    for j in range(m):
        pivot=next((i for i in range(r,n) if a[i][j]),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];v=a[r][j]
        a[r]=[x/v for x in a[r]]
        for i in range(r+1,n):
            v=a[i][j]
            if v:a[i]=[x-v*y for x,y in zip(a[i],a[r])]
        r+=1
        if r==n:break
    return r


def local_stalk(r):
    T=np.eye(r,dtype=object)
    if r%2:T[0,0]=-1; allowed=[(i,j) for i in range(1,r) for j in range(1,r)]
    else:T[0,1]=1;allowed=[(i,j) for i in range(r) if i!=1 for j in range(r) if j!=0]
    equations=[]
    basis=[]
    for i in range(r):
        for j in range(r):
            E=np.zeros((r,r),dtype=object);E[i,j]=1
            equations.append(list((T@E-E@T).flat))
            if (i,j) in allowed:
                assert np.all(T@E==E) and np.all(E@T==E)
                basis.append(list(E.flat))
    invariant_dimension=r*r-rational_rank(np.array(equations,dtype=object).T.tolist())
    actual_dimension=rational_rank(basis)
    with_identity=rational_rank(basis+[list(np.eye(r,dtype=object).flat)])
    assert (invariant_dimension,actual_dimension,with_identity)==((r-1)**2+1,(r-1)**2,(r-1)**2+1)
    return dict(rank=r,inertia_type='reflection' if r%2 else 'transvection',
                invariant_tensor_dimension=invariant_dimension,
                actual_tensor_stalk_dimension=actual_dimension,
                identity_supplies_quotient=with_identity==invariant_dimension)


def one_field(p,R):
    chi,diag,cs=matrices(p);v=diag['A'][1:];n=p-1
    S=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    ks=[None,S]
    for r in range(2,R+1):ks.append(ks[-1]@cs['A'])
    conv=convolution_kernels(p,chi,R)
    for r in range(1,R+1):
        assert all(ks[r][1,t]==conv[r][t] for t in range(1,p))
    certs=[];entry_checks=0;exceptional=[]
    for r in range(1,R+1):
        eigen=(p-3)*p**r+2
        assert eigen%(p-1)==0;eigen//=p-1
        assert sum(conv[r][t]**2 for t in range(1,p))==eigen
        for s in range(r,R+1):
            V=v[:,None]*v[None,:]*ks[r][1:,1:]*ks[s][1:,1:]
            for x in range(1,p):
                for y in range(1,p):
                    t=y*pow(x,-1,p)%p
                    assert V[x-1,y-1]==chi[t]*conv[r][t]*conv[s][t]
                    entry_checks+=1
            if r==s:
                assert np.array_equal(V@v,eigen*v)
                H=V+p**(r-1)*np.eye(n,dtype=object)-p**(r-1)*np.outer(v,v)
                small=-2*sum(p**a for a in range(r-1))
                assert np.array_equal(H@v,small*v)
                exceptional.append(dict(rank=r,raw_eigenvalue=eigen,corrected_eigenvalue=small))
                if r==1:assert np.all(H==0);continue
                operator=H;bound=(2*r-2)**2*p**(2*r-1);kind='corrected_equal_rank'
            else:
                operator=V;bound=(r+s-1)**2*p**(r+s-1);kind='unequal_rank'
            pivots=positive_principal_minors(bound*np.eye(n,dtype=object)-operator@operator)
            certs.append(dict(r=r,s=s,kind=kind,positive_principal_minors=pivots))
    # Integer Gaussian coefficients b are converted to a_r=p^((r-1)/2)b_r.
    # Consequently every tested sum and cleared bound is an integer.
    aggregate_cases=[]
    for depth in range(1,R+1):
        vectors=[([(1,0)]*depth),([((-1)**r,0) for r in range(depth)]),
                 ([(r+1,(-1)**r) for r in range(depth)])]
        C=(3*depth*depth-depth-2)//2
        for coeff in vectors:
            norm=sum((a*a+b*b)*p**r for r,(a,b) in enumerate(coeff))
            abs_square=[]
            for t in range(1,p):
                real=sum(a*conv[r+1][t] for r,(a,b) in enumerate(coeff))
                imag=sum(b*conv[r+1][t] for r,(a,b) in enumerate(coeff))
                abs_square.append(real*real+imag*imag)
            for twist in ('trivial','quadratic'):
                weights=[1]*n if twist=='trivial' else [chi[t] for t in range(1,p)]
                total=sum(w*x for w,x in zip(weights,abs_square))
                main=p*norm if twist=='trivial' else 0
                excess=max(abs(total-main)-2*norm,0)
                assert excess*excess<=C*C*p*norm*norm
                # A separately expanded Gram quadratic form, using the real
                # parts of the Gaussian coefficient products.
                expanded=sum((a*c+b*d)*sum(weights[t-1]*conv[r+1][t]*conv[s+1][t]
                             for t in range(1,p))
                             for r,(a,b) in enumerate(coeff) for s,(c,d) in enumerate(coeff))
                assert expanded==total
                aggregate_cases.append(dict(depth=depth,coefficients=coeff,twist=twist,
                                            squared_coefficient_norm=norm,total=total))
    return dict(p=p,maximum_rank=R,kernel_entry_checks=entry_checks,
                norm_certificates=certs,exceptional_modes=exceptional,
                aggregate_cases=aggregate_cases)


def main():
    local=[local_stalk(r) for r in range(1,9)]
    cases=[one_field(p,7 if p==5 else 5) for p in (5,13,17,29,41)]
    inputs=['research/parallel4-kernel-aggregate-2026-09-04.md',
            'experiments/parallel4_kernel_aggregate_2026_09_04.py',
            'research/parallel3-necklace-2026-09-04.md',
            'sources/katz-g2-hypergeometric.pdf',
            'sources/katz-gauss-kloosterman-monodromy.pdf']
    out=ROOT/'results/parallel4_kernel_aggregate_2026_09_04.json'
    record=dict(status='Pairwise rank bounds and signed quadratic aggregate; full necklace aggregate unproved.',
                arithmetic='Integer matrices, exact Bareiss certificates, rational local ranks and cleared aggregate inequalities.',
                local_stalk_checks=local,cases=cases,
                input_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs})
    out.write_text(json.dumps(record,indent=2)+'\n')
    print(json.dumps(dict(fields=len(cases),kernel_entries=sum(c['kernel_entry_checks'] for c in cases),
                          norm_certificates=sum(len(c['norm_certificates']) for c in cases),
                          positive_minors=sum(len(d['positive_principal_minors']) for c in cases for d in c['norm_certificates']),
                          signed_aggregates=sum(len(c['aggregate_cases']) for c in cases),result=str(out))))


if __name__=='__main__':main()
