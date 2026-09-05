#!/usr/bin/env python3
"""Exact tests of the twisted Legendre operator reduction at growing length."""
from hashlib import sha256
from itertools import permutations
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, direct_trace, literal, trace

ROOT=Path(__file__).resolve().parents[1]


def positive_principal_minors(matrix):
    """Bareiss elimination: pivots are exact leading principal determinants."""
    a=[[int(x) for x in row] for row in matrix]
    n=len(a)
    previous=1
    determinants=[]
    for k in range(n):
        pivot=a[k][k]
        assert pivot>0, (k,pivot)
        determinants.append(pivot)
        for i in range(k+1,n):
            for j in range(k+1,n):
                numerator=a[i][j]*pivot-a[i][k]*a[k][j]
                assert numerator%previous==0
                a[i][j]=numerator//previous
        for i in range(k+1,n):
            a[i][k]=0
        previous=pivot
    return determinants


def field_case(p):
    chi,diag,cs=matrices(p)
    s=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    d=np.diag(diag['A'])
    k2=s@d@s
    w=s*k2
    cstar=cs['A'][1:,1:]
    wstar=w[1:,1:]
    ones=np.ones(p-1,dtype=object)
    assert np.array_equal(cstar@ones,-ones)
    assert np.array_equal(wstar@ones,2*ones)
    assert np.array_equal(cs['A'][0],np.zeros(p,dtype=object))
    assert np.array_equal(cs['A'][1:,0],ones)
    assert np.array_equal(w[0,1:],-ones)
    assert np.array_equal(w[1:,0],-ones)
    assert w[0,0]==0
    # Direct complete weighted anchor average of matrices.
    average=np.zeros((p,p),dtype=object)
    for b in range(p):
        cb=np.array([chi[(x-b)%p] for x in range(p)],dtype=object)[:,None]*s
        average+=chi[b]*(cb@d@cb)
    assert np.array_equal(average,w@cs['A'])
    elliptic=[sum(chi[z]*chi[(z-1)%p]*chi[(z-t)%p] for z in range(p)) for t in range(p)]
    h=[chi[(1-t)%p]*elliptic[t] for t in range(p)]
    for x in range(1,p):
        for y in range(1,p):
            assert w[x,y]==h[(y*pow(x,-1,p))%p]
    certs={}
    for sign in (-1,1):
        certs[str(sign)]=positive_principal_minors(2*p*np.eye(p-1,dtype=object)+sign*wstar)
    quadratic=diag['A'][1:]
    v=quadratic[:,None]*quadratic[None,:]*(k2[1:,1:]**2)
    hoperator=v+p*np.eye(p-1,dtype=object)-p*np.outer(quadratic,quadratic)
    assert np.array_equal(v@quadratic,(p*p-2*p-2)*quadratic)
    assert np.array_equal(hoperator@quadratic,-2*quadratic)
    h_certificate=positive_principal_minors(4*p**3*np.eye(p-1,dtype=object)-hoperator@hoperator)
    powers=[np.eye(p-1,dtype=object)]
    for _ in range(17):
        powers.append(powers[-1]@cstar)
    values={}
    for r in range(1,17):
        word='A'*r+'BC'
        value=direct_trace(word,cs)
        rhs=trace(powers[r+1]@wstar)
        assert rhs==(p-1)*(value+(-1)**r)
        # |value| <= 2 p^((r+3)/2)+1, tested with integers.
        assert max(abs(value)-1,0)**2<=4*p**(r+3)
        values[str(r)]=value
    separated_values={}
    for r in range(1,17):
        m=r+1
        k=r+3
        value=direct_trace('A'*r+'BAC',cs)
        rhs=trace(powers[m]@hoperator)-p*trace(powers[m])+(p-1)**2*(-1)**m
        assert rhs==(p-1)*value
        if k%2==0:
            residual=max(abs(value)-p**(k//2)-(p-1),0)
            assert residual**2<=4*p**(k+1)
        else:
            residual=max(abs(value)-2*p**((k+1)//2)-(p-1),0)
            assert residual**2<=p**k
        separated_values[str(r)]=value
    kernels=[None,s]
    for j in range(2,7):
        kernels.append(kernels[-1]@cs['A'])
    gap_checks=0
    for ell in range(1,7):
        for m in range(1,7):
            # Start at B, then ell-1 copies of A, then C, then m-1 copies of A.
            word='B'+'A'*(ell-1)+'C'+'A'*(m-1)
            value=direct_trace(word,cs)
            rhs=sum(int(diag['A'][y])*int(k2[x,y])*int(kernels[ell][x,y])*int(kernels[m][x,y])
                    for x in range(p) for y in range(p))
            assert (p-1)*value==rhs,(p,ell,m)
            gap_checks+=1
    permutation_checks=0
    for r in range(1,9):
        k=r+2
        for a,b,c in permutations('ABC'):
            value=direct_trace(a*r+b+c,cs)
            if k%2==0:
                residual=max(abs(value)-1-k*p**(k//2),0)
                assert residual**2<=4*p**(k+1)
            else:
                residual=max(abs(value)-1-2*p**((k+1)//2),0)
                assert residual**2<=k*k*p**k
            permutation_checks+=1
    literal_checks=0
    if p<=13:
        for r in range(1,4):
            word='A'*r+'BC'
            assert literal(p,chi,diag,word)==values[str(r)]
            literal_checks+=1
        for r in range(1,3):
            assert literal(p,chi,diag,'A'*r+'BAC')==separated_values[str(r)]
            literal_checks+=1
    return dict(p=p, weighted_anchor_matrix_checks=1, block_and_row_checks=7,
                homogeneous_entry_checks=(p-1)**2, positive_definite_matrix_certificates=3,
                bareiss_leading_principal_minors=certs, growing_family_checks=len(values),
                sym_square_norm_principal_minors=h_certificate,
                separated_growing_family_checks=len(separated_values),
                exceptional_quadratic_mode_checks=2,
                two_gap_identity_checks=gap_checks, label_permutation_checks=permutation_checks,
                literal_coordinate_sums=literal_checks, N_A_r_B_C=values,
                N_A_r_B_A_C=separated_values,
                twisted_legendre_values=h)


def main():
    cases=[]
    for p in (5,13,17,29,41,61):
        cases.append(field_case(p))
        print(json.dumps({'p':p,'status':'exact checks passed'}),flush=True)
    files=['research/parallel2-necklace-2026-09-04.md',
           'experiments/parallel2_necklace_2026_09_04.py',
           'experiments/parallel_necklace_2026_09_04.py',
           'research/parallel-necklace-2026-09-04.md',
           'research/localized-necklace-identities.md']
    result=dict(status='Uniform growing-family proof uses Katz TwLeg Mellin theorem; no general-word claim.',
                arithmetic='Exact Python integers and Bareiss positivity certificates.',
                primary_source='https://web.math.princeton.edu/~nmk/mellin186.pdf',
                primary_source_bytes=661647,
                primary_source_sha256='08e36a0fdeb8a1dc3d8c95734f44d7569fab23cd45fb1a5da8bea5576e3e1b73',
                source_sections=['Section 15, Theorems 15.1 and 15.3, printed pages 52-53',
                                 'Corollary 4.2, printed page 22'],
                cases=cases,input_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    path=ROOT/'results/parallel2_necklace_2026_09_04.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'written':str(path)}),flush=True)


if __name__=='__main__':
    main()
