#!/usr/bin/env python3
"""Exact checks of arbitrary two-gap necklaces and hypergeometric trace identification."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, literal, trace
from parallel2_necklace_2026_09_04 import positive_principal_minors
from localized_necklace_identities import paired_trace

ROOT=Path(__file__).resolve().parents[1]


def cyclic_product(a,b):
    p=len(a)
    result=[0]*p
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            result[(i+j)%p]+=x*y
    return result


def convolution_kernels(p,chi,depth):
    # Multiplicative product-fiber convolution, independent of matrix powers.
    f=[chi[(1-t)%p] for t in range(p)]
    values=[None,f]
    for _ in range(2,depth+1):
        previous=values[-1]
        values.append([0]+[sum(f[u]*previous[t*pow(u,-1,p)%p] for u in range(1,p))
                          for t in range(1,p)])
    return values


def literal_hypergeometric(p,chi,j,kj):
    # Enumerate all nonzero x/y tuples, grouping only by their product and sum.
    tuples=Counter()
    for xs in product(range(1,p),repeat=j):
        multiplicative=1
        for x in xs:
            multiplicative=multiplicative*x%p
        tuples[(multiplicative,sum(xs)%p)]+=1
    raw={t:[0]*p for t in range(1,p)}
    for (px,sx),nx in tuples.items():
        for (py,sy),ny in tuples.items():
            t=px*pow(py,-1,p)%p
            raw[t][(sx-sy)%p]+=nx*ny*chi[py]
    gauss=[1]+[0]*(p-1)
    for _ in range(j):
        gauss=cyclic_product(gauss,chi)
    for t in range(1,p):
        difference=[raw[t][a]-kj[t]*gauss[a] for a in range(p)]
        # A polynomial of degree <p vanishes at zeta_p exactly when all
        # its coefficients are equal (a multiple of 1+X+...+X^(p-1)).
        assert len(set(difference))==1,(p,j,t)
    return dict(j=j,product_fiber_checks=p-1,boundary_t1_coefficients=raw[1],
                gauss_power_coefficients=gauss,k_j_at_1=kj[1],
                original_tuple_pairs=(p-1)**(2*j))


def field_case(p):
    depth=10 if p<=17 else 8
    chi,diag,cs=matrices(p)
    s=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    powers=[np.eye(p,dtype=object)]
    for _ in range(depth):
        powers.append(powers[-1]@cs['A'])
    kernels=[None]+[s@powers[j-1] for j in range(1,depth+1)]
    conv=convolution_kernels(p,chi,depth)
    for j in range(1,depth+1):
        for t in range(1,p):
            assert kernels[j][1,t]==conv[j][t]
    v=diag['A'][1:]
    operators={}
    certificates={}
    for j in range(1,depth+1):
        op=v[:,None]*v[None,:]*kernels[2][1:,1:]*kernels[j][1:,1:]
        operators[j]=op
        for x in range(1,p):
            for y in range(1,p):
                t=y*pow(x,-1,p)%p
                assert op[x-1,y-1]==chi[t]*conv[2][t]*conv[j][t]
        if j!=2:
            bound=(j+1)**2*p**(j+1)
            certificates[str(j)]=positive_principal_minors(bound*np.eye(p-1,dtype=object)-op@op)
    corrected=operators[2]+p*np.eye(p-1,dtype=object)-p*np.outer(v,v)
    assert np.array_equal(operators[2]@v,(p*p-2*p-2)*v)
    assert np.array_equal(corrected@v,-2*v)
    certificates['2_corrected']=positive_principal_minors(4*p**3*np.eye(p-1,dtype=object)-corrected@corrected)
    arc_depth=8
    left=[None]+[cs['B']@powers[j-1] for j in range(1,arc_depth+1)]
    right=[None]+[cs['C']@powers[m-1] for m in range(1,arc_depth+1)]
    gap_values={}
    average_checks=0
    for j in range(1,arc_depth+1):
        for m in range(1,arc_depth+1):
            value=paired_trace(left[j],right[m])
            rhs=trace(powers[m][1:,1:]@operators[j])
            k=j+m
            assert rhs==(p-1)*(value-(-1)**(k-1)),(p,j,m,value,rhs)
            bound_coefficient=min(j,m)+1
            assert max(abs(value)-1,0)**2<=bound_coefficient**2*p**(k+1)
            if p<=17:
                # Independently retain every pair including zero coordinates.
                average=sum(int(diag['A'][y])*int(kernels[2][x,y])*int(kernels[j][x,y])*int(kernels[m][x,y])
                            for x in range(p) for y in range(p))
                assert average==(p-1)*value
                average_checks+=1
            gap_values[f'{j},{m}']=value
    for j in range(1,arc_depth+1):
        for m in range(1,arc_depth+1):
            assert gap_values[f'{j},{m}']==gap_values[f'{m},{j}']
    literal_cases=[]
    if p<=13:
        pairs=[(1,1),(1,2),(2,2),(2,3)]
        if p==5:
            pairs+=[(3,3),(3,4)]
        for j,m in pairs:
            word='B'+'A'*(j-1)+'C'+'A'*(m-1)
            actual=literal(p,chi,diag,word)
            assert actual==gap_values[f'{j},{m}']
            literal_cases.append(dict(j=j,m=m,value=actual))
    hypergeom=[]
    if p<=13:
        for j in (1,2,3):
            hypergeom.append(literal_hypergeometric(p,chi,j,conv[j]))
    return dict(p=p,maximum_hypergeometric_rank=depth,
                convolution_entry_checks=depth*(p-1),homogeneous_operator_entry_checks=depth*(p-1)**2,
                gap_identity_and_bound_checks=len(gap_values),full_zero_retaining_average_checks=average_checks,
                reversal_checks=arc_depth**2,operator_norm_certificates=len(certificates),
                ranks_at_least_characteristic=[j for j in range(1,depth+1) if j>=p],
                bareiss_leading_principal_minors=certificates,
                gap_values=gap_values,literal_coordinate_cases=literal_cases,
                literal_hypergeometric_cases=hypergeom)


def main():
    cases=[]
    for p in (5,13,17,29,41):
        cases.append(field_case(p))
        print(json.dumps({'p':p,'status':'all exact checks passed'}),flush=True)
    files=['research/parallel3-necklace-2026-09-04.md',
           'experiments/parallel3_necklace_2026_09_04.py',
           'experiments/parallel2_necklace_2026_09_04.py',
           'experiments/parallel_necklace_2026_09_04.py',
           'research/parallel2-necklace-2026-09-04.md',
           'research/parallel-necklace-2026-09-04.md',
           'sources/katz-g2-hypergeometric.pdf']
    sources=[dict(url='https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf',
                  archive='sources/katz-g2-hypergeometric.pdf',bytes=479161,
                  sha256='0bf485bcc9dde2afebc8268af679566bcff65aa6d6ed8c9f4b0387d7e9e954c1',
                  locations='Section 2, printed pp. 3-5; zero-based PDF pp. 2-4'),
             dict(url='https://web.math.princeton.edu/~nmk/Katz-GKM.pdf',
                  bytes=6307344,
                  sha256='8711424f8edb14f38c0e61606ef49d6efc67d5a93ac9f5bd9fa3e13b7c7d9ffe',
                  locations='Sections 2.3.1-2.3.2 and 3.6; zero-based PDF pp. 20 and 24')]
    result=dict(status='General two-gap family proved using imported hypergeometric and weight theorems.',
                arithmetic='Python integers, cyclotomic coefficient identities, and Bareiss positivity.',
                primary_sources=sources,cases=cases,
                input_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    target=ROOT/'results/parallel3_necklace_2026_09_04.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'written':str(target)}),flush=True)


if __name__=='__main__':
    main()
