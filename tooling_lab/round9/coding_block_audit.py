#!/usr/bin/env python3
"""Exact evaluation-kernel incidence audit; known algebra, new local audit artifact."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def evaluate(f,x,p):
    y=0
    for coefficient in reversed(f):y=(y*x+coefficient)%p
    return y


def rank(rows,p):
    a=[list(row) for row in rows];r=0
    for j in range(len(a[0]) if a else 0):
        pivot=next((i for i in range(r,len(a)) if a[i][j]%p),None)
        if pivot is None:continue
        a[r],a[pivot]=a[pivot],a[r];inv=pow(a[r][j],-1,p)
        a[r]=[(x*inv)%p for x in a[r]]
        for i in range(len(a)):
            if i!=r:
                v=a[i][j];a[i]=[(x-v*y)%p for x,y in zip(a[i],a[r])]
        r+=1
    return r


def audit(p,basis,blocks):
    d=len(basis)
    if rank(basis,p)!=d:raise ValueError('dependent polynomial basis')
    dims=[d-rank([[evaluate(f,x,p) for f in basis] for x in block],p) for block in blocks]
    multiplicity=Counter(x for block in blocks for x in block)
    return {'dimension':d,'kernel_dimensions':dims,'intersection_dimension_sum':sum(dims),
            'normalized_dimension_loss':[F(sum(dims),d*len(blocks)).numerator,F(sum(dims),d*len(blocks)).denominator],
            'maximum_coordinate_reuse':max(multiplicity.values()),'coordinate_multiplicities':sorted(multiplicity.items())}


def roots_polynomial(roots,p):
    f=[1]
    for root in roots:
        g=[0]*(len(f)+1)
        for j,x in enumerate(f):g[j]=(g[j]-root*x)%p;g[j+1]=(g[j+1]+x)%p
        f=g
    return f


def main():
    p=17;gamma=3;k=4;m=2
    assert len({pow(gamma,j,p) for j in range(p-1)})==p-1
    full=[[(a*pow(gamma,j,p))%p for j in range(m)] for a in range(1,p)]
    disjoint=[[pow(gamma,m*i+j,p) for j in range(m)] for i in range((p-1)//m)]
    f=roots_polynomial([1,3,9],p);other=roots_polynomial([1,2,4],p)
    all_a=audit(p,[f],full);dis=audit(p,[f],disjoint);alternate=audit(p,[other],full)
    bound=F(k-1,m)
    assert all_a['intersection_dimension_sum']==2>bound and dis['intersection_dimension_sum']==1<=bound
    assert alternate['intersection_dimension_sum']==0
    histogram=Counter();count=0;max_all=max_dis=0
    for degree in range(k):
        for coefficients in product(range(p),repeat=degree):
            polynomial=[*coefficients,1]
            # Independent d=1 evaluation by powers, no rank or Horner routine.
            roots={x for x in range(1,p) if sum(c*pow(x,j,p) for j,c in enumerate(polynomial))%p==0}
            a=sum(set(block)<=roots for block in full);b=sum(set(block)<=roots for block in disjoint)
            assert a<=max(0,len(roots)-m+1) and b<=len(roots)//m
            histogram[a,b]+=1;count+=1;max_all=max(max_all,a);max_dis=max(max_dis,b)
    assert max_all==2 and max_dis==1
    singleton=audit(p,[[1,0,0,0],[0,1,0,0]],[[x] for x in range(1,p)])
    assert singleton['intersection_dimension_sum']==16 and F(*singleton['normalized_dimension_loss'])==F(1,2)>F(k,16)
    out={'status':'passed','scope':'Exact finite block-incidence and scalar-alphabet transfer audit. No claim that the cited paper main theorem is false.',
         'source_under_audit':'https://arxiv.org/pdf/2601.10047v1#page=10',
         'parameters':{'p':p,'gamma':gamma,'k':k,'m':m,'d':1},'polynomial_coefficients':f,
         'roots':[1,3,9],'all_nonzero_blocks':full,'disjoint_blocks':disjoint,
         'all_basepoints_audit':all_a,'disjoint_audit':dis,'displayed_lemma_4_1_rhs':[bound.numerator,bound.denominator],
         'same_root_count_control':{'polynomial_coefficients':other,'roots':[1,2,4],'audit':alternate},
         'complete_monic_polynomials_tested':count,'complete_count_histogram':[[a,b,n] for (a,b),n in sorted(histogram.items())],
         'max_all_basepoints':max_all,'max_disjoint_basepoints':max_dis,
         'ordinary_RS_two_dimensional_audit':singleton,
         'source_sha256':{'coding_block_audit.py':sha256(Path(__file__).read_bytes()).hexdigest()}}
    (HERE/'coding_block_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','monic_polynomials':count,'all_basepoint_max':max_all,'disjoint_max':max_dis,
                      'counterexample_polynomial':f,'ordinary_RS_dimension_loss':singleton['normalized_dimension_loss']}),flush=True)


if __name__=='__main__':main()
