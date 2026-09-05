#!/usr/bin/env python3
"""Exact tests of the full three-gap theorem, retaining every boundary."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, literal, all_traces
from localized_necklace_identities import paired_trace

ROOT=Path(__file__).resolve().parents[1]


def final_bound(p,k,d,value):
    """Check |value| <= d*p^((k+1)/2)+5*p^(k/2)+2 exactly."""
    left=abs(value)-2
    if k%2:
        residual=left-d*p**((k+1)//2)
        return residual<=0 or residual**2<=25*p**k
    residual=left-5*p**(k//2)
    return residual<=0 or residual**2<=d*d*p**(k+1)


def catalog():
    out={}
    for k in (5,6):
        classes=Counter();missing=[]
        for letters in product('ABC',repeat=k):
            counts=sorted(Counter(letters).values(),reverse=True)
            if len(counts)<=2:category='at_most_two_labels'
            elif counts[1:]==[1,1]:category='two_singleton_exceptions'
            elif counts[-2:]==[2,1]:category='three_gap_family'
            else:
                category='remaining'
                missing.append(''.join(letters))
            classes[category]+=1
        expected={'at_most_two_labels':93,'two_singleton_exceptions':60,'three_gap_family':90} if k==5 else {
            'at_most_two_labels':189,'two_singleton_exceptions':90,'three_gap_family':360,'remaining':90}
        assert dict(classes)==expected
        assert sum(classes.values())==3**k
        if k==6:assert all(sorted(Counter(w).values())==[2,2,2] for w in missing)
        out[str(k)]={'counts':dict(classes),'total':3**k,'remaining_words':missing}
    return out


def field_case(p):
    depth=6 if p==5 else 4
    chi,diags,cs=matrices(p)
    S=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    powers=[np.eye(p,dtype=object)]
    for _ in range(depth):powers.append(powers[-1]@cs['A'])
    K=[None]+[S@powers[j-1] for j in range(1,depth+1)]
    # Full-field rows are deliberately separate from their restrictions.
    kappa=[None]+[K[j][1,:].copy() for j in range(1,depth+1)]
    multiplicative=[None]+[kappa[j][1:].copy() for j in range(1,depth+1)]
    t=[None]+[int(kappa[j][1]) for j in range(1,depth+1)]
    for j in range(1,depth+1):
        assert kappa[j][0]==(-1)**(j-1)
        assert K[j][0,0]==0
        assert np.array_equal(K[j][0,:],(-1)**(j-1)*np.array(chi,dtype=object))
        assert int(sum(x*x for x in multiplicative[j]))*(p-1)==(p-3)*p**j+2
        assert int(sum(x*x for x in kappa[j]))<=p**j
        assert t[j]**2<=p**j
    H=np.array([[sum(chi[v]*chi[(v-x)%p]*chi[(v-y)%p]*chi[(v-1)%p]
                       for v in range(p)) for y in range(p)] for x in range(p)],dtype=object)
    qmaps={}
    for x in range(2,p):
        invx=pow(x,-1,p)
        q={z:(1-x*z)*pow(z*(1-x)%p,-1,p)%p for z in range(1,p)}
        assert q[1]==1 and q[invx]==0
        assert x*pow(x-1,-1,p)%p not in (0,1)
        assert len(set(q.values()))==p-1
        assert (1-x*0)%p!=0 and (1-x)%p!=0  # pole at zero is simple
        qmaps[x]=q
    W={}
    for a in range(1,depth+1):
        for b in range(1,depth+1):
            W[a,b]=sum(chi[(1-x)%p]*int(kappa[a][x])*int(kappa[b][x]) for x in range(1,p))
            assert W[a,b]**2<=p**(a+b)
    inner={};pointwise_checks=0
    for a in range(1,depth+1):
        for b in range(1,depth+1):
            for x in range(2,p):
                invx=pow(x,-1,p);q=qmaps[x]
                val=sum(chi[z]*int(kappa[a][z])*int(kappa[b][x*z%p])*int(kappa[2][q[z]])
                        for z in range(1,p) if z!=invx)
                no_one=sum(chi[z]*int(kappa[a][z])*int(kappa[b][x*z%p])*int(kappa[2][q[z]])
                           for z in range(1,p) if z not in (1,invx))
                # The actual included t=1 stalk contribution.
                assert val-no_one==-t[a]*int(kappa[b][x])
                full=sum(chi[z]*int(kappa[a][z])*int(kappa[b][x*z%p])*int(kappa[2][q[z]])
                         for z in range(1,p))
                # q=0 is evaluated with kappa_2(0)=-1, not a zero extension.
                assert full-val==-t[b]*int(kappa[a][x])
                direct=sum(int(K[a][x,y])*int(kappa[b][y])*(int(H[x,y])+1)
                           for y in range(2,p))
                assert direct==chi[(1-x)%p]*val,(p,a,b,x)
                assert val**2<=((3*a+1)*b)**2*p**(a+b),(p,a,b,x,val)
                inner[a,b,x]=val
                pointwise_checks+=1
    left=[None]+[cs['B']@powers[a-1] for a in range(1,depth+1)]
    right=[None]+[cs['C']@powers[c-1] for c in range(1,depth+1)]
    rows={};values={};full_quartic=0
    for a in range(1,depth+1):
        for b in range(1,depth+1):
            first=left[a]@left[b]
            for c in range(1,depth+1):
                k=a+b+c
                actual=paired_trace(first,right[c])
                U=int(kappa[c]@K[a]@kappa[b])
                Ug=int(kappa[c][1:]@K[a][1:,1:]@kappa[b][1:])
                assert U==Ug+2*(-1)**k,(p,a,b,c,U,Ug)
                assert U*U<=p**k
                M=sum(chi[(1-x)%p]*int(kappa[c][x])*inner[a,b,x] for x in range(2,p))
                rhs=M-U-t[c]*W[a,b]-t[b]*W[a,c]+p*((-1)**(a+c)*t[b]+(-1)**(a+b)*t[c])-2*(-1)**k
                assert actual==rhs,(p,a,b,c,actual,rhs)
                assert M*M<=((3*a+1)*b)**2*p**(k+1)
                assert (t[c]*W[a,b])**2<=p**k
                assert (t[b]*W[a,c])**2<=p**k
                assert (p*t[b])**2<=p**k and (p*t[c])**2<=p**k
                assert final_bound(p,k,(3*a+1)*min(b,c),actual)
                # Direct quartic sum includes zero, one, and coincident coordinates.
                literal_quartic=sum(int(K[a][x,y])*int(kappa[b][y])*int(kappa[c][x])*int(H[x,y])
                                    for x in range(p) for y in range(p))
                assert actual==literal_quartic
                full_quartic+=1
                values[a,b,c]=actual
                rows[f'{a},{b},{c}']={'necklace':actual,'M_open':M,'U_full':U,'U_nonzero':Ug,
                    'W_ab':W[a,b],'W_ac':W[a,c],
                    'restored_constant':-2*(-1)**k}
    for (a,b,c),value in values.items():assert value==values[a,c,b]
    fives={w:v for w,v in all_traces(cs,5).items() if len(w)==5}
    assert len(fives)==243
    assert all(abs(v)<=15*p**3 for v in fives.values())
    literals={}
    if p<=13:
        words=['AABCB','ABABC','AABBC']
        if p==5:words+=['BAABAC','BAABAAC','BAABAAAC']
        for word in words:
            v=literal(p,chi,diags,word)
            r=np.eye(p,dtype=object)
            for label in word:r=r@cs[label]
            assert v==sum(int(x) for x in r.diagonal())
            literals[word]=v
    return dict(p=p,maximum_rank=depth,ranks_at_least_characteristic=list(range(p,depth+1)),
                full_vs_multiplicative_kernel_checks=depth,parseval_checks=depth,
                mobius_maps_checked=p-2,weighted_pair_Cauchy_Schwarz_checks=depth**2,
                inner_factorization_and_bound_checks=pointwise_checks,
                included_t1_stalk_trace_checks=pointwise_checks,
                excluded_t_inverse_x_corrections=pointwise_checks,
                master_identity_and_bound_checks=len(values),full_U_zero_corrections=len(values),
                full_quartic_checks=full_quartic,reversal_checks=len(values),
                three_gap_values=rows,length_five_necklace_checks=len(fives),
                length_five_values=fives,literal_original_necklaces=literals)


def main():
    coverage=catalog()
    cases=[]
    for p in (5,13,17,29,41):
        cases.append(field_case(p))
        print(json.dumps({'p':p,'status':'all exact three-gap and length-five checks passed'}),flush=True)
    files=['research/parallel5-necklace-2026-09-04.md','experiments/parallel5_necklace_2026_09_04.py',
           'research/parallel4-necklace-2026-09-04.md','research/parallel3-necklace-2026-09-04.md',
           'research/parallel-necklace-2026-09-04.md','research/planar-necklace-reductions.md',
           'experiments/parallel_necklace_2026_09_04.py','experiments/localized_necklace_identities.py',
           'sources/katz-g2-hypergeometric.pdf']
    result=dict(status='Full three-gap family proved; all length-five individual words covered; no full spectral aggregate.',
                arithmetic='Python integers and exact radical inequalities; no floating-point acceptance.',
                zero_convention='kappa_j(0)=(-1)^(j-1) is the full row, distinct from a zero-extended multiplicative convolution.',
                primary_sources=[dict(url='https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf',
                    archive='sources/katz-g2-hypergeometric.pdf',locations='Section 2, printed pp. 3-5'),
                    dict(url='https://web.math.princeton.edu/~nmk/Katz-GKM.pdf',locations='Sections 2.3 and 3.6')],
                finite_word_coverage=coverage,cases=cases,
                input_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    output=ROOT/'results/parallel5_necklace_2026_09_04.json'
    output.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'written':str(output)}),flush=True)


if __name__=='__main__':main()
