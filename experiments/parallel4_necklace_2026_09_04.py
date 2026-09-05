#!/usr/bin/env python3
"""Exact checks of two arbitrary blocks and the three-gap obstruction."""
from hashlib import sha256
from itertools import product, permutations
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, literal, trace, direct_trace
from parallel2_necklace_2026_09_04 import positive_principal_minors
from parallel3_conic_2026_09_04 import reduce_poly, cyclic_product
from localized_necklace_identities import paired_trace

ROOT=Path(__file__).resolve().parents[1]


def character_data(p,chi):
    n=p-1
    g=next(a for a in range(2,p) if len({pow(a,j,p) for j in range(n)})==n)
    logs={pow(g,j,p):j for j in range(n)}
    return n,g,logs


def i_polynomials(p,chi,logs,q):
    n=p-1
    out={}
    for u in range(1,p):
        row=[0]*n
        for b in range(p):
            z=(u+(1-u)*b)%p
            if z:row[q*logs[z]%n]+=chi[b]*chi[(1-b)%p]
        out[u]=row
    return out


def mellin_checks(p,chi,qops,kernels):
    n,g,logs=character_data(p,chi)
    checks=0
    for q in range(n):
        ips=i_polynomials(p,chi,logs,q)
        if q:
            assert reduce_poly(ips[1],n)==[-1]
            # The two lists are disjoint, and zero-monodromy cannot match
            # the two quadratic characters required by an h=2 invariant.
            assert not ({(n//2+q)%n,0}&{n//2,q})
            assert sorted([(n//2+q)%n,0])!=[n//2,n//2]
        else:
            for u in range(1,p):
                assert reduce_poly(ips[u],n)==[-1-chi[u]+int(u==1)]
        if q==n//2:
            for u in range(1,p):assert reduce_poly(ips[u],n)==[int(kernels[2][1,u])]
        jacobi=[0]*n
        for y in range(1,p):jacobi[q*logs[y]%n]+=chi[(1-y)%p]
        for h,Q in qops.items():
            lhs=[0]*n
            for y in range(1,p):lhs[q*logs[y]%n]+=int(Q[1,y])
            corr=[0]*n
            for u in range(1,p):
                coeff=chi[u]*int(kernels[h][1,u])
                corr=[a+coeff*b for a,b in zip(corr,ips[u])]
            rhs=cyclic_product(jacobi,corr,n)
            assert reduce_poly([a-b for a,b in zip(lhs,rhs)],n)==[0],(p,q,h)
            checks+=1
    return dict(generator=g,mellin_factorization_checks=checks,
                hypergeometric_character_list_checks=n-1,
                actual_u1_stalk_checks=n-1)


def raw_hypergeometric_checks(p,chi):
    n,g,logs=character_data(p,chi)
    order=p*n
    cases=[]
    for q in range(1,n):
        raw={u:[0]*order for u in range(1,p)}
        # Literally enumerate the raw (2,2) definition, with lower
        # characters conjugated. No Gauss-sum factorization is used here.
        for x1,x2,y1,y2 in product(range(1,p),repeat=4):
            u=x1*x2*pow(y1*y2%p,-1,p)%p
            exponent=(n*(x1+x2-y1-y2)+p*((n//2+q)*logs[x1]
                       +n//2*logs[y1]-q*logs[y2]))%order
            raw[u][exponent]+=1
        ips=i_polynomials(p,chi,logs,q)
        for u in range(1,p):
            difference=raw[u].copy()
            for j,c in enumerate(ips[u]):difference[p*j]-=p*c
            assert reduce_poly(difference,order)==[0],(p,q,u)
        cases.append(dict(character_exponent=q,product_fibers=n,
                          raw_boundary_u1=reduce_poly(raw[1],order)))
    return dict(p=p,raw_tuple_assignments=(n-1)*n**4,
                raw_hypergeometric_fiber_checks=(n-1)*n,cases=cases)


def word_orbits():
    representatives=('AABBC','AABCB','ABABC')
    def orbit(w,full=False):
        labelings=list(permutations('ABC')) if full else [('A','B','C'),('B','A','C')]
        out=set()
        for labels in labelings:
            u=w.translate(str.maketrans('ABC',''.join(labels)))
            for v in (u,u[::-1]):
                for j in range(5):out.add(v[j:]+v[:j])
        return out
    small={w:orbit(w) for w in representatives}
    full={w:orbit(w,True) for w in representatives}
    all_words={''.join(w) for w in product('ABC',repeat=5)
               if sorted(w.count(x) for x in 'ABC')==[1,2,2]}
    assert len(all_words)==90
    assert all(len(small[w])==10 and len(full[w])==30 for w in representatives)
    assert set.union(*full.values())==all_words
    assert sum(map(len,full.values()))==len(all_words)
    return small,full


def graph_minor_checks():
    certificates={
        'AABCB': ([[0,1],[3],[6]],[[2],[4],[5]]),
        'ABABC': ([[0,1],[3],[5]],[[2],[4],[6]])}
    for word,(left,right) in certificates.items():
        edges={tuple(sorted((i,(i+1)%5))) for i in range(5)}|{(5,6)}
        for i,L in enumerate(word):
            if L in 'AC':edges.add((i,5))
            if L in 'BC':edges.add((i,6))
        assert (0,1) in edges
        assert sorted(x for group in left+right for x in group)==list(range(7))
        for a in left:
            for b in right:
                assert any(tuple(sorted((x,y))) in edges for x in a for y in b)
    return certificates


def field_case(p,covered):
    depth=8 if p==5 else 6
    chi,diag,cs=matrices(p)
    S=np.array([[chi[(x-y)%p] for y in range(p)] for x in range(p)],dtype=object)
    powers={L:[np.eye(p,dtype=object)] for L in ('A','B')}
    for L in powers:
        for _ in range(depth+2):powers[L].append(powers[L][-1]@cs[L])
    t=[None]+[trace(powers['A'][j])//(p-1) for j in range(1,depth+3)]
    kernels=[None]+[S@powers['A'][j-1] for j in range(1,depth+3)]
    qops={h:np.zeros((p,p),dtype=object) for h in range(1,depth+1)}
    for b in range(p):
        db=np.array([chi[(x-b)%p] for x in range(p)],dtype=object)
        cb=db[:,None]*S
        end=diag['A'][:,None]*cb
        power=np.eye(p,dtype=object)
        for h in range(1,depth+1):
            power=power@cb
            qops[h]+=chi[b]*(power@end)
    certificates={};blocks={};normality=0;homogeneity=0
    for h,Q in qops.items():
        op=Q[1:,1:]
        assert np.all(Q@np.ones(p,dtype=object)==0)
        assert Q[0,0]==(p-1)*((-1)**h-t[h])
        assert np.all(Q[0,1:]==t[h]-(-1)**h)
        assert np.all(op@np.ones(p-1,dtype=object)==2*(-1)**h-t[h])
        assert np.array_equal(op@diag['A'][1:],-t[h+2]*diag['A'][1:])
        assert np.array_equal(op@op.T,op.T@op)
        normality+=1
        for x in range(1,p):
            for y in range(1,p):
                assert op[x-1,y-1]==Q[1,y*pow(x,-1,p)%p]
                homogeneity+=1
        B=(h+1)**2*p**(h+2)
        certificates[str(h)]=positive_principal_minors(B*np.eye(p-1,dtype=object)-op@op.T)
        right=powers['B'][h]@cs['C']
        for r in range(1,7):
            observed=paired_trace(powers['A'][r],right)
            numerator=paired_trace(powers['A'][r][1:,1:],op)
            correction=(-1)**(r-1)*(t[h]-(-1)**h)
            assert (p-1)*observed==numerator+(p-1)*correction,(p,r,h)
            k=r+h+1;s=min(r,h)
            assert max(abs(observed)-abs(t[s]-(-1)**s),0)**2<=(s+1)**2*p**(k+1)
            reverse=paired_trace(powers['A'][h],powers['B'][r]@cs['C'])
            assert observed==reverse
            blocks[f'{r},{h}']=observed
    mellin=mellin_checks(p,chi,qops,kernels) if p<=17 else None
    H=np.array([[sum(chi[v]*chi[(v-x)%p]*chi[(v-y)%p]*chi[(v-1)%p]
                       for v in range(p)) for y in range(p)] for x in range(p)],dtype=object)
    for x in range(1,p):
        for y in range(1,p):
            assert H[x,y]==chi[x*y%p]*kernels[2][(pow(x,-1,p)-1)%p,(pow(y,-1,p)-1)%p]-1
    for y in range(p):assert H[0,y]==p*int(y==1)-1-chi[y]
    triples=list(product(range(1,4),repeat=3)) if p<=17 else [(1,1,3),(3,1,1),(2,1,2),(2,2,2)]
    gaps={}
    for a,b,c in triples:
        word='B'+'A'*(a-1)+'B'+'A'*(b-1)+'C'+'A'*(c-1)
        actual=direct_trace(word,cs)
        full=sum(int(kernels[a][x,y])*int(kernels[b][1,y])*int(kernels[c][1,x])*int(H[x,y])
                 for x in range(p) for y in range(p))
        interior=sum(int(kernels[a][x,y])*int(kernels[b][1,y])*int(kernels[c][1,x])*int(H[x,y])
                     for x in range(1,p) for y in range(1,p))
        correction=(-1)**(a+c-2)*(p*t[b]-2*(-1)**b)+(-1)**(a+b-2)*(p*t[c]-2*(-1)**c)
        assert actual==full==interior+correction,(p,a,b,c)
        gaps[f'{a},{b},{c}']=dict(value=actual,zero_correction=correction)
    covered_values={w:direct_trace(w,cs) for w in covered}
    assert all(abs(v)<=5*p**3 for v in covered_values.values())
    assert abs(covered_values['AABBC'])<=3*p**3+2
    literals={}
    if p<=13:
        words=['AABBC','AABCB','ABABC']
        if p==5:words+=['AAABBC','AABBBC','AAABBBC']
        for w in words:
            actual=literal(p,chi,diag,w)
            assert actual==direct_trace(w,cs)
            literals[w]=actual
    witness=None
    if p==13:
        witness={'ratio':2,'H(2,4,1)':int(H[2,4]),'H(3,6,1)':int(H[3,6])}
        assert H[2,4]==-3 and H[3,6]==1
    return dict(p=p,maximum_rank=depth,ranks_at_least_characteristic=list(range(p,depth+1)),
                weighted_average_operators=depth,zero_and_exceptional_mode_checks=5*depth,
                normality_checks=normality,homogeneity_entry_checks=homogeneity,
                operator_norm_certificates=len(certificates),bareiss_leading_principal_minors=certificates,
                block_identity_bound_and_reversal_checks=len(blocks),block_values=blocks,
                mellin=mellin,three_gap_full_and_boundary_checks=len(gaps),three_gap_values=gaps,
                quartic_kernel_entry_checks=(p-1)**2,quartic_zero_row_checks=p,
                covered_length_five_checks=len(covered_values),covered_length_five_values=covered_values,
                literal_coordinate_cases=literals,nonconvolution_witness=witness)


def main():
    small,full=word_orbits()
    cases=[]
    for p in (5,13,17,29,41):
        cases.append(field_case(p,sorted(full['AABBC'])))
        print(json.dumps({'p':p,'status':'all exact matrix and necklace checks passed'}),flush=True)
    hyper=[]
    for p in (5,13):
        hyper.append(raw_hypergeometric_checks(p,matrices(p)[0]))
        print(json.dumps({'p':p,'status':'literal general-character hypergeometric checks passed'}),flush=True)
    files=['research/parallel4-necklace-2026-09-04.md','experiments/parallel4_necklace_2026_09_04.py',
           'experiments/parallel_necklace_2026_09_04.py','experiments/parallel2_necklace_2026_09_04.py',
           'experiments/parallel3_conic_2026_09_04.py','experiments/localized_necklace_identities.py',
           'research/parallel3-necklace-2026-09-04.md','sources/katz-g2-hypergeometric.pdf']
    result=dict(status='Growing block family proved; two of three (2,2,1) bracelet classes remain unresolved.',
                arithmetic='Python integers and cyclotomic reduction; exact Bareiss positivity certificates.',
                primary_source=dict(url='https://web.math.princeton.edu/~nmk/g2hyper62finalcorrected.pdf',
                    archive='sources/katz-g2-hypergeometric.pdf',locations='Section 2, printed pp. 3-5'),
                covered_words=sorted(full['AABBC']),remaining_words=sorted(full['AABCB']|full['ABABC']),
                three_bracelet_orbits={w:sorted(v) for w,v in small.items()},
                nonplanar_minor_certificates=graph_minor_checks(),cases=cases,
                raw_hypergeometric_cases=hyper,
                input_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    target=ROOT/'results/parallel4_necklace_2026_09_04.json'
    target.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'written':str(target)}),flush=True)


if __name__=='__main__':main()
