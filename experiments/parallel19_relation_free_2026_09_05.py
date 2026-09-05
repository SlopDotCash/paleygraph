#!/usr/bin/env python3
"""Exact B_h extension, partition and restricted-character transfer checks.

No assertion of the unproved restricted moment estimate or Paley cancellation.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial, isqrt
from pathlib import Path
from random import Random
import json

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def matrix_product(g,h,p):
    a,b,c,d=g;e,f,k,l=h
    return ((a*e+b*k)%p,(a*f+b*l)%p,(c*e+d*k)%p,(c*f+d*l)%p)


def multiset_sums(S,h,p):
    return [sum(c)%p for c in combinations_with_replacement(sorted(S),h)]


def is_bh(S,h,p):
    sums=multiset_sums(S,h,p)
    return len(sums)==len(set(sums))


def forbidden(S,h,p):
    assert p>h and prime(p) and is_bh(S,h,p)
    if not S:return set()
    sums=[set(multiset_sums(S,j,p)) for j in range(h+1)]
    out=set(S)
    for j in range(1,h+1):
        ji=pow(j,-1,p)
        out.update((v-u)*ji%p for u in sums[h-j] for v in sums[h])
    return out


def f_bound(s,h):return s+sum(s**(2*h-j) for j in range(1,h+1))


def L(k,h):return f_bound(k-1,h)


def extension_check(S,h,p):
    F=forbidden(S,h,p)
    expected={x for x in range(p) if x in S or not is_bh(set(S)|{x},h,p)}
    assert F==expected
    assert len(F)<=f_bound(len(S),h)
    COUNTS['complete_forbidden_sets']+=1
    COUNTS['independent_candidate_extensions']+=p


def minimum_energy(k,j):
    # (j!)^2 [z^j] (sum_(a=0)^j z^a/(a!)^2)^k, exactly.
    coeff=[Fraction(1)]+[Fraction(0)]*j
    for _ in range(k):
        coeff=[sum((coeff[t-a]/factorial(a)**2 for a in range(t+1)),Fraction(0)) for t in range(j+1)]
    result=coeff[j]*factorial(j)**2
    assert result.denominator==1
    return result.numerator


def energy_check(S,h,p):
    for j in range(1,h+1):
        counts=Counter(sum(t)%p for t in product(S,repeat=j))
        actual=sum(v*v for v in counts.values())
        assert actual==minimum_energy(len(S),j)
        COUNTS['minimal_additive_energies']+=1


def block(D,k,h,p):
    assert len(D)>L(k,h)
    S=set()
    while len(S)<k:
        F=forbidden(S,h,p)
        available=set(D)-F
        assert available
        S.add(min(available))
        assert is_bh(S,h,p)
        COUNTS['greedy_extensions']+=1
    energy_check(S,h,p)
    return sorted(S)


def partition(B,k,h,p):
    assert len(set(B))==len(B) and k>=1
    remaining=set(B);blocks=[]
    while len(remaining)>L(k,h):
        C=block(remaining,k,h,p)
        blocks.append(C);remaining.difference_update(C)
    R=sorted(remaining)
    union=[x for C in blocks for x in C]+R
    assert len(union)==len(set(union))==len(B) and set(union)==set(B)
    assert len(R)<=L(k,h)<=(h+1)*k**(2*h-1)
    COUNTS['partitions']+=1
    return blocks,R


def chi_table(p):return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]


def moment(C,r,p,chi):
    return sum(sum(chi[(x-c)%p] for c in C)**(2*r) for x in range(p))


def transfer_check(A,B,k,h,p):
    chi=chi_table(p);m=len(A);n=len(B)
    ac,ar=partition(A,k,h,p);bc,br=partition(B,k,h,p)
    total=sum(chi[(a-b)%p] for a in A for b in B)
    block_sums=[sum(chi[(a-b)%p] for a in C for b in D) for C in ac for D in bc]
    signed_complete=sum(block_sums)
    remainder_pairs=m*len(br)+n*len(ar)-len(ar)*len(br)
    ars=set(ar);brs=set(br)
    exact_remainder=sum(chi[(a-b)%p] for a in A for b in B if a in ars or b in brs)
    assert total==signed_complete+exact_remainder
    assert abs(exact_remainder)<=remainder_pairs
    beta=Fraction(max(map(abs,block_sums),default=0),k*k)
    bound=beta*(m-len(ar))*(n-len(br))+remainder_pairs
    assert abs(total)<=bound
    assert bound<=m*n*(beta+Fraction(L(k,h),m)+Fraction(L(k,h),n))
    COUNTS['exact_two_sided_transfers']+=1
    for r in (1,2,3,4):
        U=max((moment(C,r,p,chi) for C in bc),default=0)
        lhs=max(0,abs(total)-m*len(br))**(2*r)
        rhs=len(bc)**(2*r)*m**(2*r-1)*U
        assert lhs<=rhs
        COUNTS['integer_moment_transfers']+=1
    return dict(p=p,h=h,k=k,m=m,n=n,blocks_A=len(ac),blocks_B=len(bc),
                remainder_A=len(ar),remainder_B=len(br),remainder_bound=L(k,h),
                sum=total,maximum_normalized_block_bias=str(beta))


def arbitrary_kernel_check(p,h,k,rng):
    A=list(range(p));B=list(range(p))
    ac,ar=partition(A,k,h,p);bc,br=partition(B,k,h,p)
    # An arbitrary signed integer kernel, bounded by one; unrelated to chi.
    K=[[rng.randrange(-1,2) for _ in range(p)] for _ in range(p)]
    total=sum(map(sum,K))
    maxblock=max((abs(sum(K[a][b] for a in C for b in D)) for C in ac for D in bc),default=0)
    remainder=p*len(ar)+p*len(br)-len(ar)*len(br)
    assert abs(total)<=len(ac)*len(bc)*maxblock+remainder
    COUNTS['arbitrary_bounded_kernel_transfers']+=1


def exponent_check():
    rows=[]
    for eps in (Fraction(1,2),Fraction(1,3),Fraction(1,10),Fraction(1,50)):
        r=4
        while True:
            h=isqrt(r);theta=Fraction(1,r+1);beta=Fraction(1,r)
            if theta+beta<eps/2 and (2*h-1)*theta<eps/2:break
            r+=1
        first=(eps-theta-beta)/(2*r)
        second=eps-(2*h-1)*theta
        delta=eps/(8*r)
        assert first>eps/(4*r)>delta and second>eps/2>delta
        rows.append(dict(epsilon=str(eps),r=r,h=h,beta=str(beta),theta=str(theta),
                         moment_saving=str(first),remainder_saving=str(second),
                         eventual_saving=str(delta)))
        COUNTS['sufficient_sequence_exponent_checks']+=1
    for h in range(2,9):
        for r in range(2,31):
            theta=Fraction(1,r+1)
            assert max(theta,(2*h-1)*theta)==(2*h-1)*theta
            COUNTS['fixed_order_threshold_checks']+=1
    assert Fraction(2*2-1,2+1)==1
    return rows


def main():
    rng=Random(19)
    exhaustive=[]
    for p,h,cap in ((5,2,5),(7,2,7),(13,2,5),(17,2,5),(7,3,4),(13,3,4),(17,3,4),(17,4,3)):
        count=0
        for size in range(cap+1):
            for S in combinations(range(p),size):
                if is_bh(S,h,p):
                    extension_check(S,h,p);count+=1
        exhaustive.append(dict(p=p,h=h,all_Bh_sets_through_size=cap,count=count))
    print(json.dumps(dict(extension_checks=exhaustive)),flush=True)
    examples=[]
    for p,h,k in ((13,2,2),(17,3,2),(257,2,3),(257,2,4),(257,3,3),(1093,2,5),(1093,3,4),(1093,4,3)):
        assert prime(p)
        examples.append(transfer_check(list(range(p)),list(range(p)),k,h,p))
        for _ in range(8):
            amin=min(p,max(1,L(k,h)+1))
            m=rng.randrange(amin,p+1);n=rng.randrange(amin,p+1)
            A=sorted(rng.sample(range(p),m));B=sorted(rng.sample(range(p),n))
            examples.append(transfer_check(A,B,k,h,p))
    # Empty block families and the exact k=1 endpoint must also work.
    examples.append(transfer_check([0,1],[2,3],4,2,13))
    examples.append(transfer_check([0,1,4],[2,3,7],1,2,13))
    for _ in range(12):arbitrary_kernel_check(17,2,2,rng)
    # Interval obstruction: enumerate every B_h subset of these small intervals.
    for p,N,h in ((17,5,2),(29,6,3),(41,7,4)):
        assert h*(N-1)<p
        for size in range(N+1):
            for C in combinations(range(N),size):
                if not is_bh(C,h,p):continue
                assert comb(size+h-1,h)<=h*(N-1)+1
                assert size**h<=factorial(h)*(h*(N-1)+1)
                COUNTS['interval_size_obstruction_checks']+=1
    exponents=exponent_check()
    # Mixed Weil quotient identity, including both exceptional cases.
    for p in (5,13,17):
        chi=chi_table(p)
        for s,t,u,v in product(range(p),repeat=4):
            d=(u-s)%p;e=(v-t)%p
            a=(1+t*d)%p;b=(v-t+t*v*d)%p;c=-d%p;dd=(1-d*v)%p
            inverse_g=(t,(1-s*t)%p,p-1,s)
            second_g=(u,(u*v-1)%p,1,v)
            assert matrix_product(inverse_g,second_g,p)==(a,b,c,dd)
            assert (a*dd-b*c)%p==1 and (a+dd-2)%p==-d*e%p
            if d and e:assert chi[(a+dd-2)%p]==chi[d]*chi[e]
            elif d:assert c and chi[c]==chi[d]
            elif e:assert (a,b,c,dd)==(1,e,0,1)
            else:assert (a,b,c,dd)==(1,0,0,1)
            COUNTS['mixed_quotient_case_checks']+=1
    inputs=['research/parallel19-relation-free-reduction-2026-09-05.md',
            'experiments/parallel19_relation_free_2026_09_05.py',
            'research/parallel18-mobius-trace-2026-09-05.md',
            'research/parallel5-classical-2026-09-04.md']
    out=dict(status='Exact extension, partition and conditional transfer identities passed; restricted cancellation and SS-B remain unproved.',
             arithmetic='Integer prime-field and multiset counts, integer moment inequalities, and rational exponent comparisons; no floating-point acceptance.',
             counts=dict(COUNTS),exhaustive_extensions=exhaustive,transfer_examples=examples,
             exponent_examples=exponents,imported_proof_theorems=[],
             input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    path=ROOT/'results/parallel19_relation_free_2026_09_05.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(written=str(path),counts=dict(COUNTS))),flush=True)


if __name__=='__main__':main()
