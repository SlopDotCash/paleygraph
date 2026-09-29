#!/usr/bin/env python3
"""Exact inverse-set relation removal and signed moment completion.

The verifier tests transfer identities, not the unproved SS-B* estimate.
"""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial
from pathlib import Path
from random import Random
import json
import parallel19_relation_free_2026_09_05 as prior

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()
INVERSION_SEEN=set()


def Q(k,h):return k+(2*h-2)*comb(comb(k+h-1,h),2)
def K(h):return max(2,4*h*(h-1),(h+1)*2**(2*h-1)+1)


@lru_cache(None)
def chi(p):return tuple(prior.chi_table(p))


def poly_mul(a,b,p=None):
    c=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):c[i+j]+=x*y
    if p is not None:c=[v%p for v in c]
    return c


def poly_value(a,z,p):
    out=0
    for v in reversed(a):out=(out*z+v)%p
    return out


def numerator(nu,p):
    T=list(nu);out=[0]*len(T)
    for c,n in nu.items():
        term=[1]
        for d in T:
            if d!=c:term=poly_mul(term,[d,-1],p)
        for i,v in enumerate(term):out[i]=(out[i]+n*v)%p
    while len(out)>1 and out[-1]==0:out.pop()
    assert any(out) and len(out)-1<=len(T)-2
    for c,n in nu.items():
        expected=n
        for d in T:
            if d!=c:expected=expected*(d-c)%p
        assert expected%p and poly_value(out,c,p)==expected%p
    return out


def field_rows(C,p,weights=None):
    if weights is None:weights=[1]*len(C)
    ch=chi(p)
    return [sum(w*ch[(x-c)%p] for c,w in zip(C,weights)) for x in range(p)]


def inversion_check(C,z,p):
    key=(p,tuple(C),z)
    D=[pow((c-z)%p,-1,p) for c in C]
    w=[chi(p)[(c-z)%p] for c in C]
    assert len(set(D))==len(C)
    if key in INVERSION_SEEN:return D,w
    INVERSION_SEEN.add(key)
    old=field_rows(C,p);new=field_rows(D,p,w);k=len(C)
    for x in range(p):
        if x==z:continue
        t=pow((x-z)%p,-1,p)
        assert new[t]==chi(p)[p-1]*chi(p)[(x-z)%p]*old[x]
        COUNTS['pointwise_inversion_rows']+=1
    assert new[0]==chi(p)[p-1]*k
    for r in range(1,5):
        oldmoment=sum(v**(2*r) for v in old)
        newmoment=sum(v**(2*r) for v in new)
        assert newmoment==oldmoment-old[z]**(2*r)+k**(2*r)
        assert oldmoment<=newmoment
        COUNTS['corrected_moment_identities']+=1
    COUNTS['inversion_cases']+=1
    return D,w


def pole_check(C,h,p):
    assert p>h and prior.prime(p)
    C=sorted(C);base=set(C)
    multisets=list(combinations_with_replacement(C,h))
    roots=set()
    for left,right in combinations(multisets,2):
        nu=Counter(left);nu.subtract(right);nu={c:n for c,n in nu.items() if n}
        assert sum(nu.values())==0 and all(0<abs(n)<=h for n in nu.values())
        coeff=numerator(nu,p)
        proots={z for z in range(p) if poly_value(coeff,z,p)==0}
        assert len(proots)<=len(coeff)-1<=2*h-2
        roots.update(proots-base)
        for z in range(p):
            if z in base:continue
            reciprocal=sum(n*pow((c-z)%p,-1,p) for c,n in nu.items())%p
            assert (z in proots)==(reciprocal==0)
            COUNTS['rational_polynomial_equivalences']+=1
        COUNTS['nonzero_relation_polynomials']+=1
    direct_bad=set()
    for z in range(p):
        if z in base:continue
        D,w=inversion_check(C,z,p)
        if not prior.is_bh(D,h,p):direct_bad.add(z)
        COUNTS['complete_pole_candidates']+=1
    assert roots==direct_bad
    assert len(base|roots)<=Q(len(C),h)
    if p>Q(len(C),h):assert len(base|roots)<p
    COUNTS['complete_pole_sets']+=1
    return dict(p=p,h=h,C=C,bad_poles=sorted(base|roots),bound=Q(len(C),h))


def complete(D,h,p):
    k=len(D);assert prior.is_bh(D,h,p) and p>prior.f_bound(2*k-1,h)
    S=set(D)
    while len(S)<2*k:
        F=prior.forbidden(S,h,p)
        x=next(x for x in range(p) if x not in F)
        S.add(x);assert prior.is_bh(S,h,p)
        COUNTS['completion_extensions']+=1
    E=sorted(S-set(D));assert len(E)==k
    for j in range(1,h+1):
        sums=Counter(sum(t)%p for t in product(S,repeat=j))
        assert sum(v*v for v in sums.values())==prior.minimum_energy(2*k,j)
        COUNTS['completion_minimal_energies']+=1
    return E


def signed_check(D,E,w,h,p,cache):
    k=len(D);P=[d for d,s in zip(D,w) if s==1];N=[d for d,s in zip(D,w) if s==-1]
    a=len(P);s=2*a-k
    plus=[tuple(sorted(P+list(J))) for J in combinations(E,k-a)]
    minus=[tuple(sorted(N+list(J))) for J in combinations(E,a)]
    q=comb(k,a);assert len(plus)==len(minus)==q
    counts_plus=Counter(x for C in plus for x in C)
    counts_minus=Counter(x for C in minus for x in C)
    weights=dict(zip(D,w))
    for x in set(D)|set(E):
        assert k*q*weights.get(x,0)==k*(counts_plus[x]-counts_minus[x])+s*q*int(x in E)
        COUNTS['signed_indicator_coordinates']+=1
    sets=set(plus+minus+[tuple(E)])
    for C in sets:
        if C not in cache:
            assert len(C)==k and prior.is_bh(C,h,p)
            rows=field_rows(C,p)
            cache[C]=[sum(v**(2*r) for v in rows) for r in range(1,5)]
            COUNTS['unsigned_completion_set_moments']+=4
    rows=field_rows(D,p,w)
    output=[]
    for r in range(1,5):
        value=sum(v**(2*r) for v in rows)
        bound=max(cache[C][r-1] for C in sets)
        assert k**(2*r)*value<=(2*k+abs(s))**(2*r)*bound
        assert value<=3**(2*r)*bound
        COUNTS['signed_moment_transfer_checks']+=1
        output.append(dict(r=r,signed_moment=value,maximum_unsigned_completion_moment=bound))
    COUNTS['signed_patterns']+=1
    return output


def transfer_example(C,h,p,all_patterns=False):
    k=len(C);assert p>h and p>Q(k,h) and p>prior.f_bound(2*k-1,h)
    z=next(z for z in range(p) if z not in C and prior.is_bh([pow((c-z)%p,-1,p) for c in C],h,p))
    D,w=inversion_check(C,z,p);E=complete(D,h,p)
    patterns=list(product((-1,1),repeat=k)) if all_patterns else [tuple(w)]
    cache={};actual=None
    for pattern in patterns:
        out=signed_check(D,E,pattern,h,p,cache)
        if list(pattern)==w:actual=out
    assert actual is not None
    old=field_rows(C,p)
    for row in actual:
        value=sum(v**(2*row['r']) for v in old)
        assert value<=row['signed_moment']<=3**(2*row['r'])*row['maximum_unsigned_completion_moment']
        row['original_unsigned_moment']=value
        COUNTS['complete_original_to_relation_free_transfers']+=1
    return dict(p=p,h=h,k=k,C=C,pole=z,inverted_support=D,signs=w,filler=E,
                bad_pole_bound=Q(k,h),extension_bound=prior.f_bound(2*k-1,h),moments=actual)


def polynomial_thresholds(h):
    d=2*h;K0=K(h)
    N=[Fraction(1)]
    for j in range(h):N=poly_mul(N,[j,1])
    N=[x/factorial(h) for x in N]
    N2=poly_mul(N,N)
    Qpoly=[Fraction(0)]*(d+1);Qpoly[1]=1
    for j in range(d+1):Qpoly[j]+=(h-1)*(N2[j]-(N[j] if j<len(N) else 0))
    Fpoly=[Fraction(-1),Fraction(2)]+[Fraction(0)]*(d-1)
    for j in range(1,h+1):
        term=[Fraction(1)]
        for _ in range(2*h-j):term=poly_mul(term,[-1,2])
        for t,v in enumerate(term):Fpoly[t]+=v
    leading=Fraction(h-1,factorial(h)**2)
    assert Qpoly[-1]==leading<=Fraction(1,4)
    certificate=[]
    for name,polynomial in [('bad_pole',Qpoly),('completion',Fpoly)]:
        diff=[-v for v in polynomial];diff[-1]+=1
        shifted=[sum(diff[j]*comb(j,i)*K0**(j-i) for j in range(i,d+1)) for i in range(d+1)]
        assert shifted[0]>0 and all(v>=0 for v in shifted)
        certificate.append(dict(name=name,nonnegative_shifted_coefficients=True,
                                strictly_positive_constant=True,degree=d))
        COUNTS['all_k_shifted_polynomial_certificates']+=1
    return dict(h=h,K=K0,critical_bad_pole_leading_coefficient=str(leading),certificates=certificate)


def main():
    rng=Random(20);families=[]
    for p,h,cap in ((5,2,3),(7,2,3),(13,2,3),(7,3,3),(13,3,3),(13,4,3)):
        count=0
        for k in range(1,cap+1):
            for C in combinations(range(p),k):pole_check(C,h,p);count+=1
        families.append(dict(p=p,h=h,all_nonempty_sets_through_size=cap,count=count))
    for p,h,k in ((17,2,4),(29,3,4),(41,3,4)):
        for _ in range(12):pole_check(sorted(rng.sample(range(p),k)),h,p)
    witness=pole_check([0,1,2],2,5)
    assert witness['bad_poles']==list(range(5))
    print(json.dumps(dict(pole_families=families,unrestricted_size_counterexample=witness)),flush=True)
    examples=[]
    for h,k in ((2,1),(2,2),(2,3),(2,4),(2,5),(3,2),(3,3),(4,2)):
        p=max(h,Q(k,h),prior.f_bound(2*k-1,h))+1
        while not prior.prime(p):p+=1
        examples.append(transfer_example(list(range(k)),h,p,all_patterns=True))
        for _ in range(3):examples.append(transfer_example(sorted(rng.sample(range(p),k)),h,p))
    thresholds=[polynomial_thresholds(h) for h in range(2,13)]
    ledger=[]
    for eps in (Fraction(1,2),Fraction(1,3),Fraction(1,10),Fraction(1,50)):
        r=3
        while Fraction(1,r+1)+Fraction(1,r)>=eps/2:r+=1
        h=(r+1)//2;theta=Fraction(1,r+1);beta=Fraction(1,r)
        assert 2*h<=r+1 and (2*h-1)*theta<1
        saving=(eps-theta-beta)/(2*r);delta=eps/(8*r)
        assert saving>eps/(4*r)>delta
        ledger.append(dict(epsilon=str(eps),r=r,h=h,beta=str(beta),theta=str(theta),
                           saving_before_constants=str(saving),eventual_saving=str(delta)))
        COUNTS['SS_B_star_exponent_checks']+=1
    for r in range(3,33):
        for h in range(2,(r+1)//2+1):
            assert 2*h<=r+1 and 2*h-1<r+1
            COUNTS['admissible_relation_order_pairs']+=1
    inputs=['research/parallel20-inversion-moments-2026-09-05.md',
            'experiments/parallel20_inversion_moments_2026_09_05.py',
            'experiments/parallel19_relation_free_2026_09_05.py',
            'research/parallel19-relation-free-reduction-2026-09-05.md',
            'research/parallel5-classical-2026-09-04.md']
    out=dict(status='Exact inversion, relation removal, signed completion and moment transfer passed. Uniform SS-B* and Paley remain unproved.',
             arithmetic='Integer finite-field, polynomial and moment calculations, rational threshold polynomials and exponent comparisons. No floating-point acceptance.',
             counts=dict(COUNTS),exhaustive_pole_families=families,
             all_poles_bad_without_size_condition=witness,transfer_examples=examples,
             all_k_threshold_certificates=thresholds,conditional_exponent_ledger=ledger,
             imported_new_proof_theorems=[],input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    path=ROOT/'results/parallel20_inversion_moments_2026_09_05.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(written=str(path),counts=dict(COUNTS))),flush=True)


if __name__=='__main__':main()
