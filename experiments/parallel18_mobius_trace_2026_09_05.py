#!/usr/bin/env python3
"""Exact Cartesian SL2 energies, escape, Weil traces and class averages.

No floating-point acceptance; no assertion of Paley cancellation.
"""
from collections import Counter,deque
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
from math import isqrt
from pathlib import Path
from random import Random
import json

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()


def mul(g,h,p):
    a,b,c,d=g;e,f,k,l=h
    return ((a*e+b*k)%p,(a*f+b*l)%p,(c*e+d*k)%p,(c*f+d*l)%p)


def inv(g,p):a,b,c,d=g;return d,-b%p,-c%p,a


def cell(s,t,p):return s%p,(s*t-1)%p,1,t%p


def sl2(p):
    out=[]
    for c in range(1,p):
        ci=pow(c,-1,p)
        out.extend((a,(a*d-1)*ci%p,c,d) for a in range(p) for d in range(p))
    out.extend((a,b,0,pow(a,-1,p)) for a in range(1,p) for b in range(p))
    assert len(out)==len(set(out))==p*(p*p-1)
    assert all((a*d-b*c)%p==1 for a,b,c,d in out)
    return out


class Cyclotomic:
    def __init__(self,p):
        assert p%4==1 and all(p%d for d in range(2,isqrt(p)+1))
        self.p=p;self.zero=(0,)*(p-1)
        self.chi=[0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
        self.monos=[self.reduce([int(i==t) for i in range(p)]) for t in range(p)]
        gauss=[0]*p
        for x in range(p):gauss[x*x%p]+=1
        self.gamma=self.reduce(gauss)
        assert self.gamma==self.reduce(self.chi)
        assert self.multiply(self.gamma,self.gamma)==self.constant(p)
        self.gs=[self.shift(self.gamma,t) for t in range(p)]
        self.qs={}
        for q,l in product(range(p),repeat=2):
            coeff=[0]*p
            for z in range(p):coeff[(q*z*z+l*z)%p]+=1
            value=self.reduce(coeff);self.qs[q,l]=value
            if q:
                predicted=self.scale(self.gs[-l*l*pow(4*q,-1,p)%p],self.chi[q])
            else:predicted=self.constant(p if l==0 else 0)
            assert value==predicted
            COUNTS['quadratic_gauss_identities']+=1
    def constant(self,a):return (a,)+(0,)*(self.p-2)
    def reduce(self,a):
        b=[0]*self.p
        for j,x in enumerate(a):b[j%self.p]+=x
        return tuple(x-b[-1] for x in b[:-1])
    def scale(self,a,k):return tuple(k*x for x in a)
    def add(self,a,b):return tuple(x+y for x,y in zip(a,b))
    def shift(self,a,k):
        b=[0]*self.p
        for j,x in enumerate(a):b[(j+k)%self.p]+=x
        return self.reduce(b)
    def multiply(self,a,b):
        out=[0]*(len(a)+len(b)-1)
        for i,x in enumerate(a):
            if x:
                for j,y in enumerate(b):out[i+j]+=x*y
        return self.reduce(out)
    def entry(self,g,x,y):
        """The exact cyclotomic integer p*rho(g)[x,y]."""
        a,b,c,d=g;p=self.p
        if c:
            exponent=(a*x*x-2*x*y+d*y*y)*pow(c,-1,p)%p
            return self.scale(self.gs[exponent],self.chi[c])
        if y==a*x%p:return self.scale(self.monos[a*b*x*x%p],p*self.chi[a])
        return self.zero


def representation_check(C):
    p=C.p;G=sl2(p);u=(1,1,0,1);w=(0,p-1,1,0)
    for g in G:
        a,b,c,d=g;trace=C.zero
        gu=mul(g,u,p);gw=mul(g,w,p)
        for x in range(p):
            trace=C.add(trace,C.entry(g,x,x))
            for y in range(p):
                assert C.shift(C.entry(g,x,y),y*y)==C.entry(gu,x,y)
                if c:
                    ci=pow(c,-1,p)
                    # Expanded rho(g)rho(w), after using Gamma^2=p.
                    actual=C.scale(C.shift(C.qs[d*ci%p,-2*(x*ci+y)%p],a*x*x*ci),C.chi[c])
                else:actual=C.scale(C.gs[(a*b*x*x-2*a*x*y)%p],C.chi[a])
                assert actual==C.entry(gw,x,y)
                COUNTS['right_generator_kernel_entries']+=2
        if (a+d)%p!=2:expected=C.constant(p*C.chi[(a+d-2)%p])
        elif c:expected=C.scale(C.gamma,p*C.chi[c])
        elif b:expected=C.scale(C.gamma,p*C.chi[b])
        else:expected=C.constant(p*p)
        assert trace==expected
        COUNTS['all_group_character_checks']+=1
    # The elementary generators really reach the whole finite group.
    seen={(1,0,0,1)};todo=deque(seen)
    while todo:
        g=todo.popleft()
        for h in (u,w):
            z=mul(g,h,p)
            if z not in seen:seen.add(z);todo.append(z)
    assert seen==set(G)
    return dict(p=p,group_order=len(G),complete_right_generator_checks=2*len(G)*p*p,
                complete_character_checks=len(G),generation_checked=True)


def additive_energy(A,p):
    return sum(x*x for x in Counter((b-a)%p for a in A for b in A).values())


def rectangular_case(C,A,B):
    p=C.p;m=len(A);n=len(B);U=[(a+2)%p for a in A];V=[-b%p for b in B]
    G=[cell(s,t,p) for s in U for t in V]
    assert len(set(G))==m*n
    products=Counter(mul(inv(g,p),h,p) for g in G for h in G)
    energy=sum(v*v for v in products.values())
    expected=n*n*additive_energy(A,p)+m*m*additive_energy(B,p)-m*m*n*n
    assert energy==expected<=m*m*n*n*(m+n-1)
    COUNTS['cartesian_energy_checks']+=1
    # All rational point transporters, including infinity (encoded p).
    for x in range(p+1):
        counts=Counter()
        for s,t in product(U,V):
            y=s if x==p else p if (x+t)%p==0 else (s-pow((x+t)%p,-1,p))%p
            counts[y]+=1
        assert max(counts.values())<=max(m,n)
        if x!=p:assert all(v<=min(m,n) for y,v in counts.items() if y!=p)
        COUNTS['rational_transporter_rows']+=1
    # All nonrational source points in F_p(sqrt(nu)); targets are counted exactly.
    nu=next(x for x in range(2,p) if C.chi[x]==-1)
    for x,y in product(range(p),range(1,p)):
        counts=Counter()
        for s,t in product(U,V):
            den=((x+t)**2-nu*y*y)%p;di=pow(den,-1,p)
            z=((s-(x+t)*di)%p,y*di%p)
            assert z[1]!=0
            counts[z]+=1
        assert max(counts.values())<=min(m,n)
        COUNTS['nonrational_transporter_rows']+=1
    bias=sum(C.chi[(a-b)%p] for a in A for b in B);overlap=len(set(A)&set(B))
    raw=[0]*p
    for a,b,x in product(A,B,range(p)):raw[(a-b)*x*x%p]+=1
    scaled_trace=C.multiply(C.gamma,raw)
    assert scaled_trace==C.add(C.constant(p*bias),C.scale(C.gamma,p*overlap))
    COUNTS['centered_trace_checks']+=1
    # Exact numerator of the factorized Hilbert-Schmidt identity.
    sA=sum(C.chi[(a-b)%p] for a in A for b in A)
    sB=sum(C.chi[(a-b)%p] for a in B for b in B)
    factors=[]
    for X,s in ((A,sA),(B,sB)):
        coeff=[0]*p
        for a,b,x in product(X,X,range(p)):coeff[(a-b)*x*x%p]+=1
        factor=C.reduce(coeff)
        assert factor==C.add(C.constant(p*len(X)),C.scale(C.gamma,s))
        factors.append(factor)
    lhs=C.multiply(*factors)
    rhs=C.add(C.constant(p*p*m*n+p*sA*sB),C.scale(C.gamma,p*(m*sB+n*sA)))
    assert lhs==rhs
    COUNTS['factorized_frobenius_checks']+=1
    return dict(p=p,A=A,B=B,m=m,n=n,bilinear_sum=bias,intersection=overlap,
                normalized_trace_rational=str(Fraction(bias,m*n)),
                normalized_trace_sqrt_p_coefficient=str(Fraction(overlap,m*n)),energy=energy)


def conjugate_average(C):
    p=C.p;assert p>=13
    G=sl2(p);g0=(3,p-1,1,0)
    klass={mul(mul(h,g0,p),inv(h,p),p) for h in G}
    eps=C.chi[5]
    assert len(klass)==p*(p+eps) and all((g[0]+g[3])%p==3 for g in klass)
    coefficients=[[[0]*p for _ in range(p)] for _ in range(p)]
    direct=[[[0]*p for _ in range(p)] for _ in range(p)]
    for a,b,c,d in klass:
        if c:
            ci=pow(c,-1,p);sign=C.chi[c]
            for x,y in product(range(p),repeat=2):
                coefficients[x][y][(a*x*x-2*x*y+d*y*y)*ci%p]+=sign
        else:
            for x in range(p):direct[x][a*x%p][a*b*x*x%p]+=p*C.chi[a]
    for x,y in product(range(p),repeat=2):
        actual=C.add(C.multiply(C.gamma,coefficients[x][y]),C.reduce(direct[x][y]))
        expected=C.constant(p*p*(int(x==y)+eps*int(x==-y%p)))
        assert actual==expected
        COUNTS['complete_class_average_entries']+=1
    d=(p+eps)//2
    assert Fraction(2,p+eps)*d==1
    return dict(p=p,quadratic_character_of_5=eps,class_size=len(klass),
                averaged_matrix='(I+chi(5)*reflection)/(p+chi(5))',
                active_parity_dimension=d,exact_trace=1,exact_operator_norm=str(Fraction(1,d)),
                exact_power_traces={j:str(Fraction(1,d**(j-1))) for j in range(1,5)},
                matrix_entries_checked=p*p,representation_entries_summed=len(klass)*p*p)


def main():
    cyclo={p:Cyclotomic(p) for p in (5,13,17,29)}
    reps=[representation_check(cyclo[p]) for p in (5,13)]
    print(json.dumps(dict(representation_checks=reps)),flush=True)
    Asets=[list(c) for n in range(1,6) for c in combinations(range(5),n)]
    for A,B in product(Asets,repeat=2):rectangular_case(cyclo[5],A,B)
    print(json.dumps(dict(p=5,complete_nonempty_subset_rectangles=len(Asets)**2)),flush=True)
    examples=[];rng=Random(18)
    families=[([0,1,2],[0,1,2]),([0,1,2],[4,5,6]),(list(range(13)),list(range(13)))]
    for _ in range(64):families.append((sorted(rng.sample(range(13),rng.randrange(1,8))),sorted(rng.sample(range(13),rng.randrange(1,8)))))
    for A,B in families:examples.append(rectangular_case(cyclo[13],A,B))
    averages=[]
    for p in (13,17,29):
        averages.append(conjugate_average(cyclo[p]))
        print(json.dumps(dict(class_average=averages[-1])),flush=True)
    inputs=['research/parallel18-mobius-trace-2026-09-05.md',
            'experiments/parallel18_mobius_trace_2026_09_05.py',
            'sources/thomas-weil-math0610644v3.html']
    out=dict(status='Exact Cartesian group identities and Weil trace obstruction passed; classical cancellation remains unproved.',
             arithmetic='Exact prime and quadratic-extension field arithmetic, integer group counts, rational norms and powers, and cyclotomic polynomial identities. No floating-point acceptance.',
             counts=dict(COUNTS),representations=reps,exhaustive_p5_rectangles=961,
             p13_rectangles=len(families),p13_examples=examples[:3],class_averages=averages,
             imported_expansion='Lyamkin primary live HTML Section1.5 Theorem11 applied after the proved coset bound. The asymptotic theorem and its unknown constants are not established by these finite tests; direct source archival requests returned403.',
             input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest=ROOT/'results/parallel18_mobius_trace_2026_09_05.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),counts=dict(COUNTS))),flush=True)


if __name__=='__main__':main()
