#!/usr/bin/env python3
"""Exact GF(p^e) adapter and independently checked tiny agreement censuses.

Polynomial-basis integers encode coefficients in base p. Irreducibility is
validated by trial division by every monic polynomial of degrees <= e/2.
This bounded validation is deliberately simple, not a scalable field factory.
"""
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path


def remainder(a, b, p):
    a = a[:]
    while len(a) >= len(b):
        c = a[-1] * pow(b[-1], -1, p) % p
        offset = len(a) - len(b)
        for j, bj in enumerate(b):
            a[offset+j] = (a[offset+j] - c*bj) % p
        while a and a[-1] == 0:
            a.pop()
    return a


def irreducible(f, p):
    e = len(f)-1
    if f[-1] != 1:
        return False
    for d in range(1, e//2+1):
        for lower in product(range(p), repeat=d):
            if not remainder(f, list(lower)+[1], p):
                return False
    return True


class Field:
    def __init__(self, p, modulus):
        assert p >= 2 and all(p%d for d in range(2, int(p**0.5)+1))
        assert irreducible(modulus, p)
        self.p, self.f, self.e = p, modulus, len(modulus)-1
        self.q = p**self.e
        self.coeffs = [self.decode(i) for i in range(self.q)]

    def decode(self, x):
        return [(x//self.p**j)%self.p for j in range(self.e)]

    def encode(self, xs):
        return sum((x%self.p)*self.p**j for j,x in enumerate(xs))

    def add(self, a, b):
        return self.encode([x+y for x,y in zip(self.coeffs[a],self.coeffs[b])])

    def neg(self, a):
        return self.encode([-x for x in self.coeffs[a]])

    def sub(self, a, b):
        return self.add(a,self.neg(b))

    def mul(self, a, b):
        c=[0]*(2*self.e-1)
        for i,x in enumerate(self.coeffs[a]):
            for j,y in enumerate(self.coeffs[b]):
                c[i+j]=(c[i+j]+x*y)%self.p
        for i in range(len(c)-1,self.e-1,-1):
            z=c[i]
            for j in range(self.e):
                c[i-self.e+j]=(c[i-self.e+j]-z*self.f[j])%self.p
        return self.encode(c[:self.e])

    def power(self, a, n):
        out=1
        while n:
            if n&1:
                out=self.mul(out,a)
            a=self.mul(a,a)
            n//=2
        return out

    def inv(self, a):
        assert a
        return self.power(a,self.q-2)

    def sum(self, xs):
        out=0
        for x in xs:
            out=self.add(out,x)
        return out

    def dot(self, a, b):
        return self.sum(self.mul(x,y) for x,y in zip(a,b))


def subgroup(F,n):
    assert (F.q-1)%n==0
    for a in range(1,F.q):
        if F.power(a,n)==1 and all(F.power(a,d)!=1 for d in range(1,n)):
            return [F.power(a,i) for i in range(n)]
    raise ValueError("No root")


def checks(F):
    for a in range(1,F.q):
        assert F.mul(a,F.inv(a))==1
    # Exhaustive field identities at q=9, deterministic probes at q=729.
    triples=product(range(F.q),repeat=3) if F.q<=9 else ((i,7*i%F.q,17*i%F.q) for i in range(F.q))
    number=0
    for a,b,c in triples:
        assert F.mul(a,F.add(b,c))==F.add(F.mul(a,b),F.mul(a,c))
        assert F.mul(F.mul(a,b),c)==F.mul(a,F.mul(b,c))
        number+=1
    return {"irreducibility":"all monic factors of degree at most e/2 excluded",
            "nonzero_inverses_checked":F.q-1,"associativity_distributivity_triples":number}


def parity(F, points, k):
    weights=[]
    for i,x in enumerate(points):
        v=1
        for j,y in enumerate(points):
            if i!=j:
                v=F.mul(v,F.sub(x,y))
        weights.append(F.inv(v))
    return [[F.mul(F.power(x,ell),w) for x,w in zip(points,weights)]
            for ell in range(len(points)-k)]


def census(F,dom,k,s,u0,u1):
    candidates,common={},[]
    for S in combinations(range(len(dom)),s):
        H=parity(F,[dom[j] for j in S],k)
        a0=[F.dot(h,[u0[j] for j in S]) for h in H]
        a1=[F.dot(h,[u1[j] for j in S]) for h in H]
        idx=next((i for i,x in enumerate(a1) if x),None)
        if idx is None:
            if not any(a0):
                common.append(list(S))
            continue
        z=F.mul(F.neg(a0[idx]),F.inv(a1[idx]))
        if all(F.add(a,F.mul(z,b))==0 for a,b in zip(a0,a1)):
            candidates.setdefault(z,[]).append(list(S))
    return {"bad_scalars":list(range(F.q)) if common else sorted(candidates),
            "common_subsets":common,"subset_provenance":candidates}


def brute(F,dom,k,s,u0,u1):
    # Full codeword enumeration for k=2; for k=1 use its equivalent exhaustive
    # constant-codeword agreement histogram to avoid q^2 repeated comparisons.
    code=[]
    if k>1:
        for cs in product(range(F.q),repeat=k):
            code.append((cs,[F.sum(F.mul(c,F.power(x,j)) for j,c in enumerate(cs)) for x in dom]))
    nodes=[]
    for z in range(F.q):
        w=[F.add(a,F.mul(z,b)) for a,b in zip(u0,u1)]
        entries=[((c,),[c]*len(dom)) for c,count in Counter(w).items() if count>=s] if k==1 else code
        for cs,cw in entries:
            support=[i for i,(a,b) in enumerate(zip(w,cw)) if a==b]
            if len(support)>=s:
                nodes.append({"scalar":z,"codeword_coefficients":list(cs),
                              "codeword":cw,"agreement_support":support})
    return nodes


def experiment(name,F,dom,k,s,u0,u1):
    result=census(F,dom,k,s,u0,u1)
    nodes=brute(F,dom,k,s,u0,u1)
    assert result['bad_scalars']==sorted({x['scalar'] for x in nodes})
    # Every subset provenance must be witnessed by exactly one decoded word.
    for z,subs in result['subset_provenance'].items():
        for S in subs:
            assert sum(x['scalar']==z and set(S)<=set(x['agreement_support']) for x in nodes)==1
    return {"name":name,"p":F.p,"extension_degree":F.e,"q":F.q,"modulus":F.f,
            "domain":dom,"k":k,"s":s,"u0":u0,"u1":u1,
            "bad_scalar_count":len(result['bad_scalars']),**result,"nodes":nodes,
            "independent_scalar_codeword_check":"passed"}


def main():
    F9=Field(3,[1,0,1])
    # Degree-six adapter is actual field arithmetic, not prime arithmetic mod 729.
    modulus=next(list(cs)+[1] for cs in product(range(3),repeat=6)
                 if cs[0] and irreducible(list(cs)+[1],3))
    F729=Field(3,modulus)
    F3=Field(3,[0,1])
    rows=[]
    for F in [F3,F9,F729]:
        rows.append(experiment(f"embedded_base_stack_q{F.q}",F,[0,1,2],1,2,[0,1,2],[0,1,1]))
    assert all(r['bad_scalars']==[1,2] for r in rows)
    dom=subgroup(F9,4)
    rows.append(experiment("F9_mu4_nonbase_stack",F9,dom,2,3,[0,1,3,5],[1,3,4,8]))
    rows.append(experiment("F9_mu4_correlated_control",F9,dom,2,3,dom,[F9.mul(3,x) for x in dom]))
    assert rows[-1]['bad_scalar_count']==9 and rows[-1]['common_subsets']
    dom6=subgroup(F729,4)
    rows.append(experiment("F729_mu4_nonbase_stack",F729,dom6,1,3,[3,3,3,4],[0,0,0,3]))
    assert rows[-1]['bad_scalar_count']==729 and rows[-1]['common_subsets']
    # A non-correlated degree-six example: precisely one scalar makes the first
    # three coordinates equal, while the fourth is arbitrary.
    rows.append(experiment("F729_mu4_noncorrelated_stack",F729,dom6,1,3,[0,3,6,9],[0,1,2,3]))
    data={"schema":"extension-field-proximity-probe/v1","encoding":"base-p polynomial coefficients, least significant first",
          "checks":{"F9":checks(F9),"F729":checks(F729)},"experiments":rows,
          "scope":"tiny exact extension-field adapter tests; not full prize parameters or a worst-case bound"}
    dest=Path(__file__).with_name('extension_results.json')
    dest.write_text(json.dumps(data,indent=2)+'\n')
    for row in rows:
        print(row['name'], 'bad scalars',row['bad_scalar_count'],'nodes',len(row['nodes']),
              'common',bool(row['common_subsets']))
    print('GF729 modulus',modulus,'Saved',dest)


if __name__=='__main__':
    main()
