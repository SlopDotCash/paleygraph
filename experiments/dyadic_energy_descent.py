#!/usr/bin/env python3
"""Exact checks of dyadic energy splitting and prime support.

No uniform bound on the mixed energy is assumed or proved by these checks.
"""
from collections import Counter
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from kernel_discriminant import discriminant, evaluate, valuation
from quadruple_orbits import Ring, find_generator, orbit_counts

ROOT = Path(__file__).resolve().parents[1]


def order_mod(p, n):
    result, x = 1, p % n
    while x != 1:
        x = x*p % n
        result += 1
        assert result <= n
    return result


class QuadraticTower:
    """One unramified quadratic step X^2=g over an existing residue ring."""
    def __init__(self, base, g, smaller_order):
        self.base, self.g = base, g
        self.p, self.r, self.modulus, self.q = base.p, base.r, base.modulus, base.q**2
        self.zero, self.one = (base.zero, base.zero), (base.one, base.zero)
        # The exact residue-field degrees prove this binomial irreducible.
        assert base.q == base.p**order_mod(base.p, smaller_order)
        assert order_mod(base.p, 2*smaller_order) == 2*order_mod(base.p, smaller_order)
        assert base.power(g, smaller_order) == base.one
        assert base.power(g, smaller_order//2) == base.scalar(-1)

    def scalar(self, a):
        return (self.base.scalar(a), self.base.zero)

    def add(self, a, b):
        return (self.base.add(a[0],b[0]),self.base.add(a[1],b[1]))

    def neg(self, a):
        return (self.base.neg(a[0]),self.base.neg(a[1]))

    def sub(self, a, b):
        return self.add(a,self.neg(b))

    def mul(self, a, b):
        c=self.base
        return (c.add(c.mul(a[0],b[0]),c.mul(self.g,c.mul(a[1],b[1]))),
                c.add(c.mul(a[0],b[1]),c.mul(a[1],b[0])))

    def power(self, a, n):
        result=self.one
        while n:
            if n%2:
                result=self.mul(result,a)
            a=self.mul(a,a)
            n//=2
        return result


def partition_orbits(counts):
    result={kind:{"old":[],"balanced":[],"unbalanced":[]} for kind in ["3","2","1"]}
    for kind, representatives in counts["orbit_representatives"].items():
        for indices in representatives:
            odd=sum(i%2 for i in indices)
            category="old" if odd in [0,4] else "balanced" if odd==2 else "unbalanced"
            result[kind][category].append(indices)
    return result


def audit(ring, g, n):
    assert n >= 8 and not n & (n-1)
    assert ring.power(g,n)==ring.one and ring.power(g,n//2)==ring.scalar(-1)
    h=[ring.power(g,j) for j in range(n)]
    k,l=h[::2],h[1::2]
    size=n//2
    assert len(set(h))==n and set(k).isdisjoint(l)
    pairs=lambda a,b:Counter(ring.add(x,y) for x in a for y in b)
    kk,ll,kl,hh=pairs(k,k),pairs(l,l),pairs(k,l),pairs(h,h)
    energy_k=sum(x*x for x in kk.values())
    energy_h=sum(x*x for x in hh.values())
    b=sum(c*ll.get(ring.neg(x),0) for x,c in kk.items())
    t=sum(c*kl.get(ring.neg(x),0) for x,c in kk.items())
    assert sum(x*x for x in kl.values())==b
    assert b >= size*size and b <= energy_k
    assert t*t <= (energy_k-size*size)*b
    assert energy_h==2*energy_k+6*b+8*t
    old,new=orbit_counts(ring,k),orbit_counts(ring,h)
    types=partition_orbits(new)
    for kind in types:
        assert len(types[kind]["old"])==len(old["orbit_representatives"][kind])
    eps=len(types["3"]["unbalanced"])
    u22,u31=len(types["2"]["balanced"]),len(types["2"]["unbalanced"])
    v22,v31=len(types["1"]["balanced"]),len(types["1"]["unbalanced"])
    assert not types["3"]["balanced"]
    assert b==size*size+4*size*u22+8*size*v22
    assert t==size*(eps+3*u31+6*v31)
    # Independent polynomial-root count of balanced orbits.
    roots=[ring.power(ring.add(ring.one,k[j]),size) for j in range(1,size//2)]
    minus_boundary=ring.scalar(-2**size)
    root_u=sum(x==minus_boundary for x in roots)
    root_v=sum(ring.add(roots[i],roots[j])==ring.zero for i in range(len(roots)) for j in range(i))
    assert (root_u,root_v)==(u22,v22)
    doubled=order_mod(ring.p,n)==2*order_mod(ring.p,size)
    if doubled:
        assert (b,t)==(size*size,0)
        assert new["D_from_orbits"]==old["D_from_orbits"]
    if ring.p % n==n-1:
        assert new["D_from_orbits"]==0
    return {"p":ring.p,"precision":ring.r,"n":n,"residue_field_size":ring.q,
            "field_degree_doubles":doubled,"energy_smaller":energy_k,"energy_larger":energy_h,
            "mixed_B":b,"mixed_T":t,"new_epsilon":eps,
            "new_u_balanced":u22,"new_u_unbalanced":u31,
            "new_v_balanced":v22,"new_v_unbalanced":v31,
            "old_D":old["D_from_orbits"],"new_D":new["D_from_orbits"],
            "new_orbit_representatives":types}


def square_root_polynomial(coefficients):
    # Q(Y)=product(Y-lambda_i^2)=(-1)^d P(sqrt(Y))P(-sqrt(Y)).
    degree=len(coefficients)-1
    product=[0]*(2*degree+1)
    for i,a in enumerate(coefficients):
        for j,b in enumerate(coefficients):
            product[i+j]+=(-1)**(degree+j)*a*b
    assert all(product[j]==0 for j in range(1,len(product),2))
    return product[::2]


def main():
    old_path=ROOT/"results/quadruple_orbits.json"
    old=json.loads(old_path.read_text())
    assert old["source_sha256"]==sha256((ROOT/"experiments/quadruple_orbits.py").read_bytes()).hexdigest()
    orders={x["n"]:x for x in old["cube_identities"]}
    factored=json.loads((ROOT/"results/kernel_discriminant.json").read_text())
    assert factored["source_sha256"]==sha256((ROOT/"experiments/kernel_discriminant.py").read_bytes()).hexdigest()
    support=[]
    for n in [8,16,32,64]:
        small,large=orders[n//2],orders[n]
        u0,v0=(int(small[key],16) for key in ["U_hex","V_hex"])
        u1,v1=(int(large[key],16) for key in ["U_hex","V_hex"])
        assert u1%u0==0 and v1%v0==0
        pcoef=small["P_coefficients_ascending"]
        qcoef=square_root_polynomial(pcoef)
        disc=discriminant(pcoef)
        disc_squared=discriminant(qcoef)
        assert disc_squared % disc == 0
        j=isqrt(disc_squared//disc)
        assert j*j==disc_squared//disc
        odd_j=j>>valuation(j,2)
        c=abs(evaluate(pcoef,-2**(n//2)))
        odd_c=c>>valuation(c,2)
        assert (v1//v0)%odd_j==0
        assert (u1//u0)%odd_c==0
        record={"n":n,"primitive_U_ratio_hex":hex(u1//u0),"primitive_V_ratio_hex":hex(v1//v0),
                "balanced_U_factor_hex":hex(odd_c),"balanced_V_factor_hex":hex(odd_j),
                "unbalanced_U_factor_hex":hex(u1//u0//odd_c),
                "unbalanced_V_factor_hex":hex(v1//v0//odd_j),"complete_prime_support_checked":n<=32}
        if n<=32:
            row=next(x for x in factored["small_orders"] if x["n"]==n)
            primitive=[]
            for p in map(int,row["complete_A_prime_factorization"]):
                if p==2:continue
                u,v=valuation(u1//u0,p),valuation(v1//v0,p)
                if u or v:
                    assert p % n in [1,n//2-1]
                    primitive.append({"p":p,"new_U_valuation":u,"new_V_valuation":v})
            record["all_primitive_primes"]=primitive
            for p in map(int,row["complete_A_prime_factorization"]):
                if p>2 and v1%p==0:
                    assert p**order_mod(p,n) <= 2**(n//2)
        support.append(record)
    # Finite valuations that the descent theorem propagates to all higher orders.
    stable=[]
    for n in [32,64]:
        v=int(orders[n]["V_hex"],16)
        fixed={p:valuation(v,p) for p in [3,5,7,17,47,79,97]}
        assert fixed=={3:0,5:0,7:1,17:5,47:1,79:1,97:9}
        quotient=v//(7*17**5*47*79)
        assert quotient % 32 in [1,31]
        stable.append({"n":n,"selected_stable_V_valuations":fixed,
                       "remaining_factor_residue_mod_32":quotient%32})
    cases=[]
    for p,n,d in [(17,8,None),(17,16,None),(97,16,None),(97,32,None),
                   (193,32,None),(193,64,None),(6700417,64,None),(6700417,128,None),
                   (67403009,128,None),(67403009,256,None),
                   (1073748737,256,None),(17179869697,512,None),
                   (5,8,2),(17,32,3),(97,64,5),(7,8,3),(7,16,3),(31,32,3),(31,64,3),(47,32,5)]:
        residue=Ring(p,1,d)
        initial=find_generator(residue,n)
        for r in [1,2]:
            ring=Ring(p,r,d)
            g=ring.power(initial,ring.q**(r-1))
            cases.append(audit(ring,g,n))
    for p,m,d in [(3,8,2),(7,16,3)]:
        residue=Ring(p,1,d)
        initial=find_generator(residue,m)
        for r in [1,2]:
            base=Ring(p,r,d)
            g=base.power(initial,base.q**(r-1))
            extension=QuadraticTower(base,g,m)
            cases.append(audit(extension,(base.zero,base.one),2*m))
    # Exact constant check for the conditional mixed-energy induction:
    # 2*22+6+8*sqrt(22) <= 4*22, since sqrt(22)<19/4.
    assert 16*22 < 19*19 and 2*22+6+2*19==4*22
    source_hashes={name:sha256((ROOT/"experiments"/name).read_bytes()).hexdigest()
                   for name in ["dyadic_energy_descent.py","quadruple_orbits.py","kernel_discriminant.py","cyclotomic_norm_audit.py"]}
    result={"status":"passed; mixed energy bound and Paley conjectures unproved",
            "source_sha256":source_hashes,"prior_orbit_result_sha256":sha256(old_path.read_bytes()).hexdigest(),
            "prior_kernel_result_sha256":sha256((ROOT/"results/kernel_discriminant.json").read_bytes()).hexdigest(),
            "integer_factor_checks":support,"stable_valuation_checks":stable,
            "energy_and_descent_checks":cases}
    path=ROOT/"results/dyadic_energy_descent.json"
    path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"energy_cases":len(cases),
                      "field_doubling_cases":sum(c["field_degree_doubles"] for c in cases),
                      "factor_orders":[c["n"] for c in support],"output":str(path)},indent=2))


if __name__=="__main__":
    main()
