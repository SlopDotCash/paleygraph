#!/usr/bin/env python3
"""Exact centered product moments for two halves of a dyadic subgroup.

NumPy is used only for explicitly bounded int64 convolution steps.
All moment products and principal-frequency subtractions use Python integers.
The uniform mixed-moment hypothesis is not proved by these finite checks.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import comb, factorial
from pathlib import Path
import json

import numpy as np

from kernel_discriminant import generator, prime
from quotient_moments import quotient

ROOT=Path(__file__).resolve().parents[1]


def intrinsic_even_count(k, q):
    if q % 2:
        return 0
    # q! [z^(q/2)] (sum_j z^j/(j!)^2)^(k/2).
    degree=q//2
    coefficients=[Fraction(1)]+[Fraction(0)]*degree
    for _ in range(k//2):
        following=[Fraction(0)]*(degree+1)
        for i,c in enumerate(coefficients):
            for j in range(degree-i+1):
                following[i+j]+=c/Fraction(factorial(j)**2)
        coefficients=following
    value=factorial(q)*coefficients[degree]
    assert value.denominator==1
    return int(value)


def exp_bounds(x, terms=80):
    assert x>=0 and x<terms+2
    term=Fraction(1)
    total=term
    for j in range(1,terms+1):
        term*=Fraction(x,j)
        total+=term
    following=term*Fraction(x,terms+1)
    return total,total+following/(1-Fraction(x,terms+2))


def coset_pairing(p, subgroup, representatives, g):
    pairing=np.empty(len(representatives),dtype=np.int64)
    assert p*p<=np.iinfo(np.int64).max
    for start in range(0,len(representatives),4096):
        rr=(representatives[start:start+4096]*g)%p
        partner=((rr[:,None]*subgroup)%p).min(axis=1)
        indices=np.searchsorted(representatives,partner)
        assert np.array_equal(representatives[indices],partner)
        pairing[start:start+len(rr)]=indices
    assert pairing[0]==0
    assert np.array_equal(pairing[pairing],np.arange(len(pairing)))
    assert np.all(pairing[1:]!=np.arange(1,len(pairing)))
    return pairing


def level(p,k,g,q,dense_check=False):
    assert q>=4 and q%2==0
    assert pow(g,2*k,p)==1 and pow(g,k,p)==p-1
    subgroup,representatives,neighbors=quotient(p,k,pow(g,2,p))
    pairing=coset_pairing(p,subgroup,representatives,g)
    counts=np.zeros(len(representatives),dtype=np.int64)
    counts[0]=1
    dense=[1]+[0]*(p-1) if dense_check else None
    moments=[]
    for exponent in range(1,q+1):
        bound=k*int(counts.max())
        assert bound<=np.iinfo(np.int64).max
        counts=counts[neighbors].sum(axis=1,dtype=np.int64)
        assert np.all(counts>=0)
        assert int(counts[0])+k*sum(map(int,counts[1:]))==k**exponent
        mixed=int(counts[0])**2+k*sum(int(a)*int(b) for a,b in zip(counts[1:],counts[pairing[1:]]))
        pure=int(counts[0])**2+k*sum(int(a)**2 for a in counts[1:])
        centered=p*mixed-k**(2*exponent)
        if exponent%2==0:
            assert centered>=0
        intrinsic=intrinsic_even_count(k,exponent)**2
        assert mixed>=intrinsic
        if dense is not None:
            following=[0]*p
            for x,c in enumerate(dense):
                for h in subgroup:
                    following[(x+int(h))%p]+=c
            dense=following
            assert all(int(counts[i])==dense[int(a)] for i,a in enumerate(representatives))
            dense_mixed=sum(dense[x]*dense[(g*x)%p] for x in range(p))
            assert mixed==dense_mixed
        moments.append({"q":exponent,"balanced_zero_count":mixed,
                        "pure_additive_energy_at_order_q":pure,
                        "centered_numerator":centered,"nonprincipal_denominator":p-1,
                        "intrinsic_complex_count":intrinsic,"extra_over_complex":mixed-intrinsic,
                        "int64_pre_step_bound":bound,
                        "CM_C1_at_even_order": centered<=(p-1)*(exponent*k)**exponent if exponent%2==0 else None})
    # Independent low-order counts from the two literal cosets.
    kk=[int(x) for x in subgroup]
    ll=[g*x%p for x in kk]
    sums_k=Counter((a+b)%p for a in kk for b in kk)
    sums_l=Counter((a+b)%p for a in ll for b in ll)
    b=sum(c*sums_l.get((-s)%p,0) for s,c in sums_k.items())
    assert moments[0]["balanced_zero_count"]==0
    assert moments[1]["balanced_zero_count"]==b
    critical_pass=moments[-1]["centered_numerator"] <= (p-1)*(q*k)**q
    assert critical_pass
    negative_constant=None
    mean=Fraction(moments[0]["centered_numerator"],p-1)
    second=Fraction(moments[1]["centered_numerator"],p-1)
    if second==mean*mean:
        assert mean<0
        assert all(Fraction(row["centered_numerator"],p-1)==mean**row["q"] for row in moments)
        negative_constant={"value":str(mean),"variance":0,"positive_part_identically_zero":True}
    print(f"p={p}, child order={k}, q={q}: centered mixed condition C=1 passed",flush=True)
    return {"p":p,"child_order":k,"parent_order":2*k,"coset_generator":g,
            "q_checked":q,"dense_crosscheck":dense_check,"moments":moments,
            "CM_C1_at_requested_order":critical_pass,"constant_negative_product_certificate":negative_constant}


def tower(p,n,g,q,compute_small_levels=False):
    assert pow(g,n,p)==1 and pow(g,n//2,p)==p-1
    results=[]
    k=2
    while k<n:
        child_generator=pow(g,n//(2*k),p)
        if k<=q and not compute_small_levels:
            assert k*k<=q*k
            results.append({"child_order":k,"parent_order":2*k,"q":q,
                            "certificate":"Pointwise |eta_K eta_gK| <= k^2 <= q*k",
                            "CM_C1_at_requested_order":True})
        else:
            results.append(level(p,k,child_generator,q,dense_check=compute_small_levels))
        k*=2
    return {"p":p,"n":n,"generator":g,"q":q,"levels":results,
            "whole_tower_CM_C1_checked":all(x["CM_C1_at_requested_order"] for x in results)}


def main():
    small=[tower(17,16,generator(17,16),8,True),tower(97,32,generator(97,32),8,True)]
    p,n,q=6700417,64,12
    m=(p-1)//n
    assert exp_bounds(q)[0]>m and exp_bounds(q-2)[1]<m
    target=tower(p,n,2,q)
    target["q_is_smallest_even_integer_at_least_log_index"]=True
    # Existing independent whole-subgroup moments check the conclusion,
    # without re-running or refining the completed spectral maximum certificate.
    whole_path=ROOT/"results/quotient_moments.json"
    whole=json.loads(whole_path.read_text())["target"]
    assert (whole["p"],whole["n"],whole["depth"])==(p,n,q)
    target["existing_q_moment_conclusion_checked"]=whole["moments"][-1]["centered_moment"] <= (p-1)*(4*q*n)**q
    assert target["existing_q_moment_conclusion_checked"]
    # At every dyadic child order k>=64 and p<=(2k)^4, the q=6
    # balanced count exceeds the complex count, even from the principal term.
    k=64
    intrinsic6=intrinsic_even_count(k,6)
    assert intrinsic6==15*k**3-45*k**2+40*k
    assert k**12 > (2*k)**4*intrinsic6**2
    assert k*k > 16*15**2  # Monotone condition proving the whole family.
    forcing={"q":6,"minimum_dyadic_child_order":64,"minimum_parent_order":128,
             "base_intrinsic_single_color_count":intrinsic6,
             "base_intrinsic_balanced_count":intrinsic6**2,
             "base_principal_lower_bound_at_p_equal_parent_fourth_power":k**12//(2*k)**4,
             "base_forced_extra_count_at_least":k**12//(2*k)**4-intrinsic6**2,
             "strict_integer_margin":k**12-(2*k)**4*intrinsic6**2,
             "scope":"All dyadic k>=64 and all eligible primes p<=(2k)^4"}
    witness_p=67403009
    assert prime(witness_p) and (witness_p-1)%(2*k)==0 and (2*k)**4<=4*witness_p<=4*(2*k)**4
    forcing["known_quartic_prime_witness"]={"p":witness_p,"parent_order":2*k,
        "balanced_count_lower_bound":(k**12+witness_p-1)//witness_p,
        "complex_count":intrinsic6**2,
        "forced_extra_count_at_least":(k**12+witness_p-1)//witness_p-intrinsic6**2}
    assert Fraction(5,4)**4>2 and 4*Fraction(5,4)+3==8
    paths=["experiments/mixed_high_moments.py","experiments/quotient_moments.py",
           "experiments/kernel_discriminant.py","experiments/paley_exact.py",
           "experiments/cyclotomic_norm_audit.py","results/quotient_moments.json"]
    result={"status":"passed; uniform product-moment and Paley bounds unproved",
            "source_and_input_sha256":{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in paths},
            "numpy_version":np.__version__,"dense_validation_towers":small,
            "critical_depth_target_tower":target,"forced_extra_balanced_relations":forcing}
    path=ROOT/"results/mixed_high_moments.json"
    path.write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps({"status":result["status"],"critical_depth":q,
                      "target_nontrivial_computed_levels":sum("moments" in x for x in target["levels"]),
                      "target_trivial_certified_levels":sum("moments" not in x for x in target["levels"]),
                      "output":str(path)},indent=2))


if __name__=="__main__":
    main()
