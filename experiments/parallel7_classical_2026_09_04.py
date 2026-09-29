#!/usr/bin/env python3
"""Exact conditioned-moment and affine quartic countermodel checks.

Finite checks support identities; the asymptotic countermodel is proved
in the note. This does not construct a character-moment counterexample.
"""
from collections import Counter
from fractions import Fraction as Q
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def falling(n, k):
    out = 1
    for i in range(k):
        out *= n-i
    return out


def inclusion(n, p, k):
    return Q(falling(n, k), falling(p, k))


def char_table(p):
    return [0] + [1 if pow(x, (p-1)//2, p) == 1 else -1
                  for x in range(1, p)]


def row_sums(p, b, chi):
    return tuple(sum(chi[(x-y) % p] for y in b) for x in range(p))


def moments(p, b, chi):
    rows = row_sums(p, b, chi)
    return sum(s**4 for s in rows), sum(rows[x]**2 for x in b)


def mean_variance(p, n):
    v = n*(p-n)
    f = falling(n, 3)*falling(p-n, 3)
    mean = Q(v*(3*v-2*p+1), p-2)
    if p % 4 == 1:
        variance = Q(24*f*(p-1)*((p-5)*(p+3)*v
                       - 3*(p**3-6*p*p+4*p+3)),
                       (p-2)**2*(p-3)*(p-4)*(p-6)*(p-7))
    else:
        variance = Q(24*f*(p+1)*((p+1)*v-3*(p*p-3*p+3)),
                       (p-2)**2*(p-4)*(p-5)*(p-6))
    return mean, variance


def conditional_coefficients(p, n):
    a1, a2, a3, a4 = [inclusion(n-4, p-4, j)
                      for j in range(1, 5)]
    c4 = 1-4*a1+6*a2-4*a3+a4
    ca = 6*a1-18*a2+18*a3-6*a4
    c2 = -4*a1+16*a2-20*a3+8*a4
    ctt = 3*a2-6*a3+3*a4
    ct = a1-7*a2+12*a3-6*a4
    constant = ((ca*(p-5)+c2)*(4*p-16)
                + ctt*(p*(p-5)**2+8*(p-5)+4)
                + ct*(p*(p-5)+4))
    assert c4 == inclusion(p-n, p-4, 4)
    return c4, ca, constant


def conditioned_average(p, n, core, chi):
    m4, aq = moments(p, core, chi)
    c4, ca, constant = conditional_coefficients(p, n)
    return c4*m4+ca*aq+constant


def completion_checks():
    cases = contributions = 0
    by_field = []
    for p in (11, 13):
        chi = char_table(p)
        cores = list(combinations(range(p), 4))
        if p == 13:
            cores = cores[:20]
        cache = {}
        count0 = cases
        contribution0 = contributions
        for n in range(4, 8):
            for core in cores:
                outside = [x for x in range(p) if x not in core]
                total = 0
                count = 0
                for rest in combinations(outside, n-4):
                    b = tuple(sorted(core+rest))
                    if b not in cache:
                        cache[b] = moments(p, b, chi)[0]
                    total += cache[b]
                    count += 1
                assert Q(total, count) == conditioned_average(p, n, core, chi)
                cases += 1
                contributions += count
        for core in cores:
            assert conditioned_average(p, p, core, chi) == 0
            cases += 1
        by_field.append({"p": p, "conditioned_identities": cases-count0,
                         "completion_contributions": contributions-contribution0,
                         "distinct_moments_evaluated": len(cache)})
    return {"identities": cases, "completion_contributions": contributions,
            "fields": by_field}


def local_checks():
    checked = 0
    fields = []
    for p in (11, 13):
        chi = char_table(p)
        sets = list(combinations(range(p), 3))
        rows = [row_sums(p, b, chi) for b in sets]
        count0 = checked
        for i, b in enumerate(sets):
            for j in range(i):
                k = len(set(b)-set(sets[j]))
                diff = [a-c for a, c in zip(rows[i], rows[j])]
                assert sum(x*x for x in diff) == 2*p*k
                assert sum(x**4 for x in diff) <= 8*p*k**3
                checked += 1
        fields.append({"p": p, "n": 3, "unordered_pairs": checked-count0})
    return {"unordered_pairs": checked, "fields": fields}


def affine_orbit(p, core):
    return {tuple(sorted((a*x+b) % p for x in core))
            for a in range(1, p) for b in range(p)}


def primitive_root(p):
    for a in range(2, p):
        if len({pow(a, j, p) for j in range(p-1)}) == p-1:
            return a
    raise AssertionError("no generator")


def cube_floor(p):
    n = 0
    while (n+1)**3 <= p:
        n += 1
    return n


def fraction_json(value):
    value = Q(value)
    return {"numerator": value.numerator, "denominator": value.denominator}


def hypergraph_check(p):
    n = cube_floor(p)
    chi = char_table(p)
    edges = set()
    orbits = []
    for core in combinations(range(n), 4):
        if core in edges:
            continue
        orbit = affine_orbit(p, core)
        assert not edges.intersection(orbit)
        edges.update(orbit)
        orbits.append((core, len(orbit)))
    h = len(edges)
    assert h <= p*(p-1)*comb(n, 4)
    generator = primitive_root(p)
    assert {tuple(sorted((x+1) % p for x in e)) for e in edges} == edges
    assert {tuple(sorted((generator*x) % p for x in e)) for e in edges} == edges
    d1, d2, d3 = Counter(), Counter(), Counter()
    for e in edges:
        d1.update(e)
        d2.update(combinations(e, 2))
        d3.update(combinations(e, 3))
    assert len(d1) == p and len(set(d1.values())) == 1
    assert len(d2) == comb(p, 2) and len(set(d2.values())) == 1
    sums = [h*h, sum(v*v for v in d1.values()),
            sum(v*v for v in d2.values()),
            sum(v*v for v in d3.values()), h]
    assert sums[1] == Q(16*h*h, p)
    assert sums[2] == Q(36*h*h, comb(p, 2))
    assert sums[3] <= 4*(p-3)*h
    pairs = [sum((-1)**(j-k)*comb(j, k)*sums[j]
                 for j in range(k, 5)) for k in range(5)]
    assert min(pairs) >= 0 and sum(pairs) == h*h and pairs[4] == h
    eu = h*inclusion(n, p, 4)
    eu2 = sum(pairs[k]*inclusion(n, p, 8-k) for k in range(5))
    vu = eu2-eu*eu
    assert vu >= 0
    upper = (inclusion(n, p, 7)*sums[1]
             + inclusion(n, p, 6)*sums[2]
             + inclusion(n, p, 5)*4*(p-3)*h
             + inclusion(n, p, 4)*h)
    assert vu <= upper
    monomial_upper = (Q(n**15, 36*p**4)+Q(n**14, 8*p**4)
                      +Q(n**9, 6*p*p)+Q(n**8, 24*p*p))
    assert upper <= monomial_upper
    vz = 144*p*vu
    assert vz <= 28*n**6+24*n**5
    mean, variance = mean_variance(p, n)
    emu = inclusion(n, p, 4)*sum(size*conditioned_average(p, n, core, chi)
                                  for core, size in orbits)
    covariance = emu-mean*eu
    assert covariance*covariance <= variance*vu
    assert comb(n, 4) > eu
    # tau^2 = variance + 144p Var(U) + 24 sqrt(p) Cov(M,U).
    # Store its exact rational and sqrt(p) coefficients; positivity follows
    # from this being a variance and is certified here without float roots.
    tau_rational = variance+vz
    tau_radical = 24*covariance
    if tau_radical < 0:
        assert tau_rational*tau_rational > p*tau_radical*tau_radical
    else:
        assert tau_rational > 0
    # Exact local example: retaining m interval points plants binom(m,4)
    # edges regardless of the remaining membership choices.
    planted = tuple(range(n))
    assert sum(e in edges for e in combinations(planted, 4)) == comb(n, 4)
    min_real_moment = Q(n*n*(p-n)**2, p)
    # M+Z is nonnegative if min_real_moment >= 12 sqrt(p) E U.
    assert min_real_moment**2 >= 144*p*eu*eu
    return {"p": p, "n": n, "edges": h,
            "affine_orbits": [{"representative": core, "size": size}
                               for core, size in orbits],
            "ordered_pair_intersection_counts": pairs,
            "mean_u": fraction_json(eu), "variance_u": fraction_json(vu),
            "variance_z": fraction_json(vz),
            "covariance_m_u": fraction_json(covariance),
            "true_moment_mean": fraction_json(mean),
            "true_moment_variance": fraction_json(variance),
            "calibration_tau_squared": {
                "rational_part": fraction_json(tau_rational),
                "sqrt_p_coefficient": fraction_json(tau_radical)},
            "planted_z_sqrt_p_coefficient": fraction_json(12*(comb(n, 4)-eu)),
            "distinct_three_point_degrees": len(d3),
            "status": "all_exact_checks_passed"}


def main():
    completion = completion_checks()
    print("conditional completion checks complete", flush=True)
    local = local_checks()
    print("local perturbation checks complete", flush=True)
    hypergraphs = []
    for p in (67, 71, 127, 131):
        hypergraphs.append(hypergraph_check(p))
        print(f"affine hypergraph checks complete: p={p}", flush=True)
    sources = ["research/parallel7-classical-2026-09-04.md",
               "experiments/parallel7_classical_2026_09_04.py",
               "research/parallel6-classical-2026-09-04.md"]
    result = {"status": "all_exact_checks_passed",
              "scope": "artificial affine quartic countermodel; no SS/CS/LM counterexample",
              "conditional_completion_checks": completion,
              "local_perturbation_checks": local,
              "hypergraph_checks": hypergraphs,
              "source_sha256": {f: sha256((ROOT/f).read_bytes()).hexdigest()
                                for f in sources}}
    path = ROOT/"results/parallel7_classical_2026_09_04.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "conditioned_identities": completion["identities"],
                      "completion_contributions": completion["completion_contributions"],
                      "local_pairs": local["unordered_pairs"],
                      "hypergraph_fields": len(hypergraphs)}, indent=2))


if __name__ == "__main__":
    main()
