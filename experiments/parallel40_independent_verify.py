#!/usr/bin/env python3
"""Independent exact checks of pass40 witnesses; no Lean or field-sized scan."""
from collections import Counter
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import isqrt
from pathlib import Path
import json


def determinant(matrix):
    """Fraction-free Bareiss, independent of the norm-descent certificate."""
    a = [row[:] for row in matrix]
    sign, previous = 1, 1
    for k in range(len(a) - 1):
        if not a[k][k]:
            r = next(r for r in range(k + 1, len(a)) if a[r][k])
            a[k], a[r] = a[r], a[k]
            sign = -sign
        pivot = a[k][k]
        for i in range(k + 1, len(a)):
            for j in range(k + 1, len(a)):
                numerator = a[i][j] * pivot - a[i][k] * a[k][j]
                assert numerator % previous == 0
                a[i][j] = numerator // previous
            a[i][k] = 0
        previous = pivot
    return sign * a[-1][-1]


def prime(p):
    assert p > 1 and all(p % d for d in range(2, isqrt(p) + 1))


def shell():
    p, q, n = 6700417, 641, 64
    for modulus in (p, q):
        prime(modulus)
        assert pow(2, n // 2, modulus) == modulus - 1
    N = n // 2
    F = [(pow(2, j, p) + p // 2) % p - p // 2 for j in range(N)]
    # Matrix columns are F*X^j modulo X^N+1.
    matrix = [[F[(i - j) % N] * (-1 if i < j else 1)
               for j in range(N)] for i in range(N)]
    norm = determinant(matrix)
    assert norm == p**31 * 1217
    prime(1217)
    certificate_path = Path(__file__).resolve().parents[1]/"results/parallel40_shell_converse_2026_09_06.json"
    reciprocal = json.loads(certificate_path.read_text())["scalar_reciprocal_T"]
    assert [sum(matrix[i][j]*reciprocal[j] for j in range(N))
            for i in range(N)] == [1217*p] + [0]*(N-1)
    reciprocal_matrix = [[reciprocal[(i-j) % N]*(-1 if i < j else 1)
                          for j in range(N)] for i in range(N)]
    assert determinant(reciprocal_matrix) == p*1217**31
    V = Fraction(sum(c*c for c in F), p*p)
    assert V == Fraction(6120237, 6700417) < 1
    assert 2**N + 1 == p*q
    lower = Fraction(n) - Fraction(1936, 49)
    upper = Fraction(1 + 9*26, 10)
    assert q < 26**2 and lower - upper == Fraction(97, 98)
    return {"p": p, "n": n, "V": str(V), "norm_F": norm,
            "norm_method": "32-by-32 Bareiss determinant",
            "auxiliary_contradiction_margin": str(lower-upper),
            "scope": "Finite arithmetic; nonprincipality follows from the companion ordinary proof."}


def distinct():
    p, n, g = 215535361, 128, 25525303
    prime(p)
    H = {pow(g, j, p) for j in range(n)}
    assert len(H) == n and p-1 in H and n**4//4 <= p <= n**4
    f = Counter((x-y) % p for x in H for y in H)
    E = sum(v*v for v in f.values())
    d = Fraction(E, 16*n*n)
    assert E == 61056 and d == Fraction(477, 2048)
    label = {x: pow(x, n, p) for x in f if x}
    levels = {j: {x for x in label if 2**(j-1)*d < f[x] <= 2**j*d}
              for j in (4, 5)}
    assert (len(levels[4]), len(levels[5])) == (6144, 768)
    reps = {}
    for x in label:
        reps.setdefault(label[x], x)
    classes = {j: {label[x] for x in levels[j]} for j in levels}
    m = Counter()
    for a in classes[4]:
        for b in classes[5]:
            for c in classes[5]:
                if len({a, b, c}) != 3:
                    continue
                ci = pow(c, -1, p)
                m[a*ci % p, b*ci % p] += 1
    rho = {}
    for u, v in m:
        # A representative of the quotient label u need not lie in f's support.
        # Recover it from one source triple, using field representatives.
        for c in classes[5]:
            a = u*c % p
            if a in classes[4]:
                representative = reps[a] * pow(reps[c], -1, p) % p
                break
        rho[u, v] = sum(1 for h in H
                         if (y := (representative*h + 1) % p)
                         and pow(y, n, p) == v)
    weighted = n * sum(m[k]*rho[k] for k in m)
    unweighted = n * sum(rho.values())
    direct = direct_weighted = 0
    for y in levels[5]:
        for z in levels[5]:
            x = (y-z) % p
            if x in levels[4] and len({label[x], label[y], label[z]}) == 3:
                direct += 1
                direct_weighted += f[x]**2*f[y]**2*f[z]**2
    assert direct == weighted == 52736 and unweighted == 51968
    assert direct_weighted == 144850944
    collisions = sum(value*(value-1) for value in m.values())
    ratios = {a*pow(b, -1, p) % p for a in classes[5] for b in classes[5]} - {1}
    intersection_sum = 0
    for q in ratios:
        t4 = len(classes[4] & {q*x % p for x in classes[4]})
        t5 = len(classes[5] & {q*x % p for x in classes[5]})
        intersection_sum += t4*t5*(t5-1)
    assert collisions == intersection_sum == 4
    repeated_rho = max(rho[key] for key, value in m.items() if value >= 2)
    assert repeated_rho == 3 and weighted-unweighted == n*repeated_rho*collisions//2
    key = (pow(674430, n, p), pow(2246831, n, p))
    assert m[key] == 2 and rho[key] == 3
    solutions = [(19584002, 19584003), (108122689, 108122690),
                 (162237922, 162237923)]
    for x, y in solutions:
        assert y-x == 1 and (pow(x, n, p), pow(y, n, p)) == key
    assert len(m) == 1438 and sum(m.values()) == 1440
    return {"p": p, "n": n, "E": E, "d": str(d),
            "direct_incidence": direct, "normalized_weighted": weighted,
            "discarded_multiplicity": unweighted, "missing_incidence": direct-unweighted,
            "witness_multiplicity": m[key], "witness_rho": rho[key],
            "normalized_image_size": len(m), "normalized_mass": sum(m.values()),
            "direct_weighted_block": direct_weighted,
            "collision_sum": collisions, "intersection_sum": intersection_sum,
            "scope": "Actual dyadic levels and three distinct edge cosets; not an asymptotic energy counterexample."}


def sidon():
    sets = translations = pairs = 0
    for p in (7, 11, 13):
        chi = [0] + [1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]
        for k in (2, 3, 4):
            for Vtuple in combinations(range(p), k):
                if len({(a+b) % p for a, b in combinations_with_replacement(Vtuple, 2)}) != k*(k+1)//2:
                    continue
                V = set(Vtuple)
                sets += 1
                F = [sum(chi[(x-v) % p] for v in V) for x in range(p)]
                assert sum(y*y for y in F) == p*k-k*k
                overlaps = [len(V & {(v+s) % p for v in V}) for s in range(p)]
                for s in range(1, p):
                    error = sum((F[(x+s) % p]-F[x])**2 for x in range(p))
                    assert overlaps[s] <= 1 and error == 2*p*(k-overlaps[s])
                    assert error <= 2*p*k
                    translations += 1
                for size in (1, 2):
                    for T in combinations(range(p), size):
                        assert len({(v+t) % p for v in V for t in T})*(k+size-1) >= k*k*size
                        zero = int(0 in T)
                        counts = [sum((x-t) % p in V for t in T) for x in range(p)]
                        l1_scaled = sum(abs(size*int(x in V)-counts[x]) for x in range(p))
                        assert l1_scaled == 2*(k*size-sum(overlaps[t] for t in T))
                        assert l1_scaled >= 2*(k-1)*(size-zero)
                        err = sum((size*F[x]-sum(F[(x-t) % p] for t in T))**2 for x in range(p))
                        assert k*err >= p*(k-1)**2*(size-zero)**2
                        pairs += 1
    assert (sets, translations, pairs) == (538, 5916, 42628)
    return {"sidon_sets": sets, "nonzero_translations": translations,
            "set_pair_and_smoothing_checks": pairs,
            "scope": "Exact identities and relative accuracy only; larger Croot-Sisask tolerance is not refuted."}


def main():
    result = {}
    for name, check in (("shell", shell), ("distinct", distinct), ("sidon", sidon)):
        result[name] = check()
        print(name, "passed", flush=True)
    out = Path(__file__).resolve().parents[1]/"results/parallel40_independent_verification_2026_09_06.json"
    out.write_text(json.dumps(result, indent=2)+"\n")
    print(out)


if __name__ == "__main__":
    main()
