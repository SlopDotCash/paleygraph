#!/usr/bin/env python3
"""Exact orbit/fiber ledger and finite quartic dyadic norm classification.

No asymptotic bound on the general unbalanced remainder is asserted.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations_with_replacement
from math import factorial, gcd, isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(ok, name):
    assert ok, name
    COUNTS[name] += 1


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def factor(v):
    out = []
    d = 2
    while d * d <= v:
        if v % d == 0:
            exponent = 0
            while v % d == 0:
                v //= d
                exponent += 1
            out.append([d, exponent])
        d += 1 if d == 2 else 2
    if v > 1:
        out.append([v, 1])
    return out


def norm_recursive(a):
    # Norm in Z[X]/(X^d+1), d a power of two. Pair roots X and -X.
    a = list(a)
    while len(a) > 1:
        d = len(a) // 2
        e, o = a[::2], a[1::2]
        b = [0] * d
        for part, shift, sign in [(e, 0, 1), (o, 1, -1)]:
            for i, x in enumerate(part):
                for j, y in enumerate(part):
                    q, r = divmod(i + j + shift, d)
                    b[r] += sign * x * y * (-1 if q % 2 else 1)
        a = b
    return abs(a[0])


def norm_determinant(a):
    # Independent exact determinant of multiplication by a modulo X^d+1.
    d = len(a)
    M = [[Fraction(0) for _ in range(d)] for _ in range(d)]
    for j in range(d):
        for k, value in enumerate(a):
            q, r = divmod(j + k, d)
            M[r][j] += value * (-1 if q % 2 else 1)
    det = Fraction(1)
    for j in range(d):
        pivot = next((k for k in range(j, d) if M[k][j]), None)
        if pivot is None:
            return 0
        if pivot != j:
            M[j], M[pivot] = M[pivot], M[j]
            det = -det
        value = M[j][j]
        det *= value
        for k in range(j + 1, d):
            ratio = M[k][j] / value
            for ell in range(j + 1, d):
                M[k][ell] -= ratio * M[j][ell]
    assert det.denominator == 1
    return abs(det.numerator)


def norm_classification(n):
    d = n // 2
    seen = set()
    records = []
    patterns = Counter()
    exceptional = defaultdict(list)
    for word in combinations_with_replacement(range(n), 6):
        support = set(word)
        if any((a + d) % n in support for a in support):
            continue
        representative = min(tuple(sorted((a - b) % n for a in word))
                             for b in support)
        if representative in seen:
            continue
        seen.add(representative)
        coeff = [representative.count(j) - representative.count(j + d)
                 for j in range(d)]
        content = 0
        for c in coeff:
            content = gcd(content, c)
        coeff = [c // content for c in coeff]
        Q = sum(c * c for c in coeff)
        norm = norm_recursive(coeff)
        check(norm > 0 and norm == norm_determinant(coeff),
              'independent_integer_norm_determinants')
        check(norm * norm <= Q ** d, 'norm_AM_GM_bound')
        fac = factor(norm)
        reconstructed = 1
        for p, e in fac:
            check(prime(p), 'factor_primality')
            reconstructed *= p ** e
        check(reconstructed == norm, 'exact_norm_factorization')
        pattern = ','.join(map(str, sorted(Counter(representative).values(),
                                          reverse=True)))
        patterns[pattern] += 1
        eligible = [p for p, _ in fac
                    if p % n == 1 and n ** 4 // 4 <= p <= n ** 4]
        for p in eligible:
            exceptional[p].append({'exponents': representative,
                                   'pattern': pattern, 'norm': norm})
        records.append({'exponents': representative, 'content': content,
                        'primitive_coefficients': coeff, 'squared_norm_mean': Q,
                        'norm': norm, 'factorization': fac,
                        'eligible_prime_factors': eligible})
    expected = [] if n < 16 else [33713, 37201, 41521]
    check(sorted(exceptional) == expected, 'complete_exceptional_prime_classification')
    check(all(row['pattern'] == '4,1,1'
              for rows in exceptional.values() for row in rows),
          'exceptional_relation_multiplicity')
    check(all(sum(row['exponents']) % 2 == 1
              for rows in exceptional.values() for row in rows),
          'exceptional_relation_product_is_nonsquare')
    return {'n': n, 'scaling_representatives': len(seen),
            'multiplicity_pattern_counts': dict(sorted(patterns.items())),
            'exceptional_primes': dict(sorted(exceptional.items())),
            'norm_ledger': records}


def subgroup(p, n):
    for a in range(2, p):
        z = pow(a, (p - 1) // n, p)
        if pow(z, n // 2, p) != 1:
            H = sorted(pow(z, j, p) for j in range(n))
            check(len(set(H)) == n and pow(z, n, p) == 1,
                  'exact_dyadic_subgroup_order')
            return z, H
    raise AssertionError('No generator')


def weight(word):
    value = factorial(len(word))
    for multiplicity in Counter(word).values():
        value //= factorial(multiplicity)
    return value


def normalized_fibers(H, p):
    n = len(H)
    pairs = defaultdict(list)
    for a, b in combinations_with_replacement(H, 2):
        pairs[(1 + a + b) % p].append((a, b, 1 if a == b else 2))
    total, intrinsic, opposite_free = Counter(), Counter(), Counter()
    patterns = Counter()
    witness = None
    for tail in combinations_with_replacement(H, 3):
        candidates = pairs.get(-sum(tail) % p, ())
        if not candidates:
            continue
        inv = pow(tail[0] * tail[1] * tail[2] % p, -1, p)
        for a, b, wp in candidates:
            word = (1, a, b, *tail)
            mass = n * wp * weight(tail)
            rho = -a * b * inv % p
            total[rho] += mass
            c = Counter(word)
            is_intrinsic = all(c[h] == c[-h % p] for h in c)
            is_free = all(-h % p not in c for h in c)
            if is_intrinsic:
                intrinsic[rho] += mass
            if is_free:
                opposite_free[rho] += mass
                pattern = ','.join(map(str, sorted(c.values(), reverse=True)))
                patterns[pattern] += mass
                if witness is None:
                    witness = word
            check(is_intrinsic or is_free, 'no_opposite_containing_nonintrinsic_word')
    return total, intrinsic, opposite_free, patterns, witness


def inspect_field(p, n):
    check(prime(p) and (p - 1) % n == 0
          and n ** 4 // 4 <= p <= n ** 4, 'eligible_actual_prime')
    g, H = subgroup(p, n)
    Hset = set(H)
    check(all(a * b % p in Hset for a in H for b in H),
          'actual_multiplicative_closure')
    r2 = Counter((a + b) % p for a in H for b in H)
    r3 = Counter()
    for s, multiplicity in r2.items():
        for a in H:
            r3[(s + a) % p] += multiplicity
    E2 = sum(c * c for c in r2.values())
    E3 = sum(c * c for c in r3.values())
    T6 = 15 * n ** 3 - 45 * n * n + 40 * n
    check(E2 == 3 * n * n - 3 * n, 'intrinsic_fourth_energy')
    total, intrinsic, free, patterns, witness = normalized_fibers(H, p)
    check(sum(total.values()) == E3, 'independent_sixth_energy_count')
    squares = {h * h % p for h in H}
    for rho in H:
        I = (6 * n ** 3 - 9 * n * n + 4 * n if rho == 1 else
             18 * n * (n - 2) if rho in squares else 0)
        check(intrinsic[rho] == I, 'exact_intrinsic_product_ratio_fiber')
        check(total[rho] == intrinsic[rho] + free[rho],
              'exact_opposite_free_residual_fiber')
    check(sum(intrinsic.values()) == T6, 'intrinsic_sixth_total')
    R6 = sum(free.values())
    expected = 480 if n == 16 and p in [33713, 37201, 41521] else 0
    check(R6 == expected and E3 == T6 + R6, 'finite_class_total_upper_and_equality')
    check(not R6 or patterns == {'4,1,1': 480}, 'actual_exception_multiplicity')
    w = Counter((a * b - a - b) % p for a in H for b in H)
    orbit_mass = Counter()
    representatives = {}
    for z, c in w.items():
        if z:
            orbit = pow(z, n, p)
            orbit_mass[orbit] += c
            representatives[orbit] = z
    for orbit, W in orbit_mass.items():
        check(W == r3[representatives[orbit]], 'orbit_mass_is_actual_triple_count')
    check(r3[0] == n * w[0], 'zero_orbit_correction')
    check(E3 == n * n * w[0] ** 2 + n * sum(W * W for W in orbit_mass.values()),
          'full_orbit_mass_identity')
    for rho in H:
        inv = pow(rho, -1, p)
        N = sum(c * w[z * inv % p] for z, c in w.items())
        check(total[rho] == n * N, 'actual_product_ratio_correlation')
    record = {'p': p, 'n': n, 'generator': g, 'E2': E2, 'E3': E3,
              'R6': R6, 'sum_rho_N': E3 // n,
              'balanced_fixed_R6': free[1], 'w_zero': w[0]}
    if R6:
        check(free[1] == 0, 'exception_is_product_unbalanced_at_fixed_split')
        check(all(rho not in squares for rho, c in free.items() if c),
              'all_exceptional_residual_fibers_are_nonsquare')
        # A 4,1,1 word could balance at some different split; check all ten.
        from itertools import combinations
        for I_tail in combinations(range(1, 6), 2):
            I = {0, *I_tail}
            left = right = 1
            for i, h in enumerate(witness):
                if i in I:
                    left = left * h % p
                else:
                    right = right * h % p
            check((left + right) % p != 0, 'exception_witness_unbalanced_all_splits')
        check(r2[-4 % p] == 2 and R6 == 15 * n * r2[-4 % p],
              'exception_exact_411_count')
        record.update({'H': H, 'normalized_opposite_free_witness': witness,
                       'ordered_pairs_summing_to_minus_four':
                       [[a, b] for a in H for b in H if (a + b + 4) % p == 0],
                       'opposite_free_product_ratio_fibers': dict(sorted(free.items()))})
    return record


def main():
    norms = [norm_classification(n) for n in [4, 8, 16]]
    fields = [inspect_field(p, n) for n in [4, 8, 16]
              for p in range(n ** 4 // 4 + 1, n ** 4 + 1, n) if prime(p)]
    inputs = ['experiments/parallel24_subgroup_unbalanced_2026_09_05.py',
              'research/parallel23-subgroup-upper-2026-09-05.md',
              'research/parallel21-subgroup-next-input-2026-09-05.md',
              'research/mixed-periods-and-shifted-energy.md']
    result = {'status': 'PASS',
              'scope': 'Exact finite-class subgroup upper bound for n=4,8,16; no uniform asymptotic aggregate saving.',
              'counts': dict(sorted(COUNTS.items())),
              'prime_counts_by_order': dict(sorted(Counter(r['n'] for r in fields).items())),
              'norm_classifications': norms, 'actual_quartic_fields': fields,
              'input_sha256': {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in inputs}}
    path = ROOT / 'results/parallel24_subgroup_unbalanced_2026_09_05.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'counts': result['counts'],
                      'prime_counts_by_order': result['prime_counts_by_order'],
                      'norm_representatives': {r['n']: r['scaling_representatives'] for r in norms},
                      'exceptions': [r for r in fields if r['R6']]}, indent=2))


if __name__ == '__main__':
    main()
