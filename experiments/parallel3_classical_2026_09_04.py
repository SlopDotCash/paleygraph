#!/usr/bin/env python3
"""Exact cyclotomic Fourier identities and interval cross-ratio counts.

No floating-point character sums or numerical plausibility enter these
checks. The asymptotic failed-input theorem is proved in the companion note.
"""
from collections import Counter
from hashlib import sha256
from itertools import permutations
from math import gcd, isqrt
from pathlib import Path
from random import Random
import json

from paley_exact import character

ROOT = Path(__file__).resolve().parents[1]


def conv(a, b):
    p = len(a)
    result = [0] * p
    for i, av in enumerate(a):
        if av:
            for j, bv in enumerate(b):
                if bv:
                    result[(i + j) % p] += av * bv
    return result


def shift(a, s):
    p = len(a)
    return [a[(i - s) % p] for i in range(p)]


def equal_at_root(a, b):
    # A polynomial of degree at most p-1 vanishes at a primitive pth root
    # iff all of its p coefficients are equal (a multiple of Phi_p).
    assert len(set(x - y for x, y in zip(a, b))) == 1


def scalar_poly(p, value):
    return [value] + [0] * (p - 1)


def fourier_poly(values, t):
    p = len(values)
    result = [0] * p
    for x, value in enumerate(values):
        result[-t * x % p] += value
    return result


def vector(p, chi, subset):
    answer = [0] * p
    for a, b, c, d in permutations(subset, 4):
        denominator = (c - a) * (d - b) % p
        lam = (c - b) * (d - a) * pow(denominator, -1, p) % p
        answer[lam] += chi[denominator]
    return answer


def totient(n):
    return sum(gcd(a, n) == 1 for a in range(1, n + 1))


def pattern_count(n, h):
    margin = n // 4
    q = margin // (h * (h + 1))
    patterns = set()
    for s in range(2, isqrt(q) + 1):
        for r in range(1, s):
            if gcd(r, s) != 1 or s == h * r:
                continue
            for ell in range(1, q // (s * s) + 1):
                quad = (0, (h * h + 1) * ell * r * s,
                        ell * r * (s - h * r), h * ell * s * (h * r - s))
                assert len(set(quad)) == 4
                assert max(abs(x) for x in quad) <= margin
                a, b, c, d = quad
                assert (c - b) * (d - a) == -h * h * (c - a) * (d - b)
                assert quad not in patterns
                patterns.add(quad)
    exact = sum(totient(s) * (q // (s * s)) for s in range(2, isqrt(q) + 1))
    if h >= 2:
        exact -= q // (h * h)
    assert exact == len(patterns)
    if n <= 128:
        all_quads = {tuple(a + v for v in quad)
                     for a in range(margin + 1, n - margin + 1)
                     for quad in patterns}
        assert len(all_quads) == (n - 2 * margin) * exact
        assert all(min(quad) >= 1 and max(quad) <= n for quad in all_quads)
    return (n - 2 * margin) * exact, exact


def literal_integer_count(n, h):
    # Each a,b,c determines at most one d. Solve the cross-ratio equation
    # over the integers, independently of the coprime parametrization.
    answer = 0
    lam = -h * h
    for a in range(1, n + 1):
        for b in range(1, n + 1):
            if a == b:
                continue
            s = b - a
            for c in range(1, n + 1):
                if c in (a, b):
                    continue
                u = c - a
                denominator = (lam - 1) * u + s
                numerator = lam * u * s
                if not denominator or numerator % denominator:
                    continue
                d = a + numerator // denominator
                if 1 <= d <= n and d not in (a, b, c):
                    assert (c - b) * (d - a) == lam * (c - a) * (d - b)
                    answer += 1
    return answer


def main():
    transform_checks = symmetry_checks = directional_checks = 0
    fields = []
    rng = Random(309042026)
    for p in (5, 7, 13, 17, 29, 41):
        chi = character(p)
        e = [sum(chi[x * (x - 1) * (x - lam) % p] for x in range(p))
             for lam in range(p)]
        f = [chi[x * (x - 1) % p] for x in range(p)]
        assert sum(e) == 0
        assert sum(v * v for v in e) == p * p - 2 * p - 1
        tau = chi[:]
        tau_bar = [tau[-x % p] for x in range(p)]
        equal_at_root(conv(tau, tau_bar), scalar_poly(p, p))
        inv2 = pow(2, -1, p)
        inv16 = pow(16, -1, p)
        kl = {}
        kl_second = [0] * p
        for t in range(1, p):
            a = t * t * inv16 % p
            counts = [0] * p
            for u in range(1, p):
                counts[(u + a * pow(u, -1, p)) % p] += 1
            kl[t] = counts
            equal_at_root(fourier_poly(f, t), shift(counts, -t * inv2))
            rhs = [chi[t] * v for v in shift(conv(tau, counts), -t * inv2)]
            equal_at_root(fourier_poly(e, t), rhs)
            square = conv(counts, counts)
            kl_second = [x + y for x, y in zip(kl_second, square)]
            transform_checks += 1
        equal_at_root(kl_second, scalar_poly(p, p * p - 2 * p - 1))

        subsets = [(), (0,), tuple(range(min(p, 4)))]
        subsets += [tuple(sorted(rng.sample(range(p), rng.randrange(2, min(p, 9) + 1))))
                    for _ in range(6)]
        denominator = (p - 2) * (p - 3)
        for subset in subsets:
            n = len(subset)
            v = vector(p, chi, subset)
            falling = n * (n - 1) * (n - 2) * (n - 3)
            # Clear the exact centering denominator before polynomial tests.
            w = [0, 0] + [denominator * v[lam] - falling * e[lam]
                           for lam in range(2, p)]
            for lam in range(2, p):
                other = (1 - lam) % p
                inverse = pow(lam, -1, p)
                assert e[other] == chi[-1 % p] * e[lam]
                assert v[other] == chi[-1 % p] * v[lam]
                assert w[other] * e[other] == w[lam] * e[lam]
                assert e[inverse] == chi[lam] * e[lam]
                assert v[inverse] == chi[lam] * v[lam]
                assert w[inverse] * e[inverse] == w[lam] * e[lam]
                symmetry_checks += 1
            u, a = 2, 3
            transformed = tuple((u * x + a) % p for x in subset)
            assert vector(p, chi, transformed) == v
            pairing = sum(x * y for x, y in zip(w, e))
            inside = [0] * p
            for t in range(1, p):
                term = shift(conv(kl[t], fourier_poly(w, t)), t * inv2)
                inside = [x + chi[t] * y for x, y in zip(inside, term)]
            equal_at_root(conv(tau_bar, inside), scalar_poly(p, p * pairing))
            directional_checks += 1
        fields.append({'p': p, 'full_elliptic_norm_squared': p * p - 2 * p - 1,
                       'subsets': len(subsets)})

    configurations = []
    for n in (32, 64, 128, 1024, 4096):
        for h in (1, 2, 3, 4):
            lower, offset_count = pattern_count(n, h)
            record = {'n': n, 'h': h, 'rational_cross_ratio': -h * h,
                      'distinct_offset_patterns': offset_count,
                      'constructed_quadruple_lower_bound': lower}
            if n <= 128:
                actual = literal_integer_count(n, h)
                assert actual >= lower
                record['independent_full_integer_count'] = actual
            configurations.append(record)

    report = {
        'claim': 'Exact Fourier/Kloosterman identities and integer configurations; no full Paley or LM proof',
        'exact_nonzero_frequency_transform_checks': transform_checks,
        'exact_signed_symmetry_checks': symmetry_checks,
        'exact_directional_pairing_and_affine_invariance_checks': directional_checks,
        'fields': fields,
        'integer_configuration_records': configurations,
        'source_sha256': {path: sha256((ROOT / path).read_bytes()).hexdigest()
                          for path in ['research/parallel3-classical-2026-09-04.md',
                                       'experiments/parallel3_classical_2026_09_04.py']},
    }
    (ROOT / 'results/parallel3_classical_2026_09_04.json').write_text(json.dumps(report, indent=2) + '\n')
    print(json.dumps(report, indent=2))


if __name__ == '__main__':
    main()
