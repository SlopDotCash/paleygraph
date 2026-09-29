#!/usr/bin/env python3
"""Exact anchored-row and arithmetic tests; no general proof is asserted."""
from collections import Counter
from hashlib import sha256
from math import comb, gcd
from pathlib import Path
import json

from parallel46_row_capacity import threshold_check
from mixed_period_collisions import prime, primitive_root

ROOT = Path(__file__).resolve().parents[1]


def prime_pair():
    n, m, p = 1024, 1094909953, 1121187791873
    c, s, witness = 16707, 16, 5
    assert m == c * 2 ** s + 1 and c % 2 and c < 2 ** s
    assert pow(witness, (m - 1) // 2, m) == m - 1
    # This forces every prime divisor of m to be 1 mod 2^s,
    # hence larger than sqrt(m); m is prime.
    assert p == n * m + 1 and m * m > p
    assert pow(2, p - 1, p) == 1 and gcd(pow(2, n, p) - 1, p) == 1
    # Every prime divisor of p then has multiplicative order divisible
    # by the certified prime m, and is larger than sqrt(p).
    assert n ** 4 < p < 2 * n ** 4
    return dict(n=n, m=m, p=p, m_proth_c=c, m_proth_s=s,
                m_witness=witness, p_witness=2, exact_prime_certificates_passed=True)


def check_anchor(a, m, n):
    corr = Counter()
    for x, vx in a.items():
        for y, vy in a.items():
            corr[(x - y) % m] += vx * vy
    h = dict(corr)
    h[0] -= n - 1
    a0 = a.get(0, 0)
    assert h[0] >= a0 * (a0 - 1)
    for i in set(a) | {(-x) % m for x in a}:
        if i:
            x, y = a.get(i, 0), a.get((-i) % m, 0)
            assert x + y <= n
            assert h.get(i, 0) >= x * (x - 1) + y * (y - 1)
    return h, sum(v * v for v in corr.values())


def centered_model(t, m, R):
    n, A, nu = t ** 5, t ** 3, (t ** 5).bit_length() - 1
    ell = (n - 2 - 2 * A * t) // 2
    I = {i % m for i in range(1 - t, t + 1)}
    gamma = ((A * t + 2 * ell * R + ell * (ell - 1)) * pow(nu - 1, -1, m)) % m
    if gamma in I or R <= gamma < R + ell:
        return None
    a = {i: A for i in I} | {j: 2 for j in range(R, R + ell)} | {gamma: 1}
    assert len(a) == 2 * t + ell + 1 and a[0] == A
    assert sum(a.values()) == n - 1 and sum(v % 2 for v in a.values()) == 1
    assert sum(i * v for i, v in a.items()) % m == nu * gamma % m
    return a, gamma, ell


def moments_actual(cert):
    n, m, p = cert['n'], cert['m'], cert['p']
    base = next(b for b in range(2, 100) if pow(pow(b, m, p), n // 2, p) == p - 1)
    generator = pow(base, m, p)
    H = {pow(generator, i, p) for i in range(n)}
    assert len(H) == n and pow(generator, n, p) == 1
    a = Counter(pow(h - 1, n, p) for h in H if h != 1)
    special = pow(2, n, p)
    assert special != 1 and pow(special, m, p) == 1
    assert all(v % 2 == (x == special) for x, v in a.items())
    product = 1
    for x, v in a.items():
        product = product * pow(x, v, p) % p
    assert product == pow(special, n.bit_length() - 1, p)
    moment_records = []
    for k in [1, 2, 3]:
        lhs = sum(v * pow(x, k, p) for x, v in a.items()) % p
        rhs = n * sum(comb(n * k, n * j) for j in range(k + 1)) % p
        assert lhs == rhs
        moment_records.append(dict(k=k, residue=lhs))
    a0 = a.get(1, 0)
    A2 = sum(v * v for v in a.values())
    assert A2 - n + 1 >= a0 * (a0 - 1)
    tested = 0
    for u in set(a) | {pow(x, -1, p) for x in a}:
        if u != 1:
            x, y = a.get(u, 0), a.get(pow(u, -1, p), 0)
            r = sum(v * a.get(z * u % p, 0) for z, v in a.items())
            assert r >= x * (x - 1) + y * (y - 1)
            tested += 1
    return dict(generator_base=base, subgroup_generator=generator, support=len(a),
                a0=a0, A2=A2, max_a=max(a.values()), anchored_indices_checked=tested,
                product_residue=product, moments=moment_records,
                counts_sha256=sha256(json.dumps(sorted(a.items())).encode()).hexdigest(),
                all_passed=True)


def geometric(q, start, length, p, m):
    assert q != 1
    return pow(q, start % m, p) * (pow(q, length, p) - 1) * pow(q - 1, -1, p) % p


def newton_case(p, n):
    assert prime(p) and n % 2 == 0 and (p - 1) % n == 0
    g = pow(primitive_root(p), (p - 1) // n, p)
    roots = [pow(pow(g, j, p) - 1, n, p) for j in range(1, n)]
    sums = [n * sum(comb(n * k, n * j) for j in range(k + 1)) % p
            for k in range(1, n)]
    assert all(sum(pow(x, k, p) for x in roots) % p == sums[k - 1] for k in range(1, n))
    elementary = [1]
    for r in range(1, n):
        value = sum((-1) ** (j - 1) * elementary[r - j] * sums[j - 1]
                    for j in range(1, r + 1))
        elementary.append(value * pow(r, -1, p) % p)
    recovered = [((-1) ** r * elementary[r]) % p for r in range(n)][::-1]
    direct = [1]
    for x in roots:
        following = [0] * (len(direct) + 1)
        for j, c in enumerate(direct):
            following[j] = (following[j] - x * c) % p
            following[j + 1] = (following[j + 1] + c) % p
        direct = following
    assert recovered == direct
    return dict(p=p, n=n, moments_checked=n-1, polynomial_coefficients_ascending=recovered,
                direct_factor_product_matches=True)


def main():
    cert = prime_pair()
    n, m, p, t = cert['n'], cert['m'], cert['p'], 4
    A, nu = t ** 3, n.bit_length() - 1
    special = pow(2, n, p)
    residues, first = [], None
    skipped = 0
    for R in range(10 * n + 1, 12 * n + 1):
        model = centered_model(t, m, R)
        if model is None:
            skipped += 1
            continue
        a, gamma, ell = model
        if first is None:
            first = (R, a, gamma, ell)
        q = pow(special, pow(gamma, -1, m), p)
        assert pow(q, gamma, p) == special and q != 1
        value = (A * geometric(q, 1 - t, 2 * t, p, m)
                 + 2 * geometric(q, R, ell, p, m) + special) % p
        if len(residues) < 3:
            assert value == sum(v * pow(q, i, p) for i, v in a.items()) % p
        residues.append([R, gamma, value])
    assert first is not None
    R, a, gamma, ell = first
    h, K = check_anchor(a, m, n)
    A2 = sum(v * v for v in a.values())
    X = 384 * n * n
    assert K <= 315 * n ** 3 and A2 <= 4 * t ** 7
    assert max(v ** 3 * count for v, count in Counter(a.values()).items()) == 2 * n * n
    assert K <= n * (2 * (n - 1) ** 2 - (n - 1) + X)
    assert X >= 3 * sum(v * (v - 1) * (v - 2) for v in a.values()) and X % 6 == 0
    thresholds = threshold_check(h, n, X)
    width = t // 2
    triples = [sum(v ** 2 * a.get((i + u) % m, 0) ** 2 * a.get((i + w) % m, 0) ** 2
                   for i, v in a.items())
               for u in range(-width, width + 1) for w in range(-width, width + 1)]
    norm_cube_lower = min(triples) ** 3 * len(triples) ** 2
    assert norm_cube_lower >= t ** 61
    actual = moments_actual(cert)
    newton = [newton_case(p0, n0) for p0, n0 in [(97, 8), (353, 16), (278177, 32)]]
    survivors = [r for r in residues if r[2] == 2 * n % p]
    assert not survivors
    old_checks = []
    for s in [2, 4, 8, 16]:
        old_n, old_A, old_d = s ** 5, s ** 3, s
        old_l = (old_n - 2 - old_A * old_d) // 2
        old_row = 4 * (old_l - old_d) + old_A
        assert old_row < old_A * (old_A - 1)
        old_checks.append(dict(t=s, required=old_A * (old_A - 1), actual_model_row_mass=old_row))
    result = dict(scope='Necessary anchored-row and cyclotomic arithmetic constraints; repaired abstract model and one finite prime-field sieve, not an improved uniform estimate.',
        prime_pair=cert, old_interval_anchor_failures=old_checks,
        repaired_model=dict(t=t, n=n, m=m, R=R, odd_label=gamma, filler_length=ell,
            a0=A, A2=A2, K=K, Q=2*n*n, formal_X_budget=X,
            anchored_rows_passed=True, formal_product_congruence_passed=True,
            total_X_threshold_check=thresholds, norm_cube_lower_bound=norm_cube_lower,
            actual_incidence_matrix_supplied=False),
        actual_subgroup=actual,
        newton_reconstruction_checks=newton,
        finite_sieve=dict(R_range=[10*n+1,12*n], candidate_positions=2*n,
            skipped_overlapping_odd_label=skipped, admissible_candidates=len(residues),
            first_power_moment_survivors=len(survivors),
            residues_sha256=sha256(json.dumps(residues).encode()).hexdigest(), first_samples=residues[:3],
            scope='Exhaustive only over the stated R range at this one n,p; no uniform exclusion.'),
        all_passed=True)
    out = ROOT/'results/parallel47_anchored_arithmetic_2026_09_06.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'n':n,'p':p,'admissible_candidates':len(residues),
        'moment_survivors':len(survivors),'actual_subgroup_max_a':actual['max_a'],'output':str(out)}))


if __name__ == '__main__':
    main()
