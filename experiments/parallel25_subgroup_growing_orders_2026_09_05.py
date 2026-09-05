#!/usr/bin/env python3
"""Repeated-coordinate residual upper bound and actual cyclic norm witnesses.

Finite identities are checked exactly; imported energy bounds and the
unresolved six-distinct unbalanced aggregate are not proved by this run.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations, combinations_with_replacement
from math import comb, gcd, isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(ok, name):
    assert ok, name
    COUNTS[name] += 1


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def subgroup(p, n):
    check(prime(p) and (p - 1) % n == 0
          and n ** 4 // 4 <= p <= n ** 4, 'actual_quartic_prime')
    for a in range(2, p):
        g = pow(a, (p - 1) // n, p)
        if pow(g, n // 2, p) != 1:
            H = sorted(pow(g, j, p) for j in range(n))
            check(len(set(H)) == n and pow(g, n, p) == 1,
                  'actual_exact_dyadic_order')
            return g, H
    raise AssertionError('No generator')


def balanced(word, p):
    for tail in combinations(range(1, 6), 2):
        I = {0, *tail}
        a = b = 1
        for j, h in enumerate(word):
            if j in I:
                a = a * h % p
            else:
                b = b * h % p
        if (a + b) % p == 0:
            return True
    return False


def norm_recursive(a):
    a = list(a)
    while len(a) > 1:
        d = len(a) // 2
        b = [0] * d
        for part, shift, sign in [(a[::2], 0, 1), (a[1::2], 1, -1)]:
            nonzero = [(j, x) for j, x in enumerate(part) if x]
            for i, x in nonzero:
                for j, y in nonzero:
                    q, r = divmod(i + j + shift, d)
                    b[r] += sign * x * y * (-1 if q % 2 else 1)
        a = b
    return abs(a[0])


def modular_rank(M, p):
    M = [[x % p for x in row] for row in M]
    rank = 0
    for col in range(len(M[0])):
        pivot = next((i for i in range(rank, len(M)) if M[i][col]), None)
        if pivot is None:
            continue
        M[rank], M[pivot] = M[pivot], M[rank]
        inv = pow(M[rank][col], -1, p)
        for j in range(col, len(M[0])):
            M[rank][j] = M[rank][j] * inv % p
        for i in range(rank + 1, len(M)):
            value = M[i][col]
            if value:
                for j in range(col, len(M[0])):
                    M[i][j] = (M[i][j] - value * M[rank][j]) % p
        rank += 1
        if rank == len(M):
            break
    return rank


def cyclic_witness(p, n, g, H, word):
    Hset = set(H)
    check(all(h in Hset for h in word) and sum(word) % p == 0
          and all(-h % p not in word for h in word)
          and not balanced(word, p), 'actual_opposite_free_unbalanced_norm_witness')
    exponents = {pow(g, j, p): j for j in range(n)}
    e = [exponents[h] for h in word]
    d = n // 2
    f = [e.count(j) - e.count(j + d) for j in range(d)]
    content = 0
    for value in f:
        content = gcd(content, value)
    f = [value // content for value in f]
    Q = sum(x * x for x in f)
    norm = norm_recursive(f)
    check(norm > 0 and norm * norm <= Q ** d, 'exact_norm_and_AM_GM')
    valuation = 0
    remainder = norm
    while remainder % p == 0:
        valuation += 1
        remainder //= p
    zeros = [u for u in range(1, n, 2)
             if sum(value * pow(g, u * j, p) for j, value in enumerate(f)) % p == 0]
    M = [[0] * d for _ in range(d)]
    for j in range(d):
        for k, value in enumerate(f):
            q, r = divmod(k + j, d)
            M[r][j] += value * (-1 if q % 2 else 1)
    rank = modular_rank(M, p)
    check(rank == d - len(zeros), 'cyclic_rank_equals_primitive_root_defect')
    check(valuation >= len(zeros) and p ** (2 * len(zeros)) <= Q ** d,
          'root_count_norm_divisibility_bound')
    for col in range(d):
        check(sum(M[row][col] * pow(g, row, p) for row in range(d)) % p == 0,
              'all_cyclic_shifts_vanish_at_generator')
    check(valuation == 1 and zeros == [1] and rank == d - 1,
          'many_rational_shifts_only_one_prime_norm_factor')
    # A nonzero norm certifies full rational rank of this multiplication matrix.
    return {'word': word, 'generator': g, 'exponents': e,
            'primitive_coefficients': f, 'Q': Q, 'norm': norm,
            'p_adic_valuation': valuation, 'norm_after_removing_p': remainder,
            'primitive_zero_exponents': zeros, 'rational_rank': d,
            'rank_mod_p': rank}


def inspect(p, n, norm_word):
    g, H = subgroup(p, n)
    check(all(a * b % p in H for a in H for b in H), 'multiplicative_closure')
    r2 = Counter((a + b) % p for a in H for b in H)
    E2 = sum(c * c for c in r2.values())
    r4_at_2 = sum(c * r2[(2 - s) % p] for s, c in r2.items())
    check(r4_at_2 == sum(c * r2[(-2 - s) % p] for s, c in r2.items()),
          'symmetric_four_sum_target')
    check(6 * n - 8 <= r4_at_2 <= E2, 'four_sum_intrinsic_baseline_and_Cauchy')
    pair_buckets = defaultdict(list)
    for a, b in combinations_with_replacement(H, 2):
        pair_buckets[(1 + a + b) % p].append((a, b, 1 if a == b else 2))
    counts = Counter()
    patterns = Counter()
    fibers = Counter()
    rep_fibers = Counter()
    distinct_unbalanced_fibers = Counter()
    for tail in combinations_with_replacement(H, 3):
        candidates = pair_buckets.get(-sum(tail) % p, ())
        if not candidates:
            continue
        wt = 1 if tail[0] == tail[2] else 3 if tail[0] == tail[1] or tail[1] == tail[2] else 6
        inv = pow(tail[0] * tail[1] * tail[2] % p, -1, p)
        for a, b, wp in candidates:
            word = (1, a, b, *tail)
            mass = n * wp * wt
            c = Counter(word)
            equal_pairs = sum(comb(v, 2) for v in c.values())
            intrinsic = all(c[h] == c[-h % p] for h in c)
            free = all(-h % p not in c for h in c)
            rho = -a * b * inv % p
            fibers[rho] += mass
            counts['E3'] += mass
            counts['all_equal_pair_weight'] += mass * equal_pairs
            if intrinsic:
                counts['intrinsic'] += mass
                counts['intrinsic_equal_pair_weight'] += mass * equal_pairs
            else:
                counts['nonintrinsic_equal_pair_weight'] += mass * equal_pairs
                if not free:
                    counts['opposite_containing_nonintrinsic'] += mass
            if free:
                counts['R6'] += mass
                has_balance = balanced(word, p)
                if has_balance:
                    counts['balanced_R6'] += mass
                if len(c) < 6:
                    counts['repeated_R6'] += mass
                    counts['repeated_R6_equal_pair_weight'] += mass * equal_pairs
                    rep_fibers[rho] += mass
                    patterns[','.join(map(str, sorted(c.values(), reverse=True)))] += mass
                elif has_balance:
                    counts['distinct_balanced_R6'] += mass
                else:
                    counts['distinct_unbalanced_R6'] += mass
                    distinct_unbalanced_fibers[rho] += mass
            COUNTS['normalized_six_multiset_matches'] += 1
    T4 = 3 * n * n - 3 * n
    T6 = 15 * n ** 3 - 45 * n * n + 40 * n
    check(counts['intrinsic'] == T6, 'exact_intrinsic_sixth_count')
    check(counts['all_equal_pair_weight'] == 15 * n * r4_at_2,
          'exact_equal_coordinate_pair_incidence')
    check(counts['intrinsic_equal_pair_weight'] == 15 * n * (6 * n - 8),
          'exact_intrinsic_equal_pair_subtraction')
    rep_upper = 15 * n * (r4_at_2 - 6 * n + 8)
    check(counts['nonintrinsic_equal_pair_weight'] == rep_upper,
          'exact_nonintrinsic_repetition_weight')
    check(counts['repeated_R6'] <= counts['repeated_R6_equal_pair_weight'] <= rep_upper,
          'aggregate_repeated_residual_upper')
    check(counts['opposite_containing_nonintrinsic'] <= 15 * n * (E2 - T4),
          'opposite_containing_upper_from_fourth_energy')
    w = Counter((a * b - a - b) % p for a in H for b in H)
    shifted_energy = sum(v * v for v in w.values())
    X = shifted_energy - (6 * n * n - 9 * n + 4)
    check(X >= 0 and counts['balanced_R6'] <= 10 * n * X,
          'existing_balanced_upper_retained')
    for rho in H:
        inv = pow(rho, -1, p)
        N = sum(c * w[z * inv % p] for z, c in w.items())
        check(fibers[rho] == n * N, 'independent_actual_product_ratio_fibers')
    D = counts['distinct_unbalanced_R6']
    error = counts['E3'] - T6 - D
    check(error == counts['opposite_containing_nonintrinsic'] + counts['repeated_R6']
          + counts['distinct_balanced_R6'], 'exact_six_distinct_remainder_partition')
    upper = 15 * n * (E2 - T4) + rep_upper + 10 * n * X
    check(0 <= error <= upper, 'uniform_structural_remainder_upper')
    return {'p': p, 'n': n, 'E2': E2, 'r4_at_2': r4_at_2, 'X': X,
            'counts': dict(sorted(counts.items())), 'repeated_patterns': dict(sorted(patterns.items())),
            'repeated_upper_after_intrinsic_subtraction': rep_upper,
            'six_distinct_remainder_error': error, 'six_distinct_remainder_upper': upper,
            'repeated_residual_fibers': dict(sorted(rep_fibers.items())),
            'distinct_unbalanced_residual_fibers': dict(sorted(distinct_unbalanced_fibers.items())),
            'cyclic_norm_witness': cyclic_witness(p, n, g, H, norm_word)}


def main():
    cases = [inspect(p, n, word) for p, n, word in [
        (6700417, 64, [1, 4, 6700409, 1, 1, 1]),
        (7204033, 64, [1, 3202259, 3519792, 1, 140968, 341012]),
        (67403009, 128, [1, 37818017, 61889035, 710836, 9585000, 24803129]),
        (1073748737, 256, [1, 914267366, 972187974, 9468345, 10646993, 240926795]),
    ]]
    inputs = ['experiments/parallel25_subgroup_growing_orders_2026_09_05.py',
              'research/parallel24-subgroup-unbalanced-2026-09-05.md',
              'research/parallel23-subgroup-upper-2026-09-05.md',
              'research/parallel21-subgroup-next-input-2026-09-05.md',
              'results/parallel23_subgroup_upper_2026_09_05.json',
              'sources/sigma-subgroup-2026-09-05/mrss-1712.00410.pdf',
              'sources/sigma-subgroup-2026-09-05/mrss-1712.00410.txt',
              'sources/mixed-periods-2026-09-04/shkredov-1504.04522.html']
    result = {'status': 'PASS',
              'scope': 'Uniform repeated-coordinate residual bound plus cyclic norm-rank restriction; no uniform total sixth-energy bound.',
              'counts': dict(sorted(COUNTS.items())), 'cases': cases,
              'input_sha256': {p: sha256((ROOT / p).read_bytes()).hexdigest() for p in inputs}}
    path = ROOT / 'results/parallel25_subgroup_growing_orders_2026_09_05.json'
    path.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps({'status': 'PASS', 'counts': result['counts'], 'cases': [
        {k: c[k] for k in ['p', 'n', 'E2', 'r4_at_2', 'counts', 'repeated_patterns',
                           'repeated_upper_after_intrinsic_subtraction',
                           'six_distinct_remainder_error', 'six_distinct_remainder_upper']}
        for c in cases]}, indent=2))


if __name__ == '__main__':
    main()
