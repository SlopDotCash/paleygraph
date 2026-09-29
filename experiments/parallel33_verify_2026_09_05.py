#!/usr/bin/env python3
"""Exact finite checks for centered distinct-coordinate moment reduction."""
from collections import Counter
from fractions import Fraction
from functools import lru_cache
from hashlib import sha256
from math import factorial, isqrt, prod
from pathlib import Path
import json
import random
import time

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()


def check(value, label):
    assert value, label
    CHECKS[label] += 1


def rational(x):
    x = Fraction(x)
    return {'numerator': x.numerator, 'denominator': x.denominator}


def partitions(total, least=1):
    if total == 0:
        yield ()
    for first in range(least, total+1):
        for rest in partitions(total-first, first):
            yield (first,)+rest


def cycle_coefficient(parts):
    r, k = sum(parts), len(parts)
    denominator = prod(parts)*prod(factorial(v) for v in Counter(parts).values())
    assert factorial(r) % denominator == 0
    return (-1 if (r-k) % 2 else 1)*(factorial(r)//denominator)


def falling(n, r):
    return prod(n-j for j in range(r)) if r <= n else 0


def weighted_convolutions(a, p):
    @lru_cache(None)
    def convolution(parts):
        if not parts:
            return (1,)+(0,)*(p-1)
        previous = convolution(parts[:-1])
        out = [0]*p
        shifts = [parts[-1]*x % p for x in a]
        for i, count in enumerate(previous):
            if count:
                for shift in shifts:
                    out[(i+shift) % p] += count
        return tuple(out)
    return convolution


def subset_zero_counts(a, p, maximum):
    # Coefficients of product_(a in A)(1+z X^a), not set-partition inversion.
    dp = [[0]*p for _ in range(maximum+1)]
    dp[0][0] = 1
    for number, value in enumerate(a, start=1):
        for k in range(min(number, maximum), 0, -1):
            for x, count in enumerate(dp[k-1]):
                if count:
                    dp[k][(x+value) % p] += count
    return [factorial(k)*dp[k][0] for k in range(maximum+1)]


def subgroup(p, n):
    assert (p-1) % n == 0
    g = next(pow(x, (p-1)//n, p) for x in range(2, p)
             if pow(pow(x, (p-1)//n, p), n//2, p) != 1)
    return sorted(pow(g, j, p) for j in range(n))


def main():
    started = time.monotonic()
    cases = [(p, subgroup(p, n), True) for p, n in
             [(17, 4), (17, 8), (41, 8), (73, 8), (97, 16),
              (193, 16), (193, 32), (257, 16), (257, 32)]]
    cases += [(101, sorted({x % 101 for x in range(-6, 7)}), False)]
    rng = random.Random(330905)
    positive = rng.sample(range(1, 129), 10)
    cases += [(257, sorted({x % 257 for t in positive for x in [t, -t]}), False)]
    rows = []
    for p, a, is_subgroup in cases:
        n = len(a)
        check(all(p % d for d in range(2, isqrt(p)+1)), 'field_primality')
        check(len(set(a)) == n < p and {-x % p for x in a} == set(a), 'proper_symmetric_set')
        if is_subgroup:
            check(all(x*y % p in a for x in a for y in a), 'multiplicative_closure')
        conv = weighted_convolutions(a, p)
        injective = subset_zero_counts(a, p, 8)
        alpha = Fraction(n*(p-n), p-1)
        lower_sqrt = isqrt(alpha.numerator//alpha.denominator)
        assert lower_sqrt >= 1
        for r in [4, 6, 8]:
            assert p > r
            s = r//2
            energy = conv((1,)*r)[0]
            centered = Fraction(energy)-Fraction(n**r, p)
            squarefree_centered = Fraction(injective[r])-Fraction(falling(n, r), p)
            check(centered > 0, 'positive_nonprincipal_even_moment')
            check(Fraction(p, p-1)*centered >= alpha**s, 'nonprincipal_jensen_lower_bound')
            inversion = 0
            principal_inversion = 0
            block_rows = []
            stirling = Counter()
            for parts in partitions(r):
                coefficient = cycle_coefficient(parts)
                k = len(parts)
                count = conv(parts)[0]
                mixed_centered = Fraction(count)-Fraction(n**k, p)
                inversion += coefficient*count
                principal_inversion += coefficient*n**k
                stirling[k] += abs(coefficient)
                check(mixed_centered*mixed_centered*alpha**(r-k) <= centered*centered,
                      'centered_mixed_partition_moment_bound')
                block_rows.append({'parts': parts, 'coefficient': coefficient,
                                   'zero_count': count, 'centered': rational(mixed_centered)})
            check(inversion == injective[r], 'mobius_matches_independent_subset_dp')
            check(principal_inversion == falling(n, r), 'exact_falling_principal_term')
            # Verify Stirling coefficients independently using cycle insertion.
            coefficients = [1]
            for j in range(r):
                new = [0]*(len(coefficients)+1)
                for k, value in enumerate(coefficients):
                    new[k] += j*value
                    new[k+1] += value
                coefficients = new
            check(all(stirling[k] == coefficients[k] for k in range(r+1)),
                  'unsigned_stirling_cycle_identity')
            epsilon = prod(Fraction(lower_sqrt+j, lower_sqrt) for j in range(1, r))-1
            check(abs(centered-squarefree_centered) <= epsilon*centered,
                  'rational_centered_distinct_coordinate_bound')
            repeated = energy-injective[r]
            if is_subgroup:
                marked = conv((1,)*(r-2))[(-2) % p]*n
                check(marked == conv((1,)*(r-2)+(2,))[0], 'normalized_equal_pair_count')
                check(repeated <= (r*(r-1)//2)*marked, 'repeated_coordinate_union_bound')
                # Eliminate fractional exponents by raising to r-2.
                choose_two = r*(r-1)//2
                check(repeated**(r-2) <= choose_two**(r-2)*n*energy**(r-3),
                      'uncentered_relative_collision_bound')
                if r == 6:
                    e2 = conv((1,)*4)[0]
                    check(repeated*repeated <= 225*e2*energy, 'sharper_six_word_collision_bound')
            rows.append({'p': p, 'n': n, 'subgroup': is_subgroup, 'r': r,
                         'elements': a, 'alpha': rational(alpha),
                         'energy': energy, 'distinct_zero_count': injective[r],
                         'centered_energy': rational(centered),
                         'centered_distinct_count': rational(squarefree_centered),
                         'relative_error': rational(abs(centered-squarefree_centered)/centered),
                         'rational_epsilon_upper': rational(epsilon), 'partition_rows': block_rows})
    # The coefficient-size/characteristic condition is substantive.
    bad_conv = weighted_convolutions([1, 2], 3)
    bad_centered = Fraction(bad_conv((3, 3))[0])-Fraction(2**2, 3)
    bad_t = Fraction(bad_conv((1,)*6)[0])-Fraction(2**6, 3)
    check(bad_centered > bad_t > 0, 'zero_block_coefficient_counterexample_outside_hypotheses')

    parameter_rows = []
    for exponent in [8, 16, 24, 28, 30, 32, 40]:
        n = 2**exponent
        s = 3*exponent
        r = 2*s
        lower = isqrt(n-1)
        epsilon = prod(Fraction(lower+j, lower) for j in range(1, r))-1
        parameter_rows.append({'log2_n': exponent, 'n': n, 's': s,
                               'sqrt_alpha_lower': lower,
                               'epsilon_upper': rational(epsilon),
                               'epsilon_upper_decimal': float(epsilon),
                               'two_sided_bound_nontrivial': epsilon < 1})
    check(Fraction(49, 20)/2+Fraction(4, 1)/2 == Fraction(129, 40), 'improved_repetition_exponent')
    check(Fraction(1, 5)/2+Fraction(1, 1)/2 == Fraction(3, 5), 'improved_repetition_log_exponent')
    result = {'status': 'exact centered collision identities and finite bounds passed; full goal unproved',
              'check_counts': dict(CHECKS), 'cases': rows,
              'logarithmic_depth_parameters': parameter_rows,
              'excluded_characteristic_example': {'p': 3, 'r': 6, 'parts': [3, 3],
                    'mixed_centered': rational(bad_centered), 'centered_even_moment': rational(bad_t)},
              'elapsed_seconds': time.monotonic()-started,
              'input_sha256': {'experiments/parallel33_verify_2026_09_05.py':
                               sha256(Path(__file__).read_bytes()).hexdigest()},
              'scope': 'Exact finite checks support an ordinary uniform proof; no upper bound on the remaining distinct-coordinate aggregate.'}
    (ROOT/'results/parallel33_verification_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'status': result['status'], 'checks': sum(CHECKS.values()),
                      'check_counts': dict(CHECKS), 'cases': len(rows),
                      'depth_parameters': [{k: row[k] for k in ['log2_n', 's', 'epsilon_upper_decimal',
                                            'two_sided_bound_nontrivial']} for row in parameter_rows],
                      'elapsed_seconds': result['elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
