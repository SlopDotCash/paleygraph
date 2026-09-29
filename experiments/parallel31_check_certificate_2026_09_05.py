#!/usr/bin/env python3
"""Read-only independent check of the finite short-multiple certificate.

Uses the row equations for multiplication by 1+X+X^19, rather than the
classifier's norm descent or generic multiplication routine.
"""
from collections import Counter
from hashlib import sha256
from itertools import combinations
from math import isqrt, prod
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()


def check(ok, label):
    assert ok, label
    COUNTS[label] += 1


def row_product(q):
    assert len(q) == 64
    return [sum(q[(k-e) % 64]*(1 if k >= e else -1) for e in [0, 1, 19])
            for k in range(64)]


def main():
    path = ROOT/'results/parallel31_short_multiples_2026_09_05.json'
    data = json.loads(path.read_text())
    p, n, g = data['p'], data['n'], data['generator']
    check((p, n, g) == (215535361, 128, 25525303), 'fixed_parameters')
    check(all(p % i for i in range(2, isqrt(p)+1)), 'independent_trial_primality')
    check(pow(g, 64, p) == p-1 and pow(g, 128, p) == 1, 'subgroup_order')
    check((1+g+pow(g, 19, p)) % p == 0, 'triangle_equation')
    h = {pow(g, e, p) for e in range(n)}
    check(len(h) == n and 2 not in h and sum((1-x) % p in h for x in h) == 6,
          'one_distinct_triangle_orbit')
    check(row_product(data['adjugate']) == [p]+[0]*63,
          'independent_integer_inverse_certificate')
    seen = set()
    profiles = Counter()
    exemplars = {}
    for row in data['orbits']:
        e = row['exponents']
        check(len(e) == 6 and len(set(e)) == 6 and all(0 <= x < n for x in e), 'six_distinct_exponents')
        check(not any((x+64) % 128 in e for x in e), 'opposite_free_exponents')
        canonical = min(tuple(sorted((x-y) % 128 for x in e)) for y in e)
        check(tuple(e) == canonical and canonical not in seen, 'distinct_canonical_orbit')
        seen.add(canonical)
        word = [pow(g, x, p) for x in e]
        check(sum(word) % p == 0, 'direct_zero_sum')
        total_product = prod(word) % p
        # Product balance for I|Ic is equivalent to product(I)^2=-product(S).
        check(all((prod(t)**2+total_product) % p != 0 for t in combinations(word, 3)),
              'independent_all_product_partitions_unbalanced')
        zero_subsets = [size for size in range(1, 6)
                        for t in combinations(word, size) if sum(t) % p == 0]
        kind = 'primitive' if not zero_subsets else 'triangular'
        check(kind == row['kind'] and (not zero_subsets or zero_subsets == [3, 3]),
              'all_proper_subsets_checked')
        target = [int(k in e)-int(k+64 in e) for k in range(64)]
        q = row['quotient']
        check(all(type(x) is int for x in q) and row_product(q) == target,
              'independent_integral_quotient_certificate')
        length = sum(abs(x) for x in q)
        check(length == row['quotient_l1'] and sum(x*x for x in q) == row['quotient_l2_squared'],
              'minimum_length_and_norm_arithmetic')
        profiles[(kind, length)] += 1
        exemplars.setdefault((kind, length), e)
        if kind == 'primitive':
            check(length % 2 == 0 and row['cancellation_graph']['cycle_rank'] == length//2-2,
                  'minimum_cycle_rank_arithmetic')
    expected = {('triangular', 2): 54, ('primitive', 4): 20,
                ('primitive', 38): 1, ('primitive', 40): 9,
                ('primitive', 42): 17, ('primitive', 44): 18}
    check(dict(profiles) == expected, 'complete_minimum_length_profile')
    prior = json.loads((ROOT/'results/parallel30_verification_2026_09_05.json').read_text())
    old = prior['independent_direct_counts'][-1]
    check(len(seen)*720*n == old['distinct_unbalanced_R6'], 'prior_independent_total_count')
    check(profiles[('triangular', 2)]*720*n == old['triangle_D6'], 'prior_independent_triangle_count')
    sources = ['experiments/parallel31_check_certificate_2026_09_05.py',
               'results/parallel31_short_multiples_2026_09_05.json',
               'results/parallel30_verification_2026_09_05.json']
    report = {'status': 'independent finite certificate checks passed; full goal unproved',
              'check_counts': dict(COUNTS),
              'profile': [{'kind': k, 'minimum_length': l, 'orbits': c, 'example': exemplars[(k, l)]}
                          for (k, l), c in sorted(profiles.items())],
              'certificate_sha256': sha256(path.read_bytes()).hexdigest(),
              'input_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in sources},
              'scope': 'Exact certificates in one field; ordinary minimum-length proof separate; not Lean verified.'}
    (ROOT/'results/parallel31_certificate_check_2026_09_05.json').write_text(json.dumps(report, indent=2)+'\n')
    print(json.dumps({'checks': sum(COUNTS.values()), 'orbits': len(seen), 'profile': report['profile']}, indent=2))


if __name__ == '__main__':
    main()
