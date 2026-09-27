#!/usr/bin/env python3
"""Bounded exact checks of actual difference levels after excluding equal edge cosets.

No asymptotic estimate is inferred from these examples.  All arithmetic is integral.
"""

from collections import Counter, defaultdict
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path


def check(p, generator, n):
    assert p > 2 and all(p % k for k in range(2, isqrt(p) + 1))
    H = sorted({pow(generator, i, p) for i in range(n)})
    assert len(H) == n and pow(generator, n, p) == 1 and p - 1 in H
    f = Counter((a - b) % p for a in H for b in H)
    assert f.pop(0) == n
    E = n * n + sum(v * v for v in f.values())
    d = Fraction(E, 16 * n * n)
    reps = {}
    cosets = {}

    def rep(x):
        x %= p
        assert x
        if x not in reps:
            elements = {x * h % p for h in H}
            r = min(elements)
            cosets[r] = elements
            reps.update((v, r) for v in elements)
        return reps[x]

    values = {rep(x): v for x, v in f.items()}
    assert all(all(f[x] == v for x in cosets[r]) for r, v in values.items())
    level = {}
    for r, v in values.items():
        i = 1
        while v > 2**i * d:
            i += 1
        assert 2 ** (i - 1) * d < v <= 2**i * d
        level[r] = i
    level_reps = defaultdict(list)
    for r, i in level.items():
        level_reps[i].append(r)

    inverse = {r: pow(r, -1, p) for r in values}
    fibers = defaultdict(list)
    for a in values:
        for b in values:
            if b == a:
                continue
            for c in values:
                if c == a or c == b:
                    continue
                u, v = rep(a * inverse[c]), rep(b * inverse[c])
                assert u != rep(1) and v != rep(1) and u != v
                fibers[((level[a], level[b], level[c]), u, v)].append((a, b, c))

    rho_cache = {}

    def rho(u, v):
        if (u, v) not in rho_cache:
            rho_cache[u, v] = sum((x + 1) % p in cosets[v] for x in cosets[u])
        return rho_cache[u, v]

    stats = defaultdict(Counter)
    witness = None
    same_level_witness = None
    weighted_W = 0
    unweighted_level_incidence = 0
    weighted_level_incidence = 0
    for (levels, u, v), triples in sorted(fibers.items()):
        m = len(triples)
        r = rho(u, v)
        stats[levels]['normalized_support'] += 1
        stats[levels]['normalized_mass'] += m
        stats[levels]['ordered_normalization_collisions'] += m * (m - 1)
        stats[levels]['direct_incidence'] += n * m * r
        stats[levels]['discard_m_incidence'] += n * r
        stats[levels]['max_m'] = max(stats[levels]['max_m'], m)
        weighted_level_incidence += n * m * r
        unweighted_level_incidence += n * r
        fiber_weight = n * r * sum(values[a]**2 * values[b]**2 * values[c]**2
                                    for a, b, c in triples)
        weighted_W += fiber_weight
        stats[levels]['weighted_triangle_sum'] += fiber_weight
        if m > 1 and r:
            candidate = {
                'levels': levels, 'u': u, 'v': v, 'multiplicity': m, 'rho': r,
                'coset_triples': triples,
                'f_values': [[values[a], values[b], values[c]] for a, b, c in triples],
                'normalized_solutions': [[x, (x + 1) % p] for x in sorted(cosets[u])
                                         if (x + 1) % p in cosets[v]],
            }
            if witness is None:
                witness = candidate
            if same_level_witness is None and len(set(levels)) == 1:
                same_level_witness = candidate

    # Independent field-level count for each selected witness's complete level triple.
    direct_checks = []
    collision_checks = []
    for chosen in [witness, same_level_witness]:
        if chosen is None or chosen['levels'] in [x['levels'] for x in direct_checks]:
            continue
        levels = chosen['levels']
        A, B, C = [{x for r in level_reps[i] for x in cosets[r]} for i in levels]
        direct = 0
        direct_weighted = 0
        for x in A:
            for y in B:
                z = (y - x) % p
                if z in C and len({rep(x), rep(y), rep(z)}) == 3:
                    direct += 1
                    direct_weighted += f[x]**2 * f[y]**2 * f[z]**2
        assert direct == stats[levels]['direct_incidence']
        assert direct_weighted == stats[levels]['weighted_triangle_sum']
        direct_checks.append({'levels': levels, 'direct_incidence': direct,
                              'direct_weighted_triangle_sum': direct_weighted,
                              **dict(stats[levels])})
        multiplicities = Counter(levels)
        quotient_levels = {i: set(level_reps[i]) for i in multiplicities}
        possible_q = set()
        for D in quotient_levels.values():
            possible_q.update(rep(a * inverse[b]) for a in D for b in D)
        collision_sum = 0
        nonzero_terms = []
        for q in sorted(possible_q - {rep(1)}):
            counts = {i: sum(rep(q * a) in D for a in D)
                      for i, D in quotient_levels.items()}
            term = 1
            for i, exponent in multiplicities.items():
                for j in range(exponent):
                    term *= counts[i] - j
            assert term >= 0
            collision_sum += term
            if term:
                nonzero_terms.append({'q': q, 'intersection_counts': counts, 'term': term})
        assert collision_sum == stats[levels]['ordered_normalization_collisions']
        collision_checks.append({'levels': levels, 'ordered_collisions': collision_sum,
                                 'nonzero_quotient_intersection_terms': nonzero_terms})

    return {
        'p': p, 'generator': generator, 'n': n, 'H': H,
        'E': E, 'd': str(d),
        'f_histogram_on_cosets': dict(sorted(Counter(values.values()).items())),
        'levels': {i: {'representatives': sorted(rs),
                       'f_values': [values[r] for r in sorted(rs)]}
                   for i, rs in sorted(level_reps.items())},
        'three_distinct_incidence_over_all_level_triples': weighted_level_incidence,
        'discard_m_incidence_over_all_level_triples': unweighted_level_incidence,
        'W_distinct_from_weighted_normalization': weighted_W,
        'first_positive_repeated_fiber': witness,
        'first_positive_repeated_same_level_fiber': same_level_witness,
        'independent_direct_checks': direct_checks,
        'independent_collision_identity_checks': collision_checks,
        'level_triple_stats': [{'levels': key, **dict(value)} for key, value in sorted(stats.items())],
    }


if __name__ == '__main__':
    cases = [(97, 64, 8), (1153, 75, 8), (6700417, 2, 64),
             (193, 64, 16), (257, 249, 16), (769, 712, 16), (65537, 64, 16),
             (215535361, 25525303, 128)]
    rows = [check(*case) for case in cases]
    previous = json.loads(Path('results/parallel39_edge_cosets_2026_09_06.json').read_text())
    for row in rows:
        old = next((r for r in previous['rows'] if (r['p'], r['n']) == (row['p'], row['n'])), None)
        if old is not None:
            assert row['W_distinct_from_weighted_normalization'] == old[
                'triangle_weight_by_number_of_distinct_edge_cosets']['3']
    output = {
        'scope': 'Eight fixed examples; actual f_H dyadic levels; all edges in distinct H cosets. No field scan.',
        'full_goal_proved': False,
        'rows': rows,
    }
    path = Path('results/parallel40_distinct_levels_2026_09_06.json')
    path.write_text(json.dumps(output, indent=2) + '\n')
    for row in rows:
        print(json.dumps({k: row[k] for k in ['p', 'n', 'f_histogram_on_cosets',
                          'first_positive_repeated_fiber',
                          'first_positive_repeated_same_level_fiber',
                          'independent_direct_checks']}, sort_keys=True))
