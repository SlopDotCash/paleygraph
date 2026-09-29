#!/usr/bin/env python3
"""Exact classification of the known finite witness's short ideal multiples.

No uniform-in-order theorem is certified by this finite computation.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import combinations
from pathlib import Path
import importlib.util
import json
import time

ROOT = Path(__file__).resolve().parents[1]
CHECKS = Counter()
SPEC = importlib.util.spec_from_file_location('triangle_check', ROOT/'experiments/parallel30_verify_2026_09_05.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)


def check(ok, label):
    assert ok, label
    CHECKS[label] += 1


def multiply(a, b):
    d = len(a)
    assert len(b) == d
    c = [0]*d
    for i, x in enumerate(a):
        if not x:
            continue
        for j, y in enumerate(b):
            if not y:
                continue
            q, r = divmod(i+j, d)
            c[r] += x*y*(-1 if q else 1)
    return c


def norm_adjugate(a):
    d = len(a)
    if d == 1:
        return a[0], [1]
    e, o = a[::2], a[1::2]
    ee, oo = multiply(e, e), multiply(o, o)
    # f(X) f(-X) = e(Y)^2-Y o(Y)^2, Y=X^2.
    descend = ee[:]
    for i, x in enumerate(oo):
        descend[(i+1) % (d//2)] -= x*(-1 if i+1 == d//2 else 1)
    value, b = norm_adjugate(descend)
    lifted = [0]*d
    lifted[::2] = b
    opposite = [x*(-1 if i % 2 else 1) for i, x in enumerate(a)]
    out = multiply(opposite, lifted)
    check(multiply(a, out) == [value]+[0]*(d-1), 'recursive_adjugate_identity')
    return value, out


def main():
    started = time.monotonic()
    p, n, g = 215535361, 128, 25525303
    h = sorted(pow(g, j, p) for j in range(n))
    exponent = {pow(g, j, p): j for j in range(n)}
    check(len(exponent) == n and pow(g, n//2, p) == p-1, 'exact_generator_order')
    f = [int(i in [0, 1, 19]) for i in range(n//2)]
    norm, adjugate = norm_adjugate(f)
    check(norm == p, 'known_prime_norm')
    bucket = defaultdict(list)
    without_one = [x for x in h if x != 1]
    for a, b in combinations(without_one, 2):
        bucket[(1+a+b) % p].append((a, b))
    normalized = set()
    for tail in combinations(without_one, 3):
        for a, b in bucket.get(-sum(tail) % p, ()):
            word = (1, a, b, *tail)
            if len(set(word)) != 6 or any(-x % p in word for x in word):
                continue
            if MOD.balanced(word, p):
                continue
            normalized.add(tuple(sorted(exponent[x] for x in word)))
    canonical = lambda e: min(tuple(sorted((x-y) % n for x in e)) for y in e)
    representatives = sorted({canonical(e) for e in normalized})
    check(len(normalized) == 6*len(representatives), 'six_normalizations_per_free_orbit')
    check(len(representatives) == 119, 'all_known_D6_orbits_recovered')
    rows = []
    profiles = defaultdict(Counter)
    for e in representatives:
        target = [0]*(n//2)
        for j in e:
            target[j % (n//2)] += 1 if j < n//2 else -1
        check(sum(x*x for x in target) == 6, 'distinct_opposite_free_coefficient_norm')
        numerator = multiply(target, adjugate)
        check(all(x % p == 0 for x in numerator), 'integral_ideal_quotient')
        quotient = [x//p for x in numerator]
        check(multiply(f, quotient) == target, 'exact_short_multiple_identity')
        word = [pow(g, j, p) for j in e]
        zero_triples = sum(sum(t) % p == 0 for t in combinations(word, 3))
        check(zero_triples in [0, 2], 'unique_zero_triple_partition')
        kind = 'triangular' if zero_triples else 'primitive'
        l1 = sum(abs(x) for x in quotient)
        l2 = sum(x*x for x in quotient)
        profiles[kind][(l1, l2)] += 1
        # Each signed monomial in q represents one zero triangle. Reduce its
        # three entries to exponents modulo128, retaining their multiplicity.
        triangles = []
        total = Counter()
        for j, c in enumerate(quotient):
            for _ in range(abs(c)):
                triangle = sorted((j+(64 if c < 0 else 0)+k) % n for k in [0, 1, 19])
                triangles.append(triangle)
                total.update(triangle)
        reduced = Counter({j: total[j]-total[j+64] for j in range(64) if total[j] != total[j+64]})
        check(all(reduced[j] == target[j] for j in range(64)), 'triangle_sum_integer_reduction')
        # Build one explicit opposite-cancellation graph. All complete pairings
        # have the same edge count; connectivity for primitive targets also
        # follows from the ordinary injectivity argument in the companion note.
        occurrences = defaultdict(list)
        for vertex, tri in enumerate(triangles):
            for entry in tri:
                occurrences[entry].append(vertex)
        edges = []
        survivors = []
        for j in range(64):
            left, right = occurrences[j], occurrences[j+64]
            paired = min(len(left), len(right))
            edges.extend(zip(left[:paired], right[:paired]))
            survivors.extend([j]*(len(left)-paired))
            survivors.extend([j+64]*(len(right)-paired))
        check(sorted(survivors) == list(e), 'graph_survivors_are_target')
        check(all(u != v for u, v in edges), 'no_internal_triangle_cancellation')
        check(2*len(edges) == 3*l1-6, 'graph_edge_count')
        parent = list(range(l1))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        for u, v in edges:
            parent[find(u)] = find(v)
        components = len({find(i) for i in range(l1)})
        check(components == (2 if kind == 'triangular' else 1), 'graph_component_count')
        cycles = len(edges)-l1+components
        if kind == 'primitive':
            check(cycles == l1//2-2, 'primitive_graph_cycle_rank')
        rows.append({'exponents': e, 'kind': kind, 'quotient': quotient,
                     'quotient_l1': l1, 'quotient_l2_squared': l2,
                     'quotient_support': sum(x != 0 for x in quotient),
                     'signed_triangle_terms': triangles,
                     'cancellation_graph': {'edges': edges, 'components': components,
                                            'cycle_rank': cycles}})
    for kind, expected in [('triangular', 54), ('primitive', 65)]:
        check(sum(profiles[kind].values()) == expected, 'known_triangle_split_recovered')
    result = {'status': 'finite short-ideal-multiple classification passed; full goal unproved',
              'p': p, 'n': n, 'generator': g, 'norm': norm, 'adjugate': adjugate,
              'normalized_six_sets': len(normalized), 'orbits': rows,
              'profiles': {kind: [{'l1': a, 'l2_squared': b, 'count': c} for (a, b), c in sorted(profile.items())]
                           for kind, profile in profiles.items()},
              'check_counts': dict(CHECKS), 'elapsed_seconds': time.monotonic()-started,
              'input_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in [
                  'experiments/parallel31_short_multiples_2026_09_05.py',
                  'experiments/parallel30_verify_2026_09_05.py',
                  'results/parallel30_verification_2026_09_05.json']},
              'scope': 'Complete orbit classification for one fixed prime and subgroup; no bound uniform in p or n.'}
    (ROOT/'results/parallel31_short_multiples_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['profiles', 'normalized_six_sets', 'check_counts', 'elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
