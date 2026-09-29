#!/usr/bin/env python3
"""Tiny separate overlap-cover discovery interface; known set-cover/RS ingredients.

This module never calls or relaxes the production partition-cap compiler.
All output completeness claims depend on exhaustive s-set coverage, not greedy
optimality or the completeness of the input candidate catalog.
"""
from collections import Counter, defaultdict
from itertools import combinations
from math import comb, isqrt


def validate_input(p, domain, k, s, u0, u1):
    n = len(domain)
    if type(p) is not int or not 2 <= p <= 5000:
        raise ValueError('This toy interface supports prime fields <=5000.')
    if any(p % d == 0 for d in range(2, isqrt(p) + 1)):
        raise ValueError('Not prime.')
    if not 1 <= k <= s <= n or comb(n, s) > 200000:
        raise ValueError('Invalid dimensions or exhaustive-cover budget exceeded.')
    if len(set(domain)) != n or len(u0) != n or len(u1) != n:
        raise ValueError('Input lengths or distinct domain.')
    if any(type(a) is not int or not 0 <= a < p for a in domain + u0 + u1):
        raise ValueError('Noncanonical field element.')


def evaluate(coefficients, x, p):
    value = 0
    for a in reversed(coefficients):
        value = (value * x + a) % p
    return value


def interpolate(xs, ys, p):
    """Direct coefficient Lagrange interpolation, independent of prior engines."""
    out = [0] * len(xs)
    for i, x in enumerate(xs):
        poly, denominator = [1], 1
        for j, y in enumerate(xs):
            if i == j:
                continue
            nxt = [0] * (len(poly) + 1)
            for d, a in enumerate(poly):
                nxt[d] = (nxt[d] - y * a) % p
                nxt[d + 1] = (nxt[d + 1] + a) % p
            poly = nxt
            denominator = denominator * (x - y) % p
        weight = ys[i] * pow(denominator, -1, p) % p
        out = [(a + weight * b) % p for a, b in zip(out, poly)]
    return out


def size_volume(n, k, s, m):
    return sum(comb(m, i) * comb(n - m, s - i)
               for i in range(k, min(s, m) + 1) if 0 <= s - i <= n - m)


def greedy_cover(n, k, s, support_masks):
    universe = [sum(1 << i for i in a) for a in combinations(range(n), s)]
    masks = sorted(set(support_masks), key=lambda m: (-m.bit_count(), m))
    if any(m < 0 or m >> n or m.bit_count() < k for m in masks):
        raise ValueError('Invalid candidate support.')
    coverage = [sum(1 << i for i, a in enumerate(universe)
                    if (a & m).bit_count() >= k) for m in masks]
    for m, covered in zip(masks, coverage):
        assert covered.bit_count() == size_volume(n, k, s, m.bit_count())
    uncovered = (1 << len(universe)) - 1
    selected, marginal_gains = [], []
    while uncovered:
        j = max(range(len(masks)), key=lambda i: (coverage[i] & uncovered).bit_count())
        gain = (coverage[j] & uncovered).bit_count()
        if gain == 0:
            raise ValueError('Candidate catalog does not cover all s-sets.')
        selected.append(j)
        marginal_gains.append(gain)
        uncovered &= ~coverage[j]
    initial_size = len(selected)
    # A bounded deterministic deletion pass; neither minimality nor optimality.
    for j in list(reversed(selected)):
        other_union = 0
        for ell in selected:
            if ell != j:
                other_union |= coverage[ell]
        if other_union.bit_count() == len(universe):
            selected.remove(j)
    answer = [masks[j] for j in selected]
    return answer, {'candidate_count': len(masks), 'universe_size': len(universe),
                    'greedy_size': initial_size, 'after_redundancy_deletion': len(answer),
                    'greedy_marginal_gains': marginal_gains,
                    'candidate_support_sizes': dict(Counter(m.bit_count() for m in masks)),
                    'selected_support_sizes': dict(Counter(m.bit_count() for m in answer))}


def track_from_base(p, domain, k, u0, u1, base):
    if len(base) != k or len(set(base)) != k:
        raise ValueError('Interpolation base.')
    aa, bb = [interpolate([domain[i] for i in base], [w[i] for i in base], p)
              for w in [u0, u1]]
    support = [i for i, x in enumerate(domain)
               if evaluate(aa, x, p) == u0[i] and evaluate(bb, x, p) == u1[i]]
    return {'base': list(base), 'intercept_coefficients': aa,
            'slope_coefficients': bb, 'joint_support': support}


def complete_nodes_by_buckets(p, domain, k, s, u0, u1, tracks):
    nodes, whole_field_ids = {}, []
    for j, tr in enumerate(tracks):
        aa, bb = tr['intercept_coefficients'], tr['slope_coefficients']
        av, bv = [[evaluate(h, x, p) for x in domain] for h in [aa, bb]]
        common, buckets = [], defaultdict(list)
        for i, (a, b, x, y) in enumerate(zip(av, bv, u0, u1)):
            da, db = (x - a) % p, (y - b) % p
            if db:
                buckets[-da * pow(db, -1, p) % p].append(i)
            elif da == 0:
                common.append(i)
        if len(common) >= s:
            zs = range(p)
            whole_field_ids.append(j)
        else:
            zs = [z for z, coords in buckets.items() if len(common) + len(coords) >= s]
        for z in zs:
            h = tuple((a + z * b) % p for a, b in zip(aa, bb))
            cw = [(a + z * b) % p for a, b in zip(av, bv)]
            support = sorted(common + buckets.get(z, []))
            assert support == [i for i, c in enumerate(cw) if c == (u0[i] + z * u1[i]) % p]
            assert len(support) >= s
            key = z, h
            if key not in nodes:
                nodes[key] = {'scalar': z, 'coefficients': list(h), 'codeword': cw,
                              'agreement_support': support, 'represented_by_tracks': []}
            nodes[key]['represented_by_tracks'].append(j)
    return [nodes[key] for key in sorted(nodes)], whole_field_ids


def compile_certificate(p, domain, k, s, u0, u1, tracks):
    validate_input(p, domain, k, s, u0, u1)
    n = len(domain)
    keys, supports = set(), []
    for tr in tracks:
        aa, bb = tr['intercept_coefficients'], tr['slope_coefficients']
        if len(aa) != k or len(bb) != k or any(type(x) is not int or not 0 <= x < p for x in aa + bb):
            raise ValueError('Degree or coefficient encoding.')
        key = tuple(aa), tuple(bb)
        if key in keys:
            raise ValueError('Duplicate affine track.')
        keys.add(key)
        support = [i for i, x in enumerate(domain)
                   if evaluate(aa, x, p) == u0[i] and evaluate(bb, x, p) == u1[i]]
        if support != tr['joint_support'] or len(support) < k:
            raise ValueError('Incorrect maximal joint support.')
        supports.append(set(support))
    multiplicities, witnesses = Counter(), []
    for a in combinations(range(n), s):
        hits = [j for j, support in enumerate(supports) if len(set(a) & support) >= k]
        if not hits:
            raise ValueError('Uncovered agreement support: ' + str(a))
        multiplicities[len(hits)] += 1
        witnesses.append(hits[0])
    nodes, whole = complete_nodes_by_buckets(p, domain, k, s, u0, u1, tracks)
    return {'schema': 'tiny_overlap_track_cover_v1', 'p': p, 'n': n, 'k': k, 's': s,
            'domain': domain, 'u0': u0, 'u1': u1, 'tracks': tracks,
            'coverage_count': comb(n, s), 'coverage_multiplicity_histogram': dict(multiplicities),
            'lexicographic_s_set_track_witnesses': witnesses,
            'nodes': nodes, 'node_count': len(nodes), 'whole_field_track_ids': whole,
            'event_definition': 'plain at-least-s agreement with degree-less-than-k polynomials',
            'completeness': 'All s-sets independently covered; coordinate affine roots enumerated over the entire field.'}
