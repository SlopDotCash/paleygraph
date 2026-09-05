#!/usr/bin/env python3
"""Exact two-anchor necklace identities; no Paley proof is claimed.

Matrix arithmetic uses NumPy object arrays containing Python integers.
Independent multiplicative convolution and small literal tuple sums check
the trace formulas. Literature-derived bounds are tested, not inferred
from these finite checks.
"""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from paley_exact import character, prime

ROOT = Path(__file__).resolve().parents[1]


def trace(a):
    return sum(int(x) for x in a.diagonal())


def paired_trace(a, b):
    return sum(int(x) for x in (a * b.T).flat)


def literal_necklace(p, chi, labels):
    """Definition, with each label a singleton anchor; includes zero entries."""
    answer = 0
    k = len(labels)
    for xs in product(range(p), repeat=k):
        term = 1
        for i, anchor in enumerate(labels):
            term *= chi[(xs[(i + 1) % k] - xs[i]) % p]
            term *= chi[(xs[i] - anchor) % p]
        answer += term
    return answer


def convolution_moments(p, chi, depth):
    """Coefficient at 1 in powers of t -> chi(1-t) on F_p^*."""
    counts = [0] * p
    counts[1] = 1
    moments = [1]
    for _ in range(depth):
        following = [0] * p
        for x in range(1, p):
            for y in range(1, p):
                following[x * y % p] += counts[x] * chi[(1 - y) % p]
        counts = following
        moments.append(counts[1])
    return moments


def case(p, depth):
    assert prime(p) and p % 4 == 1
    chi = character(p)
    s = np.array([[chi[(x - y) % p] for y in range(p)]
                  for x in range(p)], dtype=object)
    d = np.diag(np.array(chi, dtype=object))
    eye = np.eye(p, dtype=object)
    one = np.ones((p, p), dtype=object)
    assert np.array_equal(s @ s, p * eye - one)
    c = d @ s
    restricted = c[1:, 1:]
    ww = np.array(chi[1:], dtype=object)
    rr = np.ones(p - 1, dtype=object)
    gram = p * np.eye(p - 1, dtype=object) - np.outer(rr, rr) - np.outer(ww, ww)
    assert np.array_equal(restricted @ restricted.T, gram)
    assert np.array_equal(restricted.T @ restricted, gram)
    assert np.array_equal(restricted @ rr, -rr)
    assert np.array_equal(restricted @ ww, -ww)
    dd = d[1:, 1:]
    assert np.array_equal(dd @ restricted @ dd, restricted.T)

    # K_j = S(DS)^(j-1), with full zero row/column retained.
    kernels = [None, s]
    for j in range(2, depth + 1):
        kernels.append(kernels[-1] @ c)
    ts = convolution_moments(p, chi, depth)
    for j in range(1, depth + 1):
        kernel = kernels[j]
        assert np.array_equal(kernel, kernel.T)
        assert trace(d @ kernel) == (p - 1) * ts[j]
        assert kernel[0, 0] == 0
        for x in range(1, p):
            assert kernel[x, x] == chi[x] * ts[j]
            assert kernel[0, x] == (-1) ** (j - 1) * chi[x]
        # Lu-Zheng-Zheng (2.3), with tensor powers 1,1, implies this bound.
        assert ts[j] ** 2 <= (j - 1) ** 2 * p ** (j - 1)

    one_anchor_checks = 0
    two_anchor_checks = 0
    cross_trace_checks = 0
    direct_checks = 0
    checksums = sha256()
    anchors = list(range(1, p)) if p <= 17 else [1, next(x for x in range(1, p) if chi[x] == -1), p - 1]
    examples = []
    for k in range(1, depth + 1):
        for b in anchors:
            db = np.array([chi[(x - b) % p] for x in range(p)], dtype=object)
            observed = trace(db[:, None] * kernels[k])
            assert observed == -ts[k], (p, k, b, 'one')
            one_anchor_checks += 1
            if (p == 5 and k <= 5) or (p == 13 and k <= 3 and b <= 2):
                assert literal_necklace(p, chi, [b] + [0] * (k - 1)) == observed
                direct_checks += 1
        for ell in range(1, k):
            r = k - ell
            ss = min(ell, r)
            difference = abs(ell - r)
            cross_expected = (p ** ss * (p - 1) * ts[difference]
                              + 2 * (-1) ** k * (p - p ** ss))
            cross = paired_trace(kernels[ell], kernels[r])
            assert cross == cross_expected, (p, k, ell, 'cross')
            cross_trace_checks += 1
            predicted = (p * ts[ell] * ts[r] - p ** ss * ts[difference] - ts[k]
                         + 2 * (-1) ** k * sum(p ** j for j in range(1, ss)))
            assert predicted ** 2 <= k ** 4 * p ** k
            if k in (4, 6, 8) and ell <= r:
                examples.append({'k': k, 'ell': ell, 'value': predicted})
            for b in anchors:
                db = np.array([chi[(x - b) % p] for x in range(p)], dtype=object)
                observed = paired_trace(db[:, None] * kernels[ell],
                                        db[:, None] * kernels[r])
                assert observed == predicted, (p, k, ell, b, 'two')
                two_anchor_checks += 1
                checksums.update(f'{k},{ell},{b},{observed}\n'.encode())
                if (p == 5 and k <= 5) or (p == 13 and k <= 3 and b <= 2):
                    labels = [0] * k
                    labels[0] = labels[ell] = b
                    assert literal_necklace(p, chi, labels) == observed
                    direct_checks += 1
            if p <= 17:
                # Includes the coincident-anchor case b=0 and verifies the
                # averaging step directly, not just the closed formula.
                average = (p - 1) * ts[k]
                for b in range(1, p):
                    db = np.array([chi[(x - b) % p] for x in range(p)], dtype=object)
                    average += paired_trace(db[:, None] * kernels[ell],
                                            db[:, None] * kernels[r])
                assert average == p * (p - 1) * ts[ell] * ts[r] - cross

    # Binary patterns outside this note's <=2 occurrence family. Their length-six
    # formulas are subsequently proved in planar_necklace_reductions.py's note.
    triple_patterns = []
    triple_averaging_checks = 0
    for labels in ([0, 1, 0, 1, 0, 1], [0, 0, 0, 1, 1, 1]):
        power = eye.copy()
        for b in labels:
            db = np.array([chi[(x - b) % p] for x in range(p)], dtype=object)
            power = power @ (db[:, None] * s)
        value = trace(power)
        item = {'labels': labels, 'value': value,
                'status': 'Exact finite value; length-six formula proved in planar-necklace-reductions.md'}
        if p <= 17:
            marked = [i for i, b in enumerate(labels) if b == 1]
            gaps = [(marked[(i + 1) % 3] - marked[i]) % len(labels) for i in range(3)]
            repeated = distinct = 0
            for x, y, z in product(range(p), repeat=3):
                cubic = sum(chi[(x-b) % p] * chi[(y-b) % p] * chi[(z-b) % p]
                            for b in range(p))
                if len({x, y, z}) == 3:
                    assert cubic ** 2 <= 4 * p
                elif x == y == z:
                    assert cubic == 0
                else:
                    duplicate = x if x == y or x == z else y
                    singleton = next(a for a in (x, y, z) if a != duplicate)
                    assert cubic == -chi[(singleton - duplicate) % p]
                term = (cubic * kernels[gaps[0]][x, y]
                        * kernels[gaps[1]][y, z] * kernels[gaps[2]][z, x])
                if len({x, y, z}) == 3:
                    distinct += int(term)
                else:
                    repeated += int(term)
            assert distinct + repeated == (p - 1) * (value + ts[len(labels)])
            triple_averaging_checks += 1
            item['averaged_distinct_contribution'] = distinct
            item['averaged_repeated_contribution'] = repeated
        triple_patterns.append(item)
    return {'p': p, 'depth': depth, 'tested_second_anchors': anchors,
            'monochromatic_normalized_traces': ts,
            'one_anchor_checks': one_anchor_checks,
            'two_anchor_checks': two_anchor_checks,
            'cross_trace_checks': cross_trace_checks,
            'literal_tuple_checks': direct_checks,
            'triple_averaging_checks': triple_averaging_checks,
            'two_anchor_values_sha256': checksums.hexdigest(),
            'examples': examples, 'three_occurrence_examples': triple_patterns}


def main():
    cases = []
    for p, depth in [(5, 12), (13, 12), (17, 12), (29, 12), (41, 12), (61, 10), (97, 10)]:
        result = case(p, depth)
        cases.append(result)
        print(json.dumps({'p': p, 'depth': depth, 'status': 'passed'}), flush=True)
    keys = ['one_anchor_checks', 'two_anchor_checks', 'cross_trace_checks',
            'literal_tuple_checks', 'triple_averaging_checks']
    sources = ['experiments/localized_necklace_identities.py', 'experiments/paley_exact.py',
               'sources/kunisky-2303.16475v1.html', 'sources/lu-zheng-zheng-1305.3405v3.html']
    result = {'status': 'Exact identities and literature-bound specializations passed',
              'scope': 'Singleton two-anchor necklaces with one or two occurrences of one anchor; not the full necklace or Paley conjectures',
              'arithmetic': 'Python integers in NumPy object arrays',
              'numpy_version': np.__version__, 'cases': cases,
              'totals': {key: sum(c[key] for c in cases) for key in keys},
              'source_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in sources}}
    target = ROOT / 'results/localized_necklace_identities.json'
    target.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(result['totals']), flush=True)


if __name__ == '__main__':
    main()
