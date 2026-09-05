#!/usr/bin/env python3
"""Exact conditional elementary-symmetric moments of a conference sign matrix.

The input contract is a symmetric sign matrix S with diagonal zero, row sums
zero, and S S^T = q I-J, represented by cells at marked columns and ordered
cell-pair sign counts. Validation here checks the compressed consistency
conditions; a provenance check of the underlying matrix is still required.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from functools import lru_cache
import json
from math import comb
from pathlib import Path
import time


def choose(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def elementary(signs, degree):
    out = [1] + [0] * degree
    for s in signs:
        for j in range(degree, 0, -1):
            out[j] += s * out[j-1]
    return tuple(out)


def signed_binomial(positive, negative, degree):
    return sum((-1)**j * choose(negative, j) * choose(positive, degree-j)
               for j in range(degree+1))


def weights(population, sample, limit):
    out = [Fraction(1)]
    for k in range(1, limit+1):
        out.append(out[-1] * Fraction(sample-k+1, population-k+1))
    return out


def validate_counts(record):
    q = record.get('q', record.get('p'))
    marks = record['marks']
    assert isinstance(q, int) and q >= 5 and q % 4 == 1
    assert len(set(marks)) == len(marks) and all(0 <= x < q for x in marks)
    m = len(marks)
    cells = {}
    for cell in record['cells']:
        u, size = tuple(cell['pattern']), cell['size']
        assert len(u) == m and all(s in (-1, 0, 1) for s in u) and u.count(0) <= 1
        assert u not in cells and isinstance(size, int) and size > 0
        cells[u] = size
    assert sum(cells.values()) == q
    relations = Counter()
    for rel in record['relations']:
        u, v = tuple(rel['left_pattern']), tuple(rel['right_pattern'])
        s, count = rel['sign'], rel['count']
        assert u in cells and v in cells and s in (-1, 0, 1)
        assert isinstance(count, int) and count > 0 and (u, v, s) not in relations
        relations[u, v, s] = count
    assert sum(relations.values()) == q*q
    for u, nu in cells.items():
        assert sum(relations[u, v, s] for v in cells for s in (-1, 0, 1)) == nu*q
        assert sum(s * relations[u, v, s] for v in cells for s in (-1, 1)) == 0
        for v, nv in cells.items():
            assert sum(relations[u, v, s] for s in (-1, 0, 1)) == nu*nv
            assert relations[u, v, 0] == (nu if u == v else 0)
            for s in (-1, 0, 1):
                assert relations[u, v, s] == relations[v, u, s]
    for i in range(m):
        zero_cells = [u for u in cells if u[i] == 0]
        assert len(zero_cells) == 1 and cells[zero_cells[0]] == 1
        z = zero_cells[0]
        for u, nu in cells.items():
            assert relations[u, z, u[i]] == nu
        for sign, size in ((-1, (q-1)//2), (0, 1), (1, (q-1)//2)):
            assert sum(nu for u, nu in cells.items() if u[i] == sign) == size
        for j in range(m):
            assert sum(nu*u[i]*u[j] for u, nu in cells.items()) == (q-1 if i == j else -1)
    return q, tuple(marks), cells, relations


def full_pair_histogram(q, s):
    if s == 0:
        return Counter({(0, 0): 1, (1, 1): (q-1)//2, (-1, -1): (q-1)//2})
    return Counter({(0, s): 1, (s, 0): 1,
                    (1, 1): (q-3-2*s)//4,
                    (-1, -1): (q-3+2*s)//4,
                    (1, -1): (q-1)//4, (-1, 1): (q-1)//4})


def canonical_pair(u, v, s):
    # Simultaneous mark permutations, row exchange, and replacing S by -S
    # preserve a product of two degree-d coefficients, also for odd d.
    choices = []
    for flip in (1, -1):
        pairs = tuple(sorted((flip*a, flip*b) for a, b in zip(u, v)))
        choices.extend(((flip*s, pairs), (flip*s, tuple(sorted((b, a) for a, b in pairs)))))
    return min(choices)


@lru_cache(maxsize=2048)
def type_power(a, b, count, degree, limit):
    """(1+z*(a*t+b*w+a*b*t*w))**count, bounded degrees."""
    result = {}
    for i in range(degree+1):
        for j in range(degree+1):
            if (a == 0 and i) or (b == 0 and j):
                continue
            for k in range(max(i, j), min(i+j, limit, count)+1):
                lx, ly, lz = k-j, k-i, i+j-k
                value = comb(count, lx) * comb(count-lx, ly) * comb(count-lx-ly, lz)
                value *= a**i * b**j
                if value:
                    result[i, j, k] = value
    return tuple(result.items())


def multiply(left, right, degree, limit):
    out = defaultdict(int)
    for (i, j, k), a in left.items():
        for (ii, jj, kk), b in right:
            if i+ii <= degree and j+jj <= degree and k+kk <= limit:
                out[i+ii, j+jj, k+kk] += a*b
    return {key: value for key, value in out.items() if value}


@lru_cache(maxsize=16384)
def pair_coefficients(q, degree, limit, canonical):
    """Coefficients indexed by number of distinct outside columns used."""
    s, pairs = canonical
    hist = full_pair_histogram(q, s)
    hist.subtract(pairs)
    assert all(count >= 0 for count in hist.values())
    assert sum(hist.values()) == q-len(pairs)
    pol = {(0, 0, 0): 1}
    for (a, b), count in sorted(hist.items()):
        if count and (a or b):
            pol = multiply(pol, type_power(a, b, count, degree, limit), degree, limit)
    marked_left = elementary((a for a, b in pairs), degree)
    marked_right = elementary((b for a, b in pairs), degree)
    out = [0]*(limit+1)
    for (i, j, k), value in pol.items():
        out[k] += value * marked_left[degree-i] * marked_right[degree-j]
    return tuple(out)


def row_mean(q, n, degree, u):
    remaining = n-len(u)
    positive = (q-1)//2-u.count(1)
    negative = (q-1)//2-u.count(-1)
    marked = elementary(u, degree)
    limit = min(degree, remaining, q-len(u))
    prob = weights(q-len(u), remaining, limit)
    return sum((marked[degree-k] * signed_binomial(positive, negative, k) * prob[k]
                for k in range(limit+1)), Fraction())


def kernel_second_from_cells(q, degree, cells):
    """Specialization n=degree: mutual signs cancel before aggregation."""
    m = len(next(iter(cells)))
    outside_degree = degree-m
    denominator = comb(q-m, outside_degree)
    total = 0
    for u, nu in cells.items():
        for v, nv in cells.items():
            mark_product = 1
            for a, b in zip(u, v):
                mark_product *= a*b
            if not mark_product:
                continue
            positive = (q-3)//2-sum(a*b == 1 for a, b in zip(u, v))
            negative = (q-1)//2-sum(a*b == -1 for a, b in zip(u, v))
            distinct = nu*nv-(nu if u == v else 0)
            if distinct:
                assert positive >= 0 and negative >= 0
                total += distinct * mark_product * signed_binomial(positive, negative, outside_degree)
            if u == v:
                total += nu * choose(q-1-m, outside_degree)
    return Fraction(total, denominator)


def conditional_moments(record, n, degree=6, *, validate=True, audit_sign_dependence=False):
    started = time.perf_counter()
    if validate:
        q, marks, cells, relations = validate_counts(record)
    else:
        q = record.get('q', record.get('p'))
        marks = tuple(record['marks'])
        cells = {tuple(x['pattern']): x['size'] for x in record['cells']}
        relations = {(tuple(x['left_pattern']), tuple(x['right_pattern']), x['sign']): x['count']
                     for x in record['relations']}
    assert 0 <= degree <= n <= q and len(marks) <= n
    limit = min(2*degree, n-len(marks), q-len(marks))
    probability = weights(q-len(marks), n-len(marks), limit)
    grouped = Counter()
    for (u, v, s), count in relations.items():
        grouped[canonical_pair(u, v, s)] += count
    aggregate = [0]*(limit+1)
    for canonical, count in grouped.items():
        for k, coefficient in enumerate(pair_coefficients(q, degree, limit, canonical)):
            aggregate[k] += count*coefficient
    mean = sum((count*row_mean(q, n, degree, u) for u, count in cells.items()), Fraction())
    second = sum((value*p for value, p in zip(aggregate, probability)), Fraction())
    variance = second-mean*mean
    assert variance >= 0
    out = {'q': q, 'marks': list(marks), 'n': n, 'degree': degree,
           'mean': [mean.numerator, mean.denominator],
           'second_moment': [second.numerator, second.denominator],
           'variance': [variance.numerator, variance.denominator],
           'mean_float': float(mean), 'second_moment_float': float(second),
           'variance_float': float(variance),
           'cell_count': len(cells), 'relation_count': len(relations),
           'canonical_relation_classes': len(grouped),
           'outside_union_coefficients': aggregate,
           'max_aggregate_integer_bits': max(abs(x).bit_length() for x in aggregate)}
    if n == degree:
        shortcut = kernel_second_from_cells(q, degree, cells)
        assert shortcut == second
        out['n_equals_degree_cell_only_shortcut_checked'] = True
    if audit_sign_dependence:
        different = 0
        possible = 0
        max_gap = Fraction()
        for u in cells:
            for v in cells:
                if not all(relations.get((u, v, s), 0) for s in (-1, 1)):
                    continue
                values = []
                for s in (-1, 1):
                    coefs = pair_coefficients(q, degree, limit, canonical_pair(u, v, s))
                    values.append(sum((c*p for c, p in zip(coefs, probability)), Fraction()))
                possible += 1
                gap = abs(values[0]-values[1])
                different += bool(gap)
                max_gap = max(max_gap, gap)
        out['sign_dependence_audit'] = {'cell_pairs_with_both_signs': possible,
                                      'cell_pairs_with_different_coefficients': different,
                                      'max_pair_coefficient_gap': [max_gap.numerator, max_gap.denominator]}
    out['elapsed_seconds'] = time.perf_counter()-started
    return out


def paley_prime_matrix(p):
    import numpy as np
    assert p >= 5 and p % 4 == 1 and all(p % k for k in range(2, int(p**0.5)+1))
    chi = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        chi[x] = 1 if pow(x, (p-1)//2, p) == 1 else -1
    a = np.arange(p)
    return chi[(a[:, None]-a[None, :]) % p]


def matrix_counts(matrix, marks, *, validate_matrix=True):
    import numpy as np
    s = np.asarray(matrix, dtype=np.int64)
    q = len(s)
    if validate_matrix:
        assert s.shape == (q, q) and np.array_equal(s, s.T)
        assert not np.diag(s).any() and set(np.unique(s)) <= {-1, 0, 1}
        assert np.count_nonzero(s) == q*(q-1) and not s.sum(axis=1).any()
        assert np.array_equal(s @ s.T, q*np.eye(q, dtype=np.int64)-np.ones((q,q), dtype=np.int64))
    patterns = [tuple(map(int, s[x, marks])) for x in range(q)]
    cells = Counter(patterns)
    relations = Counter((patterns[x], patterns[y], int(s[x,y])) for x in range(q) for y in range(q))
    return {'q': q, 'marks': list(marks),
            'cells': [{'pattern': list(u), 'size': count} for u, count in sorted(cells.items())],
            'relations': [{'left_pattern': list(u), 'right_pattern': list(v), 'sign': sign, 'count': count}
                          for (u, v, sign), count in sorted(relations.items())],
            'provenance': 'direct matrix; symmetric balanced conference identities checked' if validate_matrix else 'matrix supplied by caller'}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('counts', type=Path)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--degree', type=int, default=6)
    parser.add_argument('--audit-sign-dependence', action='store_true')
    parser.add_argument('--output', type=Path)
    args = parser.parse_args()
    out = conditional_moments(json.loads(args.counts.read_text()), args.n, args.degree,
                              audit_sign_dependence=args.audit_sign_dependence)
    encoded = json.dumps(out, indent=2)+'\n'
    if args.output:
        args.output.write_text(encoded)
    else:
        print(encoded, end='')


if __name__ == '__main__':
    main()
