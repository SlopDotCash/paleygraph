#!/usr/bin/env python3
"""Exact first and second moments on the one-swap neighbourhood of an actual set.

The production backend streams q rows and retains an n by n internal contraction.
It never evaluates the n(q-n) neighbouring sets. Integer bounds are checked before
NumPy arithmetic; rational outputs use Python integers.
"""
from collections import Counter
from fractions import Fraction
from math import comb, isqrt
from time import perf_counter
import numpy as np

INT64 = (1 << 63) - 1


def frac(x):
    return [x.numerator, x.denominator]


def _domain(q, selected, degree):
    selected = list(selected)
    if type(q) is not int or not 2 <= q <= 10000000:
        raise ValueError('q must be an integer in [2,10000000] for this prototype')
    if any(type(x) is not int or not 0 <= x < q for x in selected):
        raise ValueError('selected points must be integer field labels')
    if len(set(selected)) != len(selected) or not 1 <= len(selected) <= min(64, q-1):
        raise ValueError('selected points must be distinct, with 1 <= n <= min(64,q-1)')
    if type(degree) is not int or not 0 <= degree <= min(6, len(selected)):
        raise ValueError('implemented degrees are 0 through min(6,n)')
    return sorted(selected)


def coefficient_table(n, degree):
    out = np.zeros((n+1, n+1), dtype=np.int64)
    if degree < 0:
        return out
    for total in range(n+1):
        for positive in range(total+1):
            out[total, positive] = sum((-1)**j * comb(total-positive, j) * comb(positive, degree-j)
                for j in range(max(0, degree-positive), min(degree, total-positive)+1))
    return out


def _compile(q, selected, degree, row_blocks, chunk_size):
    start = perf_counter()
    n = len(selected)
    target_table = coefficient_table(n, degree)
    derivative_table = coefficient_table(n, degree-1)
    bound = comb(n-1, degree-1) if degree else 0
    if q * max(1, bound) > INT64:
        raise ValueError('integer contraction bound exceeds int64')
    effective_chunk = min(chunk_size, INT64 // max(1, bound*bound))
    if effective_chunk < 1:
        raise ValueError('integer norm bound exceeds int64')
    sums = [0] * n
    norms = [0] * n
    internal = np.zeros((n, n), dtype=np.int64)
    target = 0
    row_count = 0
    histogram = Counter()
    for signs in row_blocks(effective_chunk):
        count = len(signs)
        if count > effective_chunk or signs.shape != (count, n):
            raise ValueError('invalid row block')
        positive = (signs == 1).sum(axis=1)
        nonzero = (signs != 0).sum(axis=1)
        values = derivative_table[nonzero[:, None] - (signs != 0), positive[:, None] - (signs == 1)]
        target += int(target_table[nonzero, positive].sum(dtype=np.int64))
        internal += values.T @ signs.astype(np.int64)
        for i, (s, v) in enumerate(zip(values.sum(axis=0), (values*values).sum(axis=0))):
            sums[i] += int(s)
            norms[i] += int(v)
        encoded, multiplicities = np.unique((n-nonzero)*(n+1)+positive, return_counts=True)
        histogram.update({(int(x)//(n+1), int(x)%(n+1)): int(m) for x, m in zip(encoded, multiplicities)})
        row_count += count
    assert row_count == q
    matrix = [[int(x) for x in row] for row in internal]
    denominator = n * (q-n)
    first = second = internal_energy = complete_energy = 0
    deletion_records = []
    for a, row in enumerate(matrix):
        at_deleted = row[a]
        row_sum = sum(row)
        inside_square_sum = sum(x*x for x in row)
        total_square_sum = q*norms[a] - sums[a]*sums[a]
        drift_sum = -row_sum - (q-n)*at_deleted
        square_sum = total_square_sum - inside_square_sum + 2*at_deleted*row_sum + (q-n)*at_deleted**2
        assert total_square_sum >= inside_square_sum and square_sum >= 0
        first += drift_sum
        second += square_sum
        internal_energy += inside_square_sum
        complete_energy += total_square_sum
        deletion_records.append({'point': selected[a], 'insertion_sum': sums[a], 'insertion_norm_squared': norms[a],
            'complete_transform_norm_squared': total_square_sum, 'delta_sum': drift_sum, 'delta_square_sum': square_sum})
    drift = Fraction(first, denominator)
    delta_second = Fraction(second, denominator)
    variance = delta_second - drift*drift
    assert variance >= 0
    return {'q': q, 'n': n, 'degree': degree, 'selected': selected, 'target': target,
            'neighbour_count': denominator, 'delta_mean': frac(drift), 'delta_second_moment': frac(delta_second),
            'neighbour_mean': frac(target+drift), 'neighbour_second_moment': frac(target*target+2*target*drift+delta_second),
            'neighbour_variance': frac(variance), 'internal_contraction': matrix,
            'deletion_records': deletion_records, 'complete_insertion_energy': complete_energy,
            'excluded_internal_energy': internal_energy,
            'signed_row_count_histogram': [[z, p, count] for (z, p), count in sorted(histogram.items())],
            'arithmetic': 'Exact bounded int64 block contractions, unbounded integer aggregation, rational moments.',
            'work': {'streamed_rows': q, 'internal_matrix_entries': n*n, 'neighbour_targets_evaluated': 0,
                     'effective_chunk': effective_chunk, 'seconds': perf_counter()-start}}


def prime_character(q):
    if q % 4 != 1 or any(q % a == 0 for a in range(2, isqrt(q)+1)):
        raise ValueError('prime q congruent to 1 mod 4 required')
    character = -np.ones(q, dtype=np.int8)
    x = np.arange(1, (q+1)//2, dtype=np.int64)
    character[(x*x) % q] = 1
    character[0] = 0
    assert int(character.sum()) == 0
    return character


def from_prime(q, selected, degree=6, chunk_size=8192):
    selected = _domain(q, selected, degree)
    if type(chunk_size) is not int or chunk_size <= 0:
        raise ValueError('positive integer chunk size required')
    character = prime_character(q)
    columns = np.array(selected, dtype=np.int64)
    def blocks(size):
        for start in range(0, q, size):
            rows = np.arange(start, min(q, start+size), dtype=np.int64)
            yield character[(rows[:, None]-columns) % q]
    result = _compile(q, selected, degree, blocks, chunk_size)
    result['input_contract'] = 'Prime Paley character; primality checked by trial division, conference identity follows algebraically.'
    return result


def from_matrix(matrix, selected, degree=6, chunk_size=8192):
    raw = np.asarray(matrix)
    q = len(raw)
    selected = _domain(q, selected, degree)
    if q > 257 or raw.shape != (q, q) or not np.isin(raw, (-1, 0, 1)).all():
        raise ValueError('small square sign matrix required')
    if type(chunk_size) is not int or chunk_size <= 0:
        raise ValueError('positive integer chunk size required')
    signs = raw.astype(np.int64)
    if (not np.array_equal(signs, signs.T) or np.any(np.diag(signs)) or
        np.any(signs.sum(axis=1)) or not np.array_equal(signs @ signs, q*np.eye(q, dtype=np.int64)-1)):
        raise ValueError('balanced symmetric conference identities required')
    def blocks(size):
        for start in range(0, q, size):
            yield signs[start:start+size, selected]
    result = _compile(q, selected, degree, blocks, chunk_size)
    result['input_contract'] = 'Actual sign matrix: symmetry, diagonal, balance and conference identity checked exactly.'
    return result
