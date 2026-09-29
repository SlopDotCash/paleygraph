#!/usr/bin/env python3
"""Eliminate three-mark relation data except one explicit integer contraction.

This is a locally derived specialization using established Walsh expansion and
conference identities. It is not a historical novelty certification or bound.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import combinations, product
import json
from math import prod
from pathlib import Path

from conditional_moments import (canonical_pair, conditional_moments,
                                 pair_coefficients, row_mean, validate_counts, weights)

PATTERNS = tuple(product((-1, 1), repeat=3))
SUBSETS = tuple(tuple(i for i in range(3) if mask & (1 << i)) for mask in range(8))


def char(u, indices):
    return prod(u[i] for i in indices)


def fraction_pair(x):
    return [x.numerator, x.denominator]


def coefficient(q, n, degree, u, v, sign):
    limit = min(2*degree, n-3, q-3)
    coefficients = pair_coefficients(q, degree, limit, canonical_pair(u, v, sign))
    return sum((a*b for a, b in zip(coefficients, weights(q-3, n-3, limit))), Fraction())


@lru_cache(maxsize=128)
def walsh_decomposition(q, n, degree):
    # Both hypothetical signs must admit nonnegative paired-row type counts
    # for every all-nonzero mark pattern; q >= 17 suffices for three marks.
    assert q >= 17 and q % 4 == 1 and 3 <= n <= q and 0 <= degree <= n
    gap = {(u, v): coefficient(q, n, degree, u, v, 1)-coefficient(q, n, degree, u, v, -1)
           for u in PATTERNS for v in PATTERNS}
    walsh = {(I, J): sum((char(u, I)*char(v, J)*gap[u, v]
                          for u in PATTERNS for v in PATTERNS), Fraction())/64
             for I in SUBSETS for J in SUBSETS}
    walsh = {key: value for key, value in walsh.items() if value}
    # Reconstruction is checked exactly, not inferred from expected symmetry.
    for u in PATTERNS:
        for v in PATTERNS:
            assert gap[u, v] == sum((value*char(u, I)*char(v, J)
                                     for (I, J), value in walsh.items()), Fraction())
    assert all((len(I)+len(J)) % 2 == 1 for I, J in walsh)
    for (I, J), value in walsh.items():
        assert walsh[J, I] == value
    unknown = {key: value for key, value in walsh.items() if min(map(len, key)) >= 2}
    assert all(sorted(map(len, key)) == [2, 3] for key in unknown)
    assert len(set(unknown.values())) <= 1
    if unknown:
        assert len(unknown) == 6
    c = next(iter(unknown.values()), Fraction())
    return walsh, c


def reduce_from_cells(q, marks, cells, Q, n, degree=6):
    """Compile using only a cell histogram and the one supplied contraction.

    The arithmetic backend must certify Q and the underlying conference matrix.
    This function does not need or reconstruct any cell-pair relation table.
    """
    assert len(marks) == len(set(marks)) == 3 and all(0 <= m < q for m in marks)
    assert isinstance(Q, int) and sum(cells.values()) == q
    assert all(len(u) == 3 and u.count(0) <= 1 and set(u) <= {-1,0,1}
               and isinstance(size, int) and size > 0 for u, size in cells.items())
    for i in range(3):
        zeros = [u for u in cells if u[i] == 0]
        assert len(zeros) == 1 and cells[zeros[0]] == 1
        assert sum(size*u[i] for u, size in cells.items()) == 0
        for j in range(3):
            assert sum(size*u[i]*u[j] for u, size in cells.items()) == (q-1 if i == j else -1)
    walsh, c = walsh_decomposition(q, n, degree)
    bulk = {u: size for u, size in cells.items() if 0 not in u}
    marked_rows = [next(u for u in cells if u[i] == 0) for i in range(3)]
    assert all(marked_rows[i][j] == marked_rows[j][i] for i in range(3) for j in range(3))
    h = {I: sum(size*char(u, I) for u, size in bulk.items()) for I in SUBSETS}

    def toggle(I, i):
        return tuple(j for j in range(3) if (j in I) != (j == i))

    def known_contraction(I, J):
        # f_I is its Walsh character outside M and is zero on M.
        # S f_empty = -sum_i S e_(m_i).
        # S f_{i} = q e_(m_i)-1-sum_j S_(m_j,m_i) S e_(m_j).
        if len(I) > len(J):
            I, J = J, I
        if not I:
            return -sum(h[toggle(J, i)] for i in range(3))
        assert len(I) == 1
        i = I[0]
        return -h[J]-sum(marked_rows[j][i]*h[toggle(J, j)] for j in range(3))

    cell_only = Fraction()
    for u, nu in cells.items():
        for v, nv in cells.items():
            if 0 in u or 0 in v:
                # A marked row is a singleton. Its sign to every other row
                # is already in that row's mark pattern.
                sign = v[u.index(0)] if 0 in u else u[v.index(0)]
                cell_only += nu*nv*coefficient(q, n, degree, u, v, sign)
                continue
            distinct = nu*nv-(nu if u == v else 0)
            cell_only += distinct*(coefficient(q, n, degree, u, v, 1)
                                  + coefficient(q, n, degree, u, v, -1))/2
            if u == v:
                cell_only += nu*coefficient(q, n, degree, u, u, 0)
    for (I, J), value in walsh.items():
        if min(len(I), len(J)) <= 1:
            cell_only += value*known_contraction(I, J)/2
    second = cell_only+c*Q
    mean = sum((size*row_mean(q,n,degree,u) for u,size in cells.items()), Fraction())
    assert second-mean*mean >= 0
    return {'q': q, 'marks': list(marks), 'n': n, 'degree': degree,
            'cell_only_term': fraction_pair(cell_only), 'Q': Q,
            'Q_coefficient': fraction_pair(c),
            'mean': fraction_pair(mean), 'second_moment': fraction_pair(second),
            'variance': fraction_pair(second-mean*mean),
            'identity': 'E[T_d(C)^2 | M subset C] = cell_only_term + Q_coefficient * Q',
            'Q_definition': 'h2^T S h3; h2(x)=sum_{i<j} S(x,m_i)S(x,m_j), h3(x)=product_i S(x,m_i), both set to zero on M',
            'nonzero_walsh_coefficients': len(walsh),
            'uneliminated_ordered_walsh_terms': sum(min(map(len,key)) >= 2 for key in walsh)}


def reduce_three_marks(record, n, degree=6):
    q, marks, cells, relations = validate_counts(record)
    assert len(marks) == 3
    pairs = tuple(combinations(range(3), 2))
    Q = sum(count*sign*sum(char(u, I) for I in pairs)*char(v, (0,1,2))
            for (u, v, sign), count in relations.items() if 0 not in u and 0 not in v)
    result = reduce_from_cells(q, marks, cells, Q, n, degree)
    original = conditional_moments(record, n, degree)
    for field in ('mean','second_moment','variance'):
        assert result[field] == original[field]
    result['generic_compiler_equality'] = True
    return result


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('counts', type=Path)
    parser.add_argument('--n', type=int, required=True)
    parser.add_argument('--degree', type=int, default=6)
    args = parser.parse_args()
    print(json.dumps(reduce_three_marks(json.loads(args.counts.read_text()), args.n, args.degree), indent=2))


if __name__ == '__main__':
    main()
