#!/usr/bin/env python3
"""Exact symbolic fixed-size fourth-moment mean and variance audit."""
from collections import Counter, defaultdict
from functools import lru_cache
from itertools import combinations
from fractions import Fraction
from math import comb, isqrt
from pathlib import Path
from hashlib import sha256
import json

ROOT = Path(__file__).resolve().parents[1]


def add(a, b):
    result = dict(a)
    for key, value in b.items():
        result[key] = result.get(key, 0) + value
        if not result[key]:
            del result[key]
    return result


def scale(a, c):
    return {key: value*c for key, value in a.items() if value*c}


def mul(a, b):
    result = defaultdict(int)
    for (i, j), av in a.items():
        for (k, l), bv in b.items():
            result[i+k, j+l] += av*bv
    return {key: value for key, value in result.items() if value}


def const(n):
    return {(0, 0): n} if n else {}


P = {(1, 0): 1}
N = {(0, 1): 1}


def falling(a, k):
    result = const(1)
    for j in range(k):
        result = mul(result, add(a, const(-j)))
    return result


def divide_p_linear(a, root):
    """Divide by p-root(n); return quotient and remainder."""
    degree = max((i for i, j in a), default=-1)
    carry = {}
    quotient = {}
    for i in range(degree, -1, -1):
        coefficient = {(0, j): v for (k, j), v in a.items() if k == i}
        carry = add(coefficient, mul(carry, root))
        if i:
            quotient.update({(i-1, j): v for (k, j), v in carry.items()})
    return quotient, carry


def divide_n_linear(a, root):
    swapped = {(j, i): v for (i, j), v in a.items()}
    quotient, remainder = divide_p_linear(swapped, const(root))
    return {(j, i): v for (i, j), v in quotient.items()}, remainder


def partition_types(nx, ny):
    counts = Counter()
    slots = [(1, 0)]*nx+[(0, 1)]*ny
    def recurse(i, blocks):
        if i == len(slots):
            counts[tuple(sorted(blocks))] += 1
            return
        a, b = slots[i]
        for j, (x, y) in enumerate(blocks):
            next_blocks = list(blocks)
            next_blocks[j] = x+a, y+b
            recurse(i+1, next_blocks)
        recurse(i+1, blocks+[(a, b)])
    recurse(0, [])
    return counts


def moment_coefficients(nx, ny, epsilon):
    def population(a, b):
        if not b:
            return add(P, const(-1)) if not a % 2 else {}
        if not a:
            return add(P, const(-1)) if not b % 2 else {}
        if not a % 2 and not b % 2:
            return add(P, const(-2))
        if a % 2 and b % 2:
            return const(-1)
        return const(-1 if a % 2 else -epsilon)

    @lru_cache(None)
    def distinct(types):
        if not types:
            return const(1)
        head, tail = types[0], types[1:]
        result = mul(population(*head), distinct(tail))
        for j, other in enumerate(tail):
            merged = list(tail)
            merged[j] = head[0]+other[0], head[1]+other[1]
            result = add(result, scale(distinct(tuple(sorted(merged))), -1))
        return result

    output = defaultdict(dict)
    types = partition_types(nx, ny)
    for blocks, count in types.items():
        output[len(blocks)] = add(output[len(blocks)], scale(distinct(blocks), count))
    return dict(output)


def common_numerator(coefficients, order=8):
    output = {}
    for k, coefficient in coefficients.items():
        term = mul(coefficient, falling(N, k))
        for j in range(k, order):
            term = mul(term, add(P, const(-j)))
        output = add(output, term)
    return output


def symbolic_variance(epsilon):
    single = common_numerator(moment_coefficients(8, 0, epsilon))
    joint = common_numerator(moment_coefficients(4, 4, epsilon))
    square_numerator = add(mul(P, single), mul(falling(P, 2), joint))
    v = mul(N, add(P, scale(N, -1)))
    mean_numerator = mul(v, add(scale(v, 3), add(scale(P, -2), const(1))))
    denominator = falling(P, 8)
    quotient, rem = divide_p_linear(denominator, const(2))
    assert not rem
    numerator = add(mul(square_numerator, add(P, const(-2))),
                    scale(mul(mul(mean_numerator, mean_numerator), quotient), -1))
    raw_numerator = numerator
    common_denominator = mul(denominator, add(P, const(-2)))
    v = mul(N, add(P, scale(N, -1)))
    shared = scale(mul(falling(N, 3), falling(add(P, scale(N, -1)), 3)), 24)
    if epsilon == 1:
        bracket = add(mul(mul(add(P, const(-5)), add(P, const(3))), v),
                      scale({(3, 0): 1, (2, 0): -6, (1, 0): 4, (0, 0): 3}, -3))
        expected_numerator = mul(mul(shared, add(P, const(-1))), bracket)
        expected_denominator = const(1)
        for j in (2, 2, 3, 4, 6, 7):
            expected_denominator = mul(expected_denominator, add(P, const(-j)))
    else:
        bracket = add(mul(add(P, const(1)), v),
                      scale({(2, 0): 1, (1, 0): -3, (0, 0): 3}, -3))
        expected_numerator = mul(mul(shared, add(P, const(1))), bracket)
        expected_denominator = const(1)
        for j in (2, 2, 4, 5, 6):
            expected_denominator = mul(expected_denominator, add(P, const(-j)))
    residual = add(mul(raw_numerator, expected_denominator),
                   scale(mul(expected_numerator, common_denominator), -1))
    assert not residual
    mean_raw = mul(P, common_numerator(moment_coefficients(4, 0, epsilon), 4))
    assert not add(mul(mean_raw, add(P, const(-2))),
                   scale(mul(mean_numerator, falling(P, 4)), -1))
    factors = []
    for j in [0, 1, 2, 2, 3, 4, 5, 6, 7]:
        q, rem = divide_p_linear(numerator, const(j))
        if not rem:
            numerator = q
            factors.append(('cancelled_denominator', j))
    for j in range(5):
        q, rem = divide_n_linear(numerator, j)
        if not rem:
            numerator = q
            factors.append(('n', j))
        q, rem = divide_p_linear(numerator, add(N, const(j)))
        if not rem:
            numerator = q
            factors.append(('p-n', j))
    return dict(epsilon=epsilon, factors=factors, variance_polynomial_residual={},
                mean_polynomial_residual={}, partitions_of_eight=sum(partition_types(4, 4).values()),
                raw_numerator_monomials=len(raw_numerator),
                remaining_numerator=[(i, j, v) for (i, j), v in sorted(numerator.items())])


def mean_variance(p, n):
    v = n*(p-n)
    mean = Fraction(v*(3*v-2*p+1), p-2)
    f = n*(n-1)*(n-2)*(p-n)*(p-n-1)*(p-n-2)
    if p % 4 == 1:
        numerator = 24*f*(p-1)*((p-5)*(p+3)*v-3*(p**3-6*p*p+4*p+3))
        denominator = (p-2)**2*(p-3)*(p-4)*(p-6)*(p-7)
    else:
        numerator = 24*f*(p+1)*((p+1)*v-3*(p*p-3*p+3))
        denominator = (p-2)**2*(p-4)*(p-5)*(p-6)
    return mean, Fraction(numerator, denominator)


def character(p):
    assert p >= 3 and all(p % d for d in range(2, isqrt(p)+1))
    chi = [-1]*p
    chi[0] = 0
    for x in range(1, p):
        chi[x*x % p] = 1
    return chi


def fraction_pair(q):
    return [q.numerator, q.denominator]


def exact_field(p, max_n):
    chi = character(p)
    records = []
    for n in range(max_n+1):
        total = square = bad = subsets = 0
        for B in combinations(range(p), n):
            F = [sum(chi[(a-b) % p] for b in B) for a in range(p)]
            M4 = sum(x**4 for x in F)
            assert sum(x*x for x in F) == p*n-n*n
            total += M4
            square += M4*M4
            bad += M4 > 3*p*n*n
            subsets += 1
        assert subsets == comb(p, n)
        mean, variance = mean_variance(p, n)
        assert Fraction(total, subsets) == mean
        assert Fraction(square, subsets)-mean*mean == variance
        if n:
            gap = 3*p*n*n-mean
            assert gap >= p*n
            assert Fraction(bad, subsets)*gap*gap <= variance
        records.append(dict(n=n, subsets=subsets, mean=fraction_pair(mean),
                            variance=fraction_pair(variance), above_gaussian_three=bad))
    return dict(p=p, fixed_size_classes=records)


def exceptional_slice():
    p = 1009
    B = (23, 44, 288, 319, 336, 393, 408, 485, 511, 599)
    U = (1008, 997, 982, 736, 916, 535, 826, 312)
    chi = character(p)
    n = len(B)
    assert n**3 <= p < (n+1)**3
    F = [sum(chi[(a-b) % p] for b in B) for a in range(p)]
    assert all(F[a] == n for a in U)
    M4 = sum(x**4 for x in F)
    assert M4 == 316554 > 3*p*n*n == 302700
    A = sum(F[a]**2 for a in B)
    R = p*(3*n*n-2*n)-6*n**3+14*n*n-9*n-6*A+Fraction(3*n*(n-1)*(n-2)*(n-3), p-2)
    D = M4-R
    mean, variance = mean_variance(p, n)
    expected_A = Fraction(n*(n-1)*(p-n), p-2)
    assert D == M4-mean+6*(A-expected_A)
    assert abs(D) <= abs(M4-mean)+6*n**3
    return dict(p=p, B=B, U=U, n=n, M4=M4, gaussian_three=3*p*n*n,
                internal_square_sum=A, R=fraction_pair(R), D=fraction_pair(D),
                mean=fraction_pair(mean), variance=fraction_pair(variance))


def remainder_identity():
    v = mul(N, add(P, scale(N, -1)))
    mean_numerator = mul(v, add(scale(v, 3), add(scale(P, -2), const(1))))
    expected_A_numerator = mul(falling(N, 2), add(P, scale(N, -1)))
    polynomial_R = {(1, 2): 3, (1, 1): -2, (0, 3): -6, (0, 2): 14, (0, 1): -9}
    mean_R_numerator = add(mul(polynomial_R, add(P, const(-2))),
                           add(scale(expected_A_numerator, -6), scale(falling(N, 4), 3)))
    assert mean_numerator == mean_R_numerator
    return dict(expected_remainder_polynomial_residual={})


def main():
    symbolic = [symbolic_variance(e) for e in (1, -1)]
    fields = [exact_field(p, max_n) for p, max_n in ((11, 11), (13, 13), (17, 5), (19, 4))]
    witness = exceptional_slice()
    remainder = remainder_identity()
    files = ('research/parallel6-classical-2026-09-04.md',
             'experiments/parallel6_classical_2026_09_04.py')
    result = dict(status='Exact variance and typical-set SS fourth-moment bound; no uniform SS or Paley bound.',
                  arithmetic='Integer polynomial identities, exact subset enumeration, rational mean and variance.',
                  symbolic_certificates=symbolic, fields=fields,
                  precise_slice_exception=witness, remainder_certificate=remainder,
                  source_sha256={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    destination = ROOT/'results/parallel6_classical_2026_09_04.json'
    destination.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(symbolic_variance_identities=len(symbolic),
                          fixed_size_mean_variance_checks=sum(len(f['fixed_size_classes']) for f in fields),
                          enumerated_subsets=sum(c['subsets'] for f in fields for c in f['fixed_size_classes']),
                          precise_slice_exception_M4=witness['M4'],
                          result=str(destination))))


if __name__ == '__main__':
    main()
