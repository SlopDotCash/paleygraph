#!/usr/bin/env python3
"""Exact small checks of the cyclotomic norm bound for extra zero relations."""
from collections import Counter
from hashlib import sha256
from itertools import product
from pathlib import Path
import json


def determinant(matrix):
    """Fraction-free elimination with checked exact divisions."""
    matrix = [row[:] for row in matrix]
    degree, sign, previous = len(matrix), 1, 1
    for j in range(degree-1):
        if not matrix[j][j]:
            pivot_row = next((k for k in range(j+1, degree) if matrix[k][j]), None)
            if pivot_row is None:
                return 0
            matrix[j], matrix[pivot_row] = matrix[pivot_row], matrix[j]
            sign = -sign
        pivot = matrix[j][j]
        for k in range(j+1, degree):
            for ell in range(j+1, degree):
                value = pivot*matrix[k][ell]-matrix[k][j]*matrix[j][ell]
                assert value % previous == 0
                matrix[k][ell] = value//previous
        for k in range(j+1, degree):
            matrix[k][j] = 0
        previous = pivot
    return sign*matrix[-1][-1]


def multiplication_matrix(coefficients):
    degree = len(coefficients)
    matrix = [[0]*degree for _ in range(degree)]
    for j in range(degree):
        for i, value in enumerate(coefficients):
            matrix[(i+j) % degree][j] += value*(-1 if i+j >= degree else 1)
    return matrix


def factor(value):
    factors = Counter()
    divisor = 2
    while divisor*divisor <= value:
        while value % divisor == 0:
            factors[divisor] += 1
            value //= divisor
        divisor += 1
    if value > 1:
        factors[value] += 1
    return factors


def rank_mod(matrix, p):
    matrix = [[x % p for x in row] for row in matrix]
    rank = 0
    for column in range(len(matrix[0])):
        pivot = next((i for i in range(rank, len(matrix)) if matrix[i][column]), None)
        if pivot is None:
            continue
        matrix[rank], matrix[pivot] = matrix[pivot], matrix[rank]
        inverse = pow(matrix[rank][column], -1, p)
        matrix[rank] = [x*inverse % p for x in matrix[rank]]
        for i in range(len(matrix)):
            if i != rank:
                multiple = matrix[i][column]
                matrix[i] = [(x-multiple*y) % p
                             for x, y in zip(matrix[i], matrix[rank])]
        rank += 1
    return rank


def audit(n):
    k, degree = 4, n//2
    coefficient_counts = Counter()
    for exponents in product(range(n), repeat=k-1):
        coefficients = [0]*degree
        for j in (0,)+exponents:
            coefficients[j % degree] += 1 if j < degree else -1
        coefficient_counts[tuple(coefficients)] += 1
    intrinsic_normalized = coefficient_counts[(0,)*degree]
    intrinsic = n*intrinsic_normalized
    assert intrinsic == 3*n*n-3*n
    nonzero_polynomials = []
    norm_histogram, valuations = Counter(), Counter()
    norm_product = 1
    for coefficients, multiplicity in coefficient_counts.items():
        matrix = multiplication_matrix(coefficients)
        norm = abs(determinant(matrix))
        assert (norm == 0) == (not any(coefficients))
        if not norm:
            continue
        assert norm <= k**degree
        norm_histogram[norm] += multiplicity
        norm_product *= norm**multiplicity
        for p, valuation in factor(norm).items():
            valuations[p] += multiplicity*valuation
        nonzero_polynomials.append((coefficients, matrix, multiplicity))
    nonzero_count = n**(k-1)-intrinsic_normalized
    assert sum(norm_histogram.values()) == nonzero_count
    assert norm_product <= k**(degree*nonzero_count)
    # Squared, denominator-cleared AM-GM refinement: no numerical logs.
    amgm_exponent = degree*nonzero_count
    assert norm_product**2*nonzero_count**amgm_exponent <= (k*n**(k-1))**amgm_exponent

    splitting_primes = []
    required_product = 1
    for p in sorted(valuations):
        if p % n != 1:
            continue
        seed = 2
        while pow(pow(seed, (p-1)//n, p), n//2, p) == 1:
            seed += 1
        generator = pow(seed, (p-1)//n, p)
        roots = [pow(generator, j, p) for j in range(1, n, 2)]
        assert len(set(roots)) == degree
        nullity_sum = 0
        for coefficients, matrix, multiplicity in nonzero_polynomials:
            root_count = sum(sum(c*pow(root, j, p) for j, c in enumerate(coefficients)) % p == 0
                             for root in roots)
            assert root_count == degree-rank_mod(matrix, p)
            nullity_sum += multiplicity*root_count
        subgroup = [pow(generator, j, p) for j in range(n)]
        pairs = Counter((a+b) % p for a in subgroup for b in subgroup)
        energy = sum(value**2 for value in pairs.values())
        extra = energy-intrinsic
        assert degree*extra == n*nullity_sum
        assert nullity_sum <= valuations[p]
        required_product *= p**nullity_sum
        splitting_primes.append({
            "p": p, "generator": generator, "energy": energy,
            "extra_zero_quadruples": extra,
            "nullity_sum": nullity_sum, "valuation_of_norm_product": valuations[p]})
    assert norm_product % required_product == 0
    return {
        "n": n, "k": k, "degree": degree,
        "normalized_tuples_exhausted": n**(k-1),
        "intrinsic_zero_quadruples": intrinsic,
        "nonzero_norm_histogram": dict(sorted(norm_histogram.items())),
        "norm_product_prime_factorization": dict(sorted(valuations.items())),
        "splitting_primes_with_extra_relations": splitting_primes,
        "matrix_nullities_match_all_primitive_root_counts": True,
        "required_prime_product_divides_norm_product": True,
        "norm_product_at_most_4_to_degree_times_nonzero_tuple_count": True,
        "sharper_AM_GM_bound_verified_with_integers": True,
        "scope": "Complete k=4 audit for this n; no uniform spectral bound"}


if __name__ == "__main__":
    source = Path(__file__).resolve()
    result = {"status": "passed", "source_sha256": sha256(source.read_bytes()).hexdigest(),
              "cases": [audit(n) for n in [4, 8, 16]]}
    output = source.parents[1]/"results/cyclotomic_norm_audit.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "cases": [
        {"n": case["n"], "normalized_tuples_exhausted": case["normalized_tuples_exhausted"],
         "splitting_primes_with_extra_relations": case["splitting_primes_with_extra_relations"]}
        for case in result["cases"]]}, indent=2))
