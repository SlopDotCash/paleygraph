#!/usr/bin/env python3
"""Exact MCA witnesses for RS[F_p, mu_8, 4], with no numerical arithmetic.

At radius 3/8, a witness has at least five coordinates. The proof that the
56 five-coordinate circuits give an upper bound for *all* word pairs is in
research/prize-reduction-audit.md. This script certifies an attaining pair,
and exhausts all 64 ordered monomial pairs. Polynomial division independently
checks each circuit functional and produces every witnessing cubic.
"""
from itertools import combinations
from math import comb, isqrt
from pathlib import Path
import hashlib
import json


def trim(poly):
    poly = list(poly)
    while len(poly) > 1 and poly[-1] == 0:
        poly.pop()
    return poly


def mul(left, right, p):
    result = [0] * (len(left) + len(right) - 1)
    for i, a in enumerate(left):
        for j, b in enumerate(right):
            result[i + j] = (result[i + j] + a * b) % p
    return trim(result)


def remainder(poly, monic, p):
    assert monic[-1] == 1
    result = trim([c % p for c in poly])
    while len(result) >= len(monic):
        shift, scale = len(result) - len(monic), result[-1]
        for j, c in enumerate(monic):
            result[j + shift] = (result[j + shift] - scale * c) % p
        result = trim(result)
    return result


def evaluate(poly, x, p):
    value = 0
    for coefficient in reversed(poly):
        value = (value * x + coefficient) % p
    return value


def coefficient(poly, j):
    return poly[j] if j < len(poly) else 0


def unit(j):
    return [int(i == j) for i in range(8)]


def dot(left, right, p):
    return sum(a * b for a, b in zip(left, right)) % p


def check_field(p):
    assert p > 2 and all(p % d for d in range(2, isqrt(p) + 1))
    n, k = 8, 4
    generator = next(a for a in range(2, p)
                     if pow(a, n, p) == 1 and pow(a, n // 2, p) != 1)
    domain = sorted({pow(generator, j, p) for j in range(n)})
    assert len(domain) == n
    assert all(a * b % p in domain for a in domain for b in domain)
    assert p > comb(comb(n, k + 1), 2)

    circuits, opposite_triples = [], 0
    for witness in combinations(domain, k + 1):
        weights, vanishing = [], [1]
        for x in witness:
            denominator = 1
            for y in witness:
                if y != x:
                    denominator = denominator * (x - y) % p
            weights.append(pow(denominator, -1, p))
            vanishing = mul(vanishing, [-x % p, 1], p)
        form = [sum(w * pow(x, j, p) for x, w in zip(witness, weights)) % p
                for j in range(n)]
        remainders = [remainder(unit(j), vanishing, p) for j in range(n)]
        assert form == [coefficient(poly, k) for poly in remainders]
        assert form[:k] == [0] * k and form[k] == 1
        complement = [x for x in domain if x not in witness]
        e1 = sum(complement) % p
        e2 = sum(x * y for x, y in combinations(complement, 2)) % p
        e3 = complement[0] * complement[1] * complement[2] % p
        assert form[k:] == [1, -e1 % p, e2, -e3 % p]
        if any((x + y) % p == 0 for x, y in combinations(complement, 2)):
            opposite_triples += 1
            assert all(value in domain for value in form[k:])
        circuits.append((witness, form, vanishing, remainders))
    assert len(circuits) == comb(n, k + 1) == 56
    assert opposite_triples == 24

    # For a monomial pair, each nonzero direction functional gives exactly
    # one scalar. A zero direction functional gives none: if the offset
    # functional also vanishes, the pair is jointly explainable on this set.
    counts, scalar_sets = [], {}
    for a in range(n):
        row = []
        for b in range(n):
            scalars, independently_divided = set(), set()
            for _, form, _, remainders in circuits:
                if form[b]:
                    scalars.add(-form[a] * pow(form[b], -1, p) % p)
                numerator = coefficient(remainders[a], k)
                denominator = coefficient(remainders[b], k)
                if denominator:
                    independently_divided.add(-numerator * pow(denominator, -1, p) % p)
            assert scalars == independently_divided
            row.append(len(scalars))
            scalar_sets[f"{a},{b}"] = sorted(scalars)
        counts.append(row)
    maximum = max(max(row) for row in counts)
    maximizers = [[a, b] for a in range(n) for b in range(n) if counts[a][b] == maximum]
    assert maximum == 40
    assert maximizers == [[4, 5], [5, 4], [6, 7], [7, 6]]

    offset, direction = [0, 0, 0, 0, 0, 1, 2, 3], unit(k)
    certificates, attained = [], set()
    for witness, form, vanishing, _ in circuits:
        assert dot(direction, form, p) == 1
        gamma = -dot(offset, form, p) % p
        assert gamma not in attained
        attained.add(gamma)
        folded = [(a + gamma * b) % p for a, b in zip(offset, direction)]
        cubic = remainder(folded, vanishing, p)
        assert len(cubic) <= k
        assert all(evaluate(folded, x, p) == evaluate(cubic, x, p) for x in witness)
        agreement = [x for x in domain if evaluate(folded, x, p) == evaluate(cubic, x, p)]
        assert agreement == list(witness)
        certificates.append({"witness": list(witness), "gamma": gamma,
                             "cubic_coefficients_low_to_high": cubic,
                             "vanishing_coefficients_low_to_high": vanishing})
    assert len(attained) == 56
    return {"p": p, "n": n, "k": k, "radius": "3/8", "generator": generator,
            "domain": domain, "quartic_window_n4_over4_to_n4": n**4 // 4 <= p <= n**4,
            "offset_coefficients_low_to_high": offset,
            "direction_coefficients_low_to_high": direction,
            "exact_global_mca_error": f"56/{p}",
            "structural_upper_bound": {"opposite_pair_triples": opposite_triples,
                                       "their_possible_scalars_per_monomial_pair": 8,
                                       "remaining_circuits": 32, "total_upper_bound": 40},
            "monomial_maximum": maximum, "monomial_maximizers": maximizers,
            "monomial_count_matrix_rows_offset_columns_direction": counts,
            "monomial_bad_scalar_sets": scalar_sets,
            "attaining_bad_scalars": sorted(attained), "certificates": certificates}


def main():
    source = Path(__file__).resolve()
    result = {"source_sha256": hashlib.sha256(source.read_bytes()).hexdigest(),
              "proof_scope": "Exact field computations plus the circuit theorem in the research note; not Lean.",
              "fields": [check_field(p) for p in [2017, 65537]]}
    destination = source.parents[1] / "results" / "mca_circuit_certificate.json"
    destination.write_text(json.dumps(result, indent=2) + "\n")
    for field in result["fields"]:
        print(f"p={field['p']}: all 56 circuits and 64 monomial pairs checked; "
              f"global MCA=56/{field['p']}, monomial maximum=40/{field['p']}")
    print(destination)


if __name__ == "__main__":
    main()
