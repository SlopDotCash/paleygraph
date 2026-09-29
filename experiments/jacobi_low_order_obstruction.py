#!/usr/bin/env python3
"""Exact quartic-window obstruction to beta_j <= n*j with constant one.

Standard-library integer/rational arithmetic only. This does not refute a
bound with an unspecified absolute constant, or any Paley conjecture.
"""
from collections import Counter, defaultdict
from fractions import Fraction as F
from hashlib import sha256
from math import factorial, isqrt
from pathlib import Path
import json


def inner(left, right, moments):
    return sum((a*b*moments[i+j] for i, a in enumerate(left)
                for j, b in enumerate(right)), F(0))


def certificate():
    p, n, generator = 67403009, 128, 64701253
    assert all(p % d for d in range(2, isqrt(p)+1))
    assert (p-1) % n == 0 and n**4 <= 4*p <= 4*n**4
    m = (p-1)//n
    assert pow(generator, n, p) == 1
    assert pow(generator, n//2, p) == p-1
    subgroup = sorted({pow(generator, j, p) for j in range(n)})
    assert len(subgroup) == n and {p-x for x in subgroup} == set(subgroup)

    # x^n labels a nonzero multiplicative coset of H without enumerating F_p.
    kernel = Counter(pow((1+x) % p, n, p)
                     for x in subgroup if x != p-1)
    assert sum(kernel.values()) == n-1
    kernel_square_sum = sum(value**2 for value in kernel.values())
    pair_sums = Counter((a+b) % p for a in subgroup for b in subgroup)
    energy = sum(value**2 for value in pair_sums.values())
    zero_three = sum(pair_sums[(-x) % p] for x in subgroup)
    zero_counts = [1, 0, n, zero_three, energy]
    moments = [F(p*count-n**j, n) for j, count in enumerate(zero_counts)]
    assert all(value.denominator == 1 for value in moments)
    assert moments == [m, -1, p-n, p*kernel[1]-n*n,
                       p*(n+kernel_square_sum)-n**3]
    assert energy == n*n+n*kernel_square_sum

    # First route: central moment identities.
    raw = [value/m for value in moments]
    mean = raw[1]
    beta_one = raw[2]-mean**2
    alpha_one = (raw[3]-2*mean*raw[2]+mean**3)/beta_one
    central_four = raw[4]-4*mean*raw[3]+6*mean**2*raw[2]-3*mean**4
    beta_two = central_four/beta_one-beta_one-(alpha_one-mean)**2

    # Independent route: direct Gram-Schmidt of 1, X, X^2, then norm ratios.
    polys, norms = [], []
    for degree in range(3):
        poly = [F(0)]*degree+[F(1)]
        for previous, norm in zip(polys, norms):
            coefficient = inner(poly, previous, moments)/norm
            for j, value in enumerate(previous):
                poly[j] -= coefficient*value
        norm = inner(poly, poly, moments)
        assert norm > 0
        polys.append(poly)
        norms.append(norm)
    for i in range(3):
        for j in range(3):
            assert inner(polys[i], polys[j], moments) == (norms[i] if i == j else 0)
    assert norms[1]/norms[0] == beta_one
    assert norms[2]/norms[1] == beta_two
    assert inner([F(0)]+polys[1], polys[1], moments)/norms[1] == alpha_one
    assert raw[4] == ((mean**2+beta_one)**2
                      + beta_one*(mean+alpha_one)**2+beta_one*beta_two)
    assert beta_two == F(76801470694776, 277291762225) > 2*n

    # Enumerate every unordered zero quadruple from pair-sum buckets.
    pairs = defaultdict(list)
    for i, a in enumerate(subgroup):
        for b in subgroup[i:]:
            pairs[(a+b) % p].append((a, b))
    quadruples = {tuple(sorted((a, b, c, d)))
                  for total, bucket in pairs.items() for a, b in bucket
                  for c, d in pairs.get((-total) % p, [])}
    nondegenerate = {q for q in quadruples if not any(
        (q[i]+q[j]) % p == 0 for i in range(4) for j in range(i))}
    degenerate = quadruples-nondegenerate

    def ordered_count(q):
        answer = factorial(4)
        for multiplicity in Counter(q).values():
            answer //= factorial(multiplicity)
        return answer

    assert sum(map(ordered_count, quadruples)) == energy
    assert sum(map(ordered_count, degenerate)) == 3*n*n-3*n
    representative = (1, 1074964, 1550267, 64777777)
    assert sum(representative) == p
    assert all(value in subgroup for value in representative)
    orbit = {tuple(sorted(h*value % p for value in representative))
             for h in subgroup}
    assert nondegenerate == orbit and len(orbit) == n
    assert all(len(set(q)) == 4 for q in orbit)
    assert energy == 3*n*n-3*n+24*n == 51840
    powers = {pow(generator, j, p): j for j in range(n)}

    return {
        "status": "passed",
        "scope": "B=1 fails at j=2; unspecified constants and Paley remain open",
        "p": p, "n": n, "generator": generator, "coset_count": m,
        "quartic_window": "n^4/4 <= p <= n^4",
        "primality_method": "trial division through floor(sqrt(p))",
        "subgroup": subgroup,
        "kernel_by_nth_power_coset_label": dict(sorted(kernel.items())),
        "nonzero_kernel_histogram": dict(sorted(Counter(kernel.values()).items())),
        "kernel_at_H": kernel[1], "kernel_square_sum": kernel_square_sum,
        "additive_energy": energy, "zero_counts_orders_0_through_4": zero_counts,
        "signed_coset_moments": [int(value) for value in moments],
        "alpha_0": str(mean), "alpha_1": str(alpha_one),
        "beta_1": str(beta_one), "beta_2": str(beta_two),
        "positive_beta_2_minus_2n": str(beta_two-2*n),
        "monic_orthogonal_polynomials_low_to_high": [[str(x) for x in poly] for poly in polys],
        "squared_norms": [str(x) for x in norms],
        "nondegenerate_zero_quadruple_representative": representative,
        "representative_generator_exponents": [powers[x] for x in representative],
        "nondegenerate_unordered_quadruples": len(orbit),
        "nondegenerate_ordered_quadruples": 24*n,
        "degenerate_ordered_quadruples": 3*n*n-3*n,
        "all_nondegenerate_quadruples_are_in_one_H_orbit": True,
        "energy_independently_verified_by_pair_sums_and_quadruples": True,
        "beta_2_independently_verified_by_Gram_Schmidt": True,
    }


if __name__ == "__main__":
    source = Path(__file__).resolve()
    root = source.parents[1]
    result = certificate()
    result["source_sha256"] = sha256(source.read_bytes()).hexdigest()
    output = root/"results/jacobi_low_order_obstruction.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({key: result[key] for key in [
        "status", "p", "n", "additive_energy", "beta_2", "positive_beta_2_minus_2n"]}, indent=2))
