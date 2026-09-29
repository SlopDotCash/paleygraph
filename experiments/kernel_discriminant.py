#!/usr/bin/env python3
"""Exact polynomial and prime-power certificates for fourth-energy excess.

The polynomial/discriminant identity is proved in the accompanying note.
This script does not prove a uniform energy or Paley bound.
"""
from collections import Counter
from hashlib import sha256
from math import comb, isqrt, prod
from pathlib import Path
import json

from cyclotomic_norm_audit import determinant

ROOT = Path(__file__).resolve().parents[1]


def prime(p):
    return p >= 2 and (p == 2 or (p % 2 and all(p % d for d in range(3, isqrt(p)+1, 2))))


def polynomial(n):
    degree = n//2-1
    powers = [0]
    for k in range(1, degree+1):
        twice = n*sum(comb(n*k, n*j) for j in range(k+1)) - 2**(n*k)
        assert twice % 2 == 0
        powers.append(twice//2)
    elementary = [1]
    for k in range(1, degree+1):
        numerator = sum((-1)**(i-1)*elementary[k-i]*powers[i] for i in range(1, k+1))
        assert numerator % k == 0
        elementary.append(numerator//k)
    return [(-1)**(degree-i)*elementary[degree-i] for i in range(degree+1)]


def evaluate(coefficients, value, modulus=None):
    result = 0
    for coefficient in reversed(coefficients):
        result = result*value+coefficient
        if modulus:
            result %= modulus
    return result


def discriminant(coefficients):
    degree = len(coefficients)-1
    # Multiplication by P' in Z[Y]/P. Its determinant is Res(P,P').
    matrix = [[0]*degree for _ in range(degree)]
    for j in range(degree):
        column = [0]*j + [i*coefficients[i] for i in range(1, degree+1)]
        for k in range(len(column)-1, degree-1, -1):
            leading = column[k]
            for i, coefficient in enumerate(coefficients):
                column[k-degree+i] -= leading*coefficient
        for i in range(degree):
            matrix[i][j] = column[i] if i < len(column) else 0
    return (-1)**(degree*(degree-1)//2)*determinant(matrix)


def factor_complete(value, trial_limit=1_000_000):
    sieve = bytearray(b'\x01')*(trial_limit+1)
    sieve[:2] = b'\x00\x00'
    for p in range(2, isqrt(trial_limit)+1):
        if sieve[p]:
            sieve[p*p::p] = b'\x00'*((trial_limit-p*p)//p+1)
    factors, remaining = Counter(), value
    for p in range(2, trial_limit+1):
        if not sieve[p]:
            continue
        while remaining % p == 0:
            factors[p] += 1
            remaining //= p
        if remaining == 1:
            break
    if remaining > 1:
        # Never label an unfactored large composite as a prime certificate.
        assert remaining < trial_limit**2 and prime(remaining)
        factors[remaining] += 1
    assert prod(p**v for p, v in factors.items()) == value
    return factors


def generator(p, n):
    return next(pow(a, (p-1)//n, p) for a in range(2, p)
                if pow(pow(a, (p-1)//n, p), n//2, p) == p-1)


def root_product(roots, modulus):
    coefficients = [1]
    for root in roots:
        following = [0]*(len(coefficients)+1)
        for i, c in enumerate(coefficients):
            following[i] = (following[i]-root*c) % modulus
            following[i+1] = (following[i+1]+c) % modulus
        coefficients = following
    return coefficients


def valuation(value, p):
    assert value != 0
    v = 0
    while value % p == 0:
        value //= p
        v += 1
    return v


def lift_audit(p, n, g=None, coefficients=None, expected_valuation=None):
    assert prime(p) and (p-1) % n == 0
    if g is None:
        g = generator(p, n)
    assert pow(g, n, p) == 1 and pow(g, n//2, p) == p-1
    original = g
    h = [pow(g, j, p) for j in range(n)]
    pair_counts = Counter((a+b) % p for a in h for b in h)
    energy = sum(c*c for c in pair_counts.values())
    excess = (energy-(3*n*n-3*n))//n
    assert energy == 3*n*n-3*n+n*excess
    levels, modulus = [], p
    for precision in range(1, 9):
        if precision > 1:
            previous = modulus
            modulus *= p
            residual = pow(g, n, modulus)-1
            assert residual % previous == 0
            correction = -(residual//previous)*pow(n*pow(g, n-1, p), -1, p) % p
            g += correction*previous
        assert pow(g, n, modulus) == 1 and pow(g, n//2, modulus) == modulus-1
        roots = [pow(1+pow(g, j, modulus), n, modulus) for j in range(1, n//2)]
        boundary = pow(2, n, modulus)
        counts = Counter({boundary: 1})
        for r in roots:
            counts[r] += 2
        assert sum(counts.values()) == n-1 and 0 not in counts
        d = sum(c*c for c in counts.values())-(2*n-3)
        assert d >= 0 and d % 4 == 0
        lifted_group = [pow(g, j, modulus) for j in range(n)]
        lifted_pairs = Counter((a+b) % modulus for a in lifted_group for b in lifted_group)
        lifted_energy = sum(c*c for c in lifted_pairs.values())
        assert lifted_energy == 3*n*n-3*n+n*d
        if precision == 1:
            assert d == excess
        if coefficients is not None:
            assert root_product(roots, modulus) == [c % modulus for c in coefficients]
        levels.append({"precision": precision, "modulus": modulus,
                       "excess_D": d, "direct_additive_energy": lifted_energy,
                       "nonzero_multiplicity_histogram": dict(sorted(Counter(counts.values()).items()))})
        if d == 0:
            break
    assert levels[-1]["excess_D"] == 0, "Increase precision before claiming an exact valuation"
    disc_val = 2*sum(valuation(roots[i]-roots[j], p)
                     for i in range(len(roots)) for j in range(i))
    boundary_val = sum(valuation(boundary-r, p) for r in roots)
    total_val = disc_val+boundary_val
    assert 4*total_val == sum(row["excess_D"] for row in levels)
    if expected_valuation is not None:
        assert total_val == expected_valuation
    return {"p": p, "n": n, "generator": original, "energy2": energy,
            "D": excess, "prime_power_levels": levels,
            "discriminant_valuation": disc_val, "boundary_value_valuation": boundary_val,
            "valuation_of_A": total_val, "valuation_identity_checked": True}


def small_order(n):
    coefficients = polynomial(n)
    disc = discriminant(coefficients)
    boundary = evaluate(coefficients, 2**n)
    a = disc*boundary
    assert disc > 0 and boundary > 0
    factors = factor_complete(a)
    exceptional = [lift_audit(p, n, coefficients=coefficients, expected_valuation=v)
                   for p, v in sorted(factors.items()) if p % n == 1]
    assert all(row["D"] > 0 for row in exceptional)
    quartic = [row["p"] for row in exceptional if n**4 <= 4*row["p"] <= 4*n**4]
    assert not quartic
    return {"n": n, "P_coefficients_ascending": coefficients,
            "discriminant_hex": hex(disc), "P_at_two_to_n_hex": hex(boundary),
            "A_hex": hex(a), "A_bits": a.bit_length(),
            "complete_A_prime_factorization": dict(sorted(factors.items())),
            "all_eligible_exceptional_primes": exceptional,
            "quartic_window_exceptional_primes": quartic,
            "scope": "Complete fixed-order fourth-energy classification; no spectral maximum bound."}


def main():
    small = [small_order(n) for n in [4, 8, 16, 32]]
    old_path = ROOT/"results/cyclotomic_norm_audit.json"
    old = json.loads(old_path.read_text())
    helper_hash = sha256((ROOT/"experiments/cyclotomic_norm_audit.py").read_bytes()).hexdigest()
    assert old["source_sha256"] == helper_hash
    for fresh, previous in zip(small[:3], old["cases"]):
        assert fresh["n"] == previous["n"]
        one = {x["p"]: x["energy2"] for x in fresh["all_eligible_exceptional_primes"]}
        two = {x["p"]: x["energy"] for x in previous["splitting_primes_with_extra_relations"]}
        assert one == two
    assert small[0]["P_coefficients_ascending"] == [4,1]
    assert small[1]["P_coefficients_ascending"] == [-256,-2160,120,1]
    assert max(x["p"] for x in small[3]["all_eligible_exceptional_primes"]) == 21523361
    larger = [lift_audit(6700417,64,2), lift_audit(67403009,128,64701253),
              lift_audit(1073748737,256,1064280392),
              lift_audit(17179869697,512,13395504394)]
    assert [x["D"] for x in larger] == [12,24,0,0]
    result = {"status": "passed; uniform conjecture unproved",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "determinant_helper_sha256": helper_hash,
              "prior_norm_result_sha256": sha256(old_path.read_bytes()).hexdigest(),
              "small_orders": small, "larger_prime_power_checks": larger}
    path = ROOT/"results/kernel_discriminant.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "classified_orders": [{"n": x["n"], "exceptional_prime_count": len(x["all_eligible_exceptional_primes"]),
                                             "quartic_exceptions": x["quartic_window_exceptional_primes"]} for x in small],
                      "larger_checks": [{"p": x["p"], "n": x["n"], "D": x["D"], "v_p_A": x["valuation_of_A"]} for x in larger],
                      "output": str(path)}, indent=2))


if __name__ == "__main__":
    main()
