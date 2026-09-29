#!/usr/bin/env python3
"""Exact checks of the fourth-energy orbit and odd-cube identities.

The general proof is in research/quadruple-orbits-and-cube.md.
This is not a uniform energy or Paley proof.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

from kernel_discriminant import discriminant, evaluate, polynomial, prime, valuation

ROOT = Path(__file__).resolve().parents[1]


class Ring:
    """Z/p^r, or its unramified quadratic extension with T^2=d."""
    def __init__(self, p, r, d=None):
        self.p, self.r, self.modulus, self.d = p, r, p**r, d
        self.zero, self.one = (0, 0), (1, 0)
        assert prime(p) and p % 2
        if d is not None:
            assert pow(d, (p-1)//2, p) == p-1
        self.q = p if d is None else p*p

    def scalar(self, a):
        return (a % self.modulus, 0)

    def add(self, x, y):
        return ((x[0]+y[0]) % self.modulus, (x[1]+y[1]) % self.modulus)

    def neg(self, x):
        return ((-x[0]) % self.modulus, (-x[1]) % self.modulus)

    def sub(self, x, y):
        return self.add(x, self.neg(y))

    def mul(self, x, y):
        d = self.d or 0
        return ((x[0]*y[0]+d*x[1]*y[1]) % self.modulus,
                (x[0]*y[1]+x[1]*y[0]) % self.modulus)

    def power(self, x, k):
        result = self.one
        while k:
            if k % 2:
                result = self.mul(result, x)
            x = self.mul(x, x)
            k //= 2
        return result

    def val(self, x):
        assert x != self.zero
        return min(valuation(a, self.p) for a in x if a)


def orbit_counts(ring, group):
    n = len(group)
    pairs = defaultdict(list)
    for i, a in enumerate(group):
        for j in range(i, n):
            pairs[ring.add(a, group[j])].append((i, j))
    multisets = set()
    for total, left in pairs.items():
        opposite = ring.neg(total)
        if total > opposite:
            continue
        for a in left:
            for b in pairs.get(opposite, []):
                indices = tuple(sorted(a+b))
                # An opposite pair forces the other pair to be opposite.
                if any((i+n//2) % n in indices for i in indices):
                    continue
                multisets.add(indices)
    by_type, representatives = Counter(), {3: set(), 2: set(), 1: set()}
    for indices in multisets:
        kind = max(Counter(indices).values())
        assert kind in representatives
        by_type[kind] += 1
        normal = min(tuple(sorted((i-a) % n for i in indices)) for a in indices)
        representatives[kind].add(normal)
    # The proof gives a free H action on each nontrivial multiset type.
    assert all(by_type[k] == n*len(representatives[k]) for k in representatives)
    epsilon = int(ring.power(ring.scalar(3), n) == ring.one)
    assert len(representatives[3]) == epsilon
    u, v = len(representatives[2]), len(representatives[1])
    pair_at_two = sum(ring.sub(ring.scalar(-2), x) in group for x in group)
    assert pair_at_two == 1+2*epsilon+2*u
    return {"epsilon_three_in_H": epsilon, "two_one_one_orbits": u,
            "four_distinct_orbits": v,
            "pair_count_at_minus_two": pair_at_two,
            "D_from_orbits": 4*epsilon+12*u+24*v,
            "orbit_representatives": {str(k): sorted(values) for k, values in representatives.items()}}


def find_generator(ring, n):
    assert (ring.q-1) % n == 0
    for b in range(ring.p if ring.d is not None else 1):
        for a in range(ring.p):
            candidate = ring.power((a, b), (ring.q-1)//n)
            if ring.power(candidate, n//2) == ring.scalar(-1):
                return candidate
    raise AssertionError("No generator found")


def audit(p, n, coefficients, a_value=None, d=None, generator=None):
    residue = Ring(p, 1, d)
    initial = (generator, 0) if generator is not None else find_generator(residue, n)
    assert residue.power(initial, n) == residue.one
    assert residue.power(initial, n//2) == residue.scalar(-1)
    levels = []
    for precision in range(1, 9):
        ring = Ring(p, precision, d)
        # Independent lift: q^(r-1) preserves the residue and kills principal units.
        g = ring.power(initial, ring.q**(precision-1))
        assert ring.power(g, n) == ring.one
        group = [ring.power(g, j) for j in range(n)]
        assert len(set(group)) == n
        roots = [ring.power(ring.add(ring.one, group[j]), n) for j in range(1, n//2)]
        boundary = ring.scalar(2**n)
        if coefficients is not None:
            actual = [ring.one]
            for root in roots:
                following = [ring.zero]*(len(actual)+1)
                for i, c in enumerate(actual):
                    following[i] = ring.sub(following[i], ring.mul(root, c))
                    following[i+1] = ring.add(following[i+1], c)
                actual = following
            assert actual == [ring.scalar(c) for c in coefficients]
        multiplicities = Counter({boundary: 1})
        for root in roots:
            multiplicities[root] += 2
        excess = sum(c*c for c in multiplicities.values())-(2*n-3)
        direct_pairs = Counter(ring.add(a, b) for a in group for b in group)
        energy = sum(c*c for c in direct_pairs.values())
        assert energy == 3*n*n-3*n+n*excess
        orbits = orbit_counts(ring, group)
        assert orbits["D_from_orbits"] == excess
        levels.append({"precision": precision, "D": excess, "energy2": energy, **orbits})
        if excess == 0:
            break
    assert levels[-1]["D"] == 0
    v_disc = 2*sum(ring.val(ring.sub(roots[i], roots[j])) for i in range(len(roots)) for j in range(i))
    v_boundary = sum(ring.val(ring.sub(boundary, x)) for x in roots)
    v_a = v_disc+v_boundary
    if a_value is not None:
        assert v_a == valuation(a_value, p)
    v_three = valuation(3**n-1, p)
    v_u = sum(x["two_one_one_orbits"] for x in levels)
    v_v = sum(x["four_distinct_orbits"] for x in levels)
    v_b = v_u+2*v_v
    assert sum(x["epsilon_three_in_H"] for x in levels) == v_three
    assert v_boundary == v_three+v_u
    assert v_disc == 2*v_u+6*v_v
    assert v_a == v_three+3*v_b
    assert sum(x["D"] for x in levels) == 4*v_a
    return {"p": p, "n": n, "quadratic_nonresidue": d,
            "initial_generator": initial, "levels": levels,
            "valuation_A": v_a, "valuation_three_to_n_minus_one": v_three,
            "valuation_U": v_u, "valuation_V": v_v,
            "valuation_boundary": v_boundary, "valuation_discriminant": v_disc,
            "valuation_odd_cube_root": v_b}


def integer_cube_root(value):
    assert value >= 1
    x = 1 << ((value.bit_length()+2)//3)
    while True:
        following = (2*x+value//(x*x))//3
        if following >= x:
            assert x**3 <= value < (x+1)**3
            return x
        x = following


def cube_record(n, coefficients, a_value):
    two = valuation(a_value, 2)
    odd_three = (3**n-1) >> valuation(3**n-1, 2)
    assert valuation(3**n-1, 2) == n.bit_length()+1
    odd_a = a_value >> two
    assert odd_a % odd_three == 0
    b_cubed = odd_a//odd_three
    b = integer_cube_root(b_cubed)
    assert b**3 == b_cubed and b % 2 == 1
    boundary = evaluate(coefficients, 2**n)
    assert a_value % boundary == 0
    disc = a_value//boundary
    odd_boundary = boundary >> valuation(boundary, 2)
    odd_disc = disc >> valuation(disc, 2)
    assert odd_boundary % odd_three == 0
    u = odd_boundary//odd_three
    assert odd_disc % (u*u) == 0
    v_sixth = odd_disc//(u*u)
    v = isqrt(integer_cube_root(v_sixth))
    assert v**6 == v_sixth and b == u*v*v
    return {"n": n, "v_2_A": two, "odd_part_three_to_n_minus_one_hex": hex(odd_three),
            "odd_cube_root_hex": hex(b), "A_bits": a_value.bit_length(),
            "odd_cube_root_bits": b.bit_length(), "P_coefficients_ascending": coefficients,
            "U_hex": hex(u), "V_hex": hex(v), "U_bits": u.bit_length(), "V_bits": v.bit_length(),
            "A_hex": hex(a_value), "factorization_claimed": False}


def main():
    previous_path = ROOT/"results/kernel_discriminant.json"
    previous = json.loads(previous_path.read_text())
    helper = ROOT/"experiments/kernel_discriminant.py"
    assert previous["source_sha256"] == sha256(helper.read_bytes()).hexdigest()
    cubes, checks = [], []
    for row in previous["small_orders"]:
        n, coefficients, a_value = row["n"], row["P_coefficients_ascending"], int(row["A_hex"], 16)
        cubes.append(cube_record(n, coefficients, a_value))
        for exceptional in row["all_eligible_exceptional_primes"]:
            checks.append(audit(exceptional["p"], n, coefficients, a_value, generator=exceptional["generator"]))
    # A new unfactored determinant checks the cube identity at a larger order.
    coefficients = polynomial(64)
    a_value = discriminant(coefficients)*evaluate(coefficients, 2**64)
    cubes.append(cube_record(64, coefficients, a_value))
    for row in previous["larger_prime_power_checks"]:
        n = row["n"]
        checks.append(audit(row["p"], n, coefficients if n == 64 else None,
                            a_value if n == 64 else None, generator=row["generator"]))
    for p, n, d in [(3,4,2), (3,8,2), (5,8,2), (7,16,3), (17,32,3)]:
        row = next(x for x in previous["small_orders"] if x["n"] == n)
        checks.append(audit(p, n, row["P_coefficients_ascending"], int(row["A_hex"],16), d=d))
    result = {"status": "passed; uniform energy and Paley conjectures unproved",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "kernel_helper_sha256": sha256(helper.read_bytes()).hexdigest(),
              "determinant_helper_sha256": sha256((ROOT/"experiments/cyclotomic_norm_audit.py").read_bytes()).hexdigest(),
              "prior_result_sha256": sha256(previous_path.read_bytes()).hexdigest(),
              "cube_identities": cubes, "orbit_and_lift_checks": checks}
    output = ROOT/"results/quadruple_orbits.json"
    output.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"], "rings_checked": len(checks),
                      "precision_levels_checked": sum(len(x["levels"]) for x in checks),
                      "cube_checks": [{k: x[k] for k in ["n", "v_2_A", "A_bits", "U_bits", "V_bits"]} for x in cubes],
                      "output": str(output)}, indent=2))


if __name__ == "__main__":
    main()
