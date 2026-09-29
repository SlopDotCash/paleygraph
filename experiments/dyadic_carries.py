#!/usr/bin/env python3
"""Exact ordered signed-power counts using carries, without a modulus-sized array.

For Q=a^N+1, count k-tuples from {+/-a^j:0<=j<N} whose sum is target
modulo Q. These are integer counts, not numerical Fourier integrals.
"""
from collections import Counter, defaultdict
from fractions import Fraction
from math import comb, factorial
import hashlib
import json
from pathlib import Path


def carry_counts(N, K, target=0, periodic=True, radix=2):
    assert N >= 1 and K >= 0 and radix >= 2
    Q = radix ** N + 1
    if periodic:
        target %= Q
        # Sign symmetry handles the only residue needing N+1 base-a digits.
        if target == Q - 1:
            target = 1
        initial_carries = range(-K // radix - 1, K // radix + 2)
    else:
        assert target == 0
        initial_carries = [0]
    answer = [0] * (K + 1)
    for initial in initial_carries:
        counts = {(0, initial): 1}
        for j in range(N):
            digit = target // radix ** j % radix
            following = defaultdict(int)
            for (used, carry), count in counts.items():
                for added in range(K - used + 1):
                    shuffle = comb(used + added, added)
                    for positive in range(added + 1):
                        difference = 2 * positive - added
                        numerator = carry + difference - digit
                        if numerator % radix:
                            continue
                        next_carry = numerator // radix
                        ways = shuffle * comb(added, positive)
                        following[used + added, next_carry] += count * ways
            counts = following
        final = -initial if periodic else 0
        for k in range(K + 1):
            answer[k] += counts.get((k, final), 0)
    return answer


def dense_counts(N, depth, radix=2):
    Q = radix ** N + 1
    H = sorted({sign * radix ** j % Q for j in range(N) for sign in [-1, 1]})
    assert len(H) == 2 * N
    counts = [1] + [0] * (Q - 1)
    answer = [counts]
    for k in range(1, depth + 1):
        following = [0] * Q
        for x, count in enumerate(counts):
            for h in H:
                following[(x + h) % Q] += count
        counts = following
        assert sum(counts) == (2 * N) ** k
        answer.append(counts)
    return answer


def audit_fourth_pattern(p, H):
    # A zero-sum quadruple is degenerate iff it can be partitioned into
    # opposite pairs. The exceptional pattern is (u,u,2u,-4u).
    from itertools import permutations
    n = len(H)
    patterns = set()
    for u in H:
        word = (u, u, 2 * u % p, -4 * u % p)
        orbit = set(permutations(word))
        assert len(orbit) == 12
        assert all(sum(w) % p == 0 for w in orbit)
        assert all(not any((w[0] + w[j]) % p == 0 for j in [1, 2, 3]) for w in orbit)
        assert patterns.isdisjoint(orbit)
        patterns.update(orbit)
    pairs = Counter((a + b) % p for a in H for b in H)
    energy = sum(c * c for c in pairs.values())
    assert len(patterns) == 12 * n
    assert energy >= 3 * n * n + 9 * n
    return {"p": p, "n": n, "extra_patterns": len(patterns), "energy": energy}


def main():
    root = Path(__file__).resolve().parents[1]
    crosschecks = []
    for radix, N in [(2, N) for N in [2, 3, 4, 5, 6]] + [(3, 3), (4, 3)]:
        depth = 8
        dense = dense_counts(N, depth, radix=radix)
        for target in range(radix ** N + 1):
            counted = carry_counts(N, depth, target, radix=radix)
            assert counted == [row[target] for row in dense]
        crosschecks.append({"radix": radix, "N": N, "modulus": radix ** N + 1,
                            "all_residues_checked": True, "depth": depth})

    N, K, p, divisor = 32, 24, 6700417, 641
    Q, n = (1 << N) + 1, 2 * N
    assert Q == p * divisor
    intrinsic = carry_counts(N, K)
    continuous = carry_counts(N, K, periodic=False)
    seen, aliases = set(), [0] * (K + 1)
    orbits = []
    for a in range(1, divisor):
        if a in seen:
            continue
        orbit = {a * pow(2, j, divisor) % divisor for j in range(n)}
        assert len(orbit) == n and seen.isdisjoint(orbit)
        seen.update(orbit)
        counts = carry_counts(N, K, p * a)
        for k, count in enumerate(counts):
            aliases[k] += len(orbit) * count
        first = next((k for k, count in enumerate(counts) if count), None)
        orbits.append({"representative_multiplier": a, "orbit_size": len(orbit),
                       "first_nonzero_order": first, "counts": counts})
        print(f"multiplier orbit {a}: first additional relation at order {first}", flush=True)
    assert seen == set(range(1, divisor))

    path = root / "results/quotient_moments.json"
    source = json.loads(path.read_text())["target"]
    assert (source["p"], source["n"]) == (p, n)
    checks = []
    for row in source["moments"]:
        r, energy = row["r"], row["energy"]
        k = 2 * r
        assert intrinsic[k] + aliases[k] == energy
        baseline_centered = Q * intrinsic[k] - n ** k
        assert 0 <= baseline_centered <= (Q - 1) * (160 * r * n) ** r
        # The following signed error is precisely the difference of the two
        # centered raw spectral moments, divided by p.
        discrepancy = Fraction(aliases[k]) - Fraction((divisor - 1) * n ** k, Q)
        assert p * discrepancy == row["centered_moment"] - Fraction(p * baseline_centered, Q)
        checks.append({"r": r, "intrinsic_count": intrinsic[k],
                       "additional_prime_count": aliases[k], "prime_energy": energy,
                       "centered_alias_discrepancy": str(discrepancy)})
    assert all(count == 0 for count in aliases[:8])
    assert aliases[8] == 12902400 == 64 * 5 * factorial(8)
    terms = [(0, 1), (7, -1), (9, -1), (14, 1),
             (17, -1), (19, 1), (21, -1), (23, 1)]
    assert sum(sign * (1 << j) for j, sign in terms) == p
    assert p % Q != 0
    H = {pow(2, j, p) for j in range(n)}
    pattern = audit_fourth_pattern(p, H)
    assert pattern["energy"] == 3 * n ** 2 + 9 * n == intrinsic[4]

    # The same prime subgroup has generator 8; an unqualified assertion that
    # the centered alias discrepancy is nonpositive is false for this lift.
    assert {pow(8, j, p) for j in range(n)} == H
    alternate = carry_counts(N, 4, radix=8)
    Q8 = 8 ** N + 1
    assert Q8 % p == 0 and alternate[4] == 3 * n ** 2 - 3 * n
    D4_8 = Fraction(pattern["energy"] - alternate[4]) - Fraction(n ** 4, p) + Fraction(n ** 4, Q8)
    assert D4_8 > 0

    result = {"main_conjecture_proved": False,
              "arithmetic": "exact integers and fractions",
              "dense_crosschecks": crosschecks,
              "quotient_moments_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
              "N": N, "n": n, "p": p, "Q": Q, "cofactor": divisor,
              "continuous_integer_zero_counts": continuous,
              "intrinsic_periodic_zero_counts": intrinsic,
              "prime_alias_counts": aliases, "alias_orbits": orbits,
              "independent_energy_comparisons": checks,
              "shortest_additional_relation_order": 8,
              "explicit_signed_power_relation": terms,
              "fourth_pattern_audit": pattern,
              "generator_8_comparison": {"ambient_modulus": Q8,
                  "ambient_fourth_count": alternate[4],
                  "centered_fourth_alias_discrepancy": str(D4_8),
                  "nonpositive_discrepancy_refuted": True}}
    output = root / "results/dyadic_carries.json"
    output.write_text(json.dumps(result, indent=2) + "\n")
    print(f"All twelve prime energies independently reconstructed; wrote {output}")


if __name__ == "__main__":
    main()
