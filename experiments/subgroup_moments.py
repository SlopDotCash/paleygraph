#!/usr/bin/env python3
"""Exact centered additive-energy moments for the prize-linked subgroup target.

No floating point Fourier values are used. Fixed-r experimental values cannot
establish the uniform logarithmic-depth estimate sought in this project.
"""
import argparse
import json
from fractions import Fraction
from pathlib import Path
from math import prod
from paley_exact import prime


def factors(n):
    out = []
    d = 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def primitive_root(p):
    divisors = factors(p - 1)
    return next(g for g in range(2, p)
                if all(pow(g, (p - 1) // d, p) != 1 for d in divisors))


def exact_moments(p, n, depth):
    if not prime(p) or p == 2 or (p - 1) % n or n < 2 or n & (n - 1):
        raise ValueError("Require an odd prime and a dyadic subgroup order dividing p-1")
    g = primitive_root(p)
    h = pow(g, (p - 1) // n, p)
    subgroup = sorted(pow(h, j, p) for j in range(n))
    assert len(set(subgroup)) == n
    assert all(pow(x, n, p) == 1 for x in subgroup)
    assert all((-x) % p in subgroup for x in subgroup)
    counts = [0] * p
    counts[0] = 1
    results = []
    for r in range(1, depth + 1):
        new_counts = [0] * p
        # Cyclic convolution by 1_H, expressed without modular division in the inner loop.
        for x in subgroup:
            cut = p - x
            for j in range(cut):
                new_counts[j + x] += counts[j]
            for j in range(cut, p):
                new_counts[j - cut] += counts[j]
        counts = new_counts
        assert sum(counts) == n ** r
        assert all(counts[x] == counts[(h * x) % p] for x in range(1, p))
        energy = sum(x * x for x in counts)
        centered = p * energy - n ** (2 * r)
        assert centered >= 0 and centered % n == 0
        if r == 1:
            assert energy == n and centered == n * (p - n)
        if r == 2:
            differences = [0] * p
            for a in subgroup:
                for b in subgroup:
                    differences[(a - b) % p] += 1
            assert energy == sum(x * x for x in differences)
        pairings = prod(range(1, 2 * r, 2))
        gaussian = (p - 1) * pairings * n ** r
        results.append({"r": r, "additive_energy": energy,
                        "centered_moment_sum": centered,
                        "coset_moment_sum": centered // n,
                        "gaussian_comparison_bound": gaussian,
                        "gaussian_ratio": str(Fraction(centered, gaussian)),
                        "gaussian_comparison_fails": centered > gaussian})
    return {"p": p, "n": n, "index": (p - 1) // n, "primitive_root": g,
            "subgroup_generator": h, "H": subgroup, "moments": results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/subgroup_moments.json"))
    args = parser.parse_args()
    cases = [(17, 4, 8), (97, 8, 8), (193, 16, 8), (257, 4, 8),
             (257, 16, 8), (769, 16, 8), (3329, 8, 8), (12289, 16, 8),
             (65537, 16, 10), (65537, 64, 8)]
    output = {"target": "additive Gauss periods of dyadic subgroups",
              "uniform_bound_proved": False, "cases": []}
    for p, n, depth in cases:
        result = exact_moments(p, n, depth)
        output["cases"].append(result)
        failures = [item["r"] for item in result["moments"] if item["gaussian_comparison_fails"]]
        print(f"p={p} n={n} depth={depth} Gaussian comparison failures at r={failures}", flush=True)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(output, indent=2) + "\n")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
