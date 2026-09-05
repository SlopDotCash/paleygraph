#!/usr/bin/env python3
"""Sparse exact check at a Fermat factor in the concrete quartic window."""
from collections import Counter
from fractions import Fraction
from math import prod
from pathlib import Path
import json

from paley_exact import prime


def main():
    p, n, generator = 6700417, 64, 2
    assert prime(p)
    subgroup = sorted({pow(generator, j, p) for j in range(n)})
    assert len(subgroup) == n and pow(generator, n, p) == 1
    assert all(pow(x, n, p) == 1 for x in subgroup)
    assert n ** 4 <= 4 * p and p <= n ** 4
    counts = Counter({0: 1})
    moments = []
    for r in range(1, 4):
        next_counts = Counter()
        for x, count in counts.items():
            for h in subgroup:
                next_counts[(x + h) % p] += count
        counts = next_counts
        assert sum(counts.values()) == n ** r
        energy = sum(c * c for c in counts.values())
        centered = p * energy - n ** (2 * r)
        bound = (p - 1) * prod(range(1, 2 * r, 2)) * n ** r
        moments.append({"r": r, "energy": energy, "centered": centered,
                        "gaussian_bound": bound, "gaussian_fails": centered > bound,
                        "ratio": str(Fraction(centered, bound))})
    hs = set(subgroup)
    normalized = sum((1 + a - b) % p in hs for a in subgroup for b in subgroup)
    assert normalized == 201
    assert moments[1]["energy"] == n * normalized == 12864
    assert moments[1]["centered"] == 86177387072
    assert moments[1]["gaussian_bound"] == 82334711808
    assert moments[1]["gaussian_fails"]
    result = {"p": p, "n": n, "generator": generator, "H": subgroup,
              "p_over_n_fourth": str(Fraction(p, n ** 4)), "moments": moments,
              "normalized_collisions": normalized, "primality_verified": True,
              "asymptotic_target_refuted": False}
    path = Path(__file__).resolve().parents[1] / "results/fermat-factor-moments.json"
    path.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({k: v for k, v in result.items() if k != "H"}, indent=2))


if __name__ == "__main__":
    main()
