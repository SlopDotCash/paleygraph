#!/usr/bin/env python3
"""Exact deep moments by quotienting additive convolution by multiplication by H.

NumPy int64 operations are used only after an explicit overflow bound. Moment
accumulation and zero-frequency subtraction use Python's arbitrary integers.
Floating point is not used in any acceptance condition.
"""
import argparse
import json
from math import prod
from pathlib import Path

import numpy as np

from paley_exact import prime


def ceil_root(value, degree):
    if value == 0:
        return 0
    low, high = 0, 1 << ((value.bit_length() + degree - 1) // degree)
    while low + 1 < high:
        mid = (low + high) // 2
        if mid ** degree >= value:
            high = mid
        else:
            low = mid
    assert (high - 1) ** degree < value <= high ** degree
    return high


def quotient(p, n, generator):
    assert prime(p) and n > 1 and (p - 1) % n == 0
    assert p * p <= np.iinfo(np.int64).max, "Residue products would overflow"
    assert (p - 1) // n + 1 <= np.iinfo(np.int32).max, "Orbit labels would overflow"
    assert pow(generator, n, p) == 1
    subgroup = np.array(sorted({pow(generator, j, p) for j in range(n)}), dtype=np.int64)
    assert len(subgroup) == n
    labels = np.full(p, -1, dtype=np.int32)
    labels[0] = 0
    representatives = [0]
    for b in range(1, p):
        if labels[b] < 0:
            representatives.append(b)
            labels[(subgroup * b) % p] = len(representatives) - 1
    representatives = np.array(representatives, dtype=np.int64)
    assert len(representatives) == 1 + (p - 1) // n
    assert np.all(labels >= 0)
    sizes = np.bincount(labels)
    assert sizes[0] == 1 and np.all(sizes[1:] == n)
    # Check that every constructed orbit has its own label, in bounded batches.
    for start in range(1, len(representatives), 4096):
        rr = representatives[start:start + 4096]
        actual = labels[(rr[:, None] * subgroup) % p]
        expected = np.arange(start, start + len(rr), dtype=np.int32)[:, None]
        assert np.all(actual == expected)
    neighbors = labels[(representatives[:, None] - subgroup) % p]
    return subgroup, representatives, neighbors


def run(p, n, generator, depth, dense_crosscheck=False):
    subgroup, representatives, neighbors = quotient(p, n, generator)
    counts = np.zeros(len(representatives), dtype=np.int64)
    counts[0] = 1
    dense = [1] + [0] * (p - 1) if dense_crosscheck else None
    output = []
    for r in range(1, depth + 1):
        max_before = int(counts.max())
        assert n * max_before <= np.iinfo(np.int64).max, "Would overflow int64"
        counts = counts[neighbors].sum(axis=1, dtype=np.int64)
        assert np.all(counts >= 0)
        weighted_mass = int(counts[0]) + n * sum(map(int, counts[1:]))
        assert weighted_mass == n ** r
        if dense is not None:
            next_dense = [0] * p
            for x, count in enumerate(dense):
                for h in subgroup:
                    next_dense[(x + int(h)) % p] += count
            dense = next_dense
            assert all(int(counts[i]) == dense[int(a)]
                       for i, a in enumerate(representatives))
        energy = int(counts[0]) ** 2 + n * sum(int(c) ** 2 for c in counts[1:])
        centered = p * energy - n ** (2 * r)
        assert centered >= 0 and centered % n == 0
        coset_moment = centered // n
        gaussian = (p - 1) * prod(range(1, 2 * r, 2)) * n ** r
        sg_K2 = (p - 1) * (2 * r * n) ** r
        output.append({"r": r, "energy": energy, "centered_moment": centered,
                       "coset_moment": coset_moment,
                       "M_squared_upper_integer": ceil_root(coset_moment, r),
                       "gaussian_coefficient_fails": centered > gaussian,
                       "SG_with_K_2_at_this_instance": centered <= sg_K2,
                       "max_count": int(counts.max()),
                       "overflow_bound_before_step": n * max_before})
        print(f"p={p} n={n} r={r}: M^2 <= {output[-1]['M_squared_upper_integer']}; "
              f"SG(K=2)={centered <= sg_K2}", flush=True)
    return {"p": p, "n": n, "generator": generator,
            "index": (p - 1) // n, "depth": depth, "moments": output}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path("results/quotient_moments.json"))
    args = parser.parse_args()
    small = [run(17, 8, 2, 5, True), run(97, 8, 33, 5, True)]
    target = run(6700417, 64, 2, 12)
    assert target["moments"][1]["energy"] == 12864
    assert target["moments"][2]["energy"] == 4816960
    data = {"asymptotic_conjecture_proved": False, "arithmetic": "exact integers",
            "dense_crosschecks": small, "target": target}
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(data, indent=2) + "\n")
    print(f"Wrote {args.output}", flush=True)


if __name__ == "__main__":
    main()
