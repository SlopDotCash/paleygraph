#!/usr/bin/env python3
"""Read back saved evidence and verify local certificates independently of generation loops."""
import json
from math import comb
from pathlib import Path

from deformation_microscope import coefficient_matrix, matvec, is_prime, domain
from stack_provenance import interpolate_values


def main():
    base = Path(__file__).parent
    deformation = json.loads((base / "deformation_results.json").read_text())
    stack = json.loads((base / "stack_results.json").read_text())
    counts = {"syzygy_vectors": 0, "decoded_codewords": 0, "affine_tracks": 0}
    for row in deformation["profiles"]:
        p = row["p"]
        assert is_prime(p)
        groups = row["groups"]
        roots = sum(groups, [])
        assert len(set(roots)) == len(roots)
        assert all(pow(x, row["n"], p) == 1 for x in roots)
        m = coefficient_matrix(groups, row["cofactor_degree"], p)
        for v in row["kernel_basis"]:
            assert any(v) and matvec(m, v, p) == [0] * len(m)
            counts["syzygy_vectors"] += 1
    for row in stack["records"]:
        p, n, k, s = row["p"], row["n"], row["k"], row["s"]
        dom = domain(n, p)
        assert row["exact_agreement_subsets_enumerated"] == comb(n, s)
        for node in row["nodes"]:
            z, cw = node["scalar"], node["codeword"]
            assert tuple(cw) == interpolate_values(dom[:k], cw[:k], dom, p)
            mask = sum(1 << j for j in range(n)
                       if (row["u0"][j] + z * row["u1"][j] - cw[j]) % p == 0)
            assert mask == node["agreement_mask"] and mask.bit_count() >= s
            counts["decoded_codewords"] += 1
        for z, number in row["raw_subset_certificate_count_by_scalar"].items():
            expected = sum(comb(node["agreement_size"], s) for node in row["nodes"]
                           if node["scalar"] == int(z))
            assert expected == number
        for track in row["tracks"]:
            for index in track["node_indices"]:
                node = row["nodes"][index]
                assert node["codeword"] == [(a + node["scalar"] * b) % p
                                              for a, b in zip(track["intercept"], track["slope"])]
            counts["affine_tracks"] += 1
    print(json.dumps({"status": "passed", **counts}, indent=2))


if __name__ == "__main__":
    main()
