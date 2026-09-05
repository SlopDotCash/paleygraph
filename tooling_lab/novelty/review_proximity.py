#!/usr/bin/env python3
"""Bounded independent brute-force oracles for the proximity prototypes."""
from hashlib import sha256
from itertools import combinations, product
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "proximity"))
import deformation_microscope as dm
import stack_provenance as sp


def mul(m, v, p):
    return tuple(sum(a*b for a, b in zip(row, v)) % p for row in m)


def check_jets():
    matrices = [[[v[0], v[1]], [v[2], v[3]]] for v in product(range(2), repeat=4)]
    vectors = list(product(range(2), repeat=2))
    checked = 0
    for m0 in matrices:
        for m1 in matrices:
            # The third coefficient is deliberately nonzero and need not commute.
            m2 = [[1, 1], [0, 1]]
            ms = [m0, m1, m2]
            expected = []
            for ell in range(3):
                leading = set()
                for cs in product(vectors, repeat=ell+1):
                    if all(all(sum(mul(ms[i-j], cs[j], 2)[r]
                                       for j in range(i+1)) % 2 == 0 for r in range(2))
                           for i in range(ell+1)):
                        leading.add(cs[0])
                assert len(leading) in (1, 2, 4)
                expected.append({1: 0, 2: 1, 4: 2}[len(leading)])
            assert dm.lift_survival(ms, 2) == expected
            checked += 1
    return checked


def check_census():
    # Exhaust every received-word pencil over this small code, not just samples.
    p, dom, k, s = 3, [0, 1, 2], 1, 2
    words = list(product(range(p), repeat=len(dom)))
    code = [tuple([a] * len(dom)) for a in range(p)]
    cat = sp.parity_catalog(dom, k, s, p)
    checked = 0
    for u0 in words:
        for u1 in words:
            expected = {}
            for z in range(p):
                w = tuple((a+z*b) % p for a, b in zip(u0, u1))
                good = {c for c in code if sum(x == y for x, y in zip(c, w)) >= s}
                if good:
                    expected[z] = good
            actual = sp.census(u0, u1, dom, k, s, p, cat)
            assert actual["finite_bad_scalar_count"] == len(expected)
            if actual["node_ledger_complete"]:
                got = {z: {tuple(node["codeword"]) for node in actual["nodes"] if node["scalar"] == z}
                       for z in actual["finite_bad_scalars"]}
                assert got == expected
            else:
                assert actual["whole_field_correlated"] and len(expected) == p
            checked += 1
    return checked


def check_saved_tracks():
    path = ROOT / "proximity" / "stack_results.json"
    data = json.loads(path.read_text())
    nodes_checked = tracks_checked = 0
    for case in data["records"]:
        p, n, s = case["p"], case["n"], case["s"]
        nodes = case["nodes"]
        expected_tracks = {}
        for a, b in combinations(nodes, 2):
            if a["scalar"] == b["scalar"]:
                continue
            dz = (b["scalar"]-a["scalar"]) % p
            slope = tuple((y-x)*pow(dz, -1, p) % p for x, y in zip(a["codeword"], b["codeword"]))
            intercept = tuple((x-a["scalar"]*y) % p for x, y in zip(a["codeword"], slope))
            members = tuple(i for i, node in enumerate(nodes)
                            if all(c == (x+node["scalar"]*y) % p
                                   for c, x, y in zip(node["codeword"], intercept, slope)))
            expected_tracks[intercept, slope] = members
        actual_keys = set()
        for track in case["tracks"]:
            key = tuple(track["intercept"]), tuple(track["slope"])
            assert key not in actual_keys
            actual_keys.add(key)
            assert tuple(track["node_indices"]) == expected_tracks[key]
            common = [i for i in range(n) if case["u0"][i] == key[0][i] and case["u1"][i] == key[1][i]]
            assert track["common_coordinates"] == common
            t = len(common)
            scalars = {nodes[i]["scalar"] for i in track["node_indices"]}
            assert len(scalars) == track["scalar_count"]
            if t < s:
                assert len(scalars)*(s-t) <= n-t
            # Independently check the exact coordinate allocation, not only its bound.
            for j in range(n):
                matches = sum((case["u0"][j]+z*case["u1"][j]-key[0][j]-z*key[1][j]) % p == 0
                              for z in scalars)
                assert matches == len(scalars) if j in common else matches <= 1
            tracks_checked += 1
        assert actual_keys == set(expected_tracks)
        for node in nodes:
            mask = sum(1 << i for i in range(n)
                       if (case["u0"][i]+node["scalar"]*case["u1"][i]-node["codeword"][i]) % p == 0)
            assert mask == node["agreement_mask"]
            nodes_checked += 1
    return {"stacks": len(data["records"]), "nodes": nodes_checked, "tracks": tracks_checked}


def main():
    result = {"date": "2026-09-05", "passed": True,
              "brute_force_three_coefficient_jet_systems": check_jets(),
              "exhaustive_F3_RS31_pencils": check_census(),
              "saved_track_replay": check_saved_tracks(),
              "scope": "Independent finite brute-force oracles and saved incidence replay, not full large-field census or Lean validation"}
    result["source_hashes"] = {str(p): sha256(p.read_bytes()).hexdigest()
                                for p in (ROOT / "proximity").glob("*.py")}
    Path(__file__).with_name("proximity_review.json").write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
