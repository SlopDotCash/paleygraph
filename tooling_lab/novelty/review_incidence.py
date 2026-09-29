#!/usr/bin/env python3
"""Independent exact replay of the remaining F61 quartic-incidence collision."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, permutations
import json
from math import prod
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent


def record(c, p=61):
    squares = {a*a % p for a in range(1, p)}
    def chi(a):
        a %= p
        return 0 if a == 0 else 1 if a in squares else -1
    f = [sum(chi(x-a) for a in c) for x in range(p)]
    graph = [[0]*6 for _ in c]
    for i, j in combinations(range(6), 2):
        quartic = [a for k, a in enumerate(c) if k not in (i, j)]
        graph[i][j] = graph[j][i] = sum(prod(chi(x-a) for a in quartic) for x in range(p))
    rooted = sorted(sorted(graph[i][j] for j in range(6) if j != i) for i in range(6))
    triangle = 6*sum(graph[i][j]*graph[j][k]*graph[k][i] for i, j, k in combinations(range(6), 3))
    return {"C": c, "m2": sum(x*x for x in f), "m4": sum(x**4 for x in f),
            "m6": sum(x**6 for x in f), "b2": sum(f[a]**2 for a in c),
            "b4": sum(f[a]**4 for a in c),
            "energy": sum(v*v for v in Counter((a+b) % p for a in c for b in c).values()),
            "quartic_deck": sorted(graph[i][j] for i, j in combinations(range(6), 2)),
            "quartic_incidence": rooted, "quartic_triangle": triangle, "graph": graph}


def main():
    left = record([0, 1, 2, 4, 38, 52])
    right = record([0, 1, 3, 8, 10, 21])
    for key in ("m2", "m4", "b2", "b4", "energy", "quartic_deck", "quartic_incidence"):
        assert left[key] == right[key]
    assert [left["m6"], right["m6"]] == [101790, 124830]
    assert [left["quartic_triangle"], right["quartic_triangle"]] == [-1992, 4152]
    matching = [perm for perm in permutations(range(6))
                if all(left["graph"][i][j] == right["graph"][perm[i]][perm[j]]
                       for i, j in combinations(range(6), 2))]
    assert not matching
    path = ROOT / "observables" / "results_incidence.json"
    result = json.loads(path.read_text())
    case = next(r for r in result["cases"] if r["p"] == 61)
    for expected, saved in zip((left, right), case["after"]["witness"]):
        for key, value in expected.items():
            if key != "graph":
                assert saved[key] == value
    out = {"date": "2026-09-05", "passed": True, "witnesses": [left, right],
           "permutations_checked": 720, "matching_weighted_graph_permutations": len(matching),
           "scope": "Exact witness replay only; full census counts inspected but not independently rerun",
           "result_sha256": sha256(path.read_bytes()).hexdigest(),
           "reviewer_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    HERE.joinpath("incidence_review.json").write_text(json.dumps(out, indent=2)+"\n")
    print(json.dumps({"passed": True, "m6": [left["m6"], right["m6"]],
                      "triangle_traces": [left["quartic_triangle"], right["quartic_triangle"]],
                      "weighted_graph_isomorphisms": 0}))


if __name__ == "__main__":
    main()
