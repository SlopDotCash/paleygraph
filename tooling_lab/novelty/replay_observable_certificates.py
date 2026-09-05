#!/usr/bin/env python3
"""Independent scalar replay; no imports from the prototypes under review."""
from collections import Counter
from hashlib import sha256
from itertools import combinations
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def evaluator(p):
    squares = {x*x % p for x in range(1, p)}
    chi = [0] + [1 if x in squares else -1 for x in range(1, p)]

    def evaluate(c, full=True):
        f = [sum(chi[(x-a) % p] for a in c) for x in range(p)]
        out = {"m2": sum(v**2 for v in f), "m4": sum(v**4 for v in f),
               "m6": sum(v**6 for v in f), "b2": sum(f[a]**2 for a in c),
               "b4": sum(f[a]**4 for a in c)}
        if full:
            out["energy"] = sum(v*v for v in Counter((a+b) % p for a in c for b in c).values())
            out["quartic_deck"] = sorted(sum(math.prod(chi[(x-a) % p] for a in q)
                                                for x in range(p)) for q in combinations(c, 4))
        out["interior_histogram"] = [sum(abs(f[x]) == level for x in range(p) if x not in c)
                                     for level in (0, 2, 4, 6)]
        return out
    return evaluate


def main():
    replayed = []
    inputs = []
    for filename in ("results_toy.json", "results_scale.json", "results_holdout.json"):
        path = ROOT / "observables" / filename
        if not path.exists():
            continue
        inputs.append({"path": str(path), "sha256": sha256(path.read_bytes()).hexdigest()})
        for case in json.loads(path.read_text())["cases"]:
            evaluate = evaluator(case["p"])
            for name, group in case["groups"].items():
                if group["witness"] is None:
                    continue
                left, right = group["witness"]
                assert all(left[f] == right[f] for f in group["feature_names"])
                assert left["m6"] != right["m6"]
                for w in (left, right):
                    expected = evaluate(w["C"])
                    for key, value in expected.items():
                        if key != "interior_histogram":
                            assert w[key] == value, (filename, case["p"], key)
                    replayed.append({"file": filename, "p": case["p"], "group": name,
                                     "C": w["C"], "m6": w["m6"]})

    # Verify the rank-one histogram move with elementary integer determinants.
    levels = (0, 2, 4, 6)
    move = (-10, 15, -6, 1)
    assert all(sum(v*x**k for v, x in zip(move, levels)) == 0 for k in (0, 2, 4))
    assert sum(v*x**6 for v, x in zip(move, levels)) == 23040
    # Vandermonde determinant in squared values 0,4,16 is nonzero.
    assert (4-0)*(16-0)*(16-4) == 768
    # Likewise 1,9,25 determines the three boundary masses.
    assert (9-1)*(25-1)*(25-9) == 3072

    # Independently enumerate all normalized six-sets over F29 and compare
    # actual ranges against a brute-force nonnegative histogram relaxation.
    p = 29
    evaluate = evaluator(p)
    fibers = {}
    count = 0
    for tail in combinations(range(2, p), 4):
        r = evaluate((0, 1) + tail, full=False)
        key = tuple(r[k] for k in ("m2", "m4", "b2", "b4"))
        y = r["m6"]
        if key not in fibers:
            fibers[key] = [y, y, r]
        f = fibers[key]
        f[0], f[1] = min(f[0], y), max(f[1], y)
        count += 1
    gaps = []
    for low, high, r in fibers.values():
        h = r["interior_histogram"]
        feasible = [t for t in range(-p, p+1) if all(x+t*v >= 0 for x, v in zip(h, move))]
        upper = r["m6"] + max(feasible)*23040
        assert r["m6"] + min(feasible)*23040 <= low <= high <= upper
        if upper > high:
            gaps.append(upper-high)
    census = {"p": p, "normalized_sets": count, "actual_fibers": len(fibers),
              "upper_realizability_gap_fibers": len(gaps), "max_upper_gap": max(gaps, default=0)}
    trade_path = ROOT / "observables" / "results_trades.json"
    if trade_path.exists():
        trade = json.loads(trade_path.read_text())
        saved = next(case for case in trade["realizability"] if case["p"] == p)
        for key, value in census.items():
            assert saved[key] == value, (key, saved[key], value)
        inputs.append({"path": str(trade_path), "sha256": sha256(trade_path.read_bytes()).hexdigest()})
    result = {"date": "2026-09-05", "scope": "Independent scalar witness replay and F29 exhaustive census",
              "passed": True, "inputs": inputs, "checks": len(replayed), "replayed": replayed,
              "primitive_trade_verified": True, "independent_census": census,
              "census_compared_to_saved_trades": trade_path.exists(),
              "reviewer_script_sha256": sha256(Path(__file__).read_bytes()).hexdigest()}
    out = Path(__file__).with_name("observable_review.json")
    out.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps({"passed": True, "witness_replays": len(replayed), "census": census,
                      "census_compared_to_saved_trades": trade_path.exists()}))


if __name__ == "__main__":
    main()
