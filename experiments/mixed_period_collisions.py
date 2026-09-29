#!/usr/bin/env python3
"""Exact mixed-period identities and shifted multiplicative collisions.

Finite certificates only. All arithmetic used in acceptance is integer.
"""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import permutations
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def prime(p):
    return p >= 2 and (p == 2 or (p % 2 and all(p % d for d in range(3, isqrt(p)+1, 2))))


def subgroup(p, n, generator):
    assert prime(p) and (p-1) % n == 0 and n & (n-1) == 0
    assert pow(generator, n, p) == 1 and pow(generator, n//2, p) == p-1
    h = sorted({pow(generator, j, p) for j in range(n)})
    assert len(h) == n
    return h


def falling3(k):
    return k*(k-1)*(k-2)


def shifted_energy(p, h):
    n = len(h)
    nonidentity = [a for a in h if a != 1]
    products = defaultdict(list)
    for a in nonidentity:
        for d in nonidentity:
            products[(a-1)*(d-1) % p].append((a, d))
    energy = sum(len(pairs)**2 for pairs in products.values())
    baseline = 2*(n-1)**2-(n-1)
    triples, cells = set(), defaultdict(set)
    inverse = lambda x: pow(x % p, -1, p)
    for pairs in products.values():
        for a, d in pairs:
            for b, c in pairs:
                if (a == b and d == c) or (a == c and d == b):
                    continue
                assert a != b and a != c and c != d and b != d
                x = (b-1)*inverse(a-b) % p
                y, z = a*x % p, c*x % p
                assert len({x, y, z}) == 3 and all(t not in (0, p-1) for t in (x, y, z))
                label = (pow(x, n, p), pow(x+1, n, p))
                assert all((pow(t, n, p), pow(t+1, n, p)) == label for t in (y, z))
                recovered = (y*inverse(x) % p, (y+1)*inverse(x+1) % p,
                             z*inverse(x) % p, (z+1)*inverse(x+1) % p)
                assert recovered == (a, b, c, d)
                assert (x, y, z) not in triples
                triples.add((x, y, z))
                cells[label].update([x, y, z])
    excess = energy-baseline
    assert len(triples) == excess and excess % 6 == 0
    complete_cells = []
    for label, members in sorted(cells.items()):
        representative = min(members)
        actual = {representative*a % p for a in h
                  if representative*a % p != p-1
                  and pow((representative*a+1) % p, n, p) == label[1]}
        assert members == actual
        assert all(t in triples for t in permutations(actual, 3))
        complete_cells.append({"coset_power_labels": list(label), "size": len(actual),
                               "members": sorted(actual)})
    assert sum(falling3(row["size"]) for row in complete_cells) == excess

    k = Counter(pow(x+1, n, p) for x in h if x != p-1)
    special = pow(2, n, p)
    assert all(c % 2 == (label == special) for label, c in k.items())
    assert k[special] % 2 == 1
    row_excess = sum(c*c for c in k.values())-(2*n-3)
    row_triples = sum(falling3(c) for c in k.values())
    assert row_excess >= 0 and 3*row_excess <= row_triples+6 <= excess+6
    additive_pairs = Counter((a+b) % p for a in h for b in h)
    additive_energy = sum(c*c for c in additive_pairs.values())
    assert additive_energy == 3*n*n-3*n+n*row_excess
    if excess == 0:
        assert row_excess == 0 and max(k.values()) <= 2
    return {"shift_energy_without_zero": energy, "trivial_shift_energy": baseline,
            "nontrivial_shift_energy": excess, "row_zero_energy_excess": row_excess,
            "row_zero_ordered_triples": row_triples,
            "additive_energy": additive_energy,
            "intersection_cell_size_histogram_above_two": dict(sorted(Counter(row["size"] for row in complete_cells).items())),
            "all_cells_above_two": complete_cells,
            "ordered_triple_sha256": sha256(json.dumps(sorted(triples)).encode()).hexdigest(),
            "circularity_certified": excess == 0}


def primitive_root(p):
    k, factors, d = p-1, [], 2
    while d*d <= k:
        if k % d == 0:
            factors.append(d)
            while k % d == 0:
                k //= d
        d += 1
    if k > 1:
        factors.append(k)
    return next(g for g in range(2, p) if all(pow(g, (p-1)//q, p) != 1 for q in factors))


def direct_matrix(p, h, expected):
    n, m, g = len(h), (p-1)//len(h), primitive_root(p)
    labels, cosets = {}, []
    for j in range(m):
        cell = sorted(pow(g, j, p)*x % p for x in h)
        assert all(x not in labels for x in cell)
        labels.update({x: j for x in cell})
        cosets.append(cell)
    assert len(labels) == p-1 and cosets[0] == h
    cells = defaultdict(list)
    for x in range(1, p-1):
        cells[labels[x], labels[x+1]].append(x)
    c = [Counter() for _ in range(m)]
    for (t, s), members in cells.items():
        c[t][s] = len(members)
    for t in range(m):
        assert sum(c[t].values()) == n-(t == 0)
        for s, value in c[t].items():
            assert value == c[s][t] == c[-t % m][(s-t) % m]
    total = sum(sum(row.values()) for row in c)
    double = sum(v*(v-1) for row in c for v in row.values())
    triple = sum(falling3(v) for row in c for v in row.values())
    assert total == p-2 and double == (n-1)*(n-2)
    assert triple == expected["nontrivial_shift_energy"]
    direct_triples = {t for members in cells.values() for t in permutations(members, 3)}
    assert sha256(json.dumps(sorted(direct_triples)).encode()).hexdigest() == expected["ordered_triple_sha256"]
    assert sum(c[t][t] for t in range(m)) == n-1
    for t in range(m):
        lhs = Counter((a+b) % p for a in h for b in cosets[t])
        rhs = Counter({0: n}) if t == 0 else Counter()
        for s, value in c[t].items():
            rhs.update({x: value for x in cosets[s]})
        assert lhs == rhs

    traces = []
    if p <= 97:
        # L=C-n*e0*1^T. Verify its traces independently using additive convolution.
        ell = [dict(row) for row in c]
        ell[0] = {s: c[0][s]-n for s in range(m)}
        power = [[int(i == j) for j in range(m)] for i in range(m)]
        counts = [1]+[0]*(p-1)
        for order in range(9):
            trace = sum(power[i][i] for i in range(m))
            assert n*trace == p*counts[0]-n**order
            traces.append(trace)
            if order == 8:
                break
            following = [[0]*m for _ in range(m)]
            for i in range(m):
                for j, value in enumerate(power[i]):
                    for k, edge in ell[j].items():
                        following[i][k] += value*edge
            power = following
            counts = [sum(counts[(x-a) % p] for a in h) for x in range(p)]
    return {"index": m, "primitive_root": g, "matrix_nonzero_cells": len(cells),
            "total_entries": total, "total_ordered_double_collisions": double,
            "total_ordered_triple_collisions": triple,
            "mixed_group_ring_identities_checked": m,
            "integer_multiplication_matrix_traces_orders_zero_to_eight": traces}


def run(p, n, generator, full_matrix=False):
    h = subgroup(p, n, generator)
    result = {"p": p, "n": n, "generator": generator,
              "quartic_window": n**4 <= 4*p <= 4*n**4,
              **shifted_energy(p, h)}
    if full_matrix:
        result["independent_full_matrix_audit"] = direct_matrix(p, h, result)
    return result


def main():
    specs = [(17,8,2,True),(97,8,33,True),(1049,8,223,True),
             (2017,8,438,True),(17393,16,2614,True),
             (6700417,64,2,False),(67403009,128,64701253,False),
             (1073748737,256,1064280392,False),
             (17179869697,512,13395504394,False)]
    cases = [run(*spec) for spec in specs]
    assert [row["nontrivial_shift_energy"] for row in cases] == [78,0,0,0,0,114,720,0,0]
    archive = ROOT/"sources/mixed-periods-2026-09-04"
    manifest = json.loads((archive/"manifest.json").read_text())
    for row in manifest["files"]:
        data = (archive/row["file"]).read_bytes()
        assert len(data) == row["bytes"] and sha256(data).hexdigest() == row["sha256"]
    result = {"status": "passed; no uniform Paley or prize bound proved",
              "source_sha256": sha256(Path(__file__).read_bytes()).hexdigest(),
              "archive_manifest_sha256": sha256((archive/"manifest.json").read_bytes()).hexdigest(),
              "cases": cases}
    path = ROOT/"results/mixed_period_collisions.json"
    path.write_text(json.dumps(result, indent=2)+"\n")
    print(json.dumps({"status": result["status"],
                      "cases": [{k: row[k] for k in ["p","n","nontrivial_shift_energy","additive_energy","circularity_certified"]} for row in cases],
                      "output": str(path)}, indent=2))


if __name__ == "__main__":
    main()
