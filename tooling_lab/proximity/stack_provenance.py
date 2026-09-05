#!/usr/bin/env python3
"""Exact bad-scalar certificates with decoded-codeword track provenance.

The subset parity criterion and affine-pencil bounds are known methods. This
prototype adds an inspectable witness ledger and joins it to deformation profiles.
It samples stacks, but censuses ALL agreement subsets for each sampled stack.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations
import json
from pathlib import Path
import random
import time

from deformation_microscope import domain, next_prime, kernel, matvec, F41_GROUPS


def parity(points, k, p):
    """GRS dual: x_j^ell / prod_{i!=j}(x_j-x_i), ell<s-k."""
    weights = []
    for j, x in enumerate(points):
        prod = 1
        for i, y in enumerate(points):
            if i != j:
                prod = prod * (x - y) % p
        weights.append(pow(prod, -1, p))
    return [[pow(x, ell, p) * w % p for x, w in zip(points, weights)]
            for ell in range(len(points) - k)]


def parity_catalog(dom, k, s, p):
    return [(sub, parity([dom[j] for j in sub], k, p))
            for sub in combinations(range(len(dom)), s)]


def interpolate_values(xs, ys, dom, p):
    """Lagrange interpolation evaluated on the entire domain (not coefficients)."""
    known = dict(zip(xs, ys))
    weights = []
    for i, x in enumerate(xs):
        d = 1
        for j, y in enumerate(xs):
            if i != j:
                d = d * (x - y) % p
        weights.append(ys[i] * pow(d, -1, p) % p)
    out = []
    for x in dom:
        if x in known:
            out.append(known[x])
            continue
        acc = 0
        for i, w in enumerate(weights):
            prod = w
            for j, y in enumerate(xs):
                if i != j:
                    prod = prod * (x - y) % p
            acc += prod
        out.append(acc % p)
    return tuple(out)


def syndrome(hs, sub, word, p):
    return [sum(h * word[j] for h, j in zip(row, sub)) % p for row in hs]


def census(u0, u1, dom, k, s, p, catalog):
    witnesses = {}
    shared = []
    raw_counts = Counter()
    for sub, hs in catalog:
        a0, a1 = syndrome(hs, sub, u0, p), syndrome(hs, sub, u1, p)
        idx = next((i for i, a in enumerate(a1) if a), None)
        if idx is None:
            if not any(a0):
                shared.append(sub)
            continue
        z = -a0[idx] * pow(a1[idx], -1, p) % p
        if any((a + z * b) % p for a, b in zip(a0, a1)):
            continue
        raw_counts[z] += 1
        word = [(a + z * b) % p for a, b in zip(u0, u1)]
        candidates = witnesses.setdefault(z, {})
        submask = sum(1 << j for j in sub)
        if any((submask & mask) == submask for mask in candidates):
            continue
        codeword = interpolate_values([dom[j] for j in sub[:k]], [word[j] for j in sub[:k]], dom, p)
        mask = sum(1 << j for j in range(len(dom)) if word[j] == codeword[j])
        assert mask.bit_count() >= s
        candidates[mask] = codeword
    nodes = [{"scalar": z, "agreement_mask": mask, "agreement_size": mask.bit_count(),
              "codeword": list(cw)} for z in sorted(witnesses)
             for mask, cw in sorted(witnesses[z].items())]
    return {"finite_bad_scalars": None if shared else sorted(witnesses),
            "finite_bad_scalar_count": p if shared else len(witnesses),
            "candidate_scalars_from_noncommon_subsets": sorted(witnesses),
            "node_ledger_complete": not bool(shared),
            "whole_field_correlated": bool(shared), "common_agreement_subset_count": len(shared),
            "first_common_agreement_subset": list(shared[0]) if shared else None,
            "raw_subset_certificate_count_by_scalar": dict(sorted(raw_counts.items())),
            "nodes": nodes}


def tracks(nodes, u0, u1, s, p):
    """Group pairs by exact affine codeword track f_z=A+z*B."""
    candidates = {}
    for i, a in enumerate(nodes):
        for j, b in enumerate(nodes[i + 1:], i + 1):
            dz = (b["scalar"] - a["scalar"]) % p
            if dz == 0:
                continue
            inv = pow(dz, -1, p)
            slope = tuple((y - x) * inv % p for x, y in zip(a["codeword"], b["codeword"]))
            intercept = tuple((x - a["scalar"] * y) % p for x, y in zip(a["codeword"], slope))
            candidates.setdefault((intercept, slope), set()).update((i, j))
    records = []
    for (intercept, slope), members in candidates.items():
        common = [i for i, (a, b, c, d) in enumerate(zip(u0, u1, intercept, slope))
                  if a == c and b == d]
        t = len(common)
        zs = sorted({nodes[i]["scalar"] for i in members})
        # Every non-common coordinate agrees at at most one scalar on this track.
        # Hence m*(s-t)<=n-t, provided t<s. Check the actual allocation as well.
        assert t >= s or len(zs) * (s - t) <= len(u0) - t
        records.append({"node_indices": sorted(members), "scalars": zs,
                        "scalar_count": len(zs), "common_coordinates": common,
                        "common_size": t,
                        "capacity": (len(u0) - t) // (s - t) if t < s else None,
                        "slack": len(u0) - t - len(zs) * (s - t),
                        "intercept": list(intercept), "slope": list(slope)})
    records.sort(key=lambda r: (-r["scalar_count"], -r["common_size"], r["scalars"]))
    return records


def make_stack(groups, dom, k, s, p, seed):
    d = len(groups[0])
    rem = [x for x in dom if x not in set(sum(groups, []))]
    t = s - 2 * d
    assert 0 <= t <= len(rem)
    common = rem[:t]
    core_values = [groups[0] + groups[1] + common,
                   groups[0] + groups[2] + common,
                   groups[1] + groups[2] + common]
    cores = [[dom.index(x) for x in core] for core in core_values]
    rows = []
    for z, core in zip((0, 1, 2), cores):
        for h in parity([dom[j] for j in core], k, p):
            hfull = [0] * len(dom)
            for j, value in zip(core, h):
                hfull[j] = value
            rows.append(hfull + [z * x % p for x in hfull])
    basis = kernel(rows, p)
    rng = random.Random(seed)
    cs = [rng.randrange(p) for _ in basis]
    u = [sum(c * b[j] for c, b in zip(cs, basis)) % p for j in range(2 * len(dom))]
    assert matvec(rows, u, p) == [0] * len(rows)
    return u[:len(dom)], u[len(dom):], {"cores": cores, "stack_kernel_dimension": len(basis)}


def brute_test():
    """Independent enumeration of all RS codewords and all scalar words on n=4,k=2."""
    p, k, s, dom = 5, 2, 3, [1, 2, 3, 4]
    u0, u1 = [2, 4, 0, 2], [0, 1, 4, 1]
    cat = parity_catalog(dom, k, s, p)
    result = census(u0, u1, dom, k, s, p, cat)
    code = [tuple((a + b * x) % p for x in dom) for a in range(p) for b in range(p)]
    brute = {}
    for z in range(p):
        word = [(a + z * b) % p for a, b in zip(u0, u1)]
        good = {cw for cw in code if sum(x == y for x, y in zip(cw, word)) >= s}
        if good:
            brute[z] = good
    observed = {z: {tuple(n["codeword"]) for n in result["nodes"] if n["scalar"] == z}
                for z in result["finite_bad_scalars"]}
    assert brute == observed
    # Explicit whole-field common-agreement control.
    correlated = census([x % p for x in dom], [2 * x % p for x in dom], dom, k, s, p, cat)
    assert correlated["whole_field_correlated"]
    return {"independent_all_codewords_check": "passed", "correlated_control": "passed"}


def run(output, samples):
    start = time.perf_counter()
    tests, results = brute_test(), []
    configs = [("F41_linear_middle", 20, 10, 14, 41),
               ("coset_small_characteristic", 16, 8, 11, 17),
               ("coset_medium_characteristic", 16, 8, 11, 1009),
               ("coset_large_characteristic", 16, 8, 11, next_prime(2**31, 16))]
    for name, n, k, s, p in configs:
        dom = domain(n, p)
        groups = F41_GROUPS if n == 20 else [dom[i::4] for i in range(3)]
        cat = parity_catalog(dom, k, s, p)
        for sample in range(samples):
            u0, u1, construction = make_stack(groups, dom, k, s, p, 41000 + sample)
            result = census(u0, u1, dom, k, s, p, cat)
            ts = tracks(result["nodes"], u0, u1, s, p)
            nontrivial = [track for track in ts if track["scalar_count"] >= 3]
            covered_nodes = set().union(*(set(track["node_indices"]) for track in nontrivial)) if nontrivial else set()
            covered_scalars = {result["nodes"][i]["scalar"] for i in covered_nodes}
            record = {"configuration": name, "n": n, "k": k, "s": s, "p": p,
                      "sample": sample, "seed": 41000 + sample, "u0": u0, "u1": u1,
                      "construction": construction,
                      "exact_agreement_subsets_enumerated": len(cat), **result,
                      "tracks": ts, "track_scalar_count_histogram": dict(Counter(t["scalar_count"] for t in ts)),
                      "max_track_scalar_count": max((t["scalar_count"] for t in ts), default=0),
                      "nontrivial_track_profile": {
                          "definition": "tracks through at least three distinct scalar-codeword points",
                          "track_count": len(nontrivial), "covered_node_count": len(covered_nodes),
                          "covered_scalar_count": len(covered_scalars),
                          "repeated_node_memberships": sum(len(t["node_indices"]) for t in nontrivial) - len(covered_nodes),
                          "uncovered_scalar_count": result["finite_bad_scalar_count"] - len(covered_scalars)}}
            results.append(record)
            print(name, sample, "bad scalars", record["finite_bad_scalar_count"],
                  "nodes", len(record["nodes"]), "common", record["whole_field_correlated"],
                  "track histogram", record["track_scalar_count_histogram"], flush=True)
    out = {"schema": "proximity-stack-provenance/v1", "tests": tests,
           "scope": "stacks sampled; finite affine scalars exactly censused; infinity excluded",
           "sampling_policy": "fixed scalar core constraints z=0,1,2; seeded uniform linear combinations of kernel basis",
           "warning": "These are per-stack exact counts, never a maximum over all stacks or a prize conclusion.",
           "records": results, "elapsed_seconds": round(time.perf_counter() - start, 3)}
    output.write_text(json.dumps(out, indent=2) + "\n")
    print("Saved", output, "in", out["elapsed_seconds"], "seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--samples", type=int, default=4)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("stack_results.json"))
    args = parser.parse_args()
    run(args.output, args.samples)
