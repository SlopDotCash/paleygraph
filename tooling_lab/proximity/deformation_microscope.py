#!/usr/bin/env python3
"""Exact, deterministic local syzygy deformation and admissible-edit diagnostics.

No third-party dependencies. This is a research microscope, not a prize solver.
All linear algebra uses Python integers modulo a verified prime below 2**64.
"""
from __future__ import annotations

import argparse
from collections import Counter
from itertools import combinations, product
import json
from pathlib import Path
import random
import time


def is_prime(p):
    if p < 2 or p >= 2**64:
        return False
    for q in (2, 3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37):
        if p % q == 0:
            return p == q
    d, s = p - 1, 0
    while d % 2 == 0:
        d //= 2
        s += 1
    # Deterministic Miller--Rabin bases for unsigned 64-bit integers.
    for a in (2, 325, 9375, 28178, 450775, 9780504, 1795265022):
        if a % p == 0:
            continue
        x = pow(a, d, p)
        if x in (1, p - 1):
            continue
        for _ in range(s - 1):
            x = x * x % p
            if x == p - 1:
                break
        else:
            return False
    return True


def next_prime(target, n):
    p = target + (1 - target) % n
    while not is_prime(p):
        p += n
    return p


def factors(n):
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


def domain(n, p):
    assert is_prime(p) and (p - 1) % n == 0
    fs = factors(n)
    for a in range(2, p):
        w = pow(a, (p - 1) // n, p)
        if pow(w, n, p) == 1 and all(pow(w, n // q, p) != 1 for q in fs):
            return [pow(w, i, p) for i in range(n)]
    raise ValueError("No root of requested order")


def rref(rows, p, ncols=None):
    a = [[x % p for x in row] for row in rows]
    ncols = ncols if ncols is not None else len(a[0])
    piv, r = [], 0
    for c in range(ncols):
        rr = next((i for i in range(r, len(a)) if a[i][c]), None)
        if rr is None:
            continue
        a[r], a[rr] = a[rr], a[r]
        inv = pow(a[r][c], -1, p)
        a[r] = [x * inv % p for x in a[r]]
        for i in range(len(a)):
            if i != r and a[i][c]:
                f = a[i][c]
                a[i] = [(x - f * y) % p for x, y in zip(a[i], a[r])]
        piv.append(c)
        r += 1
        if r == len(a):
            break
    return a, piv


def rank(rows, p, ncols=None):
    return len(rref(rows, p, ncols)[1])


def kernel(rows, p, ncols=None):
    ncols = ncols if ncols is not None else len(rows[0])
    a, piv = rref(rows, p, ncols)
    out = []
    for c in range(ncols):
        if c in piv:
            continue
        v = [0] * ncols
        v[c] = 1
        for i, pc in enumerate(piv):
            v[pc] = -a[i][c] % p
        out.append(v)
    return out


def transpose(a):
    return list(map(list, zip(*a)))


def matvec(a, v, p):
    return [sum(x * y for x, y in zip(row, v)) % p for row in a]


def poly_roots(roots, p):
    c = [1]
    for x in roots:
        out = [0] * (len(c) + 1)
        for i, v in enumerate(c):
            out[i] = (out[i] - x * v) % p
            out[i + 1] = (out[i + 1] + v) % p
        c = out
    return c


def coefficient_matrix(groups, q, p):
    """Columns X^j W_i, 0<=j<=q, for the three root-polynomials W_i."""
    d = max(map(len, groups))
    rows = [[0] * (len(groups) * (q + 1)) for _ in range(d + q + 1)]
    for i, roots in enumerate(groups):
        for j in range(q + 1):
            for z, c in enumerate(poly_roots(roots, p)):
                rows[z + j][i * (q + 1) + j] = c
    return rows


def jets(groups, velocities, q, p, order):
    """M(t) mod t^(order+1), roots x_i+t*v_i. Coefficient arithmetic, no derivatives."""
    d = max(map(len, groups))
    out = [[[0] * (len(groups) * (q + 1)) for _ in range(d + q + 1)]
           for _ in range(order + 1)]
    offset = 0
    for gi, roots in enumerate(groups):
        c = [[1] + [0] * order]
        for idx, root in enumerate(roots):
            velocity = velocities[offset + idx]
            nc = [[0] * (order + 1) for _ in range(len(c) + 1)]
            for xi, ts in enumerate(c):
                for ti, value in enumerate(ts):
                    nc[xi][ti] = (nc[xi][ti] - root * value) % p
                    nc[xi + 1][ti] = (nc[xi + 1][ti] + value) % p
                    if ti < order:
                        nc[xi][ti + 1] = (nc[xi][ti + 1] - velocity * value) % p
            c = nc
        for j in range(q + 1):
            for xi, ts in enumerate(c):
                for ti, value in enumerate(ts):
                    out[ti][xi + j][gi * (q + 1) + j] = value
        offset += len(roots)
    return out


def lift_survival(ms, p):
    """dim of c_0 that extend to M(t)c(t)=0 mod t^(ell+1), each ell.

    Let T_ell be the lower block Toeplitz coefficient system. Kernel elements
    with c_0=0 identify with ker(T_(ell-1)). Rank-nullity then proves
    survivor_dim(ell)=nullity(T_ell)-nullity(T_(ell-1)). This avoids choosing
    arbitrary lifts, which can otherwise manufacture false higher obstructions.
    """
    nr, nc = len(ms[0]), len(ms[0][0])
    previous, out = 0, []
    for ell in range(len(ms)):
        block = []
        for i in range(ell + 1):
            for row in range(nr):
                block.append([ms[i - j][row][c] if j <= i else 0
                              for j in range(ell + 1) for c in range(nc)])
        nullity = nc * (ell + 1) - rank(block, p)
        out.append(nullity - previous)
        previous = nullity
    assert all(a >= b for a, b in zip(out, out[1:]))
    return out


def tangent_map(groups, q, p):
    m = coefficient_matrix(groups, q, p)
    right, left = kernel(m, p), kernel(transpose(m), p)
    n = sum(map(len, groups))
    cols = []
    for i in range(n):
        v = [int(j == i) for j in range(n)]
        dm = jets(groups, v, q, p, 1)[1]
        cols.append([sum(a * b for a, b in zip(l, matvec(dm, r, p))) % p
                     for l in left for r in right])
    rows = transpose(cols)
    return rows, kernel(rows, p, n), right, left


def edit_spectrum(groups, full_domain, q, p):
    """Exhaustive one-root replacements by unused subgroup points; all remain disjoint."""
    used = set(sum(groups, []))
    unused = [x for x in full_domain if x not in used]
    base = len(coefficient_matrix(groups, q, p)[0]) - rank(coefficient_matrix(groups, q, p), p)
    counts, first = Counter(), None
    for gi, roots in enumerate(groups):
        for ri, old in enumerate(roots):
            for new in unused:
                gs = [g[:] for g in groups]
                gs[gi][ri] = new
                m = coefficient_matrix(gs, q, p)
                nullity = len(m[0]) - rank(m, p)
                counts[nullity] += 1
                if nullity < base and first is None:
                    first = {"group": gi, "index": ri, "old": old, "new": new,
                             "nullity_after": nullity}
    return {"operation": "one root replaced by an unused domain point",
            "trials": sum(counts.values()), "nullity_histogram": dict(sorted(counts.items())),
            "first_damage": first}


def profile(name, groups, full_domain, q, p, order=3, tangent_samples=12, edit=True):
    m = coefficient_matrix(groups, q, p)
    tm, tb, rk, lk = tangent_map(groups, q, p)
    n = sum(map(len, groups))
    dirs = {"translation": [1] * n, "scaling": sum(groups, []),
            "single_root": [1] + [0] * (n - 1)}
    rng = random.Random(20260905)
    for i in range(min(tangent_samples, len(tb))):
        dirs[f"tangent_basis_{i}"] = tb[i]
    for i in range(tangent_samples):
        cs = [rng.randrange(p) for _ in tb]
        dirs[f"tangent_random_{i}"] = [sum(c * b[j] for c, b in zip(cs, tb)) % p
                                        for j in range(n)]
    survival = {name: lift_survival(jets(groups, v, q, p, order), p)
                for name, v in dirs.items()}
    # Known changes of variable are exact syzygy-preserving controls in the ambient space.
    assert survival["translation"] == [len(rk)] * (order + 1)
    assert survival["scaling"] == [len(rk)] * (order + 1)
    if rk:
        for direction_name, curve in survival.items():
            if direction_name.startswith("tangent_"):
                assert curve[1] == len(rk)
    return {"name": name, "p": p, "n": len(full_domain), "groups": groups,
            "cofactor_degree": q, "matrix_shape": [len(m), len(m[0])],
            "rank": rank(m, p), "kernel_dimension": len(rk),
            "cokernel_dimension": len(lk), "kernel_basis": rk,
            "root_parameter_count": n, "tangent_constraint_rank": rank(tm, p, n),
            "tangent_dimension": len(tb), "jet_order": order,
            "jet_survival": survival,
            "tangent_survival_histogram": dict(Counter(str(v) for k, v in survival.items()
                                                        if k.startswith("tangent_"))),
            "example_tangent_velocity": tb[0] if tb else None,
            "domain_tangent_dimension": 0,
            "domain_tangent_reason": "n*x^(n-1) is nonzero at each n-th root because p does not divide n",
            "edits": edit_spectrum(groups, full_domain, q, p) if edit else None}


F41_GROUPS = [[2, 25, 37, 23, 39, 1], [8, 18, 16, 10, 5, 36], [4, 21, 33, 40, 31, 9]]


def self_test():
    """Independent tiny brute enumerations validate kernel and projected-jet dimension."""
    p = 3
    examples = [([[1, 2], [2, 1]], [[0, 1], [0, 0]]),
                ([[0, 0], [0, 0]], [[1, 0], [0, 0]]),
                ([[1, 0], [0, 0]], [[0, 0], [1, 0]])]
    for m0, m1 in examples:
        valid = set()
        for c0 in product(range(p), repeat=2):
            for c1 in product(range(p), repeat=2):
                if matvec(m0, c0, p) == [0, 0] and all(
                    (a + b) % p == 0 for a, b in zip(matvec(m0, c1, p), matvec(m1, c0, p))):
                    valid.add(c0)
        surv = lift_survival([m0, m1], p)
        assert len(valid) == p ** surv[1]
    # A known second-order obstruction invisible at order one.
    assert lift_survival([[[0]], [[0]], [[1]]], 5) == [1, 1, 0]
    # Matrix coefficients agree with direct root motion at enough nonzero t values.
    gs, vel, q, p = [[1, 2], [3, 4], [5, 6]], [2, 1, 0, 3, 2, 4], 1, 11
    js = jets(gs, vel, q, p, 2)
    for t in range(p):
        moved = [[(x + t * vel[2 * i + j]) % p for j, x in enumerate(g)]
                 for i, g in enumerate(gs)]
        direct = coefficient_matrix(moved, q, p)
        assert direct == [[sum(pow(t, k, p) * js[k][i][j] for k in range(3)) % p
                           for j in range(len(direct[0]))] for i in range(len(direct))]
    return {"tiny_exhaustive_projection_examples": len(examples), "jet_polynomial_values": p,
            "second_order_null_control": "passed"}


def run(output, quick=False):
    start = time.perf_counter()
    records, checks = [], self_test()
    d41 = domain(20, 41)
    records.append(profile("F41_linear_middle_witness", F41_GROUPS, d41, 1, 41))
    shuffled = d41[:]
    random.Random(99).shuffle(shuffled)
    records.append(profile("F41_random_partition_control", [shuffled[i*6:(i+1)*6] for i in range(3)],
                           d41, 1, 41, tangent_samples=3))
    idx = [[d41.index(x) for x in g] for g in F41_GROUPS]
    prime_sweep = []
    for target in (41, 61, 101, 1009, 100003, 2**31):
        p = next_prime(target, 20)
        dom = domain(20, p)
        gs = [[dom[j] for j in g] for g in idx]
        m = coefficient_matrix(gs, 1, p)
        prime_sweep.append({"p": p, "n": 20, "cofactor_degree": 1,
                            "matrix_rank": rank(m, p), "kernel_dimension": len(kernel(m, p))})
    for n in ([16, 24] if quick else [16, 24, 48, 96]):
        for target in ([1009] if quick else [n + 1, 1009, 2**31]):
            p = next_prime(target, n)
            dom = domain(n, p)
            gs = [dom[i::4] for i in range(3)]
            records.append(profile(f"coset_n{n}_p{p}", gs, dom, 0, p,
                                   tangent_samples=3 if n > 24 else 8, edit=n <= 48))
    result = {"schema": "proximity-deformation-microscope/v1", "seed": 20260905,
              "source_checkout": "25de4107388ce3e371c2f56252279ab6455eb1a9",
              "source_state": "read-only; contains untracked SW1 work not Lean-validated here",
              "scope": "exact finite-field local diagnostics, no MCA or prize conclusion",
              "checks": checks, "field_transport_of_F41_index_pattern": prime_sweep,
              "profiles": records, "elapsed_seconds": round(time.perf_counter() - start, 3)}
    output.write_text(json.dumps(result, indent=2) + "\n")
    for r in records:
        print(r["name"], "kernel", r["kernel_dimension"], "tangent", r["tangent_dimension"],
              "survival", r["tangent_survival_histogram"],
              "edits", r["edits"]["nullity_histogram"] if r["edits"] else "skipped")
    print("p-sweep", prime_sweep)
    print("Saved", output, "in", result["elapsed_seconds"], "seconds")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=Path(__file__).with_name("deformation_results.json"))
    parser.add_argument("--quick", action="store_true")
    args = parser.parse_args()
    run(args.output, args.quick)
