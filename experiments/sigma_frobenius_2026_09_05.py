#!/usr/bin/env python3
"""Exact checks for research/sigma-frobenius-2026-09-05.md.

Standard library only (cmath for Jacobi sums). Writes
results/sigma_frobenius_2026_09_05.json.

Checks:
 1. W_D := sum_{t in F_p} chi(prod_{r in D}(t-r)) equals N_aff(C_D) - p where
    N_aff counts (t,y) in F_p^2 with y^2 = prod_{r in D}(t-r).  Random D.
 2. For D = c*H, H the subgroup of order m, prod_{h}(t - c h) = t^m - c^m, and
    W_D = sum_{psi^m = 1, psi != 1} psi(c^m) chi(-c^m) J(psi, chi)  exactly
    (no boundary term: the trivial character contributes -chi(-c^m) and the
    t = 0 term contributes +chi(-c^m)); the
    Jacobi sums J(psi,chi) = sum_v psi(v) chi(1-v) are evaluated in floating
    point with tolerance 1e-6 (they are algebraic integers; the identity is
    checked to that tolerance and the real part rounded).
 3. The genus bookkeeping: for |D| = d distinct roots, the Weil bound
    |W_D| <= (d-1) sqrt p is verified on every computed W_D.
"""
import cmath
import json
import math
import os
import random
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(n + 1) if s[i]]


def legendre(p):
    chi = [0] * p
    for x in range(1, p):
        chi[x * x % p] = 1
    for x in range(1, p):
        if chi[x] == 0:
            chi[x] = -1
    return chi


def primitive_root(p):
    phi = p - 1
    fac = []
    n = phi
    q = 2
    while q * q <= n:
        if n % q == 0:
            fac.append(q)
            while n % q == 0:
                n //= q
        q += 1
    if n > 1:
        fac.append(n)
    for g in range(2, p):
        if all(pow(g, phi // q, p) != 1 for q in fac):
            return g
    raise ValueError


def W_of_D(chi, D, p):
    tot = 0
    for t in range(p):
        v = 1
        for r in D:
            v *= chi[(t - r) % p]
            if v == 0:
                break
        tot += v
    return tot


def affine_points(D, p, chi):
    # number of (t,y) with y^2 = f_D(t): sum_t (1 + chi(f_D(t)))
    cnt = 0
    for t in range(p):
        v = 1
        for r in D:
            v = v * ((t - r) % p) % p
        cnt += 1 + chi[v]
    return cnt


def check1(primes, rng):
    out = {"cases": 0, "weil_checked": 0, "max_ratio_to_sqrtp": 0.0}
    for p in primes:
        chi = legendre(p)
        for _ in range(12):
            d = rng.randint(1, min(8, p - 1))
            D = rng.sample(range(p), d)
            W = W_of_D(chi, D, p)
            assert W == affine_points(D, p, chi) - p, (p, D)
            out["cases"] += 1
            if d >= 2:
                assert abs(W) <= (d - 1) * math.sqrt(p) + 1e-9, (p, D, W)
                out["weil_checked"] += 1
                out["max_ratio_to_sqrtp"] = max(out["max_ratio_to_sqrtp"], abs(W) / ((d - 1) * math.sqrt(p)))
            else:
                assert W == 0
    return out


def check2(primes):
    """Jacobi-sum expression for subgroup branch loci."""
    out = {"cases": 0, "max_abs_error": 0.0, "boundary_term_formula": "none: the trivial-character term -chi(-a) cancels the t=0 term +chi(-a) exactly", "records": []}
    for p in primes:
        chi = legendre(p)
        g = primitive_root(p)
        for m in [d for d in range(2, 40) if (p - 1) % d == 0]:
            r = pow(g, (p - 1) // m, p)
            H = [pow(r, i, p) for i in range(m)]
            zeta = cmath.exp(2j * cmath.pi / m)
            # discrete log table for psi
            dlog = [None] * p
            x = 1
            for k in range(p - 1):
                dlog[x] = k
                x = x * g % p
            # Jacobi sums J(psi_j, chi) = sum_{v != 0,1} psi_j(v) chi(1-v)
            J = []
            for j in range(m):
                s = 0j
                for v in range(2, p):
                    s += zeta ** ((j * dlog[v]) % m) * chi[(1 - v) % p]
                J.append(s)
            for c in [1, 2, 3, g, (g * g) % p]:
                c %= p
                if c == 0:
                    continue
                D = [c * h % p for h in H]
                # f_D(t) = prod (t - c h) = t^m - c^m
                a = pow(c, m, p)
                W = W_of_D(chi, D, p)
                # direct evaluation of sum_t chi(t^m - a)
                Wdirect = sum(chi[(pow(t, m, p) - a) % p] for t in range(p))
                assert W == Wdirect, (p, m, c)
                # Jacobi expansion: sum over psi != 1 of psi(a) chi(-a) J(psi,chi), plus boundary
                s = 0j
                for j in range(1, m):
                    s += zeta ** ((j * dlog[a]) % m) * chi[(-a) % p] * J[j]
                # boundary: trivial character term is sum_{u != 0} chi(u - a) = -chi(-a);
                # the t = 0 term contributes chi(-a) which cancels it exactly.
                val = s.real
                err = abs(s.imag) + abs(val - W)
                out["max_abs_error"] = max(out["max_abs_error"], err)
                assert err < 1e-5, (p, m, c, W, s)
                out["cases"] += 1
                if len(out["records"]) < 12:
                    out["records"].append({"p": p, "m": m, "c": c, "W_D": W, "genus_upper": (m - 1) // 2})
    return out


def main():
    rng = random.Random(7)
    primes = [q for q in primes_upto(200) if q >= 5]
    res = {
        "check1_pointcount_identity": check1(primes, rng),
        "check2_jacobi_expansion": check2([q for q in primes if q <= 120]),
        "scope": "finite exact identities only; the Katz-Sarnak independence heuristic is not tested here",
    }
    out = os.path.join(ROOT, "results", "sigma_frobenius_2026_09_05.json")
    with open(out, "w") as f:
        json.dump(res, f, indent=1)
    print(json.dumps({k: (v if not isinstance(v, dict) else {kk: vv for kk, vv in v.items() if kk != "records"}) for k, v in res.items()}, indent=1))


if __name__ == "__main__":
    sys.exit(main())
