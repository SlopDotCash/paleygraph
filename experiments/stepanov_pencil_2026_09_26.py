#!/usr/bin/env python3
"""Exact checks for research/stepanov-pencil-2026-09-26.md.

Standard library only. Polynomials over F_p are coefficient lists, lowest
degree first. Writes results/stepanov_pencil_2026_09_26.json.
"""
import json
import math
import os
import random
import sys
import time
from itertools import combinations

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_pencil_2026_09_26.json")


def trim(a):
    a = list(a)
    while len(a) > 1 and a[-1] == 0:
        a.pop()
    return a


def padd(a, b, p):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def pmul(a, b, p):
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return trim(r)


def pscale(a, c, p):
    return trim([x * c % p for x in a])


def pderiv(a, p):
    return trim([(i * a[i]) % p for i in range(1, len(a))] or [0])


def peval(a, x, p):
    r = 0
    for c in reversed(a):
        r = (r * x + c) % p
    return r


def mult_at(f, b, p):
    """Multiplicity of b as a root of f (f nonzero)."""
    f = trim(f)
    m = 0
    while any(f) and peval(f, b, p) == 0:
        n = len(f) - 1
        q = [0] * n
        carry = f[n]
        for i in range(n - 1, -1, -1):
            q[i] = carry
            carry = (f[i] + carry * b) % p
        f = trim(q)
        m += 1
    return m


def lin_pow(a, D, p):
    """(x + a)^D."""
    return trim([math.comb(D, i) * pow(a, D - i, p) % p for i in range(D + 1)])


def inv(x, p):
    return pow(x % p, p - 2, p)


def hp_poly(A, p):
    """Hanson-Petridis polynomial F_A with its normalisation G = 1."""
    M = len(A)
    d = (p - 1) // 2
    D = d + M - 1
    c = [inv(math.prod((A[k] - A[l]) for l in range(M) if l != k), p) for k in range(M)]
    F = [p - 1]
    for k in range(M):
        F = padd(F, pscale(lin_pow(A[k], D, p), c[k], p), p)
    return F, D, c


def chi_table(p):
    t = [0] * p
    for x in range(1, p):
        t[x * x % p] = 1
    for x in range(1, p):
        if t[x] == 0:
            t[x] = -1
    return t


def elem_sym(vals_polys, i, p):
    """i-th elementary symmetric polynomial of a list of polynomials."""
    total = [0]
    for sub in combinations(vals_polys, i):
        prod = [1]
        for q in sub:
            prod = pmul(prod, q, p)
        total = padd(total, prod, p)
    return total


class Checks:
    def __init__(self):
        self.n = 0
        self.fail = []

    def check(self, ok, tag, info=None):
        self.n += 1
        if not ok:
            self.fail.append({"tag": tag, "info": info})


def main():
    t0 = time.time()
    rng = random.Random(20260926)
    C = Checks()
    primes = [29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149]
    records = []
    max_ratio = {"value": 0.0}
    for p in primes:
        d = (p - 1) // 2
        chi = chi_table(p)
        for trial in range(10):
            M = rng.randint(3, min(8, p // 5))
            A = rng.sample(range(p), M)
            F, D, c = hp_poly(A, p)
            C.check(len(F) - 1 == d, "deg_F_eq_d", (p, A))
            Fp = pderiv(F, p)
            Fpp = pderiv(Fp, p)
            invD = inv(D, p)
            W = pscale(Fp, invD, p)
            V = padd(F, pscale(pmul([0, 1], W, p), p - 1, p), p)
            # Lemma 1: exact pencil identity F_{A \ a_k} = V - a_k W
            for k in range(M):
                Fk, _, _ = hp_poly([A[l] for l in range(M) if l != k], p)
                C.check(Fk == padd(V, pscale(W, (-A[k]) % p, p), p), "pencil_identity", (p, A, k))
            # Lemma 1': F_{A \ E} = sum_i (-1)^i sigma_i(x + a_E)/(D)_i F^{(i)}
            derivs = [F]
            for _ in range(3):
                derivs.append(pderiv(derivs[-1], p))
            for e in (2, 3):
                if M - e < 2:
                    continue
                E = rng.sample(range(M), e)
                FE, _, _ = hp_poly([A[l] for l in range(M) if l not in E], p)
                lin = [[A[k] % p, 1] for k in E]
                acc = [0]
                for i in range(e + 1):
                    falling = 1
                    for s in range(i):
                        falling = falling * (D - s) % p
                    coef = pscale(elem_sym(lin, i, p) if i else [1], ((-1) ** i) % p * inv(falling, p) % p, p)
                    acc = padd(acc, pmul(coef, derivs[i], p), p)
                C.check(acc == FE, "linear_system_identity", (p, A, E))
            # Wronskian N = V'W - VW' = ((D-1)F'^2 - D F F'')/D^2, degree 2d-2
            N = padd(pmul(pderiv(V, p), W, p), pscale(pmul(V, pderiv(W, p), p), p - 1, p), p)
            N2 = pscale(padd(pscale(pmul(Fp, Fp, p), (D - 1) % p, p),
                             pscale(pmul(F, Fpp, p), (-D) % p, p), p), inv(D * D, p), p)
            C.check(N == N2, "wronskian_formula", (p, A))
            C.check(len(N) - 1 == 2 * d - 2, "deg_N", (p, A, len(N) - 1))
            negA = {(-a) % p for a in A}
            N0 = N1 = 0
            for b in range(p):
                if b in negA:
                    continue
                E = [k for k in range(M) if chi[(b + A[k]) % p] == -1]
                y = [(b + A[k]) % p for k in range(M)]
                T = lambda t: sum(c[k] * pow(y[k], t, p) for k in E) % p
                # value of N at a rational point: 4(D-1)(T_{M-2}^2 - T_{M-1} T_{M-3})
                pred = 4 * (D - 1) * (T(M - 2) ** 2 - T(M - 1) * T(M - 3)) % p
                C.check(peval(N, b, p) == pred, "N_value_formula", (p, A, b))
                if len(E) == 0:
                    N0 += 1
                    C.check(mult_at(N, b, p) >= 2 * M - 2, "mult_complete", (p, A, b))
                elif len(E) == 1:
                    N1 += 1
                    C.check(mult_at(N, b, p) >= M - 2, "mult_one_bad", (p, A, b))
                    C.check(peval(W, b, p) != 0, "W_nonzero_one_bad", (p, A, b))
                elif len(E) == 2:
                    C.check(peval(N, b, p) != 0, "N_nonzero_two_bad", (p, A, b))
            C.check((2 * M - 2) * N0 + (M - 2) * N1 <= 2 * d - 2, "ramification_inequality", (p, A, N0, N1))
            ratio = M * (N0 + N1) / d
            records.append({"p": p, "M": M, "A": A, "N0": N0, "N1_exact": N1,
                            "M_times_N_le1_over_d": round(ratio, 4),
                            "ramification_budget_used": round(((2 * M - 2) * N0 + (M - 2) * N1) / (2 * d - 2), 4)})
            if ratio > max_ratio["value"]:
                max_ratio = {"value": ratio, "p": p, "A": A}
    result = {
        "checks": C.n,
        "failures": C.fail,
        "n_failures": len(C.fail),
        "seconds": round(time.time() - t0, 1),
        "max_M_times_N_le1_over_d": max_ratio,
        "records": records,
        "scope": "exact identities at small primes; no asymptotic claim",
    }
    with open(OUT, "w") as f:
        json.dump(result, f, indent=1)
    print("checks", C.n, "failures", len(C.fail), "seconds", result["seconds"])
    return 1 if C.fail else 0


if __name__ == "__main__":
    sys.exit(main())
