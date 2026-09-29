#!/usr/bin/env python3
"""Verifier for paper/robust-hanson-petridis.md (2026-09-29).

Re-checks, by exact integer / modular / rational computation at small primes, every identity,
inequality and displayed constant of the write-up.  Floating point is used only where marked
("float"), always with an explicit margin, and never for a statement that the paper proves by
an exact identity.  No code is imported from the earlier verifiers.

Sections (keys of the JSON "checks" dictionary carry the same letters):
  A  Lagrange weights (Lemma 2.1) and the local form of the Hanson-Petridis polynomial (Lemma 3.1)
  B  full Hankel determinants H_{e+1}: u_s, Step 1, exact degree, leading coefficient, orders
  C  inequality (star), supersaturation (Cor. B), robust HP (Thm C) incl. the size-free eta-form
  D  constant bias (Thm D): core lemma, final bound, non-vacuity thresholds
  E  sharpness family |A| = 2 and the cap 2k/(1+2k)
  F  leading coefficient Lambda: Krattenthaler product, p does not divide Lambda, the
     elementary interpolation argument, the counterexample with D >= p
  G  constants and expansions (series, crossings, Section 9 numbers)
  H  Paley-graph density corollary (Cor. E)
  I  second-moment (Vinogradov) identity and bound

Run: /opt/miniconda3/bin/python3 experiments/paper_robust_hp_2026_09_29.py
Writes: results/paper_robust_hp_2026_09_29.json
"""
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction

import numpy as np
import sympy
import mpmath

T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "paper_robust_hp_2026_09_29.json")
RNG = random.Random(20260929)

CHECKS = {}
FAILURES = []
WITNESSES = {}
NOTES = []


def ok(key, cond, info=None):
    CHECKS[key] = CHECKS.get(key, 0) + 1
    if not cond:
        if len(FAILURES) < 200:
            FAILURES.append({"check": key, "info": info})
    return cond


def okn(key, n_cases, n_bad, info=None):
    """Record a vectorised batch of n_cases checks, n_bad of which failed."""
    CHECKS[key] = CHECKS.get(key, 0) + int(n_cases)
    if n_bad:
        FAILURES.append({"check": key, "failed": int(n_bad), "info": info})


def primes_between(a, b):
    return [q for q in range(max(a, 2), b + 1) if sympy.isprime(q)]


def legendre_table(p):
    d = (p - 1) // 2
    t = [0] * p
    for x in range(1, p):
        t[x] = 1 if pow(x, d, p) == 1 else -1
    return t


def weights(A, p):
    c = []
    for k, ak in enumerate(A):
        prod = 1
        for l, al in enumerate(A):
            if l != k:
                prod = prod * (ak - al) % p
        c.append(pow(prod, p - 2, p))
    return c


def falling(x, j):
    r = 1
    for t in range(j):
        r *= (x - t)
    return r


def poly_F(A, p):
    """F(x) = -1 + sum_k c_k (x + a_k)^D, coefficient list low -> high, mod p."""
    d = (p - 1) // 2
    m = len(A)
    D = d + m - 1
    c = weights(A, p)
    F = [0] * (D + 1)
    for l in range(D + 1):
        s = 0
        for ck, ak in zip(c, A):
            s += ck * pow(ak, D - l, p)
        F[l] = math.comb(D, l) * s % p
    F[0] = (F[0] - 1) % p
    return F, c, D


def deriv_at(P, j, b, p):
    """j-th ordinary derivative of P at b, mod p."""
    s = 0
    bp = 1
    for l in range(j, len(P)):
        if P[l]:
            s += P[l] * falling(l, j) * pow(b, l - j, p)
    return s % p


def order_at(P, b, p, cap=None):
    """Multiplicity of b as a root of P (mod p) by repeated synthetic division; P != 0."""
    Q = [x % p for x in P]
    while Q and Q[-1] == 0:
        Q.pop()
    assert Q, "zero polynomial"
    o = 0
    while len(Q) > 1:
        # divide by (x - b)
        n = len(Q) - 1
        quot = [0] * n
        acc = Q[n]
        quot[n - 1] = acc
        for i in range(n - 1, 0, -1):
            acc = (Q[i] + acc * b) % p
            quot[i - 1] = acc
        rem = (Q[0] + acc * b) % p
        if rem != 0:
            break
        o += 1
        Q = quot
        if cap is not None and o >= cap:
            break
    return o


# ---------------------------------------------------------------------------------------------
# A. Lagrange weights and the local form of F
# ---------------------------------------------------------------------------------------------

def section_A():
    for p in [5, 7, 11, 13]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        sets = []
        for m in range(1, (p + 1) // 2 + 1):
            for rest in itertools.combinations(range(1, p), m - 1):
                sets.append((0,) + rest)
        for A in sets:
            check_local_form(A, p, chi, d, "A_exhaustive")
    for p in primes_between(17, 43):
        chi = legendre_table(p)
        d = (p - 1) // 2
        for _ in range(4):
            m = RNG.randint(2, (p + 1) // 2)
            A = tuple(sorted(RNG.sample(range(p), m)))
            check_local_form(A, p, chi, d, "A_random")


def check_local_form(A, p, chi, d, tag):
    m = len(A)
    F, c, D = poly_F(A, p)
    # Lemma 2.1 (Fact 0), translated form, at a few b
    for b in ([0, 1, p - 1] if m > 1 else [0]):
        for j in range(m):
            s = sum(ck * pow(b + ak, j, p) for ck, ak in zip(c, A)) % p
            ok("A_fact0_" + tag, s == (1 if j == m - 1 else 0), (p, A, b, j))
    # Lemma 3.1(a): deg F = d, leading coefficient C(D, m-1)
    top = max(l for l in range(len(F)) if F[l] % p)
    ok("A_degF_" + tag, top == d and F[d] == math.comb(D, m - 1) % p, (p, A))
    negA = {(-a) % p for a in A}
    for b in range(p):
        E = [k for k, ak in enumerate(A) if chi[(b + ak) % p] == -1]
        delta = 1 if b in negA else 0
        for j in range(m):
            val = deriv_at(F, j, b, p)
            ff = falling(D, j) % p
            pred = -2 * ff * sum(c[k] * pow(b + A[k], m - 1 - j, p) for k in E)
            if delta and j == m - 1 and m >= 2:
                k0 = [k for k, ak in enumerate(A) if (b + ak) % p == 0][0]
                pred = ff * (-c[k0] - 2 * sum(c[k] for k in E))
            if delta and m == 1:
                pred = -c[0]
            ok("A_lemma31bc_" + tag, val == pred % p, (p, A, b, j))
        # Lemma 3.1(d): ord_b (F - 2 R_E) >= m - delta_b
        G = list(F)
        for k in E:
            ck, ak = c[k], A[k]
            for l in range(D + 1):
                G[l] = (G[l] - 2 * ck * math.comb(D, l) * pow(ak, D - l, p)) % p
        if any(G):
            o = order_at(G, b, p, cap=m)
        else:
            o = m
        ok("A_lemma31d_" + tag, o >= m - delta, (p, A, b, o))


# ---------------------------------------------------------------------------------------------
# B. Full Hankel determinants
# ---------------------------------------------------------------------------------------------

def u_poly(A, c, p, d, s):
    m = len(A)
    D = d + m - 1
    N = D - s
    coef = np.zeros(N + 1, dtype=np.int64)
    for l in range(N + 1):
        tot = 0
        for ck, ak in zip(c, A):
            tot += ck * pow(ak, N - l, p)
        coef[l] = math.comb(N, l) * tot % p
    if s == 0:
        coef[0] = (coef[0] - 1) % p
    return coef


def trim(a):
    a = np.array(a, dtype=np.int64)
    nz = np.nonzero(a)[0]
    if len(nz) == 0:
        return np.zeros(1, dtype=np.int64)
    return a[: nz[-1] + 1]


def pmul(a, b, p):
    return np.convolve(a, b) % p


def padd(a, b, p):
    n = max(len(a), len(b))
    r = np.zeros(n, dtype=np.int64)
    r[: len(a)] += a
    r[: len(b)] += b
    return r % p


def perm_sign(perm):
    s = 1
    seen = [False] * len(perm)
    for i in range(len(perm)):
        if not seen[i]:
            j = i
            L = 0
            while not seen[j]:
                seen[j] = True
                j = perm[j]
                L += 1
            if L % 2 == 0:
                s = -s
    return s


def lambda_mod(p, d, m, e):
    D = d + m - 1
    M = [[math.comb(D - i - j, m - 1) % p for j in range(e + 1)] for i in range(e + 1)]
    return det_mod(M, p)


def det_mod(M, p):
    M = [row[:] for row in M]
    n = len(M)
    det = 1
    for col in range(n):
        piv = None
        for r in range(col, n):
            if M[r][col] % p:
                piv = r
                break
        if piv is None:
            return 0
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            det = -det
        det = det * M[col][col] % p
        inv = pow(M[col][col], p - 2, p)
        for r in range(col + 1, n):
            f = M[r][col] * inv % p
            if f:
                for cc in range(col, n):
                    M[r][cc] = (M[r][cc] - f * M[col][cc]) % p
    return det % p


def coeff_rows(A, c, p, N):
    """Row k: coefficients (low -> high) of c_k (x + a_k)^N mod p."""
    V = np.zeros((len(A), N + 1), dtype=np.int64)
    binom = [math.comb(N, l) % p for l in range(N + 1)]
    for k, (ck, ak) in enumerate(zip(c, A)):
        for l in range(N + 1):
            V[k, l] = ck * binom[l] * pow(ak, N - l, p) % p
    return V


def hankel_data(A, p, smax):
    d = (p - 1) // 2
    m = len(A)
    D = d + m - 1
    c = weights(A, p)
    F, _, _ = poly_F(A, p)
    V = [coeff_rows(A, c, p, D - s) for s in range(smax + 1)]
    us = []
    for s in range(smax + 1):
        u = V[s].sum(axis=0) % p
        if s == 0:
            u[0] = (u[0] - 1) % p
        # u_s = F^{(s)} / (D)_s
        Fs = [F[l + s] * falling(l + s, s) % p for l in range(len(F) - s)]
        lhs = [(int(x) * (falling(D, s) % p)) % p for x in u]
        ok("B_us_is_normalised_derivative", lhs == Fs, (p, A, s))
        ok("B_us_top_vanish", all(int(x) == 0 for x in u[d - s + 1:]), (p, A, s))
        ok("B_us_lead", int(u[d - s]) == math.comb(D - s, m - 1) % p, (p, A, s))
        us.append(u)
    return c, V, us


def check_step1(A, p, chi, V, us, smax):
    m = len(A)
    negA = {(-a) % p for a in A}
    for b in range(p):
        E = [k for k, a in enumerate(A) if chi[(a + b) % p] == -1]
        delta = 1 if b in negA else 0
        for s in range(smax + 1):
            rho = V[s][E].sum(axis=0) % p if E else np.zeros(V[s].shape[1], dtype=np.int64)
            eps = (us[s] - 2 * rho) % p
            need = m - delta - s
            if need <= 0:
                continue
            if eps.any():
                oe = order_at([int(x) for x in eps], b, p, cap=need)
            else:
                oe = need
            ok("B_step1_eps_order", oe >= need, (p, A, b, s))


def check_hankel(A, p, e, chi, us):
    d = (p - 1) // 2
    m = len(A)
    ut = [trim(us[s][: d - s + 1]) for s in range(2 * e + 1)]
    n = e + 1
    H = np.zeros(1, dtype=np.int64)
    for perm in itertools.permutations(range(n)):
        term = np.ones(1, dtype=np.int64)
        for i in range(n):
            term = pmul(term, ut[i + perm[i]], p)
        if perm_sign(perm) < 0:
            term = (-term) % p
        H = padd(H, term, p)
    H = trim(H)
    degH = len(H) - 1
    lam = lambda_mod(p, d, m, e)
    ok("B_H_degree_exact", degH == (e + 1) * (d - e), (p, A, e, degH))
    ok("B_H_leading_coeff_is_Lambda", int(H[-1]) == lam and lam != 0, (p, A, e))
    Hl = [int(x) for x in H]
    negA = {(-a) % p for a in A}
    total_bound = 0
    total_ord = 0
    tight = 0
    for b in range(p):
        eb = sum(1 for a in A if chi[(a + b) % p] == -1)
        delta = 1 if b in negA else 0
        Tb = sum(m - e - i - delta for i in range(eb, e + 1))
        o = order_at(Hl, b, p)
        ok("B_H_order_ge_step3_bound", o >= Tb, (p, A, e, b, o, Tb))
        if o == Tb and Tb > 0 and e >= 1:
            tight += 1
        total_bound += Tb
        total_ord += o
    ok("B_sum_orders_le_degree", total_bound <= total_ord <= degH, (p, A, e))
    return tight


def section_B():
    tight = 0
    npoly = 0
    for p in primes_between(11, 151):
        chi = legendre_table(p)
        Q = sorted({x * x % p for x in range(1, p)})
        cands = []
        cands.append(("interval", tuple(range((p + 1) // 4))))
        cands.append(("Q_and_0", tuple([0] + Q)))  # m = (p+1)/2, D = p - 1
        m1 = int(math.isqrt(p)) + 1
        cands.append(("random_sqrt", tuple(sorted(RNG.sample(range(p), m1)))))
        cands.append(("random_big", tuple(sorted(RNG.sample(range(p), (p + 1) // 2 - 1)))))
        S = [0]
        for x in range(1, p):
            if all(chi[(x + y) % p] >= 0 for y in S) and chi[(2 * x) % p] >= 0:
                S.append(x)
        cands.append(("sum_clique", tuple(S)))
        emax = 4 if p <= 61 else 3
        for name, A in cands:
            m = len(A)
            ee = min(emax, (m - 1) // 2)
            c, V, us = hankel_data(A, p, 2 * ee)
            check_step1(A, p, chi, V, us, 2 * ee)
            for e in range(0, ee + 1):
                tight += check_hankel(A, p, e, chi, us)
                npoly += 1
    WITNESSES["B_polynomials"] = {"count": npoly, "primes": "11..151", "e_max": "4 (p <= 61), 3", "points_with_order_equal_bound_e_ge_1": tight}


# ---------------------------------------------------------------------------------------------
# C. (star), supersaturation, robust HP
# ---------------------------------------------------------------------------------------------

def all_sets_containing_zero(p, m):
    rows = []
    for rest in itertools.combinations(range(1, p), m - 1):
        rows.append((0,) + rest)
    if not rows:
        return np.zeros((0, p), dtype=np.int8)
    ind = np.zeros((len(rows), p), dtype=np.int8)
    for r, A in enumerate(rows):
        ind[r, list(A)] = 1
    return ind


def e_and_delta(ind, p, chi):
    NR = np.zeros((p, p), dtype=np.int32)
    for a in range(p):
        for b in range(p):
            NR[a, b] = 1 if chi[(a + b) % p] == -1 else 0
    eb = ind.astype(np.int32) @ NR
    neg = [(-b) % p for b in range(p)]
    delta = ind[:, neg].astype(np.int32)
    return eb, delta


def section_C():
    star_best = {}
    for p in [5, 7, 11, 13, 17, 19]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        best = (Fraction(0), None)
        for m in range(1, (p + 1) // 2 + 1):
            ind = all_sets_containing_zero(p, m)
            eb, delta = e_and_delta(ind, p, chi)
            nsets = ind.shape[0]
            for e in range(0, (m - 1) // 2 + 1):
                live = eb <= e
                T2 = np.where(live, (e + 1 - eb) * (2 * m - 3 * e - eb - 2 * delta), 0)
                okn("C_star_terms_nonneg", live.sum(), int((T2 < 0).sum()))
                lhs2 = T2.sum(axis=1)
                rhs2 = 2 * (e + 1) * (d - e)
                okn("C_star", nsets, int((lhs2 > rhs2).sum()), (p, m, e))
                if e >= 1:
                    i = int(np.argmax(lhs2))
                    fr = Fraction(int(lhs2[i]), rhs2)
                    if fr > best[0]:
                        best = (fr, {"p": p, "A": [int(x) for x in np.nonzero(ind[i])[0]], "e": e,
                                     "lhs": str(Fraction(int(lhs2[i]), 2)), "rhs": rhs2 // 2})
                # Corollary B, all B at once: worst B = {b : w_b > 0}
                w = (m - 2 * e) * (e + 1 - eb) - (e + 1) * delta
                okn("C_corB_allB", nsets, int((np.maximum(w, 0).sum(axis=1) > (e + 1) * (d - e)).sum()),
                    (p, m, e))
                # Theorem C first and second displays on level sets B = {e_b <= e'}
                for ep in range(0, e + 1):
                    lev = eb <= ep
                    n = lev.sum(axis=1)
                    r = (lev & (delta == 1)).sum(axis=1)
                    lhs = (e + 1 - ep) * ((2 * m - 3 * e - ep) * n - 2 * r)
                    okn("C_thmC_display1", nsets, int((lhs > 2 * (e + 1) * (d - e)).sum()), (p, m, e, ep))
                    if e == ep and 2 * ep <= m - 1:
                        okn("C_second_moment_levelset", nsets,
                            int((n * (m - 1 - 2 * ep) ** 2 > m * (p - m)).sum()), (p, m, ep))
                    lhs2b = (m * n - r) * (e + 1 - ep) * (m - 2 * e)
                    rhs2b = (e + 1) * m * d + 2 * e * r * (e + 1 - ep)
                    okn("C_thmC_display2", nsets, int((lhs2b > rhs2b).sum()), (p, m, e, ep))
            # eta-form (exact), level sets with 8 e' <= m
            for ep in range(0, m // 8 + 1):
                lev = eb <= ep
                n = lev.sum(axis=1)
                r = (lev & (delta == 1)).sum(axis=1)
                bad = 0
                for nn, rr in set(zip(n.tolist(), r.tolist())):
                    bad += 0 if eta_form_ok(m, nn, rr, ep, d) else 1
                okn("C_thmC_eta_form", nsets, bad, (p, m, ep))
        star_best[p] = best
        # size-free eta-form: m > (p+1)/2
        for m in range((p + 1) // 2 + 1, p + 1):
            ind = all_sets_containing_zero(p, m)
            eb, delta = e_and_delta(ind, p, chi)
            ep = m // 8
            lev = eb <= ep
            n = lev.sum(axis=1)
            okn("C_thmC_sizefree_n_le_1", ind.shape[0], int((n > 1).sum()), (p, m))
    WITNESSES["C_star_tightest_e_ge_1"] = {str(p): (str(v[0]), v[1]) for p, v in star_best.items()}
    # the equality witness quoted in the paper
    p, A, e = 13, (0, 2, 3, 5), 1
    chi = legendre_table(p)
    d = 6
    terms = {}
    for b in range(p):
        eb = sum(1 for a in A if chi[(a + b) % p] == -1)
        delta = 1 if (-b) % p in A else 0
        if eb <= e:
            terms[b] = Fraction((e + 1 - eb) * (2 * len(A) - 3 * e - eb - 2 * delta), 2)
    lhs = sum(terms.values())
    ok("C_equality_witness_p13", lhs == 10 == (e + 1) * (d - e)
       and terms == {1: 2, 7: 2, 9: 2, 12: 2, 10: 1, 11: 1}, terms)
    WITNESSES["C_equality_e1"] = {"p": 13, "A": list(A), "e": 1,
                                  "terms": {str(k): str(v) for k, v in terms.items()}, "lhs": str(lhs),
                                  "rhs": 10}
    # the count used in the size-free argument
    for p in primes_between(5, 97):
        chi = legendre_table(p)
        for cc in range(1, p):
            cnt = sum(1 for y in range(p) if chi[y] != -1 and chi[(y + cc) % p] != -1)
            pred = Fraction(p - 3 - chi[cc] - chi[(-cc) % p], 4) + (chi[cc] == 1) + (chi[(-cc) % p] == 1)
            ok("C_jacobsthal_count", cnt == pred and 4 * cnt <= p + 3, (p, cc))
    # larger primes: Theorem C (both displays, eta-form) on random / structured A
    for p in [101, 211, 409]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        Q = sorted({x * x % p for x in range(1, p)})
        fams = [tuple(range(m)) for m in (8, 20, (p + 1) // 4, (p + 1) // 2)]
        fams += [tuple([0] + Q)]
        fams += [tuple(sorted(RNG.sample(range(p), RNG.randint(8, (p + 1) // 2)))) for _ in range(12)]
        for A in fams:
            m = len(A)
            negA = {(-a) % p for a in A}
            ebs = [sum(1 for a in A if chi[(a + b) % p] == -1) for b in range(p)]
            for ep in range(0, (m - 1) // 2 + 1):
                lev = [b for b in range(p) if ebs[b] <= ep]
                n = len(lev)
                r = sum(1 for b in lev if b in negA)
                if 8 * ep <= m:
                    ok("C_thmC_eta_form_large_p", eta_form_ok(m, n, r, ep, d), (p, m, ep))
                for e in range(ep, (m - 1) // 2 + 1, max(1, (m - 1) // 12)):
                    ok("C_thmC_display1_large_p",
                       (e + 1 - ep) * ((2 * m - 3 * e - ep) * n - 2 * r) <= 2 * (e + 1) * (d - e), (p, m, e, ep))


def eta_form_ok(m, n, r, ep, d):
    """Exact test of  m n - r <= (1 - sqrt(2 ep/m))^{-2} d + m  for 0 <= ep/m <= 1/8."""
    X = m * n - r - m
    if X <= 0:
        return True
    # need sqrt(2ep/m) + sqrt(d/X) >= 1
    a = Fraction(2 * ep, m)
    bb = Fraction(d, X)
    if a >= 1 or bb >= 1:
        return True
    t = 1 + bb - a
    if t <= 0:
        return True
    return 4 * bb >= t * t


# ---------------------------------------------------------------------------------------------
# D. Constant bias
# ---------------------------------------------------------------------------------------------

def factor_float(kappa, p):
    u = 1.0 / np.sqrt(1.0 + 2.0 * kappa)
    return 1.0 - (1.0 - u) ** 2 + (np.sqrt((0.5 + kappa) * p) + 1.0) / (2.0 * (p - 1))


def core_bracket_float(kappa, m, d):
    u = 1.0 / np.sqrt(1.0 + 2.0 * kappa)
    return 1.0 - (1.0 - u) ** 2 + u * (1.0 - u) * m / d


def check_bias_rows(p, m_arr, n_arr, smax_arr, smin_arr, tag, grid=200):
    """Final bound of Theorem D and the core lemma, float with margin 1e-9."""
    d = (p - 1) / 2.0
    mn = m_arr * n_arr
    K = np.minimum(1.5, mn / p - 0.5)
    adm = K > 0
    absS = np.maximum(smax_arr, -smin_arr)
    # final bound: unique (mn, |S|)
    key = np.unique(np.stack([mn[adm], absS[adm]], axis=1), axis=0)
    bad = 0
    min_slack_nonvac = None
    for t in np.linspace(0, 1, grid + 1)[1:]:
        Kk = np.minimum(1.5, key[:, 0] / p - 0.5) * t
        fac = factor_float(Kk, p)
        slack = fac * key[:, 0] - key[:, 1]
        bad += int((slack < -1e-9).sum())
        nv = fac < 1
        if nv.any():
            s = float(slack[nv].min())
            min_slack_nonvac = s if min_slack_nonvac is None else min(min_slack_nonvac, s)
    # one check per admissible (A, n) pair; each is tested at `grid` values of kappa in (0, K]
    okn("D_thmD_final_" + tag, int(adm.sum()), bad, (p, "grid", grid))
    # core lemma (m <= (p+1)/2): upper bound for S at kappa = K (concavity in u), and for -S
    sel = adm & (m_arr <= (p + 1) // 2)
    for sgn, S in ((1, smax_arr), (-1, -smin_arr)):
        key2 = np.unique(np.stack([m_arr[sel], n_arr[sel], S[sel]], axis=1), axis=0)
        if len(key2) == 0:
            continue
        Kk = np.minimum(1.5, key2[:, 0] * key2[:, 1] / p - 0.5)
        br = core_bracket_float(Kk, key2[:, 0], d)
        viol = key2[:, 2] > br * key2[:, 0] * key2[:, 1] + 1e-9
        okn("D_core_lemma_" + tag, int(sel.sum()), int(viol.sum()), p)
    return min_slack_nonvac


def section_D():
    info = {}
    for p in [11, 13, 17, 19]:
        chi = np.array(legendre_table(p), dtype=np.int32)
        CH = np.zeros((p, p), dtype=np.int32)
        for a in range(p):
            for b in range(p):
                CH[a, b] = chi[(a + b) % p]
        # all A containing 0 (translation (A,B) -> (A+t, B-t) preserves S and sizes)
        idx = np.arange(2 ** (p - 1), dtype=np.int64)
        rows = ((idx[:, None] >> np.arange(p - 1)) & 1).astype(np.int8)
        ind = np.concatenate([np.ones((rows.shape[0], 1), dtype=np.int8), rows], axis=1)
        m = ind.sum(axis=1).astype(np.int64)
        rs = ind.astype(np.int32) @ CH  # row sums s_A(b)
        srt = np.sort(rs, axis=1)
        smax = np.cumsum(srt[:, ::-1], axis=1)  # n = 1..p
        smin = np.cumsum(srt, axis=1)
        n = np.arange(1, p + 1)
        M = np.repeat(m[:, None], p, axis=1).ravel()
        N = np.tile(n, ind.shape[0])
        slack = check_bias_rows(p, M, N, smax.ravel().astype(np.int64), smin.ravel().astype(np.int64),
                                "exhaustive")
        info[p] = {"sets": int(ind.shape[0]), "pairs": int(ind.shape[0] * p),
                   "min_slack_nonvacuous": slack}
    WITNESSES["D_exhaustive"] = info
    # larger primes: |A| = 2 family (closest cases) and random small A, best B for each n
    big = {}
    for p in [101, 211, 409, 1009, 2003, 10009]:
        chi = np.array(legendre_table(p), dtype=np.int64)
        fams = [(0, 1)]
        if p <= 2003:
            for _ in range(20):
                mm = RNG.randint(2, int(math.isqrt(2 * p)) + 1)
                fams.append(tuple(RNG.sample(range(p), mm)))
            fams.append(tuple(range(int(math.isqrt(p)))))
        worst = None
        for A in fams:
            bidx = np.arange(p)
            rs = np.zeros(p, dtype=np.int64)
            for a in A:
                rs += chi[(bidx + a) % p]
            srt = np.sort(rs)
            smax = np.cumsum(srt[::-1])
            smin = np.cumsum(srt)
            nn = np.arange(1, p + 1)
            mm = np.full(p, len(A))
            s = check_bias_rows(p, mm, nn, smax, smin, "large_p", grid=100)
            if s is not None:
                worst = s if worst is None else min(worst, s)
        big[p] = worst
    WITNESSES["D_large_p_min_slack_nonvacuous"] = big
    # non-vacuity thresholds (mpmath, 40 digits): factor < 1 iff eps_p < (1-u)^2; eps decreasing in p
    mpmath.mp.dps = 40

    def fac(kappa, p):
        kappa = mpmath.mpf(kappa)
        u = 1 / mpmath.sqrt(1 + 2 * kappa)
        return 1 - (1 - u) ** 2 + (mpmath.sqrt((mpmath.mpf(1) / 2 + kappa) * p) + 1) / (2 * (p - 1))

    thr = {}
    for kap, expect in ((Fraction(1, 10), 2741), (Fraction(1, 2), 47), (Fraction(3, 2), 17)):
        kk = mpmath.mpf(kap.numerator) / kap.denominator
        first = None
        for q in primes_between(11, 5000):
            if fac(kk, q) < 1:
                first = q
                break
        prev = max(x for x in primes_between(11, first) if x < first) if first > 11 else None
        ok("D_nonvacuity_threshold", first == expect and (prev is None or fac(kk, prev) >= 1),
           (str(kap), first))
        thr[str(kap)] = first
    # monotonicity of eps in p (derivative sign) checked symbolically
    P, c = sympy.symbols("P c", positive=True)
    der = sympy.diff((sympy.sqrt(c * P) + 1) / (P - 1), P)
    num = sympy.simplify(der * (P - 1) ** 2)
    ok("D_eps_decreasing_in_p", sympy.simplify(num + sympy.sqrt(c) * (P + 1) / (2 * sympy.sqrt(P)) + 1) == 0,
       str(num))
    # vacuous for p = 11, 13 at every kappa in (0, 3/2]: both eps_p(kappa) and (1-u(kappa))^2 are
    # increasing in kappa, so on [k_i, k_{i+1}] it suffices that eps_p(k_i) >= (1-u(k_{i+1}))^2
    # (monotone interval argument on a 3000-point grid, 40-digit arithmetic; k_0 = 0)

    def eps_(kappa, q):
        return (mpmath.sqrt((mpmath.mpf(1) / 2 + kappa) * q) + 1) / (2 * (q - 1))

    def gap_(kappa):
        return (1 - 1 / mpmath.sqrt(1 + 2 * kappa)) ** 2

    for q in (11, 13):
        grid_ = [mpmath.mpf(3) / 2 * t / 3000 for t in range(0, 3001)]
        worst = min(eps_(grid_[i], q) - gap_(grid_[i + 1]) for i in range(3000))
        ok("D_vacuous_p_le_13_interval", worst >= 0, (q, float(worst)))
        thr[f"p{q}_min_interval_margin"] = float(worst)
    WITNESSES["D_thresholds"] = thr
    # m0 <= sqrt(2p) + 1 <= (p+1)/2 for p >= 11, fails at p = 7
    for q in primes_between(3, 20000):
        holds = (math.sqrt(2 * q) + 1 <= (q + 1) / 2)
        ok("D_m0_bound", holds == (q >= 11), q)
        if q >= 11:
            ok("D_m0_bound_ge11", 8 * q <= (q - 1) ** 2, q)
    ok("D_m0_bound_fails_p7", not (math.sqrt(14) + 1 <= 4), 7)


# ---------------------------------------------------------------------------------------------
# E. Sharpness family
# ---------------------------------------------------------------------------------------------

def section_E():
    for p in primes_between(5, 3000):
        chi = legendre_table(p)
        A = (0, 1)
        eb = [sum(1 for a in A if chi[(a + b) % p] == -1) for b in range(p)]
        N0 = eb.count(0)
        N1 = eb.count(1)
        B0 = [b for b in range(p) if eb[b] == 0]
        r = sum(1 for b in B0 if (-b) % p in A)
        S = sum(chi[(a + b) % p] for a in A for b in B0)
        if p % 4 == 1:
            ok("E_pair_p1mod4", N0 == (p + 3) // 4 and r == 2 and S == (p - 1) // 2 and N1 == (p - 1) // 2,
               (p, N0, r, S, N1))
            # Corollary B with e = 0 for m = 2 gives S <= d for every B
        else:
            ok("E_pair_p3mod4", N0 == (p + 1) // 4 and r == 1 and S == (p - 1) // 2, (p, N0, r, S))
        # m = 2: S(A,B) <= d for all B (max over B = sum of positive row sums)
        rsum = [chi[b % p] + chi[(b + 1) % p] for b in range(p)]
        ok("E_m2_S_le_d", sum(x for x in rsum if x > 0) <= (p - 1) // 2, p)
    # padding: saving tends to 2k/(1+2k) for k <= 1
    pad = {}
    for p in [10009, 100049]:
        if not sympy.isprime(p) or p % 4 != 1:
            continue
        chi = legendre_table(p)
        eb = [(chi[b] == -1) + (chi[(b + 1) % p] == -1) for b in range(p)]
        B0 = [b for b in range(p) if eb[b] == 0]
        B1 = [b for b in range(p) if eb[b] == 1]
        for kap in (Fraction(1, 20), Fraction(1, 10), Fraction(1, 4), Fraction(1, 2), Fraction(1, 1)):
            n = math.ceil((Fraction(1, 2) + kap) * p / 2)
            extra = n - len(B0)
            feasible = 0 <= extra <= len(B1)
            B = B0 + B1[:max(extra, 0)]
            S = sum(chi[(a + b) % p] for a in (0, 1) for b in B)
            saving = 1 - Fraction(S, 2 * len(B))
            target = Fraction(2 * kap, 1 + 2 * kap)
            ok("E_padding", feasible and 2 * len(B) >= (Fraction(1, 2) + kap) * p and S == (p - 1) // 2
               and abs(float(saving - target)) < 8.0 / p, (p, str(kap), float(saving), float(target)))
            pad[f"{p}:{kap}"] = [float(saving), float(target)]
    WITNESSES["E_padding_saving_vs_2k_over_1p2k"] = pad


# ---------------------------------------------------------------------------------------------
# F. Leading coefficient Lambda
# ---------------------------------------------------------------------------------------------

def bareiss(M):
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = None
            for r in range(k + 1, n):
                if M[r][k] != 0:
                    sw = r
                    break
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


def section_F():
    fact = math.factorial
    # (i) Krattenthaler product (with sign) equals the integer determinant
    for D in range(0, 37):
        for m in range(1, D + 2):
            d = D - m + 1
            if d < 0:
                continue
            for e in range(0, (m - 1) // 2 + 1):
                if e > d:
                    continue
                M = [[math.comb(D - i - j, m - 1) for j in range(e + 1)] for i in range(e + 1)]
                lam = bareiss(M)
                Np = D - 2 * e
                c = m - 1 - e
                n = e + 1
                num = 1
                for k in range(1, n):
                    num *= fact(k)
                for i in range(1, n + 1):
                    num *= fact(Np + i - 1)
                den = 1
                for i in range(1, n + 1):
                    den *= fact(c - i + n) * fact(Np - c + i - 1)
                sgn = (-1) ** (e * (e + 1) // 2)
                ok("F_krattenthaler_product", num % den == 0 and lam == sgn * (num // den), (D, m, e))
    # (ii) p does not divide Lambda, Paley d, all m <= (p+1)/2, e <= min((m-1)/2, 5)
    for p in primes_between(3, 300):
        d = (p - 1) // 2
        for m in range(1, (p + 1) // 2 + 1):
            for e in range(0, min((m - 1) // 2, 5) + 1):
                ok("F_lambda_nonzero_mod_p_paley", lambda_mod(p, d, m, e) != 0, (p, m, e))
    # general d (as in the formalization): 2e <= d, 2e+1 <= m, d + m - 1 < p
    for p in primes_between(3, 41):
        for d in range(1, p):
            for m in range(1, p - d + 1):
                for e in range(0, min(d // 2, (m - 1) // 2, 4) + 1):
                    ok("F_lambda_nonzero_mod_p_general_d", lambda_mod(p, d, m, e) != 0, (p, d, m, e))
    # counterexample without D < p: d = 2, m = 3, e = 1 gives Lambda = -3
    lam = bareiss([[math.comb(4 - i - j, 2) for j in range(2)] for i in range(2)])
    ok("F_counterexample_D_ge_p", lam == -3, lam)
    WITNESSES["F_counterexample"] = {"d": 2, "m": 3, "e": 1, "D": 4, "Lambda": lam, "p": 3}
    # (iii) the elementary argument
    for d in range(1, 25):
        for K in range(0, 25):
            D = d + K
            for e in range(0, min(d // 2, K) + 1):
                for i in range(e + 1):
                    for j in range(e + 1):
                        lhs = math.comb(D - i - j, K) * fact(K) * fact(d - i)
                        rhs = fact(D - i - e) * falling(d - i, j) * falling(D - i - j, e - j)
                        ok("F_elem_factorisation", lhs == rhs, (d, K, e, i, j))
    for p in primes_between(3, 200):
        for _ in range(8):
            d = RNG.randint(1, p - 1)
            K = RNG.randint(0, p - 1 - d)  # D = d + K <= p - 1
            emax = min(d // 2, K, 6)
            e = RNG.randint(0, emax)
            D = d + K

            def P(j, x):
                return falling(d - x, j) * falling(D - j - x, e - j)

            # triangularity at the nodes d - k
            for k in range(e + 1):
                for j in range(e + 1):
                    v = P(j, d - k) % p
                    if j > k:
                        ok("F_elem_triangular_zero", v == 0, (p, d, K, e, k, j))
                    if j == k:
                        ok("F_elem_diagonal_unit", v == fact(k) * falling(K, e - k) % p and v != 0,
                           (p, d, K, e, k))
            # Lambda = prod alpha_i * det[P_j(i)]  (mod p)
            lam = det_mod([[math.comb(D - i - j, K) % p for j in range(e + 1)] for i in range(e + 1)], p)
            alpha = 1
            for i in range(e + 1):
                alpha = alpha * fact(D - i - e) * pow(fact(K) * fact(d - i) % p, p - 2, p) % p
            dP = det_mod([[P(j, i) % p for j in range(e + 1)] for i in range(e + 1)], p)
            ok("F_elem_lambda_factorises", lam == alpha * dP % p and dP != 0, (p, d, K, e))


# ---------------------------------------------------------------------------------------------
# G. Constants and expansions
# ---------------------------------------------------------------------------------------------

def section_G():
    k = sympy.symbols("kappa", positive=True)
    u = 1 / sympy.sqrt(1 + 2 * k)
    ser = sympy.series((1 - u) ** 2, k, 0, 4).removeO()
    ok("G_series_(1-u)^2", sympy.expand(ser - (k ** 2 - 3 * k ** 3)) == 0, str(ser))
    # u in [1/2, 1) and theta in (0, 1/4] for 0 < kappa <= 3/2
    ok("G_u_range", sympy.nsimplify(u.subs(k, sympy.Rational(3, 2))) == sympy.Rational(1, 2), None)
    # (1+2k) d < (1/2+k) p  (i.e. mn > d/u^2)
    P = sympy.symbols("P", positive=True)
    ok("G_mn_gt_d_over_u2", sympy.simplify((sympy.Rational(1, 2) + k) * P - (1 + 2 * k) * (P - 1) / 2
                                            - (sympy.Rational(1, 2) + k)) == 0, None)
    # crossing of (1 - sqrt(2 eta))^-2 and 2/(1 - 2 eta)^2 at eta = (3 - 2 sqrt 2)/2
    eta = sympy.symbols("eta", positive=True)
    root = (3 - 2 * sympy.sqrt(2)) / 2
    f1 = (1 - sympy.sqrt(2 * eta)) ** -2
    f2 = 2 / (1 - 2 * eta) ** 2
    ok("G_crossing_eta", sympy.simplify((f1 - f2).subs(eta, root)) == 0 and abs(float(root) - 0.0858) < 5e-5,
       float(root))
    ok("G_crossing_side", float((f1 - f2).subs(eta, sympy.Rational(1, 20))) < 0
       and float((f1 - f2).subs(eta, sympy.Rational(1, 9))) > 0, None)
    WITNESSES["G_crossing_eta"] = float(root)
    # Section 9 (sharpened constant): eta_star = 1 - u(w), w = 1/(1+2k)
    w = 1 / (1 + 2 * k)
    uw = (sympy.sqrt(12 * w - 3 * w ** 2) - w) / 2
    ser2 = sympy.series(1 - uw, k, 0, 4).removeO()
    ok("G_eta_star_series", sympy.expand(ser2 - (sympy.Rational(4, 3) * k ** 2 - sympy.Rational(40, 9) * k ** 3)) == 0,
       str(ser2))
    table = {Fraction(1, 100): (1.290e-4, 0.971e-4), Fraction(1, 10): (9.838e-3, 7.591e-3),
             Fraction(1, 2): (0.10436, 0.08579), Fraction(1, 1): (0.20924, 0.17863),
             Fraction(3, 2): (0.28647, 0.25)}
    vals = {}
    for kap, (es, rr) in table.items():
        kk = sympy.Rational(kap.numerator, kap.denominator)
        e1 = float((1 - uw).subs(k, kk))
        e2 = float(((1 - u) ** 2).subs(k, kk))
        ok("G_table_values", abs(e1 - es) <= 5e-4 * max(es, 1e-3) + 1e-7 and abs(e2 - rr) <= 5e-4 * max(rr, 1e-3) + 1e-7,
           (str(kap), e1, e2))
        vals[str(kap)] = [e1, e2]
    WITNESSES["G_eta_star_and_(1-u)^2"] = vals
    # beta_m(q): mean of 1 - 2j/m over the lowest mass q of Binomial(m, 1/2)

    def beta(mm, q):
        left = q
        tot = Fraction(0)
        for j in range(mm + 1):
            mass = Fraction(math.comb(mm, j), 2 ** mm)
            take = min(mass, left)
            tot += take * (1 - Fraction(2 * j, mm))
            left -= take
            if left <= 0:
                break
        return tot / q

    wv = Fraction(1, 2)  # kappa = 1/2
    b5 = beta(5, 1 / (2 * 5 * wv))
    ok("G_beta5_kappa_half", b5 == Fraction(51, 80) and 1 - b5 == Fraction(29, 80), str(b5))
    WITNESSES["G_saving_cap_kappa_half_m5"] = str(1 - b5)
    # beta_2 = w for w >= 1/3 and beta_5 = 3/5 + w/8 on [8/15, 1]
    for num in range(1, 200):
        wv = Fraction(1, 3) + Fraction(2, 3) * Fraction(num, 200)
        ok("G_beta2_equals_w", beta(2, 1 / (4 * wv)) == wv, str(wv))
        if wv >= Fraction(8, 15):
            ok("G_beta5_formula", beta(5, 1 / (10 * wv)) == Fraction(3, 5) + wv / 8, str(wv))
    ok("G_11_over_48", (Fraction(3, 5) + Fraction(24, 35) / 8 == Fraction(24, 35))
       and Fraction(1, 1 + 2 * Fraction(11, 48)) == Fraction(24, 35), None)


# ---------------------------------------------------------------------------------------------
# H. Paley-graph density corollary
# ---------------------------------------------------------------------------------------------

def corE_bounds(kappa, p, m):
    u = 1.0 / math.sqrt(1 + 2 * kappa)
    c = (1 - u) ** 2
    eps = (math.sqrt((0.5 + kappa) * p) + 1) / (2 * (p - 1))
    Phi = 1 - c + eps
    lo = c / 2 - eps / 2 - Phi / (2 * (m - 1))
    hi = 1 - c / 2 + eps / 2 + Phi / (2 * (m - 1))
    return lo, hi


def section_H():
    for p in [13, 17, 29, 37, 41]:
        chi = legendre_table(p)
        adj = np.zeros((p, p), dtype=np.int32)
        for x in range(p):
            for y in range(p):
                adj[x, y] = 1 if chi[(x - y) % p] == 1 else 0
        mmax = {13: 7, 17: 8, 29: 6, 37: 5, 41: 5}[p]
        for m in range(3, mmax + 1):
            if m * m < p / 2:
                continue
            combs = np.array(list(itertools.combinations(range(1, p), m - 1)), dtype=np.int64)
            A = np.concatenate([np.zeros((combs.shape[0], 1), dtype=np.int64), combs], axis=1)
            E = np.zeros(A.shape[0], dtype=np.int64)
            for i in range(m):
                for j in range(i + 1, m):
                    E += adj[A[:, i], A[:, j]]
            Emin, Emax = int(E.min()), int(E.max())
            # identity S(A,-A) = 4E - m(m-1) on a sample
            for t in range(0, A.shape[0], max(1, A.shape[0] // 50)):
                S = sum(chi[(a - b) % p] for a in A[t] for b in A[t])
                ok("H_identity_S_4E", S == 4 * int(E[t]) - m * (m - 1), (p, m))
            npairs = m * (m - 1) // 2
            Kmax = min(1.5, m * m / p - 0.5)
            if Kmax <= 0:
                continue
            for t in range(1, 101):
                kap = Kmax * t / 100
                lo, hi = corE_bounds(kap, p, m)
                ok("H_corE_exhaustive", lo - 1e-12 <= Emin / npairs and Emax / npairs <= hi + 1e-12,
                   (p, m, kap))
    # larger p: greedy dense and sparse sets
    dens = {}
    for p in [101, 401, 1009]:
        chi = legendre_table(p)
        for kap in (0.05, 0.25, 0.5, 1.0, 1.5):
            m = math.ceil(math.sqrt((0.5 + kap) * p))
            for sign in (1, -1):
                A = [0]
                cand = list(range(1, p))
                RNG.shuffle(cand)
                while len(A) < m:
                    best = max(cand, key=lambda x: sum(sign * chi[(x - a) % p] for a in A))
                    A.append(best)
                    cand.remove(best)
                E = sum(1 for i in range(m) for j in range(i + 1, m) if chi[(A[i] - A[j]) % p] == 1)
                rho = E / (m * (m - 1) / 2)
                lo, hi = corE_bounds(kap, p, m)
                ok("H_corE_greedy", lo - 1e-12 <= rho <= hi + 1e-12, (p, kap, sign, rho, lo, hi))
                dens[f"{p}:{kap}:{'dense' if sign > 0 else 'sparse'}"] = [round(rho, 4), round(lo, 4), round(hi, 4)]
    WITNESSES["H_greedy_density_vs_bounds"] = dens


# ---------------------------------------------------------------------------------------------
# I. Second moment
# ---------------------------------------------------------------------------------------------

def section_I():
    for p in primes_between(5, 97):
        chi = legendre_table(p)
        for c in range(1, p):
            ok("I_jacobsthal_minus1", sum(chi[x] * chi[(x + c) % p] for x in range(p)) == -1, (p, c))
        for _ in range(6):
            n = RNG.randint(1, p)
            B = RNG.sample(range(p), n)
            tot = sum(sum(chi[(x + b) % p] for b in B) ** 2 for x in range(p))
            ok("I_second_moment", tot == n * (p - n), (p, n))
            m = RNG.randint(1, p)
            A = RNG.sample(range(p), m)
            S = sum(chi[(a + b) % p] for a in A for b in B)
            ok("I_vinogradov_bound", S * S <= m * n * (p - n) <= m * n * p, (p, m, n))


def main():
    timing = {}
    for name, fn in [("A", section_A), ("B", section_B), ("C", section_C), ("D", section_D),
                     ("E", section_E), ("F", section_F), ("G", section_G), ("H", section_H),
                     ("I", section_I)]:
        t = time.time()
        fn()
        timing[name] = round(time.time() - t, 1)
        print(f"section {name}: {timing[name]} s, checks so far {sum(CHECKS.values())}, "
              f"failures {len(FAILURES)}", flush=True)
    out = {
        "description": "Verifier for paper/robust-hanson-petridis.md; exact unless the key says float "
                       "(D_* bias checks and H_* density checks use float64 with margin 1e-9 / 1e-12).",
        "checks": CHECKS,
        "total_checks": sum(CHECKS.values()),
        "failures": FAILURES,
        "n_failures": len(FAILURES),
        "witnesses": WITNESSES,
        "notes": NOTES,
        "timing_seconds": timing,
        "elapsed_seconds": round(time.time() - T0, 1),
        "cpu_seconds": round(time.process_time(), 1),
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print(json.dumps({"total_checks": out["total_checks"], "n_failures": out["n_failures"],
                      "elapsed": out["elapsed_seconds"]}))
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    sys.exit(main())
