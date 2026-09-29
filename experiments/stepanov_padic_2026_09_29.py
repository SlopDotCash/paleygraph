#!/usr/bin/env python3
"""Verifier for research/stepanov-padic-2026-09-29.md (worker `padic`).

Lifting the Hanson-Petridis (HP) polynomial and S(A,B) from F_p to Z/p^2 and Z_p.
Everything is exact integer arithmetic (Python ints, numpy int64 only for residues < p^2).

Sections (match the note):
  A  Fermat quotients: refined Euler criterion, shift rule, homomorphism, Teichmuller kernel.
  B  The Fermat-quotient cocycle, the truncated logarithm L1, decomposition of q(a+b),
     degree of the canonical-lift function u.
  C  The lifted HP polynomial: degree, Taylor data mod p^2 at every point, complete points,
     lift independence, Teichmuller decomposition, extra mod-p^2 vanishing statistics,
     Newton polygons (Eisenstein clusters), the Frobenius-exponent polynomial F_*.
  D  S(A,B) mod p^2: the identity, recovery of S from its residue, decomposition of T,
     the translation-average identity.
  E  Orbit test: T and the second-digit data along the symmetry orbit (sA+t, sB-t), s in Q.
  F  Randomness diagnostics (HEURISTIC): value distribution, ranks, joint (chi, q) counts,
     degree of the Teichmuller second-digit function.
  G  Height illustration: p-adic valuation versus archimedean size of Disc(F_A) over Q.

Run:  /opt/miniconda3/bin/python3 experiments/stepanov_padic_2026_09_29.py
Writes results/stepanov_padic_2026_09_29.json
"""
import json
import math
import os
import random
import sys
import time
from math import comb

import numpy as np

T0 = time.time()
random.seed(20260929)
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "stepanov_padic_2026_09_29.json")

CHECKS = {}
FAILS = []
RES = {}


def check(name, cond, info=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not cond:
        if len(FAILS) < 200:
            FAILS.append({"check": name, "info": info})


def primes_upto(N):
    s = bytearray([1]) * (N + 1)
    s[0:2] = b"\x00\x00"
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(N + 1) if s[i]]


PRIMES = primes_upto(4000)


def leg(x, p):
    r = pow(x % p, (p - 1) // 2, p)
    return 0 if r == 0 else (1 if r == 1 else -1)


def fq(x, p):
    """Fermat quotient (x^(p-1)-1)/p mod p of an integer x with p not dividing x."""
    M = p * p
    r = pow(x % M, p - 1, M)
    assert (r - 1) % p == 0
    return ((r - 1) // p) % p


def ucan(x, p):
    """u(x) = (x^p - x)/p mod p for the canonical representative x in [0,p)."""
    M = p * p
    return ((pow(x, p, M) - x) // p) % p


def inv(x, m):
    return pow(x % m, -1, m)


def vp(x, p, cap):
    """p-adic valuation of an integer residue known mod p^cap (returns cap if 0)."""
    if x == 0:
        return cap
    v = 0
    while x % p == 0:
        x //= p
        v += 1
    return v


_L1 = {}


def L1_table(p):
    """L1(s) = sum_{k=1}^{p-1} s^k/k mod p, computed directly (numpy)."""
    if p in _L1:
        return _L1[p]
    s = np.arange(p, dtype=np.int64)
    acc = np.zeros(p, dtype=np.int64)
    pw = np.ones(p, dtype=np.int64)
    for k in range(1, p):
        pw = (pw * s) % p
        acc = (acc + pw * pow(k, -1, p)) % p
    _L1[p] = [int(v) for v in acc]
    return _L1[p]


def chi_table(p):
    t = np.zeros(p, dtype=np.int64)
    sq = set((x * x) % p for x in range(1, p))
    for x in range(1, p):
        t[x] = 1 if x in sq else -1
    return t


def partner_set(A, p, chit=None):
    if chit is None:
        chit = chi_table(p)
    x = np.arange(p)
    ok = np.ones(p, dtype=bool)
    for a in A:
        ok &= chit[(x + a) % p] >= 0
    return [int(b) for b in np.nonzero(ok)[0]]


# ----------------------------------------------------------------------------------------
# Section A: Fermat quotients
# ----------------------------------------------------------------------------------------
def section_A():
    t = time.time()
    for p in [q for q in PRIMES if 3 <= q <= 61]:
        M = p * p
        d = (p - 1) // 2
        inv2 = (p + 1) // 2
        ker = 0
        vals = set()
        for x in range(1, M):
            if x % p == 0:
                continue
            q = fq(x, p)
            ch = leg(x, p)
            vals.add(q)
            # A1 refined Euler: x^d = chi(x)(1 + p q/2) mod p^2
            check("A1.refined_euler", pow(x, d, M) == (ch * (1 + p * (q * inv2 % p))) % M, (p, x))
            # A1b: x^(d + j(p-1)) = chi(x)(1 + p (j+1/2) q) mod p^2
            j = random.randrange(0, 3 * p)
            check("A1b.exponent_family",
                  pow(x, d + j * (p - 1), M) == (ch * (1 + p * ((2 * j + 1) * inv2 * q % p))) % M, (p, x, j))
            # A1c: x^(p d) = chi(x) mod p^2
            check("A1c.frobenius_exponent", pow(x, p * d, M) == ch % M, (p, x))
            # A2 shift rule q(x + p k) = q(x) - k/x
            k = random.randrange(-3 * p, 3 * p)
            check("A2.shift_rule", fq(x + p * k, p) == (q - k * inv(x, p)) % p, (p, x, k))
            # A3 homomorphism
            y = random.randrange(1, 5 * M)
            if y % p:
                check("A3.homomorphism", fq(x * y, p) == (q + fq(y, p)) % p, (p, x, y))
            # A4 Teichmuller: omega(x) = x^p = x(1 + p q(x)) mod p^2 and q(omega) = 0
            w = pow(x, p, M)
            check("A4.teichmuller", w == (x * (1 + p * q)) % M and fq(w, p) == 0, (p, x))
            if q == 0:
                ker += 1
        check("A4.kernel_is_mu_{p-1}", ker == p - 1, (p, ker))
        check("A5.surjective", len(vals) == p, (p, len(vals)))
    RES["A_time"] = round(time.time() - t, 2)


# ----------------------------------------------------------------------------------------
# Section B: cocycle, L1, decomposition of q(a+b), degree of u
# ----------------------------------------------------------------------------------------
def poly_degree_of_function(vals, p):
    """Degree of the interpolating polynomial (deg <= p-1) of a function F_p -> F_p.
    Coefficient of x^k (1 <= k <= p-1) is -sum_x f(x) x^(p-1-k)."""
    v = np.array(vals, dtype=np.int64) % p
    x = np.arange(p, dtype=np.int64)
    # powers x^e for e = 0..p-2
    deg = 0
    pw = np.ones(p, dtype=np.int64)
    pw[0] = 1  # 0^0 = 1
    coeffs = {}
    # e = p-1-k runs 0..p-2 as k runs p-1..1
    for e in range(0, p - 1):
        if e == 0:
            pe = np.ones(p, dtype=np.int64)
        else:
            pe = (pw * x) % p
            pw = pe
        c = int((-(v * pe).sum()) % p)
        k = p - 1 - e
        coeffs[k] = c
    for k in range(p - 1, 0, -1):
        if coeffs[k] != 0:
            return k
    return 0


def section_B():
    t = time.time()
    # B1 exact integer cocycle f(x+y) = f(x) + f(y) + C_p(x,y), f(x) = (x^p - x)/p = x q(x)
    for p in [q for q in PRIMES if 3 <= q <= 31]:
        cpk = [comb(p, k) // p for k in range(p)]
        for _ in range(60):
            x = random.randrange(-3 * p * p, 3 * p * p)
            y = random.randrange(-3 * p * p, 3 * p * p)
            f = lambda z: (z ** p - z) // p
            C = sum(cpk[k] * x ** k * y ** (p - k) for k in range(1, p))
            check("B1.exact_cocycle", f(x + y) == f(x) + f(y) + C, (p, x, y))
            if x % p and y % p and (x + y) % p:
                # consequence: (x+y) q(x+y) = x q(x) + y q(y) + C  (mod p)
                check("B1b.cocycle_mod_p",
                      ((x + y) * fq(x + y, p) - x * fq(x, p) - y * fq(y, p) - C) % p == 0, (p, x, y))
    # B2..B5 over all pairs of canonical representatives
    stats = {}
    for p in [q for q in PRIMES if 3 <= q <= 151]:
        M = p * p
        L1 = L1_table(p)
        inv2 = (p + 1) // 2
        u = [ucan(x, p) for x in range(p)]
        invt = [0] + [inv(k, p) for k in range(1, p)]
        # L1 fast formula: L1(s) = -(s^p + (1-s)^p - 1)/p mod p
        for s in range(p):
            fast = (-((pow(s, p, M) + pow((1 - s) % M, p, M) - 1) // p)) % p
            check("B5a.L1_formula", fast == L1[s], (p, s))
            check("B5b.L1_symmetry", L1[s] == L1[(1 - s) % p], (p, s))
        check("B5c.sum_L1_is_1", sum(L1) % p == 1 % p, p)
        for s in range(1, p):
            check("B5h.L1_inversion", L1[inv(s, p)] == (-L1[s] * inv(s, p)) % p, (p, s))
        check("B5d.sum_u_is_half", sum(u) % p == inv2 % p, p)
        for x in range(p - 1):
            check("B5e.u_difference", (u[x + 1] - u[x]) % p == (-L1[(-x) % p]) % p, (p, x))
        # sum_{x<p} x^p = 0 mod p^2
        check("B5f.sum_xp_mod_p2", sum(pow(x, p, M) for x in range(p)) % M == 0, p)
        # degree of u and of L1 as functions on F_p
        if p <= 101:
            du = poly_degree_of_function(u, p)
            check("B6.deg_u_is_p-1", du == p - 1, (p, du))
            dl = poly_degree_of_function(L1, p)
            check("B6b.deg_L1_is_p-1", dl == p - 1, (p, dl))
        omega = [pow(a, p, M) for a in range(p)]
        carries = 0
        for a in range(p):
            for b in range(p):
                if (a + b) % p == 0:
                    continue
                z = a + b
                Q = fq(z, p)
                s = a * invt[(a + b) % p] % p
                pred = (L1[s] + (u[a] + u[b]) * invt[(a + b) % p]) % p
                check("B3.decomposition_q(a+b)", Q == pred, (p, a, b))
                if z > p:
                    carries += 1
                    check("B4.carry_rule", Q == (fq(z - p, p) - invt[z - p]) % p, (p, a, b))
                # carry form: q(a+b) = phi(a+b mod p) - eps/(a+b), eps = [a+b >= p]
                eps = 1 if z >= p else 0
                zr = z % p
                check("B7a.carry_form", Q == (fq(zr, p) - eps * invt[zr]) % p, (p, a, b))
                # the carry bit is a polynomial function mod p: eps = u(a+b) - u(a) - u(b) - Gamma(a,b)
                check("B7b.carry_equals_u_coboundary_minus_Gamma",
                      eps % p == (u[zr] - u[a] - u[b] - zr * L1[s]) % p, (p, a, b))
                # Teichmuller version
                check("B5g.teichmuller_sum", fq(omega[a] + omega[b], p) == L1[s], (p, a, b))
                if p <= 61:
                    G = sum(((-1) ** (k - 1)) * invt[k] * pow(a, k, p) * pow(b, p - k, p) for k in range(1, p)) % p
                    check("B2.Gamma_equals_(a+b)L1", G == (a + b) * L1[s] % p, (p, a, b))
        stats[p] = {"pairs_with_carry": carries}
    RES["B_carry_stats_sample"] = {str(k): v for k, v in list(stats.items())[:5]}
    RES["B_time"] = round(time.time() - t, 2)


# ----------------------------------------------------------------------------------------
# Configurations (complete bicliques A + B subset of Q u {0})
# ----------------------------------------------------------------------------------------
def greedy_clique(p, seed):
    rng = random.Random(seed)
    chit = chi_table(p)
    A = [0]
    cand = [x for x in range(1, p) if chit[x] == 1]
    rng.shuffle(cand)
    for x in cand:
        if all(chit[(x - a) % p] == 1 for a in A):
            A.append(x)
    return sorted(A)


def greedy_balanced(p, seed):
    """Greedy: grow A keeping |B(A)| large; stop when |B(A)| <= |A| + 1."""
    rng = random.Random(seed)
    chit = chi_table(p)
    A = [0]
    while True:
        B = partner_set(A, p, chit)
        if len(B) <= len(A) + 1:
            break
        cands = [b for b in B if b not in A]
        rng.shuffle(cands)
        best, bestn = None, -1
        for x in cands[:60]:
            nb = len(partner_set(A + [x], p, chit))
            if nb > bestn:
                best, bestn = x, nb
        if best is None or bestn < len(A) + 1:
            break
        A.append(best)
    B = partner_set(A, p, chit)
    return sorted(A), B


def build_configs():
    cfg = []
    # tight / named examples
    for p in [13, 17, 29, 37, 41, 53, 61, 101]:
        if p % 4 == 1:
            A = [0, 1]
            cfg.append(("A={0,1}", p, A, partner_set(A, p)))
    cfg.append(("tight p=13 clique", 13, [0, 1, 4], [0, 9, 12]))
    cfg.append(("tight p=37", 37, [0, 1, 11], partner_set([0, 1, 11], 37)))
    cfg.append(("tight p=41 5-clique", 41, [0, 1, 2, 10, 33], [0, 40, 39, 31, 8]))
    for p in [19, 43, 103]:  # p = 3 mod 4
        A = [0, 1]
        cfg.append(("A={0,1}", p, A, partner_set(A, p)))
    # extremal witnesses from the balanced worker's data file (A maximising |B(A)|)
    path = os.path.join(ROOT, "results", "stepanov_balanced_2026_09_29_data.json")
    if os.path.exists(path):
        rows = json.load(open(path))["rows"]
        for p in [101, 197, 509, 997]:
            r = rows.get(str(p))
            if not r:
                continue
            for mstr, A in sorted(r["witness"].items(), key=lambda kv: int(kv[0])):
                if len(A) >= 3:
                    B = partner_set(A, p)
                    if len(B) >= 2:
                        cfg.append(("extremal |A|=%d" % len(A), p, list(A), B))
    # greedy Paley cliques (B = B(A) contains -A)
    for p in [101, 197, 509, 1009]:
        A = greedy_clique(p, p)
        cfg.append(("greedy clique", p, A, partner_set(A, p)))
    # greedy balanced pairs
    for p in [101, 211, 509, 1009]:
        A, B = greedy_balanced(p, 7 * p)
        if len(A) >= 2 and len(B) >= 2:
            cfg.append(("greedy balanced", p, A, B))
    # random small A, B = B(A)
    rng = random.Random(99)
    for p in [103, 211, 499, 1013]:
        for m in [3, 4, 5]:
            for _ in range(10):
                A = sorted(rng.sample(range(p), m))
                B = partner_set(A, p)
                if len(B) >= 2:
                    cfg.append(("random A, B=B(A)", p, A, B))
                    break
    return cfg


# ----------------------------------------------------------------------------------------
# Section C: the lifted HP polynomial
# ----------------------------------------------------------------------------------------
class LiftedHP:
    def __init__(self, A, p, K=4, nodes=None):
        self.p = p
        self.K = K
        self.MK = p ** K
        self.A = list(A)
        self.m = len(A)
        self.d = (p - 1) // 2
        self.D = self.d + self.m - 1
        self.nodes = list(nodes) if nodes is not None else list(A)  # integer lifts of the a_k
        MK = self.MK
        self.c = []
        for k, ak in enumerate(self.nodes):
            pr = 1
            for l, al in enumerate(self.nodes):
                if l != k:
                    pr = pr * (ak - al) % MK
            self.c.append(inv(pr, MK))
        self.cmodp = [c % p for c in self.c]

    def hasse(self, j, z, exponent=None):
        """F^[j](z) = C(N,j) sum_k c_k (z + a_k)^(N-j) - [j=0]  mod p^K, N = exponent (default D)."""
        N = self.D if exponent is None else exponent
        MK = self.MK
        s = 0
        for ck, ak in zip(self.c, self.nodes):
            s += ck * pow((z + ak) % MK, N - j, MK)
        val = comb(N, j) % MK * s - (1 if j == 0 else 0)
        return val % MK

    def coeff(self, i, exponent=None):
        N = self.D if exponent is None else exponent
        MK = self.MK
        s = 0
        for ck, ak in zip(self.c, self.nodes):
            s += ck * pow(ak % MK, N - i, MK)
        val = comb(N, i) % MK * s - (1 if i == 0 else 0)
        return val % MK


def section_C(configs):
    t = time.time()
    summary = []
    tau_examples = []
    for (kind, p, A, B) in configs:
        M = p * p
        m = len(A)
        d = (p - 1) // 2
        inv2 = (p + 1) // 2
        L1 = L1_table(p)
        u = [ucan(x, p) for x in range(p)]
        H = LiftedHP(A, p)
        D = H.D
        negA = set((-a) % p for a in A)
        B1 = [b for b in B if b not in negA]
        B0 = [b for b in B if b in negA]
        # completeness sanity
        for b in B:
            check("C0.complete", all(leg(a + b, p) >= 0 for a in A), (p, A, b))
        # C1: degree exactly d, leading coefficient C(D, m-1)
        for i in range(d + 1, D + 1):
            check("C1.coeff_above_d_vanish", H.coeff(i) == 0, (p, A, i))
        check("C1.leading_coeff", H.coeff(d) == comb(D, m - 1) % H.MK, (p, A))
        # C1b: Hasse-derivative evaluation agrees with the coefficient vector (small p)
        if p <= 101:
            coeffs = [H.coeff(i) for i in range(d + 1)]
            for z in random.sample(range(-p * p, p * p), 4):
                for j in range(min(m + 1, d + 1)):
                    direct = sum(comb(i, j) * coeffs[i] * pow(z % H.MK, i - j, H.MK) for i in range(j, d + 1)) % H.MK
                    check("C1b.hasse_vs_coeffs", direct == H.hasse(j, z), (p, A, z, j))
        # C2: general Taylor formula mod p^2 at every b not in -A (all b for p <= 509, else sample)
        blist = [b for b in range(p) if b not in negA]
        if p > 509:
            blist = random.sample(blist, 150) + B1
        for b in blist:
            ys = [b + a for a in A]
            qs = [fq(y, p) for y in ys]
            chs = [leg(y, p) for y in ys]
            for j in range(m):
                pred = comb(D, j) * sum(H.c[k] * chs[k] * (1 + p * (qs[k] * inv2 % p)) * pow(ys[k], m - 1 - j, M)
                                        for k in range(m)) - (1 if j == 0 else 0)
                check("C2.taylor_mod_p2_all_points", H.hasse(j, b) % M == pred % M, (p, A, b, j))
        # C3/C5/C6/C9 at complete points b in B1
        nzero_tau0 = 0
        nzero_tau01 = 0
        n_eis = 0
        n_full_order_lift = 0
        tau0_vals = []
        for b in B1:
            ys = [b + a for a in A]
            qs = [fq(y, p) for y in ys]
            hs = [H.hasse(j, b) for j in range(m + 1)]
            for j in range(m):
                pred = p * (comb(D, j) * inv2 % p * sum(H.cmodp[k] * qs[k] * pow(ys[k], m - 1 - j, p)
                                                       for k in range(m)) % p)
                check("C3.complete_point_formula", hs[j] % M == pred % M, (p, A, b, j))
            check("C3b.exact_order_m_mod_p", hs[m] % p != 0, (p, A, b))
            # lift independence for j <= m-2
            tt = random.randrange(1, p)
            for j in range(m - 1):
                check("C5.lift_independence", H.hasse(j, b + p * tt) % M == hs[j] % M, (p, A, b, j))
            # decomposition (canonical = Teichmuller part + u part), j <= m-2
            for j in range(m - 1):
                yk = [(b + a) % p for a in A]
                wit = sum(H.cmodp[k] * L1[A[k] * inv(yk[k], p) % p] * pow(yk[k], m - 1 - j, p) for k in range(m))
                upart = sum(H.cmodp[k] * u[A[k]] * pow(yk[k], m - 2 - j, p) for k in range(m))
                pred = comb(D, j) * inv2 * (wit + upart) % p
                check("C6.decomposition", (hs[j] // p) % p == pred, (p, A, b, j))
            tau0 = (hs[0] // p) % p
            tau0_vals.append(tau0)
            if tau0 == 0:
                nzero_tau0 += 1
                if m >= 3 and (hs[1] // p) % p == 0:
                    nzero_tau01 += 1
            # Newton polygon at b: v(hs[j]) >= 1 for j < m, v(hs[m]) = 0
            vs = [vp(h, p, H.K) for h in hs]
            check("C9.newton_polygon_disc_count", all(v >= 1 for v in vs[:m]) and vs[m] == 0, (p, A, b, vs))
            if vs[0] == 1:
                n_eis += 1
            # full order m mod p^2 at some lift  <=>  some lift b+pt with q(b+pt+a_k) = 0 for all k
            if p <= 211:
                lhs = False
                if all(hs[j] % M == 0 for j in range(m - 1)):
                    for t2 in range(p):
                        if H.hasse(m - 1, b + p * t2) % M == 0:
                            lhs = True
                            break
                rhs = any(all(fq(b + p * t2 + a, p) == 0 for a in A) for t2 in range(p))
                check("C8.full_order_lift_criterion", lhs == rhs, (p, A, b))
                if lhs:
                    n_full_order_lift += 1
        # B0 points with lift b = -a_{k0}
        for b in B0:
            k0 = [(-a) % p for a in A].index(b)
            bt = -A[k0]
            ys = [bt + a for a in A]
            for j in range(m - 1):
                pred = p * (comb(D, j) * inv2 % p * sum(H.cmodp[k] * fq(ys[k], p) * pow(ys[k], m - 1 - j, p)
                                                       for k in range(m) if k != k0) % p)
                check("C4.B0_point_formula", H.hasse(j, bt) % M == pred % M, (p, A, b, j))
            check("C4b.B0_exact_order", H.hasse(m - 1, bt) % p != 0, (p, A, b))
            tt = random.randrange(1, p)
            for j in range(m - 2):
                check("C5b.lift_independence_B0", H.hasse(j, bt + p * tt) % M == H.hasse(j, bt) % M, (p, A, b, j))
        # C7 Teichmuller-lifted HP polynomial
        if p <= 1013:
            omK = [pow(a, p ** (H.K - 1), H.MK) for a in A]
            HT = LiftedHP(A, p, nodes=omK)
            if p <= 101:
                for i in range(d + 1):
                    check("C7a.teich_same_reduction", (HT.coeff(i) - H.coeff(i)) % p == 0, (p, A, i))
            for b in B1[:40]:
                wb = pow(b, p ** (H.K - 1), H.MK)
                yk = [(b + a) % p for a in A]
                for j in range(m):
                    val = HT.hasse(j, wb)
                    wit = sum(HT.cmodp[k] * L1[A[k] * inv(yk[k], p) % p] * pow(yk[k], m - 1 - j, p) for k in range(m))
                    pred = p * (comb(D, j) * inv2 * wit % p)
                    check("C7b.teichmuller_second_digit", val % M == pred % M, (p, A, b, j))
        # Frobenius-exponent polynomial F_*: exponent N = p d + m - 1
        N = p * d + m - 1
        check("C10a.Fstar_degree_pd", comb(N, m - 1) % p == 1 % p, (p, m))
        for b in B1[:30]:
            tt = random.randrange(0, p)
            for j in range(m):
                check("C10b.Fstar_vanishes_mod_p2_order_m", H.hasse(j, b + p * tt, exponent=N) % M == 0, (p, A, b, j))
        n1 = len(B1)
        summary.append({"kind": kind, "p": p, "m": m, "n": len(B), "r": len(B0), "mn_minus_r_over_d":
                        round((m * len(B) - len(B0)) / d, 4), "complete_points_B1": n1,
                        "tau0_zero": nzero_tau0, "tau0_tau1_zero": nzero_tau01, "expected_tau0_zero": round(n1 / p, 3),
                        "eisenstein_clusters": n_eis, "full_order_m_mod_p2_some_lift": n_full_order_lift})
        if len(tau_examples) < 4 and n1 > 0:
            tau_examples.append({"kind": kind, "p": p, "A": A, "B1": B1[:20], "tau0": tau0_vals[:20]})
    RES["C_summary"] = summary
    RES["C_tau0_witnesses"] = tau_examples
    tot_pts = sum(s["complete_points_B1"] for s in summary)
    tot_z = sum(s["tau0_zero"] for s in summary)
    tot_exp = sum(s["expected_tau0_zero"] for s in summary)
    tot_e = sum(s["eisenstein_clusters"] for s in summary)
    RES["C_totals"] = {"configs": len(summary), "complete_points": tot_pts, "tau0_zero": tot_z,
                       "expected_if_uniform": round(tot_exp, 2), "eisenstein": tot_e,
                       "full_order_lift": sum(s["full_order_m_mod_p2_some_lift"] for s in summary)}
    RES["C_time"] = round(time.time() - t, 2)


def section_C_wieferich():
    """A = {0,1}: F(x) = -1 - x^D + (x+1)^D, so F(1) = 2^(d+1) - 2.  For p = +-1 mod 8 the point
    b = 1 is complete (1, 2 in Q) and F(1) = p q_p(2) mod p^2: extra mod-p^2 vanishing at b = 1
    happens exactly for Wieferich primes."""
    zeros = []
    for p in [q for q in PRIMES if 7 <= q <= 3600 and q % 8 in (1, 7)]:
        H = LiftedHP([0, 1], p, K=2)
        M = p * p
        v = H.hasse(0, 1) % M
        check("C11.A01_b1_is_p_q(2)", v == p * fq(2, p) % M, p)
        if v == 0:
            zeros.append(p)
    check("C11b.only_wieferich_3511_below_3600", zeros == [3511], zeros)
    RES["C_wieferich_zeros"] = zeros
    # Cauchy-Mirimanoff: (x^2+x+1)^2 divides (x+1)^D - x^D - 1 over Z when D = 1 mod 6
    import sympy as sp
    X = sp.symbols("X")
    for D in range(7, 200, 6):
        rem = sp.rem(sp.Poly((X + 1) ** D - X ** D - 1, X), sp.Poly((X ** 2 + X + 1) ** 2, X))
        check("C12a.cauchy_mirimanoff_divisibility", rem.is_zero, D)
    # hence for p = 1 mod 12 and A = {0,1}: b = -zeta (zeta primitive 6th root of 1) is complete and
    # F_A vanishes mod p^2 to order 2 there (exact double root -omega(zeta) over Z_p)
    rows = []
    for p in [q for q in PRIMES if q % 12 == 1 and q <= 600]:
        M = p * p
        H = LiftedHP([0, 1], p, K=2)
        zetas = [z for z in range(2, p) if (z * z - z + 1) % p == 0]
        for z in zetas:
            b = (-z) % p
            check("C12b.minus_zeta_complete", leg(b, p) == 1 and leg(b + 1, p) == 1, (p, z))
            w = pow(b, p, M)  # Teichmuller lift of b mod p^2
            check("C12c.double_root_mod_p2_at_teichmuller_lift",
                  H.hasse(0, w) % M == 0 and H.hasse(1, w) % M == 0, (p, z))
        rows.append({"p": p, "zetas": zetas})
    RES["C_cauchy_mirimanoff_sample"] = rows[:6]


# ----------------------------------------------------------------------------------------
# Section D: S(A,B) modulo p^2
# ----------------------------------------------------------------------------------------
def T_value(A, B, p, qtab, chit):
    Z = np.array(A, dtype=np.int64)[:, None] + np.array(B, dtype=np.int64)[None, :]
    ch = chit[Z % p]
    return int((ch * qtab[Z]).sum() % p)


def q_table(p):
    """Fermat quotients of the integers 0..2p-1 (0 where p | z)."""
    return np.array([fq(z, p) if z % p else 0 for z in range(2 * p)], dtype=np.int64)


def section_D(configs):
    t = time.time()
    rng = random.Random(4)
    cases = []
    for (kind, p, A, B) in configs:
        cases.append((kind, p, A, B))
    for p in [5, 7, 11, 13, 29, 101, 199, 503]:
        for _ in range(8):
            m = rng.randrange(1, min(p, 12))
            n = rng.randrange(1, min(p, 12))
            cases.append(("random rectangle", p, sorted(rng.sample(range(p), m)), sorted(rng.sample(range(p), n))))
    for (kind, p, A, B) in cases:
        if p < 5:
            continue
        M = p * p
        d = (p - 1) // 2
        inv2 = (p + 1) // 2
        L1 = L1_table(p)
        u = [ucan(x, p) for x in range(p)]
        chit = chi_table(p)
        qtab = q_table(p)
        Sint = sum(pow(a + b, d, M) for a in A for b in B) % M
        S = sum(leg(a + b, p) for a in A for b in B)
        T = sum(leg(a + b, p) * fq(a + b, p) for a in A for b in B if (a + b) % p) % p
        check("D1.S_mod_p2_identity", Sint == (S + p * (inv2 * T % p)) % M, (p, A, B))
        check("D1b.T_numpy_agrees", T == T_value(A, B, p, qtab, chit), (p,))
        m, n = len(A), len(B)
        r = sum(1 for b in B if (-b) % p in set(A))
        if m * n < p:
            r0 = Sint % p
            cand = [s for s in range(-m * n, m * n + 1) if (s - r0) % p == 0 and (s - (m * n - r)) % 2 == 0]
            check("D2.S_recovered_from_residue", cand == [S], (p, m, n, cand, S))
        TW = sum(leg(a + b, p) * L1[a * inv(a + b, p) % p] for a in A for b in B if (a + b) % p) % p
        TU = sum(leg(a + b, p) * (u[a] + u[b]) * inv(a + b, p) for a in A for b in B if (a + b) % p) % p
        check("D3.T_decomposition", T == (TW + TU) % p, (p,))
        if p <= 1013 and m * n * p <= 3_000_000:
            tot = 0
            for tt in range(p):
                A2 = [(a + tt) % p for a in A]
                B2 = [(b - tt) % p for b in B]
                tot += T_value(A2, B2, p, qtab, chit)
            pred = (S + sum(leg(a + b, p) * inv(a + b, p) for a in A for b in B if (a + b) % p)) % p
            check("D4.translation_average", tot % p == pred, (p, A, B))
    RES["D_cases"] = len(cases)
    RES["D_time"] = round(time.time() - t, 2)


# ----------------------------------------------------------------------------------------
# Section E: orbit test
# ----------------------------------------------------------------------------------------
def section_E(configs):
    """Orbit test.  The maps (A,B) -> (sA+t, sB-t), s in Q, t in F_p, preserve A+B in Q u {0}
    and every count-level (profile) quantity.  We record how the mod-p^2 data move, and check
    the structural decomposition T = T_phi - T_carry, where T_phi depends only on s and the
    Witt part T_W only on t/s."""
    t0 = time.time()
    out = []
    chosen = []
    seen = set()
    for (kind, p, A, B) in configs:
        key = (p, tuple(A))
        if key in seen:
            continue
        seen.add(key)
        if p in (13, 37, 41, 101, 197) or (p in (509, 997, 1009) and kind in ("greedy clique", "greedy balanced", "extremal |A|=3", "extremal |A|=8")):
            chosen.append((kind, p, A, B))
    for (kind, p, A, B) in chosen:
        chit = chi_table(p)
        qtab = q_table(p)
        L1 = np.array(L1_table(p), dtype=np.int64)
        invtab = np.array([0] + [inv(z, p) for z in range(1, p)], dtype=np.int64)
        m = len(A)
        negA = set((-a) % p for a in A)
        B1 = [b for b in B if b not in negA]
        squares = sorted(set((x * x) % p for x in range(1, p)))
        full = p <= 200
        svals = squares if full else squares[:12]
        Tcount = np.zeros(p, dtype=np.int64)
        byTW = {}
        byPhi = {}
        byCarry = {}
        cA = []
        for k in range(m):
            pr = 1
            for l in range(m):
                if l != k:
                    pr = pr * (A[k] - A[l]) % p
            cA.append(inv(pr, p))
        inv2 = (p + 1) // 2
        Aa = np.array(A, dtype=np.int64)
        Bb = np.array(B1, dtype=np.int64)
        Ball = np.array(B, dtype=np.int64)
        okTW = okPhi = okCarry = okSplit = True
        wit = {}
        for s in svals:
            sinv = inv(s, p)
            sinv_pow = pow(sinv, m - 1, p)
            c = np.array([ck * sinv_pow % p for ck in cA], dtype=np.int64)
            for tt in range(p):
                A2 = (s * Aa + tt) % p
                B2all = (s * Ball - tt) % p
                Z = A2[:, None] + B2all[None, :]
                Zr = Z % p
                ch = chit[Zr]
                T = int((ch * qtab[Z]).sum() % p)
                Tcount[T] += 1
                # Witt part: sum chi(z) L1(a'/z), z = a'+b' (mod p), z != 0
                nz = Zr != 0
                sarg = (A2[:, None] * invtab[Zr]) % p
                TW = int((ch * L1[sarg] * nz).sum() % p)
                tprime = tt * sinv % p
                if tprime in byTW and byTW[tprime] != TW:
                    okTW = False
                byTW[tprime] = TW
                # phi part: sum chi(z) phi(z mod p); carry part: sum chi(z) [Z >= p] / z
                Tphi = int((ch * qtab[Zr]).sum() % p)
                Tcar = int((ch * (Z >= p) * invtab[Zr]).sum() % p)
                if (Tphi - Tcar - T) % p:
                    okSplit = False
                if s in byPhi and byPhi[s] != Tphi:
                    okPhi = False
                byPhi[s] = Tphi
                if len(B1):
                    B2 = (s * Bb - tt) % p
                    Z1 = A2[:, None] + B2[None, :]
                    Q = qtab[Z1]
                    Y = Z1 % p
                    Yp = np.ones_like(Y)
                    for _ in range(m - 1):
                        Yp = (Yp * Y) % p
                    tau = ((c[:, None] * Q % p) * Yp % p).sum(axis=0) % p
                    tau = tau * inv2 % p
                    carry = (Z1 >= p)
                    nz0 = int((tau == 0).sum())
                    if nz0 and "zero" not in wit:
                        wit["zero"] = {"s": int(s), "t": int(tt), "A'": [int(v) for v in A2],
                                       "points_with_tau0=0": [int(B2[j]) for j in range(len(B1)) if tau[j] == 0]}
                    if nz0 == 0 and "nonzero" not in wit:
                        wit["nonzero"] = {"s": int(s), "t": int(tt), "A'": [int(v) for v in A2]}
                    for j, b in enumerate(B1):
                        key = (s, b, tuple(carry[:, j].tolist()))
                        v = int(tau[j])
                        if key in byCarry and byCarry[key] != v:
                            okCarry = False
                        byCarry[key] = v
        check("E2.T_split_phi_minus_carry", okSplit, (p, A))
        check("E3.Witt_part_depends_only_on_t/s", okTW, (p, A))
        check("E4.phi_part_depends_only_on_s", okPhi, (p, A))
        check("E5.tau0_depends_only_on_s_and_carries", okCarry, (p, A))
        norb = int(Tcount.sum())
        rec = {"kind": kind, "p": p, "m": m, "n": len(B), "orbit_elements": norb, "full_orbit": full,
               "distinct_T_values": int((Tcount > 0).sum()),
               "max_multiplicity_of_a_T_value": int(Tcount.max()),
               "distinct_TW_values": len(set(byTW.values())), "distinct_Tphi_values": len(set(byPhi.values())),
               "distinct_(s,b,carry)_classes": len(byCarry),
               "witness_tau0_zero_vs_nonzero_same_profile": wit if len(wit) == 2 else None}
        out.append(rec)
        check("E1.T_not_orbit_invariant", rec["distinct_T_values"] > 1, rec)
    RES["E_orbits"] = out
    RES["E_time"] = round(time.time() - t0, 2)


# ----------------------------------------------------------------------------------------
# Section F: randomness diagnostics (HEURISTIC)
# ----------------------------------------------------------------------------------------
def rank_mod_p(Mat, p):
    Mt = [list(map(lambda v: int(v) % p, row)) for row in Mat]
    rows, cols = len(Mt), len(Mt[0]) if Mt else 0
    r = 0
    for c in range(cols):
        piv = None
        for i in range(r, rows):
            if Mt[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        Mt[r], Mt[piv] = Mt[piv], Mt[r]
        iv = inv(Mt[r][c], p)
        Mt[r] = [v * iv % p for v in Mt[r]]
        for i in range(rows):
            if i != r and Mt[i][c]:
                f = Mt[i][c]
                Mt[i] = [(vi - f * vr) % p for vi, vr in zip(Mt[i], Mt[r])]
        r += 1
    return r


def expsum_max(vals, p, weights=None):
    """max over h != 0 of |sum_i w_i e(h v_i / p)|."""
    v = np.array(vals, dtype=np.float64)
    w = np.ones(len(vals)) if weights is None else np.array(weights, dtype=np.float64)
    best = 0.0
    for h in range(1, p):
        z = np.exp(2j * np.pi * h * v / p)
        best = max(best, abs((w * z).sum()))
    return best


def section_F(configs):
    """HEURISTIC diagnostics.  Note that q(a~+b~) is a function of the integer sum a~+b~, so all
    statistics are taken over DISTINCT integer sums (repeated sums repeat the value)."""
    t0 = time.time()
    rng = random.Random(11)
    rows = []
    for (kind, p, A, B) in configs:
        if p < 29:
            continue
        negA = set((-a) % p for a in A)
        B1 = [b for b in B if b not in negA]
        if len(B1) < 2:
            continue
        sums = sorted(set(a + b for a in A for b in B1))
        qv = [fq(z, p) for z in sums]
        Qm = [[fq(a + b, p) for b in B1] for a in A]
        k = min(len(A), len(B1))
        rkQ = rank_mod_p(Qm, p)
        # random rectangle of the same shape (avoiding sums = 0 mod p)
        Ar = sorted(rng.sample(range(p), len(A)))
        Br = [b for b in rng.sample(range(p), p) if all((a + b) % p for a in Ar)][:len(B1)]
        sums_r = sorted(set(a + b for a in Ar for b in Br))
        qr = [fq(z, p) for z in sums_r]
        Qr = [[fq(a + b, p) for b in Br] for a in Ar]
        pairs = [fq(a + b, p) for a in A for b in B1]
        pairs_r = [fq(a + b, p) for a in Ar for b in Br]
        # normalise by the L2 norm of the representation function r(z) of the integer sums
        from collections import Counter
        l2 = math.sqrt(sum(v * v for v in Counter(a + b for a in A for b in B1).values()))
        l2r = math.sqrt(sum(v * v for v in Counter(a + b for a in Ar for b in Br).values()))
        rows.append({"kind": kind, "p": p, "m": len(A), "n1": len(B1),
                     "distinct_integer_sums": len(sums), "zeros_among_distinct_sums": qv.count(0),
                     "expected_zeros": round(len(sums) / p, 3),
                     "rank_Q": rkQ, "min(m,n1)": k, "rank_Q_random_rectangle": rank_mod_p(Qr, p),
                     "max_h|sum_pairs e(hQ/p)|/||r||_2": round(expsum_max(pairs, p) / l2, 3),
                     "same_for_random_rectangle": round(expsum_max(pairs_r, p) / l2r, 3),
                     "zeros_random_distinct_sums": qr.count(0)})
        check("F1.rank_Q_full", rkQ == k or p < 50, (p, A, rkQ, k))
    RES["F_matrix_stats"] = rows
    # exponential sums of the canonical Fermat quotient phi on [1, p-1], with and without chi
    es = {}
    for p in [q for q in PRIMES if 29 <= q <= 400][::3]:
        phi = [fq(r, p) for r in range(1, p)]
        chs = [leg(r, p) for r in range(1, p)]
        es[p] = {"max_h|sum e(h phi/p)|/sqrt(p)": round(expsum_max(phi, p) / math.sqrt(p), 3),
                 "max_h|sum chi e(h phi/p)|/sqrt(p)": round(expsum_max(phi, p, chs) / math.sqrt(p), 3)}
    RES["F_phi_exponential_sums"] = es
    # degree of the Teichmuller second-digit function sigma_A(x)
    degs = []
    for (kind, p, A, B) in configs:
        if p > 211 or len(degs) >= 25:
            continue
        L1 = L1_table(p)
        m = len(A)
        cA = []
        for k in range(m):
            pr = 1
            for l in range(m):
                if l != k:
                    pr = pr * (A[k] - A[l]) % p
            cA.append(inv(pr, p))
        vals = []
        for x in range(p):
            s = 0
            for k, a in enumerate(A):
                y = (x + a) % p
                if y == 0:
                    continue
                s += cA[k] * leg(y, p) * L1[a * inv(y, p) % p] * pow(y, m - 1, p)
            vals.append(s % p)
        dg = poly_degree_of_function(vals, p)
        degs.append({"kind": kind, "p": p, "m": m, "deg_sigma": dg, "p-1": p - 1, "d": (p - 1) // 2})
        check("F4.second_digit_function_degree_exceeds_d", dg > (p - 1) // 2, (p, A, dg))
    RES["F_second_digit_degrees"] = degs
    RES["F_time"] = round(time.time() - t0, 2)


# ----------------------------------------------------------------------------------------
# Section G: height illustration
# ----------------------------------------------------------------------------------------
def section_G():
    t = time.time()
    try:
        import sympy as sp
    except Exception as ex:  # pragma: no cover
        RES["G"] = "sympy unavailable: %s" % ex
        return
    from fractions import Fraction
    x = sp.symbols("x")
    out = []
    for (p, A) in [(13, [0, 1, 4]), (17, [0, 1]), (29, [0, 1]), (37, [0, 1, 11]), (29, [0, 1, 5])]:
        m = len(A)
        d = (p - 1) // 2
        D = d + m - 1
        cs = []
        for k in range(m):
            pr = 1
            for l in range(m):
                if l != k:
                    pr *= (A[k] - A[l])
            cs.append(Fraction(1, pr))
        L = 1
        for c in cs:
            L = L * c.denominator // math.gcd(L, c.denominator)
        poly = -L * sp.Integer(1) + sum(sp.Integer(int(c * L)) * (x + a) ** D for c, a in zip(cs, A))
        P = sp.Poly(sp.expand(poly), x)
        check("G0.degree_d_over_Q", P.degree() == d, (p, A))
        disc = int(sp.discriminant(P))
        B = partner_set(A, p)
        negA = set((-a) % p for a in A)
        n1 = sum(1 for b in B if b not in negA)
        vdisc = vp(abs(disc), p, 10 ** 9) if disc else None
        out.append({"p": p, "A": A, "d": d, "n": len(B), "r": len(B) - n1,
                    "disc_nonzero": disc != 0, "v_p(disc)": vdisc,
                    "forced_by_clusters_approx n1(m-1)": n1 * (m - 1),
                    "log_p|disc|": round(math.log(abs(disc)) / math.log(p), 1) if disc else None})
        check("G1.disc_nonzero", disc != 0, (p, A))
    RES["G_height"] = out
    RES["G_time"] = round(time.time() - t, 2)


# ----------------------------------------------------------------------------------------
# Section H: the Witt part is field-agnostic (Galois ring W_2(F_{p^2}) = GR(p^2, 2))
# ----------------------------------------------------------------------------------------
def section_H():
    t0 = time.time()
    for p in [3, 5, 7, 11]:
        M = p * p
        q = p * p
        # monic quadratic x^2 + c1 x + c0 irreducible mod p
        c0 = c1 = None
        for u1 in range(p):
            for u0 in range(1, p):
                if all((x * x + u1 * x + u0) % p for x in range(p)):
                    c1, c0 = u1, u0
                    break
            if c1 is not None:
                break

        def mul(x, y, mod):
            a0 = x[0] * y[0]
            a1 = x[0] * y[1] + x[1] * y[0]
            a2 = x[1] * y[1]
            return ((a0 - a2 * c0) % mod, (a1 - a2 * c1) % mod)

        def pw(x, e, mod):
            r = (1, 0)
            b = (x[0] % mod, x[1] % mod)
            while e:
                if e & 1:
                    r = mul(r, b, mod)
                b = mul(b, b, mod)
                e >>= 1
            return r

        def add(x, y, mod):
            return ((x[0] + y[0]) % mod, (x[1] + y[1]) % mod)

        def fqinv(x):
            return pw(x, q - 2, p)

        elems = [(a0, a1) for a0 in range(p) for a1 in range(p)]
        nz = [e for e in elems if e != (0, 0)]
        half = (p + 1) // 2  # 1/2 mod p
        omega = {e: pw(e, q, M) for e in elems}
        dq = (q - 1) // 2
        for z in nz:
            for _ in range(3):
                zl = (z[0] + p * random.randrange(p), z[1] + p * random.randrange(p))
                r = pw(zl, q - 1, M)
                assert r[0] % p == 1 and r[1] % p == 0
                qq = (((r[0] - 1) // p) % p, (r[1] // p) % p)   # Fermat quotient in F_q
                chi = pw(z, dq, p)  # +-1 in F_q
                chi_int = 1 if chi == (1, 0) else -1
                assert chi in [(1, 0), (p - 1, 0)]
                pred = ((chi_int * (1 + p * (qq[0] * half % p))) % M, (chi_int * p * (qq[1] * half % p)) % M)
                check("H1.refined_euler_F_{p^2}", pw(zl, dq, M) == pred, (p, zl))
            check("H2.teichmuller_is_q-th_power", pw(omega[z], q - 1, M) == (1, 0), (p, z))
        invk = [0] + [inv(k, p) for k in range(1, p)]
        for a in elems:
            for b in elems:
                s = add(a, b, p)
                if s == (0, 0):
                    continue
                c = add(omega[a], omega[b], M)
                r = pw(c, q - 1, M)
                qq = (((r[0] - 1) // p) % p, (r[1] // p) % p)
                # Gamma(a,b) = sum_{k=1}^{p-1} (-1)^(k-1)/k a^k b^(p-k)  in F_q
                G = (0, 0)
                for k in range(1, p):
                    term = mul(pw(a, k, p), pw(b, p - k, p), p)
                    coef = ((-1) ** (k - 1)) * invk[k] % p
                    G = add(G, (term[0] * coef % p, term[1] * coef % p), p)
                Gtw = pw(G, p, p)  # inverse Frobenius on F_{p^2} is x -> x^p
                pred = mul(Gtw, fqinv(s), p)
                check("H3.witt_second_digit_F_{p^2}", qq == pred, (p, a, b))
    RES["H_time"] = round(time.time() - t0, 2)


def main():
    section_A()
    print("A done", round(time.time() - T0, 1), flush=True)
    section_B()
    print("B done", round(time.time() - T0, 1), flush=True)
    configs = build_configs()
    RES["configs"] = [{"kind": k, "p": p, "A": A, "n": len(B)} for (k, p, A, B) in configs]
    section_C(configs)
    section_C_wieferich()
    print("C done", round(time.time() - T0, 1), flush=True)
    section_D(configs)
    print("D done", round(time.time() - T0, 1), flush=True)
    section_E(configs)
    print("E done", round(time.time() - T0, 1), flush=True)
    section_F(configs)
    print("F done", round(time.time() - T0, 1), flush=True)
    section_G()
    print("G done", round(time.time() - T0, 1), flush=True)
    section_H()
    print("H done", round(time.time() - T0, 1), flush=True)
    out = {"description": __doc__, "checks": CHECKS, "n_checks": sum(CHECKS.values()),
           "failures": FAILS, "n_failures": len(FAILS), "results": RES,
           "elapsed_seconds": round(time.time() - T0, 1)}
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(out, f, indent=1, default=str)
    print("checks:", sum(CHECKS.values()), "failures:", len(FAILS), "elapsed:", round(time.time() - T0, 1))
    for k in sorted(CHECKS):
        print("  %-45s %d" % (k, CHECKS[k]))
    if FAILS:
        print("FAILURES (first 10):")
        for f in FAILS[:10]:
            print("  ", f)


if __name__ == "__main__":
    main()
