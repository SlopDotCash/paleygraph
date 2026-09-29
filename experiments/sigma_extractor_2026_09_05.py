#!/usr/bin/env python3
"""Sigma pass, direction `extractor` (2026-09-05): exact verifier.

Every pass/fail decision below is made by exact integer or rational
arithmetic (Python ints, fractions.Fraction) or by exact arithmetic in the
cyclotomic ring Z[zeta_p], represented as integer coefficient vectors of
length p modulo the relation 1 + zeta + ... + zeta^(p-1) = 0.  Floating point
is used only to choose set sizes (p^0.45) and to print logarithms in reports.

Sections (numbers match research/sigma-extractor-2026-09-05.md):
  S1  Gram identity M M^T = pI - J and the weighted Chung bound          (PROVED)
  S2  bias <-> statistical distance for the one-bit Paley extractor       (PROVED)
  S3  greedy flat-source decomposition of min-entropy sources (exact),
      bilinearity of the bias, data processing                            (PROVED)
  S4  level-set counting behind the strong-extractor statement            (PROVED)
  S5  Gauss-sum identity in Z[zeta_p], S^2 <= p E+(A,B), never below Chung (PROVED)
  S6  the squaring trick: Hadamard grows (Rao, Lemma 3.3), Paley is
      entropy-neutral (exact identity)                                    (PROVED)
  S7  parabola encoding: E+(Enc A) = 2|A|^2 - |A| exactly; 3-fold numbers  (illustration)
  S8  BIW-type composition chi(ab+c): identities, subgroup collapse, and the
      witness refuting `the dilation average controls S(A,B)'             (PROVED / REFUTED)
  S9  additive-condenser interval witnesses via least non-residues        (finite witnesses)

Writes results/sigma_extractor_2026_09_05.json.  Exit code 1 if any check fails.
"""
import json
import math
import random
import time
import traceback
from fractions import Fraction
from math import isqrt
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "results" / "sigma_extractor_2026_09_05.json"
RNG = random.Random(20260905)

CHECKS = {}
WITNESSES = []
NOTES = {}
FAILURES = []


def bump(sec, n=1):
    CHECKS[sec] = CHECKS.get(sec, 0) + n


# ---------------------------------------------------------------- utilities
def is_prime(n):
    return n >= 2 and all(n % d for d in range(2, isqrt(n) + 1))


def primes_upto(n):
    sieve = bytearray([1]) * (n + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, isqrt(n) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(range(i * i, n + 1, i)))
    return [i for i in range(n + 1) if sieve[i]]


def prime_factors(n):
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


class Field:
    def __init__(self, p):
        assert is_prime(p) and p > 2
        self.p = p
        sq = {x * x % p for x in range(1, p)}
        self.chi = [0] + [1 if x in sq else -1 for x in range(1, p)]
        self.inv = [0] + [pow(x, p - 2, p) for x in range(1, p)]
        self.chi_np = np.array(self.chi, dtype=np.int64)
        self._g = None

    def S(self, A, B):
        chi, p = self.chi, self.p
        return sum(chi[(a + b) % p] for a in A for b in B)

    def S_np(self, A, B):
        A = np.asarray(A, dtype=np.int64)
        B = np.asarray(B, dtype=np.int64)
        return int(self.chi_np[(A[:, None] + B[None, :]) % self.p].sum())

    def f(self, A):
        chi, p = self.chi, self.p
        return [sum(chi[(a + b) % p] for a in A) for b in range(p)]

    def generator(self):
        if self._g is None:
            p = self.p
            fac = prime_factors(p - 1)
            for g in range(2, p):
                if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
                    self._g = g
                    break
        return self._g

    def subgroup(self, d):
        p = self.p
        assert (p - 1) % d == 0
        g = pow(self.generator(), (p - 1) // d, p)
        H, h = [], 1
        for _ in range(d):
            H.append(h)
            h = h * g % p
        return H

    def gp(self, n):
        g, p = self.generator(), self.p
        out, h = [], 1
        for _ in range(n):
            out.append(h)
            h = h * g % p
        return out


def greedy_sidon(p, n):
    out, diffs = [], set()
    for x in range(1, p):
        nd = {(x - y) % p for y in out} | {(y - x) % p for y in out}
        if x not in out and not (nd & diffs) and len(nd) == 2 * len(out):
            out.append(x)
            diffs |= nd
            if len(out) == n:
                break
    return out


def families(F, n):
    p = F.p
    fam = {
        "interval": list(range(1, n + 1)),
        "random": RNG.sample(range(p), n),
        "gp": F.gp(n),
        "sidon": greedy_sidon(p, n),
    }
    if (p - 1) % n == 0:
        fam["subgroup"] = F.subgroup(n)
    return fam


def energy_sum(p, A, B):
    """E+(A,B) = #{(a,b,a',b') : a+b = a'+b'} = sum_s r_{A+B}(s)^2."""
    r = {}
    for a in A:
        for b in B:
            s = (a + b) % p
            r[s] = r.get(s, 0) + 1
    return sum(c * c for c in r.values())


def energy_mult(p, A, B):
    r = {}
    for a in A:
        for b in B:
            s = (a * b) % p
            r[s] = r.get(s, 0) + 1
    return sum(c * c for c in r.values()), r


# ------------------------------------------------ exact arithmetic in Z[zeta_p]
def cmul(u, v, p):
    w = [0] * p
    for i, ui in enumerate(u):
        if ui:
            for j, vj in enumerate(v):
                if vj:
                    w[(i + j) % p] += ui * vj
    return w


def cconj(u, p):
    return [u[(-i) % p] for i in range(p)]


def cadd(u, v):
    return [a + b for a, b in zip(u, v)]


def cscale(u, c):
    return [c * a for a in u]


def ceq(u, v):
    d0 = u[0] - v[0]
    return all(a - b == d0 for a, b in zip(u, v))


def cint(u):
    """The rational integer represented by u, if u is one."""
    assert all(u[i] == u[1] for i in range(1, len(u))), "not a rational integer"
    return u[0] - u[1]


# ------------------------------------------------------------------ S1
def section1():
    for p in [5, 7, 11, 13, 17, 19, 23, 29, 31]:
        F = Field(p)
        chi = F.chi
        for x in range(p):
            for x2 in range(p):
                g = sum(chi[(x + y) % p] * chi[(x2 + y) % p] for y in range(p))
                assert g == (p - 1 if x == x2 else -1), (p, x, x2, g)
                bump("S1_gram_identity")
        for trial in range(24):
            alpha = [RNG.randint(-3, 3) for _ in range(p)]
            beta = [RNG.randint(-3, 3) for _ in range(p)]
            if trial % 2:
                alpha = [abs(a) for a in alpha]
                beta = [abs(b) for b in beta]
            lhs = sum(alpha[x] * beta[y] * chi[(x + y) % p]
                      for x in range(p) for y in range(p))
            a2, b2 = sum(a * a for a in alpha), sum(b * b for b in beta)
            a1, b1 = sum(alpha), sum(beta)
            assert lhs * lhs <= a2 * (p * b2 - b1 * b1), (p, alpha, beta)
            assert lhs * lhs <= b2 * (p * a2 - a1 * a1), (p, alpha, beta)
            bump("S1_weighted_chung", 2)
        for trial in range(24):
            m, n = RNG.randint(1, p - 1), RNG.randint(1, p - 1)
            A, B = RNG.sample(range(p), m), RNG.sample(range(p), n)
            s = F.S(A, B)
            assert s * s <= m * n * (p - n) and s * s <= m * n * (p - m)
            bump("S1_flat_chung")


# ------------------------------------------------------------------ S2
def section2():
    for p in [7, 11, 13, 17, 19, 23]:
        F = Field(p)
        chi = F.chi
        for trial in range(30):
            m, n = RNG.randint(1, p), RNG.randint(1, p)
            A, B = RNG.sample(range(p), m), RNG.sample(range(p), n)
            ext = lambda a, b: 1 if ((a + b) % p == 0 or chi[(a + b) % p] == 1) else 0
            ones = sum(ext(a, b) for a in A for b in B)
            sd = abs(Fraction(ones, m * n) - Fraction(1, 2))
            s = F.S(A, B)
            z = sum(1 for a in A for b in B if (a + b) % p == 0)
            assert sd == Fraction(abs(s + z), 2 * m * n)
            assert z <= min(m, n)
            e = sum((-1) ** ext(a, b) for a in A for b in B)
            assert e == -(s + z)
            # Satake, Definition 13: Ext(x,y)=1 if x=y, else (chi(x-y)+1)/2
            ext2 = lambda a, b: 1 if (a == b or chi[(a - b) % p] == 1) else 0
            ones2 = sum(ext2(a, b) for a in A for b in B)
            sm = sum(chi[(a - b) % p] for a in A for b in B)
            zz = len(set(A) & set(B))
            assert abs(Fraction(ones2, m * n) - Fraction(1, 2)) == Fraction(abs(sm + zz), 2 * m * n)
            assert sm == F.S(A, [(-b) % p for b in B])
            bump("S2_bias_vs_sd", 5)


# ------------------------------------------------------------------ S3
def random_source(p, K):
    """A random distribution on F_p with max probability <= 1/K, with some
    atoms exactly at the maximum, given as {x: Fraction}."""
    m = RNG.randint(K, p)
    supp = RNG.sample(range(p), m)
    w = {x: RNG.randint(1, 10) for x in supp}
    for x in RNG.sample(supp, RNG.randint(0, K - 1)):
        w[x] = 10
    while sum(w.values()) < 10 * K:
        cand = [x for x in w if w[x] < 10]
        if not cand:
            x = RNG.choice([y for y in range(p) if y not in w])
            w[x] = 1
        else:
            w[RNG.choice(cand)] += 1
    W = sum(w.values())
    P = {x: Fraction(v, W) for x, v in w.items()}
    assert max(P.values()) <= Fraction(1, K)
    return P


def greedy_flat_decomposition(P, K):
    """P: {x: Fraction}, sum 1, max <= 1/K.  Returns [(lambda_i, S_i)] with
    |S_i| = K, lambda_i > 0, sum lambda_i = 1, sum lambda_i U_{S_i} = P."""
    P = {x: v for x, v in P.items() if v > 0}
    n0 = len(P)
    parts, remaining, steps = [], Fraction(1), 0
    while remaining > 0:
        steps += 1
        assert steps <= n0, "greedy decomposition did not terminate in n steps"
        order = sorted(P, key=lambda x: (-P[x], x))
        assert len(order) >= K
        S = order[:K]
        pK = P[order[K - 1]]
        pnext = P[order[K]] if len(order) > K else Fraction(0)
        t = min(K * pK, remaining - K * pnext)
        assert t > 0
        parts.append((t, tuple(sorted(S))))
        for x in S:
            P[x] -= t / K
        remaining -= t
        P = {x: v for x, v in P.items() if v > 0}
        assert all(v <= remaining / K for v in P.values())
        assert all(v >= 0 for v in P.values())
    return parts


def section3():
    for p in [11, 13, 17, 19, 23, 29, 31, 37]:
        F = Field(p)
        chi = F.chi
        for K in [2, 3, 4, 5]:
            for trial in range(5):
                PX, PY = random_source(p, K), random_source(p, K)
                DX, DY = greedy_flat_decomposition(dict(PX), K), greedy_flat_decomposition(dict(PY), K)
                for D, P in ((DX, PX), (DY, PY)):
                    assert sum(l for l, _ in D) == 1
                    recon = {x: Fraction(0) for x in range(p)}
                    for l, Sset in D:
                        assert len(Sset) == K
                        for x in Sset:
                            recon[x] += l / K
                    assert all(recon[x] == P.get(x, 0) for x in range(p))
                    bump("S3_flat_decomposition")
                bias = sum(PX[x] * PY[y] * chi[(x + y) % p] for x in PX for y in PY)
                terms, tot = [], Fraction(0)
                for lx, Sx in DX:
                    for ly, Sy in DY:
                        b = Fraction(F.S(Sx, Sy), K * K)
                        terms.append(b)
                        tot += lx * ly * b
                assert tot == bias
                assert abs(bias) <= max(abs(b) for b in terms)
                bump("S3_bilinearity", 2)
        for trial in range(20):
            P = random_source(p, 3)
            fmap = [RNG.randrange(p) for _ in range(p)]
            Q = {}
            for x, v in P.items():
                Q[fmap[x]] = Q.get(fmap[x], 0) + v
            assert max(Q.values()) >= max(P.values())
            assert sum(v * v for v in Q.values()) >= sum(v * v for v in P.values())
            bump("S3_data_processing", 2)


# ------------------------------------------------------------------ S4
def section4():
    for p in [31, 37, 41, 43, 47, 53]:
        F = Field(p)
        for trial in range(10):
            m = RNG.randint(2, p // 2)
            A = RNG.sample(range(p), m)
            fA = F.f(A)
            for t in range(0, m + 1):
                Bp = [b for b in range(p) if fA[b] > t]
                Bm = [b for b in range(p) if fA[b] < -t]
                Sp, Sm = F.S(A, Bp), F.S(A, Bm)
                assert Sp == sum(fA[b] for b in Bp) and Sm == sum(fA[b] for b in Bm)
                if Bp:
                    assert Sp > t * len(Bp)
                if Bm:
                    assert Sm < -t * len(Bm)
                cnt = len(Bp) + len(Bm)
                if cnt:
                    assert Sp - Sm > t * cnt
                bump("S4_levelset")


# ------------------------------------------------------------------ S5
def section5():
    min_ratio_energy_vs_chung = None
    max_ratio_bound_vs_S = None
    ledger = []
    # exact cyclotomic identities at small p
    for p in [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61]:
        F = Field(p)
        chi = F.chi
        G = list(chi)                       # G = sum_t chi(t) zeta^t
        assert cint(cmul(G, cconj(G, p), p)) == p
        assert cint(cmul(G, G, p)) == chi[p - 1] * p
        bump("S5_gauss_sum", 2)
        for trial in range(4):
            m, n = RNG.randint(1, p - 1), RNG.randint(1, p - 1)
            A, B = RNG.sample(range(p), m), RNG.sample(range(p), n)
            hatA = [[0] * p for _ in range(p)]
            hatB = [[0] * p for _ in range(p)]
            for t in range(p):
                for a in A:
                    hatA[t][(t * a) % p] += 1
                for b in B:
                    hatB[t][(t * b) % p] += 1
            rhs = [0] * p
            for t in range(1, p):
                rhs = cadd(rhs, cscale(cmul(hatA[t], hatB[t], p), chi[t]))
            lhs = cscale(G, F.S(A, B))
            assert ceq(lhs, rhs), (p, A, B)
            bump("S5_gauss_fourier_identity")
            # Parseval quartic: sum_t |hatA(t)|^4 = p E+(A,A)
            quart = [0] * p
            for t in range(p):
                sq = cmul(hatA[t], cconj(hatA[t], p), p)
                quart = cadd(quart, cmul(sq, sq, p))
            assert cint(quart) == p * energy_sum(p, A, A)
            # mixed: sum_t |hatA(t)|^2 |hatB(t)|^2 = p E+(A,B)
            mixed = [0] * p
            for t in range(p):
                mixed = cadd(mixed, cmul(cmul(hatA[t], cconj(hatA[t], p), p),
                                         cmul(hatB[t], cconj(hatB[t], p), p), p))
            assert cint(mixed) == p * energy_sum(p, A, B)
            bump("S5_parseval", 2)
    # the inequalities, on random and structured sets, larger p
    for p in [101, 211, 401, 601, 1009, 2003, 4001, 8009]:
        F = Field(p)
        for n in sorted({4, 8, 12, 16, 24, 32, 48, 64, int(p ** 0.45), int(p ** 0.5) + 1}):
            if n >= p:
                continue
            famA = families(F, n)
            famB = families(F, n)
            for na, A in famA.items():
                for nb, B in famB.items():
                    s = F.S_np(A, B)
                    EAB = energy_sum(p, A, B)
                    EA, EB = energy_sum(p, A, A), energy_sum(p, B, B)
                    assert s * s <= p * EAB, (p, na, nb)
                    assert EAB >= len(A) * len(B)
                    assert EAB * EAB <= EA * EB
                    assert s ** 4 <= p * p * EA * EB
                    bump("S5_energy_bound", 4)
                    r1 = Fraction(p * EAB, len(A) * len(B) * (p - len(B)))
                    if min_ratio_energy_vs_chung is None or r1 < min_ratio_energy_vs_chung:
                        min_ratio_energy_vs_chung = r1
                    if s != 0:
                        r2 = Fraction(p * EAB, s * s)
                        if max_ratio_bound_vs_S is None or r2 > max_ratio_bound_vs_S:
                            max_ratio_bound_vs_S = r2
                    if na == nb and n == int(p ** 0.45):
                        ledger.append({"p": p, "family": na, "n": n, "S": s,
                                       "E_plus_AB": EAB,
                                       "energy_bound_sq_over_chung_sq": str(r1)})
    NOTES["S5_min_ratio_(pE+(A,B))/(|A||B|(p-|B|))"] = str(min_ratio_energy_vs_chung)
    NOTES["S5_max_ratio_(pE+(A,B))/S^2"] = str(max_ratio_bound_vs_S)
    NOTES["S5_ledger_equal_families_at_p^0.45"] = ledger


# ------------------------------------------------------------------ S6
def section6():
    hadamard_growth = []
    paley_neutral = []
    for p in [7, 11, 13, 17, 19, 23, 29, 31]:
        F = Field(p)
        chi, inv = F.chi, F.inv
        for trial in range(6):
            X = [RNG.randint(0, 4) for _ in range(p)]
            Y = [RNG.randint(0, 4) for _ in range(p)]
            if trial == 0:  # interval-like X, sparse Y
                X = [1 if x <= p // 3 else 0 for x in range(p)]
                Y = [1 if x % 3 == 0 else 0 for x in range(p)]
            if trial == 1:  # sparse flat sets
                X, Y = [0] * p, [0] * p
                for x in RNG.sample(range(p), max(2, p // 4)):
                    X[x] = 1
                for y in RNG.sample(range(p), max(2, p // 4)):
                    Y[y] = 1
            WX, WY = sum(X), sum(Y)
            if WX == 0 or WY == 0:
                continue
            # ---- Hadamard kernel zeta^{xy}: Rao, Lemma 3.3 (bias^2 <= bias(X-X, Y)),
            # certified by exact identities in Z[zeta_p].
            z = [0] * p
            for x in range(p):
                for y in range(p):
                    z[(x * y) % p] += X[x] * Y[y]
            h = []
            for y in range(p):
                hy = [0] * p
                for x in range(p):
                    hy[(x * y) % p] += X[x]
                h.append(hy)
            XX = [0] * p
            for x in range(p):
                for x2 in range(p):
                    XX[(x - x2) % p] += X[x] * X[x2]
            w = [0] * p                       # sum_{d,y} (X*X)(d) Y(y) zeta^{dy}
            for d in range(p):
                if XX[d]:
                    for y in range(p):
                        if Y[y]:
                            w[(d * y) % p] += XX[d] * Y[y]
            hh = [cmul(h[y], cconj(h[y], p), p) for y in range(p)]
            w2 = [0] * p                      # sum_y Y(y) |h_y|^2
            for y in range(p):
                if Y[y]:
                    w2 = cadd(w2, cscale(hh[y], Y[y]))
            assert ceq(w, w2), "difference-source identity failed"
            # Lagrange identity:  W_Y * w - |z|^2 = sum_{y<y'} Y(y)Y(y') |h_y - h_{y'}|^2  (>= 0)
            lhs = cadd(cscale(w, WY), cscale(cmul(z, cconj(z, p), p), -1))
            rhs = [0] * p
            for y in range(p):
                if not Y[y]:
                    continue
                for y2 in range(y + 1, p):
                    if not Y[y2]:
                        continue
                    diff = [a - b for a, b in zip(h[y], h[y2])]
                    rhs = cadd(rhs, cscale(cmul(diff, cconj(diff, p), p), Y[y] * Y[y2]))
            assert ceq(lhs, rhs), "Lagrange identity failed"
            bump("S6_hadamard_squaring", 2)
            cpX = Fraction(sum(x * x for x in X), WX * WX)
            cpXX = Fraction(sum(v * v for v in XX), WX ** 4)
            hadamard_growth.append({"p": p, "trial": trial, "cp(X)": str(cpX), "cp(X-X)": str(cpXX),
                                    "log_p(1/cp(X))": round(math.log(1 / cpX) / math.log(p), 3),
                                    "log_p(1/cp(X-X))": round(math.log(1 / cpXX) / math.log(p), 3)})
            # ---- Paley kernel chi(x+y): squaring is entropy-neutral (exact identity)
            PX = [Fraction(x, WX) for x in X]
            PY = [Fraction(y, WY) for y in Y]
            bias = sum(PX[x] * PY[y] * chi[(x + y) % p] for x in range(p) for y in range(p))
            rhs_sum = Fraction(0)
            for x1 in range(p):
                Q = sum(X[x2] * Y[y] * chi[(x1 + y) % p] * chi[(x2 + y) % p]
                        for x2 in range(p) for y in range(p))
                diag = X[x1] * sum(Y[y] for y in range(p) if (x1 + y) % p != 0)
                T = 0
                for d in range(1, p):
                    xd = X[(x1 + d) % p]
                    if not xd:
                        continue
                    for y in range(p):
                        if (x1 + y) % p == 0 or not Y[y]:
                            continue
                        T += xd * Y[y] * chi[d] * chi[(inv[d] + inv[(x1 + y) % p]) % p]
                assert Q == diag + T, (p, x1)
                alpha = [0] * p
                beta = [0] * p
                for u in range(1, p):
                    alpha[u] = chi[u] * X[(x1 + inv[u]) % p]
                    beta[u] = Y[(inv[u] - x1) % p]
                T2 = sum(alpha[u] * beta[v] * chi[(u + v) % p] for u in range(p) for v in range(p))
                assert T2 == T
                assert sum(a * a for a in alpha) == sum(x * x for x in X) - X[x1] ** 2
                assert sum(b * b for b in beta) == sum(y * y for y in Y) - Y[(-x1) % p] ** 2
                assert max(abs(a) for a in alpha) <= max(X)
                assert max(abs(b) for b in beta) <= max(Y)
                a2, b2, b1 = sum(a * a for a in alpha), sum(b * b for b in beta), sum(beta)
                assert T * T <= a2 * (p * b2 - b1 * b1)
                bump("S6_paley_squaring_identity", 7)
                rhs_sum += PX[x1] * Fraction(abs(T), WX * WY)
            assert bias * bias <= cpX + rhs_sum, (p, X, Y)
            bump("S6_paley_squaring_inequality")
            paley_neutral.append({"p": p, "trial": trial, "bias^2": str(bias * bias),
                                  "cp(X)+E|T|": str(cpX + rhs_sum)})
    NOTES["S6_hadamard_growth_examples"] = hadamard_growth[:8]
    NOTES["S6_paley_neutral_examples"] = paley_neutral[:8]


# ------------------------------------------------------------------ S7
def section7():
    rows = []
    for p in [211, 401, 809, 1009, 2003]:
        F = Field(p)
        n = int(p ** 0.45)
        for name, A in families(F, n).items():
            enc = [(a, a * a % p) for a in A]
            cnt = {}
            for e1 in enc:
                for e2 in enc:
                    key = ((e1[0] + e2[0]) % p, (e1[1] + e2[1]) % p)
                    cnt[key] = cnt.get(key, 0) + 1
            E2 = sum(c * c for c in cnt.values())
            assert E2 == 2 * len(A) ** 2 - len(A), (p, name)
            bump("S7_parabola_sidon")
            cnt3 = {}
            for e1 in enc:
                for e2 in enc:
                    for e3 in enc:
                        key = ((e1[0] + e2[0] + e3[0]) % p, (e1[1] + e2[1] + e3[1]) % p)
                        cnt3[key] = cnt3.get(key, 0) + 1
            E3 = sum(c * c for c in cnt3.values())
            N = len(A)
            rows.append({"p": p, "family": name, "n": N,
                         "rho=log n/log p": round(math.log(N) / math.log(p), 3),
                         "E3": E3, "H2(3EncA)/log p": round(math.log(N ** 6 / E3) / math.log(p), 3),
                         "H2(2EncA)/log p": round(math.log(N ** 4 / E2) / math.log(p), 3)})
    NOTES["S7_encoding_growth"] = rows


# ------------------------------------------------------------------ S8
def section8():
    for p in [13, 17, 19, 23, 29, 31]:
        F = Field(p)
        chi, inv = F.chi, F.inv
        for trial in range(8):
            n = RNG.randint(2, p - 2)
            A, B = RNG.sample(range(1, p), n), RNG.sample(range(1, p), n)
            C = RNG.sample(range(p), n)
            direct = sum(chi[(a * b + c) % p] for a in A for b in B for c in C)
            Emul, r = energy_mult(p, A, B)
            fC = F.f(C)
            assert direct == sum(cnt * fC[u] for u, cnt in r.items())
            assert direct * direct <= Emul * n * (p - n)
            # additive condenser A + T against B: sum_{a,t,b} chi(a+t+b) = sum_u r_{A+T}(u) f_B(u)
            T = RNG.sample(range(p), n)
            add_direct = sum(chi[(a + t + b) % p] for a in A for t in T for b in B)
            rT = {}
            for a in A:
                for t in T:
                    u = (a + t) % p
                    rT[u] = rT.get(u, 0) + 1
            fB = F.f(B)
            assert add_direct == sum(cnt * fB[u] for u, cnt in rT.items())
            EAT = sum(cnt * cnt for cnt in rT.values())
            assert add_direct * add_direct <= EAT * n * (p - n)
            assert EAT >= n * n
            bump("S8_additive_condenser_identity", 3)
            direct2 = sum(chi[(c * a + b) % p] for c in C if c for a in A for b in B)
            assert direct2 == sum(chi[c] * F.S(A, [(inv[c] * b) % p for b in B]) for c in C if c)
            bump("S8_condenser_identities", 3)
        for d in [d for d in range(2, p - 1) if (p - 1) % d == 0]:
            H = F.subgroup(d)
            for trial in range(2):
                C = RNG.sample(range(p), RNG.randint(1, p - 1))
                three = sum(chi[(a * b + c) % p] for a in H for b in H for c in C)
                assert three == d * F.S(H, C)
                Emul, _ = energy_mult(p, H, H)
                assert Emul == d ** 3
                bump("S8_subgroup_collapse", 2)
    # witness search: H_avg  |S(A,B)| <= 2 * mean_c |S(A,cB)|
    refuted = []
    for p in [503, 1009, 2003, 4001, 8009]:
        F = Field(p)
        chi = F.chi
        N = int(p ** 0.45)
        A = RNG.sample(range(1, p), N)
        fA = F.f(A)
        order = sorted(range(1, p), key=lambda b: (-fA[b], b))
        B = order[:N]
        SAB = sum(fA[b] for b in B)
        assert SAB == F.S_np(A, B)
        Bn = np.array(B, dtype=np.int64)
        Sc = [F.S_np(A, (c * Bn) % p) for c in range(1, p)]
        assert Sc[0] == SAB
        avg_abs = Fraction(sum(abs(s) for s in Sc), p - 1)
        signed = sum(chi[c] * Sc[c - 1] for c in range(1, p))
        C = RNG.sample(range(1, p), N)
        threeC = sum(chi[c] * Sc[c - 1] for c in C)   # = sum_{c in C, a, b} chi(c^{-1} a + b)... see note
        # direct recomputation of the three-source sum sum_{c in C,a,b} chi(c a + b):
        direct3 = sum(chi[c] * F.S_np(A, (pow(c, p - 2, p) * Bn) % p) for c in C)
        direct3b = int(sum(F.chi_np[(np.int64(c) * np.array(A)[:, None] + Bn[None, :]) % p].sum()
                           for c in C))
        assert direct3 == direct3b
        bump("S8_three_source_direct", 1)
        rec = {"p": p, "N": N, "A": sorted(A), "B": sorted(B),
               "S(A,B)": SAB, "bias2=|S|/N^2": str(Fraction(abs(SAB), N * N)),
               "mean_c|S(A,cB)|/N^2": str(avg_abs / (N * N)),
               "max_{c!=1}|S(A,cB)|/N^2": str(Fraction(max(abs(s) for s in Sc[1:]), N * N)),
               "signed_average_|sum_c chi(c)S(A,cB)|/((p-1)N^2)": str(Fraction(abs(signed), (p - 1) * N * N)),
               "three_source_random_C_|sum chi(ca+b)|/N^3": str(Fraction(abs(direct3), N ** 3)),
               "H_avg_refuted": bool(abs(SAB) > 2 * avg_abs)}
        refuted.append(rec)
        assert abs(SAB) > 2 * avg_abs, "H_avg not refuted at this p"
        bump("S8_H_avg_witness")
    # six-source BIW-type composition chi((ab+c) + (a'b'+c')) with weighted Chung,
    # and the condensed collision probabilities cp(AB+C) as an illustration of BIW Lemma 3.1
    comp_rows = []
    for p in [101, 211, 401, 1009]:
        F = Field(p)
        chi = F.chi
        n = int(p ** 0.45)
        fams = families(F, n)
        triples = {"random": (fams["random"], RNG.sample(range(p), n), RNG.sample(range(p), n)),
                   "interval^3": (fams["interval"], fams["interval"], fams["interval"]),
                   "gp^3": (fams["gp"], fams["gp"], fams["gp"]),
                   "gp,gp,random": (fams["gp"], fams["gp"], RNG.sample(range(p), n))}
        if "subgroup" in fams:
            triples["subgroup,subgroup,random"] = (fams["subgroup"], fams["subgroup"], RNG.sample(range(p), n))
        A2, B2, C2 = (RNG.sample(range(p), n) for _ in range(3))
        r2 = {}
        for a in A2:
            for b in B2:
                ab = a * b % p
                for c in C2:
                    u = (ab + c) % p
                    r2[u] = r2.get(u, 0) + 1
        E2 = sum(v * v for v in r2.values())
        for name, (A, B, C) in triples.items():
            r = {}
            for a in A:
                for b in B:
                    ab = a * b % p
                    for c in C:
                        u = (ab + c) % p
                        r[u] = r.get(u, 0) + 1
            E = sum(v * v for v in r.values())
            tot = sum(r[u] * r2[v] * chi[(u + v) % p] for u in r for v in r2)
            R = n ** 3
            assert tot * tot <= E * (p * E2 - R * R)
            assert tot * tot <= E2 * (p * E - R * R)
            bump("S8_six_source_composition", 2)
            comp_rows.append({"p": p, "triple": name, "n": n,
                              "log_p(1/cp(A))": round(math.log(n) / math.log(p), 3),
                              "log_p(1/cp(AB+C))": round(math.log(n ** 6 / E) / math.log(p), 3),
                              "bias": str(Fraction(abs(tot), n ** 6)),
                              "chung_bound": round(math.sqrt(E * (p * E2 - R * R)) / n ** 6, 4)})
    NOTES["S8_six_source_composition"] = comp_rows
    WITNESSES.extend({"hypothesis": "H_avg: |S(A,B)| <= 2*mean_{c in F_p^*}|S(A,cB)| for all A,B",
                      "status": "REFUTED", **r} for r in refuted)


# ------------------------------------------------------------------ S9
def section9():
    P = primes_upto(2 * 10 ** 6)
    small = P[:100]
    best3, best2 = [], []
    for p in P:
        if p < 5:
            continue
        e = (p - 1) // 2
        n_p = None
        for q in small:
            if pow(q, e, p) == p - 1:
                n_p = q
                break
        assert n_p is not None
        N3, N2 = (n_p - 1) // 3, (n_p - 1) // 2
        if N3 >= 2:
            best3.append((math.log(N3) / math.log(p), p, n_p, N3))
        if N2 >= 2:
            best2.append((math.log(N2) / math.log(p), p, n_p, N2))
        bump("S9_least_nonresidue")
    best3.sort(reverse=True)
    best2.sort(reverse=True)
    out3, out2 = [], []
    for rate, p, n_p, N in best3[:5]:
        e = (p - 1) // 2
        assert all(pow(s, e, p) == 1 for s in range(1, 3 * N + 1))
        total = sum(1 for a in range(1, N + 1) for t in range(1, N + 1) for b in range(1, N + 1))
        assert total == N ** 3
        out3.append({"p": p, "n_p": n_p, "N": N, "rate": round(rate, 4),
                     "sum_{a,t,b in [1,N]} chi(a+t+b)": N ** 3})
        bump("S9_interval_witness_3")
    for rate, p, n_p, N in best2[:5]:
        e = (p - 1) // 2
        assert all(pow(s, e, p) == 1 for s in range(1, 2 * N + 1))
        out2.append({"p": p, "n_p": n_p, "N": N, "rate": round(rate, 4), "S([1,N],[1,N])": N * N})
        bump("S9_interval_witness_2")
    WITNESSES.append({"hypothesis": "chi(a+t+b) on three interval sources [1,N] has bias < 1",
                      "status": "finite witnesses (bias exactly 1 at these p; asymptotics open = Vinogradov)",
                      "cases": out3})
    WITNESSES.append({"hypothesis": "chi(a+b) on two interval sources [1,N] has bias < 1",
                      "status": "finite witnesses (bias exactly 1; asymptotics open = Vinogradov)",
                      "cases": out2})


# ------------------------------------------------------------------ main
def main():
    t0 = time.time()
    OUT.parent.mkdir(exist_ok=True)
    for name, fn in [("S1", section1), ("S2", section2), ("S3", section3), ("S4", section4),
                     ("S5", section5), ("S6", section6), ("S7", section7), ("S8", section8),
                     ("S9", section9)]:
        t1 = time.time()
        try:
            fn()
            status = "ok"
        except AssertionError as exc:
            status = "FAILED: " + repr(exc)
            FAILURES.append({"section": name, "error": repr(exc), "trace": traceback.format_exc()})
        except Exception as exc:  # noqa
            status = "ERROR: " + repr(exc)
            FAILURES.append({"section": name, "error": repr(exc), "trace": traceback.format_exc()})
        print(f"{name}: {status} ({time.time() - t1:.1f}s)")
    total = sum(CHECKS.values())
    result = {
        "status": "passed" if not FAILURES else "FAILED",
        "total_checks": total,
        "checks_by_section": CHECKS,
        "witnesses": WITNESSES,
        "notes": NOTES,
        "failures": FAILURES,
        "seed": 20260905,
        "runtime_seconds": round(time.time() - t0, 1),
        "scope": "Exact finite identities, inequalities and witnesses only; no asymptotic claim.",
    }
    OUT.write_text(json.dumps(result, indent=2))
    print(f"total checks {total}, failures {len(FAILURES)}, wrote {OUT}")
    return 0 if not FAILURES else 1


if __name__ == "__main__":
    raise SystemExit(main())
