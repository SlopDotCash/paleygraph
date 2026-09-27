#!/usr/bin/env python3
"""Stepanov wave, worker `adversary` (2026-09-26): exact verifier.

Checks, by exact integer / rational arithmetic (numpy integer arrays, fractions),
every claim and every stored witness of research/stepanov-adversary-2026-09-26.md:

  A. the fixed-m (Weil-range) limit of R_e: the binomial constants g(m,e), the
     character-sum expansion of |B_e(A)| (exact at small p), the explicit
     Weil error bound (on random sets at several primes), the finite + entropy
     proof that g(m,e) <= 7/8 whenever 1 <= e < m/6, and the step function
     F(eta) = sup_{e/m <= eta} g(m,e);
  B. the exhaustive small-prime tables of the C++ helper: every witness is
     recomputed, the maxima are re-derived independently in Python for
     p <= 47 (all m <= floor(sqrt p)) and p <= 23 (m <= floor(sqrt p) + 3),
     and Hanson-Petridis (R_0 <= 1) is checked on every table entry;
  C. F_{p^2}: the field is built explicitly (x^2 = nu), chi_q is computed by
     exponentiation and compared with chi_p(Norm); HP with r retained fails for
     A = B = F_p; R_e of subfield and near-subfield constructions; the exact
     exhaustive F_9, F_25, F_49 tables and all F_{p^2} search witnesses;
  D. every stored RHP search witness over F_p (|B_e|, r_e, obj) and the
     derived empirical f table;
  E. every stored rectangle witness (task 5): S(A,B) recomputed exactly.

Writes results/stepanov_adversary_2026_09_26.json.  Standard library + numpy.
Heavy parts use a multiprocessing pool of at most 6 workers.
"""
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SEARCH = ROOT / "results" / "stepanov_adversary_2026_09_26_search.json"
OUT = ROOT / "results" / "stepanov_adversary_2026_09_26.json"
NPROC = int(os.environ.get("ADV_NPROC", "6"))

CHECKS = {"n": 0, "fail": 0}
FAILS = []


def check(cond, label, data=None):
    CHECKS["n"] += 1
    if not cond:
        CHECKS["fail"] += 1
        if len(FAILS) < 200:
            FAILS.append({"label": label, "data": data})


def add_counts(n, fails):
    CHECKS["n"] += n
    CHECKS["fail"] += len(fails)
    FAILS.extend(fails[: max(0, 200 - len(FAILS))])


# ------------------------------------------------------------------ arithmetic

def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def chi_table(p):
    """Legendre symbol table from the set of squares (exact, no Euler criterion)."""
    chi = np.full(p, -1, dtype=np.int64)
    sq = (np.arange(1, p, dtype=np.int64) ** 2) % p
    chi[sq] = 1
    chi[0] = 0
    return chi


def least_nonresidue(p):
    chi = chi_table(p)
    return int(np.nonzero(chi == -1)[0][0])


def profile_fp(p, chi, A):
    """e_b for all b and the index set -A."""
    A = np.asarray(sorted(set(int(a) % p for a in A)), dtype=np.int64)
    b = np.arange(p, dtype=np.int64)
    M = chi[(A[:, None] + b[None, :]) % p]
    eb = (M == -1).sum(axis=0)
    return A, eb, (-A) % p


def obj_table(eb, negA, m):
    """For e = 0..m: (|B_e|, r_e, obj_e = m|B_e| - r_e) as Python ints."""
    hist = np.bincount(eb, minlength=m + 1)
    cum = np.cumsum(hist)
    ebn = eb[negA]
    out = []
    for e in range(m + 1):
        Bs = int(cum[e])
        r = int((ebn <= e).sum())
        out.append((Bs, r, m * Bs - r))
    return out


def comb(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


_CUM = {}


def cbin(m, e):
    """sum_{j <= e} C(m, j), cached per m."""
    if m not in _CUM:
        acc, row = 0, []
        c = 1
        for j in range(m + 1):
            acc += c
            row.append(acc)
            c = c * (m - j) // (j + 1)
        _CUM[m] = row
    if e < 0:
        return 0
    return _CUM[m][min(e, m)]


def g_limit(m, e):
    """g(m,e) = 2m 2^{-m} sum_{j<=e} C(m,j): the p -> infinity value of R_e at fixed m."""
    return Fraction(2 * m * cbin(m, e), 2 ** m)


# ------------------------------------------------------------------ F_{p^2}

class Fq:
    """F_{p^2} = F_p[x]/(x^2 - nu), nu the least non-residue; element u + v x has index u + p v."""

    def __init__(self, p):
        self.p = p
        self.q = p * p
        self.nu = least_nonresidue(p)
        q = self.q
        U = np.arange(q) % p
        V = np.arange(q) // p
        self.U, self.V = U, V
        # chi_q by exponentiation z^((q-1)/2) in the field
        chi = np.zeros(q, dtype=np.int64)
        for idx in range(1, q):
            w = self.pow((int(U[idx]), int(V[idx])), (q - 1) // 2)
            if w == (1, 0):
                chi[idx] = 1
            elif w == (p - 1, 0):
                chi[idx] = -1
            else:
                raise AssertionError("chi_q not +-1")
        self.chi = chi
        chip = chi_table(p)
        N = (U * U - self.nu * V * V) % p
        self.norm_ok = bool(np.all(chi[1:] == chip[N[1:]]))

    def mul(self, a, b):
        p, nu = self.p, self.nu
        return ((a[0] * b[0] + nu * a[1] * b[1]) % p, (a[0] * b[1] + a[1] * b[0]) % p)

    def pow(self, a, k):
        r = (1, 0)
        while k:
            if k & 1:
                r = self.mul(r, a)
            a = self.mul(a, a)
            k >>= 1
        return r

    def idx(self, u, v):
        return (u % self.p) + self.p * (v % self.p)

    def add_idx(self, A, B):
        """matrix of indices of a + b for a in A, b in B"""
        p = self.p
        A = np.asarray(A)
        B = np.asarray(B)
        return (self.U[A][:, None] + self.U[B][None, :]) % p + p * ((self.V[A][:, None] + self.V[B][None, :]) % p)

    def neg_idx(self, A):
        p = self.p
        A = np.asarray(A)
        return (-self.U[A]) % p + p * ((-self.V[A]) % p)

    def profile(self, A):
        A = np.asarray(sorted(set(int(a) for a in A)), dtype=np.int64)
        M = self.chi[self.add_idx(A, np.arange(self.q))]
        eb = (M == -1).sum(axis=0)
        return A, eb, self.neg_idx(A)


# ------------------------------------------------------------------ part A

def part_A(rng):
    res = {}
    # A1: g(m,e) <= 7/8 for 1 <= e < m/6, m <= 400, equality only at (7,1)
    worst = (Fraction(0), None)
    for m in range(7, 401):
        for e in range(1, (m - 1) // 6 + 1):
            if 6 * e >= m:
                continue
            g = g_limit(m, e)
            check(g <= Fraction(7, 8), "A1 g<=7/8", (m, e))
            if g == Fraction(7, 8):
                check((m, e) == (7, 1), "A1 equality only at (7,1)", (m, e))
            if g > worst[0]:
                worst = (g, (m, e))
    res["A1_max_g_below_one_sixth"] = {"value": str(worst[0]), "at": worst[1], "m_range": [7, 400]}
    # A2: entropy tail for m > 400: sum_{j<=e} C(m,j) <= 2^{m H(e/m)} (e <= m/2), H(1/6) < 0.6501
    H16 = -(1 / 6) * math.log2(1 / 6) - (5 / 6) * math.log2(5 / 6)
    check(H16 < 0.6501, "A2 H(1/6)")
    # 2m 2^{-m(1-0.6501)} < 7/8 for all m >= 401: decreasing for m >= 5, check at 401
    check(2 * 401 * 2 ** (-401 * (1 - 0.6501)) < 7 / 8, "A2 tail")
    # entropy inequality spot check (exact integers) at m = 401..420
    for m in range(401, 421):
        e = (m - 1) // 6
        check(cbin(m, e) <= 2 ** (m * H16) * (1 + 1e-12), "A2 entropy spot", m)
    res["A2_entropy"] = {"H_1_6": H16}
    # A3: table of g(m,e), and the step function F(eta) = sup_{e/m <= eta} g(m,e)
    table = []
    for m in range(1, 41):
        for e in range(0, m):
            if 2 * e >= m and not (m == 1 and e == 0):
                continue
            g = g_limit(m, e)
            table.append({"m": m, "e": e, "g": str(g), "g_float": float(g)})
    res["A3_g_table_m_le_40"] = table

    def F_step(eta, mmax=400):
        best, arg = Fraction(0), None
        for m in range(1, mmax + 1):
            e = math.floor(eta * m + 1e-12)
            if 2 * e >= m and m > 1:
                e = (m - 1) // 2
            g = g_limit(m, e)
            if g > best:
                best, arg = g, (m, e)
        return best, arg

    Fvals = {}
    for name, eta in [("1/32", 1 / 32), ("1/16", 1 / 16), ("1/8", 1 / 8), ("1/7", 1 / 7), ("0.16", 0.16),
                      ("1/6", 1 / 6), ("1/5", 1 / 5), ("2/9", 2 / 9), ("1/4", 1 / 4), ("2/7", 2 / 7),
                      ("3/10", 0.3), ("1/3", 1 / 3), ("2/5", 0.4), ("0.45", 0.45)]:
        v, arg = F_step(eta)
        ceiling = 2 / (1 - 2 * eta) ** 2
        Fvals[name] = {"eta": eta, "F": str(v), "F_float": float(v), "argmax_m_e": arg,
                       "second_moment_ceiling": ceiling, "ratio_to_ceiling": float(v) / ceiling}
    check(Fvals["1/7"]["F"] == "1", "A3 F(1/7) = 1")
    check(Fvals["1/6"]["F"] == "21/16", "A3 F(1/6) = 21/16")
    check(Fvals["1/5"]["F"] == "15/8", "A3 F(1/5) = 15/8")
    check(Fvals["1/4"]["F"] == "5/2", "A3 F(1/4) = 5/2")
    res["A3_F_step"] = Fvals

    # A4: exact character-sum expansion of #{b not in -A : e_b <= e} at small p
    n_exp = 0
    for p in (11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        chi = chi_table(p)
        for m in (2, 3, 4, 5):
            for _ in range(3):
                A = sorted(rng.sample(range(p), m))
                Aarr, eb, negA = profile_fp(p, chi, A)
                notneg = np.ones(p, dtype=bool)
                notneg[negA] = False
                # T(I) = sum_{b not in -A} chi(prod_{a in I}(a+b))
                T = {}
                for k in range(0, m + 1):
                    for I in itertools.combinations(range(m), k):
                        prod = np.ones(p, dtype=np.int64)
                        for i in I:
                            prod = prod * chi[(A[i] + np.arange(p)) % p]
                        T[I] = int(prod[notneg].sum())
                for e in range(0, m):
                    lhs = int(((eb <= e) & notneg).sum())
                    tot = 0
                    for k in range(0, e + 1):
                        for J in itertools.combinations(range(m), k):
                            Js = set(J)
                            for I, t in T.items():
                                tot += (-1) ** len(Js.intersection(I)) * t
                    check(tot == lhs * 2 ** m, "A4 expansion", (p, A, e))
                    n_exp += 1
    res["A4_expansion_checks"] = n_exp

    # A5: explicit Weil error bound on random sets
    def weil_bound(m, e, p):
        C = cbin(m, e)
        W = (m * 2 ** (m - 1) - m) + comb(m, 2) + math.sqrt(p) * ((m - 2) * 2 ** (m - 1) + 1 - comb(m, 2))
        return m + C * W / 2 ** m

    n_w, worst_ratio = 0, 0.0
    for p in (1009, 10007, 100003):
        chi = chi_table(p)
        for m in range(2, 9):
            for _ in range(4):
                A = rng.sample(range(p), m)
                Aarr, eb, negA = profile_fp(p, chi, A)
                tab = obj_table(eb, negA, m)
                for e in range(0, (m + 1) // 2):
                    dev = abs(tab[e][0] - p * cbin(m, e) / 2 ** m)
                    bnd = weil_bound(m, e, p)
                    check(dev <= bnd, "A5 Weil bound", (p, m, e))
                    worst_ratio = max(worst_ratio, dev / bnd)
                    n_w += 1
    res["A5_weil_bound_checks"] = {"n": n_w, "max_dev_over_bound": worst_ratio}

    # A6: explicit range of Theorem 3: largest M such that for all 7 <= m <= M and 1 <= e < m/6
    #     g(m,e) + err(m,e,p) <= 1, err = g/(p-1) + 2m(m+1+2^{-m} C W_m(p))/(p-1)
    def err_bound(m, e, p):
        C = cbin(m, e)
        W = (m * 2 ** (m - 1) - m) + comb(m, 2) + math.sqrt(p) * ((m - 2) * 2 ** (m - 1) + 1 - comb(m, 2))
        g = float(g_limit(m, e))
        return g / (p - 1) + 2 * m * (m + 1 + C * W / 2 ** m) / (p - 1)

    ranges = {}
    for k in (6, 9, 12, 15, 20, 30, 50, 100):
        p = 10.0 ** k
        M = 6
        for m in range(7, 400):
            ok = all(float(g_limit(m, e)) + err_bound(m, e, p) <= 1 for e in range(1, (m - 1) // 6 + 1) if 6 * e < m)
            if not ok:
                break
            M = m
        ranges[f"1e{k}"] = {"M": M, "log2_p": k * math.log2(10), "M_over_log2_p": M / (k * math.log2(10))}
    res["A6_theorem3_explicit_range"] = ranges
    # A7: Weil-limit share s_lim(m,e) = (1/m) 2^{-m} sum_{j<=e} C(m,j)(m-2j)^2 satisfies
    #     g(m,e)(1-2e/m)^2 <= 2 s_lim(m,e) <= 1 (binomial symmetry), equality in the second only for
    #     e = floor((m-1)/2) or (m,e) in {(1,0),(2,0)}; hence F(eta) <= (1-2 eta)^{-2}.
    n7, worst7 = 0, Fraction(0)
    for m in range(1, 201):
        for e in range(0, (m + 1) // 2):
            if 2 * e >= m and m > 1:
                continue
            sl = Fraction(sum(comb(m, j) * (m - 2 * j) ** 2 for j in range(e + 1)), m * 2 ** m)
            lhs = g_limit(m, e) * Fraction(m - 2 * e, m) ** 2
            check(lhs <= 2 * sl <= 1, "A7 share sandwich", (m, e))
            if 2 * sl == 1:
                check(e == (m - 1) // 2 or (m, e) in ((1, 0), (2, 0)), "A7 equality cases", (m, e))
            worst7 = max(worst7, lhs)
            n7 += 1
    res["A7_share_limit"] = {"n": n7, "max_g_times_(1-2eta)^2": str(worst7)}
    # the bound itself must hold on the random sets of A5 (it is implied by A5's bound): spot check
    for p in (1009, 10007):
        chi = chi_table(p)
        d = (p - 1) // 2
        for m in range(7, 9):
            for _ in range(3):
                A = rng.sample(range(p), m)
                Aarr, eb, negA = profile_fp(p, chi, A)
                tab = obj_table(eb, negA, m)
                for e in range(0, (m + 1) // 2):
                    R = Fraction(tab[e][2], d)
                    check(abs(float(R) - float(g_limit(m, e))) <= err_bound(m, e, p) + 1e-9, "A6 R-g bound", (p, m, e))
    return res


# ------------------------------------------------------------------ part B (exhaustive F_p)

def py_exhaustive(args):
    """Independent Python maximum of obj_e over normalised A of size m (contains {0,1} or {0,nu})."""
    p, m = args
    chi = chi_table(p)
    nu = least_nonresidue(p)
    bad = (chi[(np.arange(p)[:, None] + np.arange(p)[None, :]) % p] == -1).astype(np.int16)
    best = {}
    if m == 1:
        A, eb, negA = profile_fp(p, chi, [0])
        for e, (Bs, r, o) in enumerate(obj_table(eb, negA, 1)):
            best[e] = o
        return p, m, best
    for s in (1, nu):
        cands = [c for c in range(p) if c not in (0, s)]
        if m == 2:
            combos = np.zeros((1, 0), dtype=np.int64)
        else:
            combos = np.array(list(itertools.combinations(cands, m - 2)), dtype=np.int64).reshape(-1, m - 2)
        for lo in range(0, len(combos), 50000):
            C = combos[lo:lo + 50000]
            N = len(C)
            full = np.concatenate([np.zeros((N, 1), np.int64), np.full((N, 1), s, np.int64), C], axis=1)
            eb = bad[0][None, :] + bad[s][None, :]
            eb = np.broadcast_to(eb, (N, p)).copy()
            for j in range(m - 2):
                eb += bad[C[:, j]]
            negA = (-full) % p
            ebn = np.take_along_axis(eb, negA, axis=1)
            for e in range(0, m + 1):
                cnt = (eb <= e).sum(axis=1)
                r = (ebn <= e).sum(axis=1)
                o = int((m * cnt - r).max())
                if o > best.get(e, -1):
                    best[e] = o
    return p, m, best


def recheck_exh_record(rec):
    """Recompute one exhaustive-table witness."""
    p, m, e = rec["p"], rec["m"], rec["e"]
    chi = chi_table(p)
    A, eb, negA = profile_fp(p, chi, rec["A"])
    ok = len(A) == m
    Bs, r, o = obj_table(eb, negA, m)[e]
    return ok and (Bs, r, o) == (rec["Bsize"], rec["r"], rec["obj"])


def part_B(search, pool):
    res = {}
    recs = search.get("exh", [])
    done = {r["p"]: r for r in search.get("exh_done", [])}
    res["primes"] = sorted(done)
    res["mmax"] = {p: done[p]["mmax"] for p in sorted(done)}
    # witnesses
    oks = pool.map(recheck_exh_record, recs, chunksize=32)
    for rec, ok in zip(recs, oks):
        check(ok, "B witness recompute", (rec["p"], rec["m"], rec["e"]))
    res["n_witnesses_rechecked"] = len(recs)
    # HP: R_0 <= 1, i.e. obj_0 <= d, on every table entry; record the equality cases
    # Prop 5(ii): s_0(A) <= (d - r + r/m)/(p - m) on every e = 0 witness
    eq_cases = []
    for rec in recs:
        if rec["e"] == 0:
            d = (rec["p"] - 1) // 2
            check(rec["obj"] <= d, "B HP R_0<=1", (rec["p"], rec["m"]))
            p, m = rec["p"], rec["m"]
            if m < p:
                chi = chi_table(p)
                Aarr, eb, negA = profile_fp(p, chi, rec["A"])
                Fa = (m - 2 * eb).astype(np.int64)
                Fa[negA] -= 1
                s0 = Fraction(int((Fa[eb == 0] ** 2).sum()), m * (p - m))
                check(s0 <= Fraction(d - rec["r"], p - m) + Fraction(rec["r"], m * (p - m)), "B Prop5(ii) share",
                      (p, m))
            if rec["obj"] == d:
                eq_cases.append({"p": rec["p"], "m": rec["m"], "A": rec["A"], "Bsize": rec["Bsize"], "r": rec["r"]})
    res["HP_equality_cases"] = eq_cases
    # independent Python maxima
    tasks = []
    for p in sorted(done):
        if p > 47:
            continue
        mtop = int(math.isqrt(p)) + (3 if p <= 23 else 0)
        for m in range(1, min(mtop, done[p]["mmax"]) + 1):
            tasks.append((p, m))
    tab = {(r["p"], r["m"], r["e"]): r["obj"] for r in recs}
    n_cmp = 0
    for p, m, best in pool.map(py_exhaustive, tasks, chunksize=1):
        for e, o in best.items():
            if (p, m, e) in tab:
                check(o == tab[(p, m, e)], "B python max == C++ max", (p, m, e, o, tab[(p, m, e)]))
                n_cmp += 1
    res["n_python_maxima_compared"] = n_cmp
    res["python_exhaustive_tasks"] = len(tasks)
    # summary table: max R_e for m <= floor(sqrt p) and, separately, all m in range
    summary = []
    for rec in recs:
        p, m, e = rec["p"], rec["m"], rec["e"]
        if 2 * e >= m and m > 1:
            continue
        d = (p - 1) // 2
        share = None
        if m < p:
            chi = chi_table(p)
            Aarr, eb, negA = profile_fp(p, chi, rec["A"])
            Fa = (m - 2 * eb).astype(np.int64)
            Fa[negA] -= 1
            share = float(Fraction(int((Fa[eb <= e] ** 2).sum()), m * (p - m)))
        summary.append({"p": p, "m": m, "e": e, "share": share, "R": str(Fraction(rec["obj"], d)),
                        "R_float": rec["obj"] / d, "A": rec["A"], "Bsize": rec["Bsize"], "r": rec["r"],
                        "m_le_sqrt_p": m * m <= p, "n_normalised_attaining": rec["n_normalised_attaining"]})
    res["table"] = summary
    # sup of R_e over the exhaustive range with e/m < 1/6, e >= 1, m <= sqrt p
    small = [s for s in summary if s["m_le_sqrt_p"] and s["e"] >= 1 and 6 * s["e"] < s["m"]]
    res["max_R_e_ge1_eta_lt_1_6_m_le_sqrtp"] = max(small, key=lambda s: s["R_float"]) if small else None
    small2 = [s for s in summary if s["e"] >= 1 and 6 * s["e"] < s["m"]]
    res["max_R_e_ge1_eta_lt_1_6_all_m"] = max(small2, key=lambda s: s["R_float"]) if small2 else None
    sh = [s for s in summary if s["share"] is not None and s["m"] >= 3 and s["m"] * s["m"] <= s["p"]]
    res["max_share_m_ge3_le_sqrtp"] = max(sh, key=lambda s: s["share"]) if sh else None
    sh1 = [s for s in sh if s["e"] >= 1 and 6 * s["e"] < s["m"]]
    res["max_share_e_ge1_eta_lt_1_6"] = max(sh1, key=lambda s: s["share"]) if sh1 else None
    return res


# ------------------------------------------------------------------ part C (F_{p^2})

def part_C(search, rng):
    res = {}
    fields = {}
    constructions = []
    for p in (3, 5, 7, 11, 13, 17, 19, 23, 29, 31):
        F = Fq(p)
        fields[p] = F
        check(F.norm_ok, "C chi_q = chi_p(Norm)", p)
        q, d = F.q, (F.q - 1) // 2
        sub = [F.idx(u, 0) for u in range(p)]
        # A + B subset of Q_q u {0} for A = B = F_p
        M = F.chi[F.add_idx(sub, sub)]
        check(bool(np.all(M >= 0)), "C subfield sums are squares", p)
        r = p
        check(p * p - r > d, "C HP with r fails over F_{p^2}", p)
        # second moment of A = F_p: F_A(b) = p-1 on F_p, -1 off F_p
        Fa = F.chi[F.add_idx(sub, np.arange(q))].sum(axis=0)
        onsub = np.zeros(q, dtype=bool)
        onsub[sub] = True
        check(bool(np.all(Fa[onsub] == p - 1)) and bool(np.all(Fa[~onsub] == -1)), "C F_A values", p)
        check(int((Fa * Fa).sum()) == p * (q - p), "C second moment identity", p)

        def R_of(A, e):
            Aarr, eb, negA = F.profile(A)
            m = len(Aarr)
            Bs, rr, o = obj_table(eb, negA, m)[e]
            share = Fraction(int((F.chi[F.add_idx(Aarr, np.nonzero(eb <= e)[0])].sum(axis=0) ** 2).sum()),
                             m * (q - m))
            return m, Bs, rr, o, Fraction(o, d), share

        rows = []
        m, Bs, rr, o, R, sh = R_of(sub, 0)
        check(R == Fraction(2 * p, p + 1), "C R_0(F_p) = 2p/(p+1)", p)
        rows.append({"construction": "subfield F_p", "m": m, "e": 0, "Bsize": Bs, "r": rr, "R": str(R),
                     "R_float": float(R), "share": float(sh)})
        off = [x for x in range(q) if not onsub[x]]
        for k in (1, 2, 3):
            if k > p // 3:
                continue
            X = rng.sample(off, k)
            A = sub + X
            m, Bs, rr, o, R, sh = R_of(A, k)
            check(R >= Fraction(2 * p, p + 1), "C R_k(F_p u X) >= 2p/(p+1)", (p, k))
            rows.append({"construction": f"F_p plus {k} random points", "m": m, "e": k, "eta": k / m, "Bsize": Bs,
                         "r": rr, "R": str(R), "R_float": float(R), "share": float(sh)})
            Y = rng.sample(sub, k)
            A = [x for x in sub if x not in Y] + X
            m, Bs, rr, o, R, sh = R_of(A, k)
            rows.append({"construction": f"F_p with {k} points swapped out", "m": m, "e": k, "eta": k / m,
                         "Bsize": Bs, "r": rr, "R": str(R), "R_float": float(R), "share": float(sh)})
        for frac in (0.5, 0.75):
            k = max(2, int(round(frac * p)))
            A = rng.sample(sub, k)
            m, Bs, rr, o, R, sh = R_of(A, 0)
            rows.append({"construction": f"random {k}-subset of F_p", "m": m, "e": 0, "Bsize": Bs, "r": rr,
                         "R": str(R), "R_float": float(R), "share": float(sh)})
        # a square multiple of the subfield (another line through 0)
        cands = [x for x in off if F.chi[x] == 1]
        c = cands[0]
        cu, cv = int(F.U[c]), int(F.V[c])
        line = [F.idx(*F.mul((cu, cv), (u, 0))) for u in range(p)]
        m, Bs, rr, o, R, sh = R_of(line, 0)
        rows.append({"construction": "square multiple c*F_p", "m": m, "e": 0, "Bsize": Bs, "r": rr, "R": str(R),
                     "R_float": float(R), "share": float(sh)})
        # random baseline of the same size
        A = rng.sample(range(q), p)
        for e in (0, 1, 2):
            if 2 * e < p:
                m, Bs, rr, o, R, sh = R_of(A, e)
                rows.append({"construction": "random p-set", "m": m, "e": e, "Bsize": Bs, "r": rr, "R": str(R),
                             "R_float": float(R), "share": float(sh)})
        constructions.append({"p": p, "q": q, "rows": rows})
    res["constructions"] = constructions
    # exhaustive F_q tables
    ex = search.get("exhq", [])
    for rec in ex:
        F = fields.get(rec["p"]) or Fq(rec["p"])
        fields[rec["p"]] = F
        Aarr, eb, negA = F.profile(rec["A"])
        Bs, rr, o = obj_table(eb, negA, rec["m"])[rec["e"]]
        check(len(Aarr) == rec["m"] and (Bs, rr, o) == (rec["Bsize"], rec["r"], rec["obj"]), "C exhq witness",
              (rec["p"], rec["m"], rec["e"]))
    exq = []
    for rec in ex:
        if 2 * rec["e"] < rec["m"] or rec["m"] == 1:
            d = (rec["q"] - 1) // 2
            exq.append({"q": rec["q"], "m": rec["m"], "e": rec["e"], "R": str(Fraction(rec["obj"], d)),
                        "R_float": rec["obj"] / d, "A": rec["A"]})
    res["exhaustive_Fq"] = exq
    # search witnesses over F_q
    out_s = []
    for key in ("rhpq",):
        for rec in search.get(key, []):
            F = fields.get(rec["p"]) or Fq(rec["p"])
            fields[rec["p"]] = F
            Aarr, eb, negA = F.profile(rec["A"])
            Bs, rr, o = obj_table(eb, negA, rec["m"])[rec["e"]]
            check(len(Aarr) == rec["m"] and (Bs, rr, o) == (rec["Bsize"], rec["r"], rec["obj"]), "C rhpq witness",
                  (rec["p"], rec["m"], rec["e"]))
            # is the witness contained in a line (affine subfield translate)?
            pts = [(int(F.U[a]), int(F.V[a])) for a in Aarr]
            on_line = max_on_affine_line(pts, rec["p"])
            out_s.append({"p": rec["p"], "q": rec["q"], "m": rec["m"], "e": rec["e"],
                          "R_float": o / ((rec["q"] - 1) / 2), "max_points_on_one_line": on_line})
    res["rhpq"] = out_s
    rq = []
    for rec in search.get("rectq", []):
        F = fields.get(rec["p"]) or Fq(rec["p"])
        S = int(F.chi[F.add_idx(rec["A"], rec["B"])].sum())
        check(S == rec["S"] and len(set(rec["A"])) == rec["m"] and len(set(rec["B"])) == rec["n"], "C rectq witness",
              (rec["p"], rec["m"]))
        pts = [(int(F.U[a]), int(F.V[a])) for a in rec["A"]]
        rq.append({"p": rec["p"], "q": rec["q"], "m": rec["m"], "n": rec["n"], "S": S,
                   "bias": S / (rec["m"] * rec["n"]), "mn_over_q": rec["m"] * rec["n"] / rec["q"],
                   "max_points_of_A_on_one_line": max_on_affine_line(pts, rec["p"])})
    res["rectq"] = rq
    return res


def max_on_affine_line(pts, p):
    """Largest number of the points (u,v) in F_p^2 lying on one affine F_p-line."""
    best = min(len(pts), 1)
    n = len(pts)
    for i in range(n):
        dirs = {}
        for j in range(n):
            if i == j:
                continue
            du, dv = (pts[j][0] - pts[i][0]) % p, (pts[j][1] - pts[i][1]) % p
            if du:
                key = (1, dv * pow(du, -1, p) % p)
            else:
                key = (0, 1)
            dirs[key] = dirs.get(key, 0) + 1
        if dirs:
            best = max(best, 1 + max(dirs.values()))
    return best


# ------------------------------------------------------------------ part D (RHP search witnesses)

def recheck_rhp_group(args):
    p, recs = args
    chi = chi_table(p)
    d = (p - 1) // 2
    out = []
    fails = []
    for rec in recs:
        A, eb, negA = profile_fp(p, chi, rec["A"])
        m = rec["m"]
        tab = obj_table(eb, negA, m)
        ok = len(A) == m and tab[rec["e"]] == (rec["Bsize"], rec["r"], rec["obj"])
        if not ok:
            fails.append({"label": "D rhp witness", "data": (p, m, rec["e"], rec.get("family"))})
        # share of the second moment carried by B_e
        Fa = m - 2 * eb
        Fa[negA] -= 1  # b = -a: one zero sum, e_b bad and m-1-e_b good -> F = m - 1 - 2 e_b
        good = eb <= rec["e"]
        share = Fraction(int((Fa[good].astype(np.int64) ** 2).sum()), m * (p - m))
        e = rec["e"]
        if m > 2 * e + 1:
            bound = Fraction(2 * (p - m), p - 1) * share / (Fraction(m - 2 * e - 1, m) ** 2)
            if not Fraction(rec["obj"], d) <= bound:
                fails.append({"label": "D Prop5(i)", "data": (p, m, e)})
        allR = [str(Fraction(t[2], d)) for t in tab[: (m + 1) // 2 + 1]]
        out.append({"p": p, "m": m, "e": rec["e"], "family": rec.get("family"), "job": rec.get("job"),
                    "obj": rec["obj"], "R": str(Fraction(rec["obj"], d)), "R_float": rec["obj"] / d,
                    "share": float(share), "R_all_e": allR})
    # sanity of the F_A formula on the first record
    if recs:
        rec = recs[0]
        A, eb, negA = profile_fp(p, chi, rec["A"])
        Fa = m_minus = len(A) - 2 * eb
        Fa = Fa.copy()
        Fa[negA] -= 1
        b = np.arange(p)
        direct = chi[(A[:, None] + b[None, :]) % p].sum(axis=0)
        if not np.array_equal(direct, Fa):
            fails.append({"label": "D F_A formula", "data": p})
    return p, out, fails, 2 * len(recs) + 1


def part_D(search, pool):
    res = {}
    groups = {}
    for rec in search.get("rhp", []) + search.get("rhpcal", []):
        groups.setdefault(rec["p"], []).append(rec)
    allrows = []
    ntot = 0
    for p, out, fails, n in pool.map(recheck_rhp_group, sorted(groups.items()), chunksize=1):
        allrows.extend(out)
        add_counts(n, fails)
        ntot += n
    res["n_witnesses_rechecked"] = ntot
    # best per (p, m, e)
    best = {}
    for r in allrows:
        if r["job"] and r["job"].endswith(":cal"):
            continue
        k = (r["p"], r["m"], r["e"])
        if k not in best or r["obj"] > best[k]["obj"]:
            best[k] = r
    rows = []
    for (p, m, e), r in sorted(best.items()):
        eta = e / m
        g = g_limit(m, e)
        rows.append({"p": p, "m": m, "e": e, "eta": eta, "R_best": r["R_float"], "R_exact": r["R"],
                     "family": r["family"], "g_limit": float(g), "excess_over_limit": r["R_float"] - float(g),
                     "second_moment_ceiling": 2 / (1 - 2 * eta) ** 2,
                     "ratio_to_ceiling": r["R_float"] / (2 / (1 - 2 * eta) ** 2), "share": r["share"]})
    res["best_per_p_m_e"] = rows
    # family win counts
    wins = {}
    for r in rows:
        wins[r["family"]] = wins.get(r["family"], 0) + 1
    res["family_wins"] = wins
    # empirical f: for each p and eta in {1/16,1/8,1/4} with e = floor(eta m); and e = 1, 2 fixed
    emp = []
    for p in sorted({r["p"] for r in rows}):
        for lab, eta in (("1/16", 1 / 16), ("1/8", 1 / 8), ("1/4", 1 / 4)):
            cand = [r for r in rows if r["p"] == p and r["e"] == math.floor(eta * r["m"]) and r["m"] >= 8]
            if cand:
                top = max(cand, key=lambda r: r["R_best"])
                emp.append({"p": p, "eta": lab, "max_R": top["R_best"], "at_m": top["m"], "at_e": top["e"],
                            "per_m": {r["m"]: round(r["R_best"], 4) for r in sorted(cand, key=lambda r: r["m"])}})
        for e in (1, 2):
            cand = [r for r in rows if r["p"] == p and r["e"] == e and r["m"] >= 8]
            if cand:
                top = max(cand, key=lambda r: r["R_best"])
                emp.append({"p": p, "e_fixed": e, "max_R": top["R_best"], "at_m": top["m"],
                            "per_m": {r["m"]: round(r["R_best"], 4) for r in sorted(cand, key=lambda r: r["m"])}})
    res["empirical_f"] = emp
    # the sup over all runs with e >= 1 and e/m < 1/6
    sub = [r for r in rows if r["e"] >= 1 and 6 * r["e"] < r["m"]]
    res["max_R_e_ge1_eta_lt_1_6"] = max(sub, key=lambda r: r["R_best"]) if sub else None
    res["max_ratio_to_ceiling"] = max(rows, key=lambda r: r["ratio_to_ceiling"]) if rows else None
    res["max_share"] = max(rows, key=lambda r: r["share"]) if rows else None
    # calibration against exact maxima
    cal = []
    tab = {(r["p"], r["m"], r["e"]): r["obj"] for r in search.get("exh", [])}
    calbest = {}
    for r in allrows:
        if r["job"] and r["job"].endswith(":cal"):
            k = (r["p"], r["m"], r["e"])
            calbest[k] = max(calbest.get(k, -1), r["obj"])
    for k, o in sorted(calbest.items()):
        if k in tab:
            check(o <= tab[k], "D calibration search <= exact max", k)
            cal.append({"p": k[0], "m": k[1], "e": k[2], "search_obj": o, "exact_max": tab[k], "hit": o == tab[k]})
    res["calibration_vs_exact"] = cal
    return res


# ------------------------------------------------------------------ part E (rectangles)

def recheck_rect_group(args):
    p, recs = args
    chi = chi_table(p)
    out, fails = [], []
    for rec in recs:
        A = np.asarray(rec["A"], dtype=np.int64)
        B = np.asarray(rec["B"], dtype=np.int64)
        S = int(chi[(A[:, None] + B[None, :]) % p].sum())
        zeros = int(((A[:, None] + B[None, :]) % p == 0).sum())
        ok = S == rec["S"] and len(set(rec["A"])) == rec["m"] and len(set(rec["B"])) == rec["n"]
        if not ok:
            fails.append({"label": "E rect witness", "data": (p, rec["m"], rec["n"])})
        # column profile: e_b for b in B
        eb = (chi[(A[:, None] + B[None, :]) % p] == -1).sum(axis=0)
        out.append({"p": p, "m": rec["m"], "n": rec["n"], "S": S, "bias": S / (rec["m"] * rec["n"]),
                    "mn_over_p": rec["m"] * rec["n"] / p, "zeros": zeros, "job": rec["job"],
                    "max_e_b": int(eb.max()), "mean_e_b_over_m": float(eb.mean() / rec["m"])})
    return p, out, fails, len(recs)


def part_E(search, pool):
    res = {}
    groups = {}
    for rec in search.get("rect", []):
        groups.setdefault(rec["p"], []).append(rec)
    rows = []
    n = 0
    for p, out, fails, k in pool.map(recheck_rect_group, sorted(groups.items()), chunksize=1):
        rows.extend(out)
        add_counts(k, fails)
        n += k
    res["n_witnesses_rechecked"] = n
    best = {}
    for r in rows:
        parts = r["job"].split(":")
        k = (r["p"], parts[3], parts[4])  # (p, rho, shape)
        if k not in best or r["bias"] > best[k]["bias"]:
            best[k] = r
    table = []
    for (p, rho, shape), r in sorted(best.items(), key=lambda kv: (kv[0][2], float(kv[0][1]), kv[0][0])):
        table.append({"p": p, "rho_target": float(rho), "shape": shape, "m": r["m"], "n": r["n"],
                      "mn_over_p": r["mn_over_p"], "max_bias": r["bias"], "S": r["S"],
                      "sqrt_log_p_over_sqrt_mn": math.sqrt(math.log(p)) / (r["m"] * r["n"]) ** 0.25})
    res["max_bias_table"] = table
    return res


# ------------------------------------------------------------------ part F (random-function null model)

def null_table(p, seed):
    """Must agree with experiments/stepanov_adversary_2026_09_26_driver.py: balanced random +-1 on F_p^*."""
    rng = random.Random(1000003 * p + seed)
    d = (p - 1) // 2
    signs = [1] * d + [-1] * d
    rng.shuffle(signs)
    return np.asarray([0] + signs, dtype=np.int64)


def recheck_null_group(args):
    p, rhp_recs, rect_recs = args
    chi = null_table(p, 1)
    d = (p - 1) // 2
    out_r, out_q, fails = [], [], []
    for rec in rhp_recs:
        A, eb, negA = profile_fp(p, chi, rec["A"])
        m = rec["m"]
        tab = obj_table(eb, negA, m)
        if not (len(A) == m and tab[rec["e"]] == (rec["Bsize"], rec["r"], rec["obj"])):
            fails.append({"label": "F null rhp witness", "data": (p, m, rec["e"])})
        out_r.append({"p": p, "m": m, "e": rec["e"], "obj": rec["obj"], "R_float": rec["obj"] / d,
                      "family": rec["family"]})
    for rec in rect_recs:
        A = np.asarray(rec["A"], dtype=np.int64)
        B = np.asarray(rec["B"], dtype=np.int64)
        S = int(chi[(A[:, None] + B[None, :]) % p].sum())
        if not (S == rec["S"] and len(set(rec["A"])) == rec["m"] and len(set(rec["B"])) == rec["n"]):
            fails.append({"label": "F null rect witness", "data": (p, rec["m"], rec["n"])})
        out_q.append({"p": p, "m": rec["m"], "n": rec["n"], "bias": S / (rec["m"] * rec["n"]),
                      "mn_over_p": rec["m"] * rec["n"] / p, "rho": rec["job"].split(":")[3]})
    return p, out_r, out_q, fails, len(rhp_recs) + len(rect_recs)


def part_F(search, pool, D, E):
    res = {}
    groups = {}
    for rec in search.get("rhpnull", []):
        groups.setdefault(rec["p"], ([], []))[0].append(rec)
    for rec in search.get("rectnull", []):
        groups.setdefault(rec["p"], ([], []))[1].append(rec)
    tasks = [(p, a, b) for p, (a, b) in sorted(groups.items())]
    rows_r, rows_q = [], []
    n = 0
    for p, out_r, out_q, fails, k in pool.map(recheck_null_group, tasks, chunksize=1):
        rows_r.extend(out_r)
        rows_q.extend(out_q)
        add_counts(k, fails)
        n += k
    res["n_witnesses_rechecked"] = n
    # exact: max over |A| = 2 of R_0 for the null function (translation invariance: A = {0, c});
    # for the Legendre symbol this maximum is exactly 1 (HP equality at A = {0,1})
    hp2 = []
    for p in (2003, 8009):
        for label, chi in (("null", null_table(p, 1)), ("legendre", chi_table(p))):
            gmask = (chi != -1).astype(np.int64)
            d = (p - 1) // 2
            best, arg = -1, None
            for c in range(1, p):
                Bs = int((gmask * np.roll(gmask, -c)).sum())  # b with f(b) != -1 and f(b+c) != -1
                r = int(gmask[c] == 1) + int(gmask[(-c) % p] == 1)  # b = 0 needs f(c) != -1; b = -c needs f(-c) != -1
                o = 2 * Bs - r
                if o > best:
                    best, arg = o, c
            hp2.append({"p": p, "function": label, "max_obj_m2": best, "d": d, "R0_max_m2": str(Fraction(best, d)),
                        "R0_float": best / d, "A": [0, arg]})
            if label == "legendre":
                check(best == d, "F HP equality at m = 2 for Legendre", p)
    res["hp_m2_null_vs_legendre"] = hp2
    # the null function has no HP: record max R_0 (m = 2 and others)
    bestn = {}
    for r in rows_r:
        k = (r["p"], r["m"], r["e"])
        if k not in bestn or r["obj"] > bestn[k]["obj"]:
            bestn[k] = r
    paley = {(r["p"], r["m"], r["e"]): r["R_best"] for r in D.get("best_per_p_m_e", [])}
    cmp_rows = []
    for k, r in sorted(bestn.items()):
        cmp_rows.append({"p": k[0], "m": k[1], "e": k[2], "eta": k[2] / k[1], "R_null": r["R_float"],
                         "R_paley": paley.get(k), "g_limit": float(g_limit(k[1], k[2]))})
    res["paley_vs_null_rhp"] = cmp_rows
    bq = {}
    for r in rows_q:
        k = (r["p"], r["rho"])
        if k not in bq or r["bias"] > bq[k]["bias"]:
            bq[k] = r
    pal = {(r["p"], str(r["rho_target"])): r["max_bias"] for r in E.get("max_bias_table", []) if r["shape"] == "bal"}
    res["paley_vs_null_rect"] = [{"p": k[0], "rho": float(k[1]), "m": r["m"], "n": r["n"], "bias_null": r["bias"],
                                  "bias_paley": pal.get((k[0], str(float(k[1]))))} for k, r in sorted(bq.items())]
    return res


# ------------------------------------------------------------------ part G (multiplicative subgroups, exact scan)

def prime_factors(n):
    out, k = [], 2
    while k * k <= n:
        if n % k == 0:
            out.append(k)
            while n % k == 0:
                n //= k
        k += 1
    if n > 1:
        out.append(n)
    return out


def subgroup_scan_prime(p):
    """R_e(tH) for every subgroup H of F_p^* of order m, 4 <= m <= sqrt p (R_e(tH) = R_e(H) for t in Q;
    t non-residue is the other class).  Uses e_{bh} = e_b (chi(h) = 1) and m - e_b - z_b (chi(h) = -1)."""
    chi = chi_table(p)
    d = (p - 1) // 2
    fs = prime_factors(p - 1)
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // f, p) != 1 for f in fs))
    nu = least_nonresidue(p)
    out = []
    for m in range(4, math.isqrt(p) + 1):
        if (p - 1) % m:
            continue
        k = (p - 1) // m
        H = np.array([pow(g, k * i, p) for i in range(m)], dtype=np.int64)
        for t in (1, nu):
            A = (t * H) % p
            reps = np.array([pow(g, j, p) for j in range(k)], dtype=np.int64)
            # e_b and zero-count for b in reps (b ranges over coset reps of H)
            S = (A[:, None] + reps[None, :]) % p
            eb_rep = (chi[S] == -1).sum(axis=0)
            z_rep = (S == 0).sum(axis=0)
            chiH = chi[H]
            npos = int((chiH == 1).sum())
            nneg = m - npos
            e0 = int((chi[A] == -1).sum())
            negA = (-A) % p
            # full e_b only needed on -A (for r_e) and b = 0
            S2 = (A[:, None] + negA[None, :]) % p
            eb_negA = (chi[S2] == -1).sum(axis=0)
            row = {"p": p, "m": m, "t": "Q" if t == 1 else "nu", "H_in_Q": nneg == 0, "R": {}}
            for e in sorted({0, 1, 2, m // 16, m // 8, m // 4}):
                if 2 * e >= m:
                    continue
                Bs = int(e0 <= e) + npos * int((eb_rep <= e).sum()) + nneg * int(((m - eb_rep - z_rep) <= e).sum())
                r = int((eb_negA <= e).sum())
                row["R"][e] = Fraction(m * Bs - r, d)
            out.append(row)
    # direct cross-check of the coset formula on the first two groups
    fails = []
    for row in out[:2]:
        m = row["m"]
        k = (p - 1) // m
        H = np.array([pow(g, k * i, p) for i in range(m)], dtype=np.int64)
        A = (H * (1 if row["t"] == "Q" else nu)) % p
        Aarr, eb, negA = profile_fp(p, chi, A)
        tab = obj_table(eb, negA, m)
        for e, R in row["R"].items():
            if Fraction(tab[e][2], d) != R:
                fails.append({"label": "G coset formula", "data": (p, m, e)})
    return out, fails


def part_G(pool):
    primes = [p for p in range(1001, 60000, 2) if is_prime(p)][::6]
    res = {"n_primes": len(primes), "p_range": [primes[0], primes[-1]]}
    allrows = []
    n = 0
    for out, fails in pool.map(subgroup_scan_prime, primes, chunksize=8):
        allrows.extend(out)
        add_counts(len(out) + 1, fails)
        n += len(out)
    res["n_subgroup_cosets"] = n
    for row in allrows:
        if 0 in row["R"]:
            check(row["R"][0] <= 1, "G HP on subgroups", (row["p"], row["m"]))
    # max R_e over the scan, by eta class, for e >= 1 and e/m < 1/6, and for e = floor(eta m)
    best = {}
    for row in allrows:
        m = row["m"]
        for e, R in row["R"].items():
            for lab, eta in (("1/16", 1 / 16), ("1/8", 1 / 8), ("1/4", 1 / 4)):
                if e == math.floor(eta * m) and m >= 8:
                    key = "eta=" + lab
                    if key not in best or R > best[key]["R"]:
                        best[key] = {"R": R, "p": row["p"], "m": m, "e": e, "t": row["t"]}
            if e >= 1 and 6 * e < m:
                if "e>=1,eta<1/6" not in best or R > best["e>=1,eta<1/6"]["R"]:
                    best["e>=1,eta<1/6"] = {"R": R, "p": row["p"], "m": m, "e": e, "t": row["t"]}
            if e == 0:
                if "e=0" not in best or R > best["e=0"]["R"]:
                    best["e=0"] = {"R": R, "p": row["p"], "m": m, "e": e, "t": row["t"]}
    for v in best.values():
        v["R_float"] = float(v["R"])
        v["g_limit"] = float(g_limit(v["m"], v["e"]))
        v["R"] = str(v["R"])
    res["max_by_class"] = best
    # growth with m: max over primes of R_1 for each m bucket
    buck = {}
    for row in allrows:
        if 1 in row["R"]:
            m = row["m"]
            b = "4-7" if m < 8 else ("8-15" if m < 16 else ("16-31" if m < 32 else ("32-63" if m < 64 else "64+")))
            v = float(row["R"][1])
            if b not in buck or v > buck[b]["R_float"]:
                buck[b] = {"R_float": v, "p": row["p"], "m": m, "t": row["t"]}
    res["max_R1_by_m_bucket"] = buck
    return res


# ------------------------------------------------------------------ part H (second-moment share search)

def share_value(F_vals, thr):
    F_vals = F_vals.astype(np.int64)
    return int((F_vals[F_vals >= thr] ** 2).sum())


def recheck_share_group(args):
    key, recs = args
    kind, p = key
    fails, out = [], []
    if kind == "fq":
        F = Fq(p)
        q = F.q
    else:
        chi = chi_table(p) if kind == "fp" else null_table(p, 1)
        q = p
    for rec in recs:
        m, thr = rec["m"], rec["thr"]
        A = np.asarray(sorted(set(rec["A"])), dtype=np.int64)
        if kind == "fq":
            Fa = F.chi[F.add_idx(A, np.arange(q))].sum(axis=0)
        else:
            Fa = chi[(A[:, None] + np.arange(p)[None, :]) % p].sum(axis=0)
        obj = share_value(Fa, thr)
        if not (len(A) == m and obj == rec["obj"]):
            fails.append({"label": "H share witness", "data": (kind, p, m, thr)})
        # the second moment identity holds for the quadratic characters of F_p and F_q only;
        # the random +-1 null model ("fr") has no such identity, so its share is normalised by
        # its actual second moment (root fix 2026-09-27)
        total = int((Fa.astype(np.int64) ** 2).sum())
        if kind in ("fp", "fq") and total != m * (q - m):
            fails.append({"label": "H second moment identity", "data": (kind, p, m)})
        denom = m * (q - m) if kind in ("fp", "fq") else total
        out.append({"kind": kind, "p": p, "q": q, "m": m, "thr": thr, "family": rec["family"],
                    "share": obj / denom, "share_exact": str(Fraction(obj, denom))})
    return out, fails, 2 * len(recs)


def part_H(search, pool):
    res = {}
    groups = {}
    for rec in search.get("share", []):
        f = str(rec["field"])
        kind = "fp" if f == "fp" else ("fq" if f == "fq" else "fr")
        groups.setdefault((kind, rec["p"]), []).append(rec)
    rows = []
    n = 0
    for out, fails, k in pool.map(recheck_share_group, sorted(groups.items()), chunksize=1):
        rows.extend(out)
        add_counts(k, fails)
        n += k
    res["n_checks"] = n
    best = {}
    for r in rows:
        m = r["m"]
        if r["thr"] == 1:
            lab = "positive part"
        elif r["thr"] == m - 2 * (m // 8) - 1:
            lab = "rows with e<=m/8"
        else:
            lab = "rows with e<=m/4"
        k = (r["kind"], r["p"], m, lab)
        if k not in best or r["share"] > best[k]["share"]:
            best[k] = dict(r, cls=lab)
    res["best"] = [best[k] for k in sorted(best)]
    return res


# ------------------------------------------------------------------ main

def main():
    t0 = time.time()
    c0 = time.process_time()
    rng = random.Random(20260926)
    search = json.load(SEARCH.open()) if SEARCH.exists() else {}
    out = {"description": "Exact verifier output for research/stepanov-adversary-2026-09-26.md"}
    out["A_fixed_m_limit"] = part_A(rng)
    print("part A done", CHECKS, flush=True)
    with Pool(NPROC) as pool:
        out["B_exhaustive_Fp"] = part_B(search, pool)
        print("part B done", CHECKS, flush=True)
        out["C_Fp2"] = part_C(search, rng)
        print("part C done", CHECKS, flush=True)
        out["D_rhp_search"] = part_D(search, pool)
        print("part D done", CHECKS, flush=True)
        out["E_rectangles"] = part_E(search, pool)
        print("part E done", CHECKS, flush=True)
        out["F_null_model"] = part_F(search, pool, out["D_rhp_search"], out["E_rectangles"])
        print("part F done", CHECKS, flush=True)
        out["G_subgroups"] = part_G(pool)
        print("part G done", CHECKS, flush=True)
        out["H_share"] = part_H(search, pool)
        print("part H done", CHECKS, flush=True)
    out["checks"] = CHECKS["n"]
    out["failures"] = CHECKS["fail"]
    out["failure_samples"] = FAILS
    out["wall_seconds"] = round(time.time() - t0, 1)
    out["main_process_cpu_seconds"] = round(time.process_time() - c0, 1)
    OUT.write_text(json.dumps(out, indent=1, default=str))
    print(f"checks={CHECKS['n']} failures={CHECKS['fail']} wall={out['wall_seconds']}s -> {OUT}")
    return 0 if CHECKS["fail"] == 0 else 1


if __name__ == "__main__":
    sys.exit(main())
