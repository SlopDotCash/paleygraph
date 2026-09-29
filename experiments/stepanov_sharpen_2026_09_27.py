#!/usr/bin/env python3
"""Exact verifier for research/stepanov-sharpen-2026-09-27.md (worker `sharpen`).

Run with /opt/miniconda3/bin/python3 (needs numpy, scipy, sympy).  Writes
results/stepanov_sharpen_2026_09_27.json.

Sections (all inequalities checked in exact integer / rational arithmetic unless the
section says "float"):
  A. Re-check of (star) (k=2) and of its order-k analogue (Prop. 2.9 of the robust note)
     exhaustively over all A containing 0 for small p.
  B. Lemma 1.1 (mean / Jensen form of (star)) and Theorem 1.4 of the note (pair form,
     |S| <= u(max(w',1/4)) mn; result keys "B_thm13") on exhaustive A and many B (top-n
     sets, level sets, random sets).
  C. Theorem 1.3 of the note (the one-variable inequality min(HP, Jensen) <= u(w)): the
     analytic constants, a rigorous rational interval verification of the finite range
     m <= 10, w in [1/4, 0.55], and a dense float sanity scan.
  D. The linear programme of task 1: exact rational LPs (sympy) for many (p,m,n), the
     comparison LP <= u, float LPs for large m, and the explicit primal profiles of
     Proposition 2.1 (optimality of eta_star within (star)).
  E. Exhaustive maxima of the bias for p <= 23 compared with the exact LP.
  F. Fixed-|A| (Weil range) profiles: the counts N_j(A) against the Weil error bars, and
     the refined upper bound on the saving (beta_m).
  G. Balanced sets: biclique data, fourth-moment inequality.
  H. Characters of order k: Jensen form per value, the vertex bound, the final statement,
     sharpness example.
"""
import itertools, json, math, os, random, sys, time
from fractions import Fraction as Fr
import numpy as np

T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_sharpen_2026_09_27.json")
R = {"checks": {}, "witnesses": {}, "tables": {}, "notes": {}, "timing": {}}
FAIL = []


def bump(key, n=1):
    R["checks"][key] = R["checks"].get(key, 0) + n


def fail(key, info):
    FAIL.append((key, info))
    R["witnesses"].setdefault("FAILURES", []).append({"key": key, "info": str(info)[:600]})


def check(key, cond, info=""):
    bump(key)
    if not cond:
        fail(key, info)
    return cond


# ---------------------------------------------------------------- basics
def primes_upto(n):
    s = [True] * (n + 1); s[0] = s[1] = False
    for i in range(2, int(n ** .5) + 1):
        if s[i]:
            s[i * i::i] = [False] * len(s[i * i::i])
    return [i for i in range(n + 1) if s[i]]


def leg(x, p):
    x %= p
    if x == 0:
        return 0
    return 1 if pow(x, (p - 1) // 2, p) == 1 else -1


def legtab(p):
    return np.array([leg(x, p) for x in range(p)], dtype=np.int64)


def primroot(p):
    fac = [q for q in primes_upto(p) if (p - 1) % q == 0]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise ValueError


def dlog_table(p):
    g = primroot(p); L = [None] * p; x = 1
    for i in range(p - 1):
        L[x] = i; x = x * g % p
    return L


# exact square-root bounds: lo <= sqrt(X) <= hi for a nonnegative Fraction X
SQN = 10 ** 15


def sqrt_lo(X):
    X = Fr(X)
    if X <= 0:
        return Fr(0)
    q = (X.numerator * SQN * SQN) // X.denominator
    return Fr(math.isqrt(q), SQN)


def sqrt_hi(X):
    X = Fr(X)
    if X <= 0:
        return Fr(0)
    q = (X.numerator * SQN * SQN) // X.denominator
    return Fr(math.isqrt(q) + 1, SQN)


def u_float(w):
    return (math.sqrt(12 * w - 3 * w * w) - w) / 2


def u_lo(w):
    w = Fr(w); return (sqrt_lo(12 * w - 3 * w * w) - w) / 2


def u_hi(w):
    w = Fr(w); return (sqrt_hi(12 * w - 3 * w * w) - w) / 2


def le_u_exact(x, w):
    """exact test of x <= u(w) = (sqrt(12w-3w^2)-w)/2 for rationals x, 0<w<=4."""
    x = Fr(x); w = Fr(w)
    lhs = 2 * x + w
    if lhs < 0:
        return True
    return lhs * lhs <= 12 * w - 3 * w * w


def star_weight2(m, e, j, delta):
    """2 * (e+1-j) * (m - (3e+j)/2 - delta)  (integer)."""
    return (e + 1 - j) * (2 * m - 3 * e - j - 2 * delta)


def W_E(m, E, x):
    """W_E(x) = (E-x)(m+3/2-(3E+x)/2) for x<=E, 0 beyond (exact, x rational)."""
    x = Fr(x)
    if x >= E:
        return Fr(0)
    return (E - x) * (Fr(2 * m + 3 - 3 * E) - x) / 2


def jensen_hbound(m, E, w):
    """upper bound h_m(E,w) = -a-beta+2 sqrt(a^2-a w+(1+beta) w)  (float)."""
    beta = 1.5 / m; a = 1 - 2 * E / m + beta
    return -a - beta + 2 * math.sqrt(a * a - a * w + (1 + beta) * w)


def jensen_hbound_hi(m, E, w):
    beta = Fr(3, 2 * m); a = 1 - Fr(2 * E, m) + beta; w = Fr(w)
    return -a - beta + 2 * sqrt_hi(a * a - a * w + (1 + beta) * w)


def hp_bound(m, w):
    return 1 - Fr(2, m) * (1 - Fr(w))


def best_bound_float(m, w):
    b = 1 - 2 * (1 - w) / m
    for E in range(1, (m + 1) // 2 + 1):
        b = min(b, jensen_hbound(m, E, w))
    return b


# ---------------------------------------------------------------- set enumeration helpers
def all_sets_containing0(p):
    """0/1 matrix (N x p) of all subsets of F_p containing 0 (row = indicator)."""
    N = 1 << (p - 1)
    idx = np.arange(N, dtype=np.int64)
    X = np.zeros((N, p), dtype=np.int8)
    X[:, 0] = 1
    for i in range(1, p):
        X[:, i] = (idx >> (i - 1)) & 1
    return X


def bad_matrix(p, k):
    """Bad[a,b] = 1 iff a+b != 0 and (a+b)^{(p-1)/k} != 1 (for k=2: chi(a+b) = -1)."""
    dk = (p - 1) // k
    bad = np.array([0 if x == 0 else (1 if pow(x, dk, p) != 1 else 0) for x in range(p)], dtype=np.int64)
    M = np.zeros((p, p), dtype=np.int64)
    for a in range(p):
        M[a] = bad[(a + np.arange(p)) % p]
    return M


# ---------------------------------------------------------------- A. (star) and (star)_k
def section_A():
    t = time.time()
    R["tables"]["A_star"] = []
    for (p, klist) in [(7, [2, 3]), (11, [2, 5]), (13, [2, 3, 4, 6]), (17, [2, 4, 8]), (19, [2, 3, 6, 9])]:
        X = all_sets_containing0(p)
        m = X.sum(axis=1).astype(np.int64)
        negidx = [(-b) % p for b in range(p)]
        delta = X[:, negidx].astype(np.int64)          # delta[A,b] = [b in -A]
        for k in klist:
            dk = (p - 1) // k
            Bad = bad_matrix(p, k)
            e = X.astype(np.int64) @ Bad                # e[A,b]
            ncheck = 0; worst = -10 ** 9; wit = None
            for ee in range(0, (p + 1) // 2):
                if 2 * ee > dk:
                    break
                ok_rows = (2 * ee <= m - 1) & (m + dk - 1 <= p - 1)
                if not ok_rows.any():
                    continue
                w = np.where(e <= ee, (ee + 1 - e) * (2 * m[:, None] - 3 * ee - e - 2 * delta), 0)
                lhs = w.sum(axis=1)
                rhs = 2 * (ee + 1) * (dk - ee)
                diff = np.where(ok_rows, lhs - rhs, -10 ** 9)
                ncheck += int(ok_rows.sum())
                i = int(np.argmax(diff))
                if diff[i] > worst:
                    worst = int(diff[i]); wit = (ee, [int(v) for v in np.nonzero(X[i])[0]])
            check("A_star" if k == 2 else "A_star_k", worst <= 0, (p, k, worst, wit))
            bump("A_star_instances" if k == 2 else "A_stark_instances", ncheck - 1)
            R["tables"]["A_star"].append({"p": p, "k": k, "instances": ncheck, "max_LHS_minus_RHS_x2": worst,
                                           "argmax(e,A)": wit})
    # random larger p for order k
    rng = random.Random(20260927)
    for (p, klist) in [(31, [2, 3, 5, 6]), (37, [2, 3, 4, 6, 9]), (61, [3, 4, 5, 6])]:
        for k in klist:
            dk = (p - 1) // k
            Bad = bad_matrix(p, k)
            worst = -10 ** 9
            for trial in range(400):
                mm = rng.randint(1, min(p - dk, (p + 1) // 2 + 3))
                A = rng.sample(range(p), mm)
                x = np.zeros(p, dtype=np.int64); x[A] = 1
                e = x @ Bad
                delta = x[[(-b) % p for b in range(p)]]
                for ee in range(0, (mm - 1) // 2 + 1):
                    if 2 * ee > dk:
                        break
                    w = np.where(e <= ee, (ee + 1 - e) * (2 * mm - 3 * ee - e - 2 * delta), 0)
                    dd = int(w.sum()) - 2 * (ee + 1) * (dk - ee)
                    worst = max(worst, dd)
                    check("A_star" if k == 2 else "A_star_k", dd <= 0, (p, k, sorted(A), ee))
            R["tables"]["A_star"].append({"p": p, "k": k, "random_sets": 400, "max_LHS_minus_RHS_x2": worst})
    R["timing"]["A"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- B. Lemma 1.1 and Theorem 1.3
def jensen_and_thm13_rows(m, d, n, r, N, S, key):
    """vectorised exact checks for arrays (same shape) of pairs (A,B) with |A|=m (scalar).
    Lemma 1.1: for 1<=E<=(m+1)/2:  N < E n  =>  (En-N)(n(2m+3-3E)-N) <= 2 n E (d-E+1+r).
    Theorem 1.3: if d+r < mn then |S| <= u(max(w',1/4)) mn, w'=(d+r)/(mn)."""
    n = n.astype(object); r = r.astype(object); N = N.astype(object); S = S.astype(object)
    cnt = 0; bad = []
    for E in range(1, (m + 1) // 2 + 1):
        act = N < E * n
        lhs = (E * n - N) * (n * (2 * m + 3 - 3 * E) - N)
        rhs = 2 * n * E * (d - E + 1 + r)
        viol = act & (lhs > rhs)
        cnt += int(act.sum())
        if viol.any():
            i = int(np.nonzero(viol)[0][0]); bad.append(("jensen", m, E, int(n[i]), int(r[i]), int(N[i])))
    bump(key + "_jensen", cnt)
    mn = m * n
    dr = d + r
    act = dr < mn
    big = 4 * dr >= mn
    for sgn in (1, -1):
        SS = sgn * S
        t1 = 2 * SS + dr
        okA = (t1 < 0) | (t1 * t1 <= 12 * dr * mn - 3 * dr * dr)
        t2 = 8 * SS + mn
        okB = (t2 < 0) | (t2 * t2 <= 45 * mn * mn)
        ok = np.where(big, okA, okB)
        viol = act & ~ok
        bump(key + "_thm13", int(act.sum()))
        if viol.any():
            i = int(np.nonzero(viol)[0][0]); bad.append(("thm13", m, int(n[i]), int(r[i]), int(S[i])))
    return bad


def section_B():
    t = time.time()
    rng = np.random.default_rng(12345)
    summary = []
    for p in [7, 11, 13, 17]:
        d = (p - 1) // 2
        X = all_sets_containing0(p)
        if p == 17:
            X = X[rng.choice(X.shape[0], 20000, replace=False)]
        m_all = X.sum(axis=1)
        delta = X[:, [(-b) % p for b in range(p)]].astype(np.int64)
        e = X.astype(np.int64) @ bad_matrix(p, 2)
        F = m_all[:, None] - 2 * e - delta
        tight = 0
        for m in range(1, (p + 1) // 2 + 1):
            rows = np.nonzero(m_all == m)[0]
            if len(rows) == 0:
                continue
            ee, dd, FF = e[rows], delta[rows], F[rows]
            for order in ("top", "bottom"):
                o = np.argsort(-FF if order == "top" else FF, axis=1, kind="stable")
                es = np.take_along_axis(ee, o, 1).cumsum(1)
                ds = np.take_along_axis(dd, o, 1).cumsum(1)
                Fs = np.take_along_axis(FF, o, 1).cumsum(1)
                for n in range(1, p + 1):
                    bad = jensen_and_thm13_rows(m, d, np.full(len(rows), n), ds[:, n - 1], es[:, n - 1],
                                                Fs[:, n - 1], "B")
                    for bb in bad:
                        fail("B", (p,) + bb)
            # level sets {e_b <= e'} and random subsets
            for ep in range(0, m + 1):
                mask = ee <= ep
                n = mask.sum(1); keep = n > 0
                bad = jensen_and_thm13_rows(m, d, n[keep], (dd * mask).sum(1)[keep], (ee * mask).sum(1)[keep],
                                            (FF * mask).sum(1)[keep], "B")
                for bb in bad:
                    fail("B", (p, "level") + bb)
            for rep in range(3):
                mask = rng.random(ee.shape) < rng.random()
                n = mask.sum(1); keep = n > 0
                bad = jensen_and_thm13_rows(m, d, n[keep], (dd * mask).sum(1)[keep], (ee * mask).sum(1)[keep],
                                            (FF * mask).sum(1)[keep], "B")
                for bb in bad:
                    fail("B", (p, "random") + bb)
        summary.append({"p": p, "sets_A": int(X.shape[0])})
    R["tables"]["B_sets"] = summary
    R["timing"]["B"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- C. Theorem 1.2
def section_C():
    t = time.time()
    out = {}
    # C1: analytic constants used in the proof of Theorem 1.2
    # (i) 1-u(w) <= eps^2/2 for eps = 1-w <= sqrt(3)-1 (algebraic proof in the note; grid sanity check)
    worst = -1.0
    for i in range(1, 7321):
        eps = i / 10000.0; w = 1 - eps
        worst = max(worst, (1 - u_float(w)) - eps * eps / 2)
    check("C1_one_minus_u_le_eps2_over_2", worst <= 1e-15, worst)
    out["max_(1-u)-eps^2/2_on_grid"] = worst
    # (ii) M_J(eps) <= 4/eps for eps <= 9/20: (1+2/(3 sqrt c))(3-2eps)/(2(1-eps)) <= 4 at eps=9/20 (monotone)
    eps = Fr(9, 20); w = 1 - eps; c = w * (4 - w) / 4
    val = (1 + Fr(2, 3) / sqrt_lo(c)) * (3 - 2 * eps) / (2 * (1 - eps))
    check("C1_MJ_le_MHP_at_0.45", val <= 4, float(val)); out["MJ*eps_at_eps=0.45(upper)"] = float(val)
    # (iii) on eps in [9/20,3/4]: M_J <= (1+2/(3 sqrt c(3/4))) * max(f(9/20),f(3/4))/2 < 11
    f = lambda e: (3 - 2 * e) / (e * (1 - e))
    e1, e2 = Fr(9, 20), Fr(3, 4)
    c2 = (1 - e2) * (4 - (1 - e2)) / 4
    MJmax = (1 + Fr(2, 3) / sqrt_lo(c2)) * max(f(e1), f(e2)) / 2
    check("C1_MJ_lt_11_on_[0.45,0.75]", MJmax < 11, float(MJmax)); out["MJ_max_on_[0.45,0.75]"] = float(MJmax)
    # f has its only critical point at (6-sqrt12)/4 in (9/20,3/4) and it is a minimum (f' sign change - to +)
    fp = lambda e: -2 * e * e + 6 * e - 3
    check("C1_f_shape", fp(Fr(9, 20)) < 0 and fp(Fr(3, 4)) > 0, "")
    # (iv) the Taylor / gain inequalities are re-evaluated numerically for many (m,w) in the Jensen range
    cnt = 0; worst = -1e9
    for w in [0.25 + 0.75 * i / 400 for i in range(400)]:
        eps = 1 - w; c = w * (4 - w) / 4
        MJ = (1 + 2 / (3 * math.sqrt(c))) * (3 - 2 * eps) / (2 * w * eps)
        for m in [math.ceil(MJ), math.ceil(MJ) + 1, 2 * math.ceil(MJ), 10 * math.ceil(MJ)]:
            b = min(jensen_hbound(m, E, w) for E in range(1, (m + 1) // 2 + 1))
            worst = max(worst, b - u_float(w)); cnt += 1
    check("C1_jensen_range_float", worst <= 1e-12, worst); bump("C1_jensen_range_float_pairs", cnt - 1)
    out["C1_max(J-u)_for_m>=MJ"] = worst
    # C2: rigorous rational interval verification, m<=10, w in [1/4, 11/20]
    boxes = 0; maxdepth = 0; margins = []
    for m in range(1, 11):
        stack = [(Fr(1, 4) + Fr(3, 10) * i / 64, Fr(1, 4) + Fr(3, 10) * (i + 1) / 64, 0) for i in range(64)]
        minmargin = None
        while stack:
            w1, w2, dep = stack.pop()
            ulo = u_lo(w1)
            opts = [hp_bound(m, w2)] + [jensen_hbound_hi(m, E, w2) for E in range(1, (m + 1) // 2 + 1)]
            best = min(opts)
            if best <= ulo:
                boxes += 1; maxdepth = max(maxdepth, dep)
                mg = float(ulo - best); minmargin = mg if minmargin is None else min(minmargin, mg)
            elif dep < 30:
                mid = (w1 + w2) / 2; stack += [(w1, mid, dep + 1), (mid, w2, dep + 1)]
            else:
                fail("C2_interval", (m, float(w1), float(w2)))
        margins.append((m, minmargin))
    bump("C2_interval_boxes", boxes)
    out["C2_boxes"] = boxes; out["C2_maxdepth"] = maxdepth; out["C2_min_margin_by_m"] = margins
    # C3: float scan m <= 3000, w in [1/4, 0.999]
    worst = (-1, None)
    ws = [0.25 + 0.749 * i / 150 for i in range(151)]
    for m in list(range(1, 400)) + list(range(400, 3001, 20)):
        for w in ws:
            v = best_bound_float(m, w) - u_float(w)
            if v > worst[0] or worst[1] is None:
                worst = (v, (m, w))
    check("C3_float_scan", worst[0] <= 1e-12, worst)
    bump("C3_float_scan_pairs", (399 + 131) * 151 - 1)
    out["C3_max(min(HP,J)-u)"] = worst
    R["tables"]["C"] = out
    R["timing"]["C"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- D. the linear programme
def build_lp(p, m, n, dk=None):
    """LP of task 1 in the counts: variables n_j (b not in -A, e_b=j, 0<=j<=m) and n'_j
    (b in -A, 0<=j<=m-1).  maximise S = sum (m-2j) n_j + sum (m-1-2j) n'_j subject to
    (star)_e for 0<=e<=(m-1)/2 (and 2e<=dk), sum n'_j <= min(m,n), sum = n.
    Returned in 'minimise c.x, A x <= b, Aeq x = beq, x >= 0' form with Fractions."""
    d = (p - 1) // 2 if dk is None else dk
    nv = (m + 1) + m
    c = [Fr(-(m - 2 * j)) for j in range(m + 1)] + [Fr(-(m - 1 - 2 * j)) for j in range(m)]
    A = []; b = []
    for e in range(0, (m - 1) // 2 + 1):
        if 2 * e > d:
            break
        row = [Fr(0)] * nv
        for j in range(e + 1):
            row[j] = Fr(star_weight2(m, e, j, 0), 2)
            if j <= m - 1:
                row[m + 1 + j] = Fr(star_weight2(m, e, j, 1), 2)
        A.append(row); b.append(Fr((e + 1) * (d - e)))
    A.append([Fr(0)] * (m + 1) + [Fr(1)] * m); b.append(Fr(min(m, n)))
    Aeq = [[Fr(1)] * nv]; beq = [Fr(n)]
    return c, A, b, Aeq, beq


def solve_lp_exact(c, A, b, Aeq, beq):
    """exact primal and dual via sympy; returns (value_of_max, x, y, z) after verifying
    primal feasibility, dual feasibility and equality of objectives in exact arithmetic."""
    from sympy import Matrix, Rational
    from sympy.solvers.simplex import linprog as slp
    toR = lambda v: Rational(v.numerator, v.denominator)

    def toF(v):
        v = Rational(v); return Fr(int(v.p), int(v.q))
    cM = Matrix([[toR(v)] for v in c])
    res = slp(cM, Matrix([[toR(v) for v in row] for row in A]), Matrix([toR(v) for v in b]),
              Matrix([[toR(v) for v in row] for row in Aeq]), Matrix([toR(v) for v in beq]))
    val = -toF(res[0])
    x = [toF(v) for v in res[1]]
    # primal feasibility
    okp = all(v >= 0 for v in x)
    okp &= all(sum(a * xx for a, xx in zip(row, x)) <= bb for row, bb in zip(A, b))
    okp &= all(sum(a * xx for a, xx in zip(row, x)) == bb for row, bb in zip(Aeq, beq))
    okp &= (-sum(cc * xx for cc, xx in zip(c, x)) == val)
    # dual: min b.y + beq.z  s.t. A^T y + Aeq^T z >= -c, y >= 0, z free  (z split z+ - z-)
    ny = len(A); nz = len(Aeq); nv = len(c)
    dc = list(b) + list(beq) + [-v for v in beq]
    dA = []; db = []
    for j in range(nv):
        dA.append([-A[i][j] for i in range(ny)] + [-Aeq[i][j] for i in range(nz)] + [Aeq[i][j] for i in range(nz)])
        db.append(c[j])
    res2 = slp(Matrix([[toR(v)] for v in dc]), Matrix([[toR(v) for v in row] for row in dA]),
               Matrix([toR(v) for v in db]))
    yz = [toF(v) for v in res2[1]]
    y = yz[:ny]; z = [yz[ny + i] - yz[ny + nz + i] for i in range(nz)]
    okd = all(v >= 0 for v in y)
    for j in range(nv):
        okd &= sum(A[i][j] * y[i] for i in range(ny)) + sum(Aeq[i][j] * z[i] for i in range(nz)) >= -c[j]
    dval = sum(bb * yy for bb, yy in zip(b, y)) + sum(bb * zz for bb, zz in zip(beq, z))
    return val, x, y, z, okp, okd, dval


def section_D():
    t = time.time()
    rows = []
    maxratio = -1.0; comb_gain = 0.0
    for p in [101, 1009, 10007]:
        d = (p - 1) // 2
        for m in [2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 25, 30]:
            for kap in [Fr(1, 20), Fr(1, 10), Fr(1, 4), Fr(1, 2), Fr(1)]:
                n = -((-(Fr(1, 2) + kap) * p) // m)   # ceil
                n = int(n)
                if m > (p + 1) // 2:
                    continue
                c, A, b, Aeq, beq = build_lp(p, m, n)
                val, x, y, z, okp, okd, dval = solve_lp_exact(c, A, b, Aeq, beq)
                check("D1_lp_certificate", okp and okd and dval == val, (p, m, n, okp, okd))
                mn = m * n
                w2 = Fr(d + min(m, n), mn)
                usew = max(w2, Fr(1, 4))
                if w2 < 1:
                    check("D1_lp_le_u", le_u_exact(val / mn, usew), (p, m, n, float(val / mn), float(usew)))
                    # LP <= single-constraint bounds (HP, Jensen) evaluated at w''
                    hb = min([hp_bound(m, usew)] + [jensen_hbound_hi(m, E, usew) for E in range(1, (m + 1) // 2 + 1)])
                    check("D1_lp_le_single_E", val / mn <= hb, (p, m, n))
                    comb_gain = max(comb_gain, float(hb - val / mn))
                    maxratio = max(maxratio, float(val / mn) / u_float(float(usew)))
                rused = sum(x[m + 1:])
                rows.append({"p": p, "m": m, "n": n, "kappa": str(kap), "LP_bias": float(val / mn),
                             "u(w'')": u_float(float(usew)) if w2 < 1 else 1.0,
                             "w''": float(w2), "r_used": float(rused),
                             "support": [(j, float(v)) for j, v in enumerate(x) if v != 0]})
    R["tables"]["D1_exact_lp"] = rows
    R["tables"]["D1_summary"] = {"max LP/u": maxratio, "max gain of full LP over best single E": comb_gain,
                                 "LPs": len(rows)}
    # D2: float LPs for large m (d = infinity, r = 0), convergence to u(w)
    from scipy.optimize import linprog
    d2 = []
    for w in [0.5, 0.8, 0.9, 0.95]:
        for m in [10, 25, 50, 100, 200, 400, 800]:
            J = np.arange(m + 1)
            cc = -(1 - 2 * J / m)
            AA = [np.where(J < E, (E - J) * (m + 1.5 - (3 * E + J) / 2), 0.0) for E in range(1, (m + 1) // 2 + 1)]
            bb = [E * m * w for E in range(1, (m + 1) // 2 + 1)]
            res = linprog(cc, A_ub=np.array(AA), b_ub=np.array(bb), A_eq=np.ones((1, m + 1)), b_eq=[1],
                          bounds=[(0, None)] * (m + 1), method="highs")
            d2.append({"w": w, "m": m, "LP_bias(float)": -res.fun, "u(w)": u_float(w),
                       "m*(u-LP)": m * (u_float(w) + res.fun), "single_E_bound": best_bound_float(m, w)})
    R["tables"]["D2_float_lp_large_m"] = d2
    # D3: explicit primal profiles (Proposition 2.1): all e_b = J, r = 0, n = floor(d/(m w))
    d3 = []
    for (w, m, dd) in [(Fr(1, 2), 200, 10 ** 7), (Fr(4, 5), 400, 10 ** 8), (Fr(9, 10), 1000, 10 ** 9),
                       (Fr(9, 10), 4000, 10 ** 10), (Fr(19, 20), 4000, 10 ** 10), (Fr(2, 3), 2000, 10 ** 9)]:
        wf = float(w); ts = (1 - u_float(wf)) / 2
        need = 1.5 / m + m / (2 * dd)
        # smallest delta with delta*(1/2-1/(2m)-ts-delta) >= need, then round J up
        a0 = 0.5 - 0.5 / m - ts
        delta = (a0 - math.sqrt(a0 * a0 - 4 * need)) / 2
        J = math.ceil((ts + delta) * m + 1e-9)
        n = int(Fr(dd) / (m * w))
        okall = True; cnt = 0
        for e in range(J, (m - 1) // 2 + 1):
            lhs = n * (e + 1 - J) * (2 * m - 3 * e - J)
            rhs = 2 * (e + 1) * (dd - e)
            cnt += 1
            if lhs > rhs:
                okall = False; fail("D3_profile", (float(w), m, dd, e)); break
        bump("D3_profile_constraints", cnt)
        d3.append({"w": str(w), "m": m, "d": dd, "J": J, "n": n, "profile_bias": 1 - 2 * J / m,
                   "u(w)": u_float(wf), "w_eff=d/(mn)": dd / (m * n), "all_star_ok": okall})
    R["tables"]["D3_primal_profiles"] = d3
    R["timing"]["D"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- E. exhaustive maxima, p <= 23
def exhaustive_maxS(p, chunk=1 << 17):
    """maxS[m][n] = max over |A|=m, |B|=n of S(A,B) (= max |S| by the non-residue dilation),
    all A containing 0 (translation invariance), B = top-n values of F_A."""
    chi = legtab(p)
    Chi = np.array([[chi[(a + b) % p] for b in range(p)] for a in range(p)], dtype=np.int16)
    N = 1 << (p - 1)
    maxS = np.full((p + 1, p + 1), -10 ** 6, dtype=np.int64)
    wit = {}
    for start in range(0, N, chunk):
        idx = np.arange(start, min(N, start + chunk), dtype=np.int64)
        X = np.zeros((len(idx), p), dtype=np.int16); X[:, 0] = 1
        for i in range(1, p):
            X[:, i] = (idx >> (i - 1)) & 1
        F = X @ Chi
        m = X.sum(1)
        T = np.cumsum(-np.sort(-F, axis=1), axis=1)
        for mm in np.unique(m):
            rows = np.nonzero(m == mm)[0]
            sub = T[rows]
            j = np.argmax(sub, axis=0)
            vals = sub[j, np.arange(p)]
            for n in range(1, p + 1):
                if vals[n - 1] > maxS[mm, n]:
                    maxS[mm, n] = int(vals[n - 1])
                    wit[(int(mm), n)] = int(idx[rows[j[n - 1]]])
    return maxS, wit


def mask_to_set(p, mask):
    return [0] + [i for i in range(1, p) if (mask >> (i - 1)) & 1]


def section_E():
    t = time.time()
    from scipy.optimize import linprog
    tables = {}
    kappas = [0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.35, 0.5, 0.75, 1.0]
    for p in [7, 11, 13, 17, 19, 23]:
        d = (p - 1) // 2
        maxS, wit = exhaustive_maxS(p)
        per = []
        for m in range(1, p + 1):
            for n in range(1, p + 1):
                S = int(maxS[m, n]); mn = m * n
                # symmetric: S(A,B)=S(B,A)
                check("E_symmetry", S == int(maxS[n, m]), (p, m, n))
                mm, nn = (m, n) if m <= n else (n, m)
                if mm > (p + 1) // 2:
                    continue
                w2 = Fr(d + mm, mn)            # r <= min(m,n)
                if w2 < 1:
                    usew = max(w2, Fr(1, 4))
                    check("E_truth_le_u", le_u_exact(Fr(S, mn), usew), (p, m, n, S))
                    hb = min([hp_bound(mm, usew)] + [jensen_hbound_hi(mm, E, usew) for E in range(1, (mm + 1) // 2 + 1)])
                    check("E_truth_le_singleE", Fr(S, mn) <= hb, (p, m, n, S))
                if mn > d and m <= n:
                    # float LP (r allowed), as a comparison
                    c, A, b, Aeq, beq = build_lp(p, mm, nn)
                    res = linprog([float(v) for v in c], A_ub=[[float(v) for v in row] for row in A],
                                  b_ub=[float(v) for v in b], A_eq=[[1.0] * len(c)], b_eq=[float(nn)],
                                  bounds=[(0, None)] * len(c), method="highs")
                    lpv = -res.fun
                    check("E_truth_le_LP(float)", S <= lpv + 1e-7, (p, m, n, S, lpv))
                    per.append((m, n, S, round(lpv, 4)))
        # maximal bias as a function of kappa
        kt = []
        for kap in kappas:
            best = (-2, None)
            for m in range(1, p + 1):
                for n in range(1, p + 1):
                    if m * n >= (0.5 + kap) * p:
                        bval = maxS[m, n] / (m * n)
                        if bval > best[0] + 1e-12:
                            best = (bval, (m, n, mask_to_set(p, wit[(m, n)])))
            wv = 1 / (1 + 2 * kap)
            kt.append({"kappa": kap, "max_bias": best[0], "argmax(m,n,A)": best[1],
                       "1-2k/(1+2k)=w": wv, "u(w)": u_float(max(wv, 0.25))})
        tables[p] = {"per_mn(m<=n, mn>d)": per, "kappa_table": kt}
    R["tables"]["E_exhaustive"] = tables
    R["timing"]["E"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- F. fixed |A| (Weil range)
def beta_m(m, q):
    """limit maximal bias for |A| = m fixed, |B| = q p (0 < q <= 1): the mean of 1-2j/m over
    the lowest-j mass q of Binomial(m,1/2)."""
    q = Fr(q); rem = q; tot = Fr(0)
    for j in range(m + 1):
        mass = Fr(math.comb(m, j), 2 ** m)
        take = min(mass, rem)
        if take <= 0:
            break
        tot += take * (1 - Fr(2 * j, m)); rem -= take
    return tot / q


def counts_by_e(p, A, chi):
    x = np.zeros(p, dtype=np.int64); x[A] = 1
    e = np.zeros(p, dtype=np.int64)
    for a in A:
        e += (chi[(a + np.arange(p)) % p] == -1)
    delta = x[[(-b) % p for b in range(p)]]
    return e, delta


def section_F():
    t = time.time()
    rng = random.Random(777)
    # F1: Weil error bars for N_j(A):  |N_j - C(m,j) p/2^m| <= C(m,j)(m/2)(sqrt p + 1) + 2m
    worst = 0.0
    for p in [10007, 100003]:
        chi = legtab(p)
        for trial in range(40 if p < 50000 else 12):
            m = rng.randint(1, 8)
            A = rng.sample(range(p), m)
            e, delta = counts_by_e(p, A, chi)
            for j in range(m + 1):
                Nj = int((e == j).sum())
                Cj = math.comb(m, j)
                diff = abs(Nj * 2 ** m - Cj * p)            # 2^m * |N_j - C p/2^m|
                # need diff <= 2^m [ C (m/2)(sqrt p + 1) + 2m ]
                X = Fr(diff, 2 ** m) - Fr(Cj * m, 2) - 2 * m  # need X <= C (m/2) sqrt p
                ok = X <= 0 or X * X <= Fr(Cj * m, 2) ** 2 * p
                check("F1_weil_counts", ok, (p, sorted(A), j, Nj))
                worst = max(worst, float(Fr(diff, 2 ** m)) / (Cj * m / 2 * (math.sqrt(p) + 1) + 2 * m))
    R["tables"]["F1_max_ratio_to_error_bar"] = worst
    # F2: beta_m and the refined upper bound on the saving
    rows = []
    for kap in [Fr(1, 100), Fr(1, 20), Fr(1, 10), Fr(1, 5), Fr(1, 4), Fr(3, 10), Fr(2, 5), Fr(1, 2), Fr(3, 4), Fr(1), Fr(3, 2)]:
        w = 1 / (1 + 2 * kap)
        best = (Fr(-2), None); vals = {}
        for m in range(1, 81):
            q = 1 / (2 * m * w)
            if q > 1:
                continue
            b = beta_m(m, q); vals[m] = float(b)
            if b > best[0]:
                best = (b, m)
        rows.append({"kappa": str(kap), "w": float(w), "sup_m beta_m": float(best[0]), "argmax m": best[1],
                     "beta_2": vals.get(2), "beta_3": vals.get(3), "beta_4": vals.get(4),
                     "upper bound on saving 1-sup beta": 1 - float(best[0]), "1-w": 1 - float(w),
                     "eta_star=1-u(w)": 1 - u_float(max(float(w), 0.25))})
        # beta_2 = w exactly when w >= 1/3 (q = 1/(4w) <= 3/4)
        if w >= Fr(1, 3):
            check("F2_beta2_equals_w", beta_m(2, 1 / (4 * w)) == w, float(w))
    R["tables"]["F2_refined_upper_bound"] = rows
    # F3: finite-p realisations: best B for random A of size m at p=100003 vs beta_m
    p = 100003; chi = legtab(p); f3 = []
    for m in [2, 3, 4, 6]:
        for w in [Fr(9, 10), Fr(1, 2)]:
            n = int(math.ceil(p / (2 * m * float(w))))
            if n > p:
                continue
            A = rng.sample(range(p), m)
            e, delta = counts_by_e(p, A, chi)
            F = m - 2 * e - delta
            S = int(np.sort(F)[::-1][:n].sum())
            f3.append({"m": m, "w": float(w), "n": n, "bias": S / (m * n), "beta_m": float(beta_m(m, Fr(n, p)))})
    R["tables"]["F3_realised"] = f3
    # F4: crossover.  For 3 <= m < 200 the function w -> beta_m(1/(2mw)) - w is piecewise linear in w
    # (breakpoints where q = 1/(2mw) equals a cumulative binomial mass), so checking its sign at all
    # breakpoints and endpoints of [24/35, 1] is exact.  m >= 200: Hoeffding gives beta_m <= 1/2 + 2m e^{-m/8}.
    wlo, whi = Fr(24, 35), Fr(1)
    cnt = 0; eq = []
    for m in range(3, 200):
        pts = {wlo, whi}
        cum = Fr(0)
        for j in range(m + 1):
            cum += Fr(math.comb(m, j), 2 ** m)
            wb = 1 / (2 * m * cum)
            if wlo < wb < whi:
                pts.add(wb)
        for w in pts:
            q = 1 / (2 * m * w)
            if q > 1:
                continue
            val = beta_m(m, q) - w
            cnt += 1
            if not check("F4_crossover_m>=3_below_w", val <= 0, (m, str(w), float(val))):
                pass
            if val == 0:
                eq.append((m, str(w)))
    bump("F4_breakpoints", cnt - 1)
    check("F4_m5_formula", all(beta_m(5, 1 / (10 * w)) == Fr(3, 5) + w / 8
                               for w in [Fr(8, 15), Fr(24, 35), Fr(7, 10), Fr(9, 10), Fr(1)]), "")
    check("F4_hoeffding_tail", 2 * 200 * math.exp(-200 / 8) < Fr(24, 35) - Fr(1, 2), "")
    R["tables"]["F4_crossover"] = {"w_c": "24/35", "kappa_c": "11/48", "equality_cases": eq}
    R["timing"]["F"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- G. balanced sets
def section_G():
    t = time.time()
    rng = random.Random(4242)
    # G1: exhaustive biclique numbers with min(|A|,|B|) >= k, p <= 23
    g1 = []
    for p in [11, 13, 17, 19, 23]:
        chi = legtab(p)
        N = 1 << (p - 1)
        Mm = np.zeros(p + 1, dtype=np.int64)       # M_m(p) = max_{|A|=m} N_0(A)
        NegOne = np.array([[1 if chi[(a + b) % p] == -1 else 0 for b in range(p)] for a in range(p)], dtype=np.int16)
        for start in range(0, N, 1 << 17):
            idx = np.arange(start, min(N, start + (1 << 17)), dtype=np.int64)
            X = np.zeros((len(idx), p), dtype=np.int16); X[:, 0] = 1
            for i in range(1, p):
                X[:, i] = (idx >> (i - 1)) & 1
            N0 = ((X @ NegOne) == 0).sum(1)
            m = X.sum(1)
            np.maximum.at(Mm, m, N0)
        row = {"p": p, "M_m(p)": [int(v) for v in Mm[1:]]}
        for k in range(1, 6):
            best = 0
            for m in range(k, p + 1):
                if Mm[m] >= k:
                    best = max(best, m * int(Mm[m]))
            row["bic_k=%d" % k] = best
            row["bic/p_k=%d" % k] = best / p
            row["k/2^k"] = None
        g1.append(row)
        # HP: m * M_m <= d + m  (check)
        d = (p - 1) // 2
        for m in range(1, (p + 1) // 2 + 1):
            check("G1_HP", m * int(Mm[m]) <= d + m, (p, m, int(Mm[m])))
    R["tables"]["G1_bicliques"] = g1
    # G2: fourth moment  sum_b F_A(b)^4 <= 3 m^2 p + 3 m^4 sqrt p ; and Hoelder consequence
    worst = 0.0
    for p in [1009, 10007]:
        chi = legtab(p)
        for trial in range(30):
            m = rng.randint(2, 60)
            A = rng.sample(range(p), m)
            F = np.zeros(p, dtype=np.int64)
            for a in A:
                F += chi[(a + np.arange(p)) % p]
            M4 = int((F ** 4).sum())
            X = M4 - 3 * m * m * p
            check("G2_fourth_moment", X <= 0 or X * X <= 9 * m ** 8 * p, (p, m))
            worst = max(worst, M4 / (3 * m * m * p + 3 * m ** 4 * math.sqrt(p)))
            for n in [m, 5 * m, int(12 * math.sqrt(p))]:
                if n > p:
                    continue
                B = np.argsort(-F)[:n]
                S = int(F[B].sum())
                bound = (3 * p / (m * m * n) + 3 * math.sqrt(p) / n) ** 0.25
                check("G2_holder", abs(S) <= bound * m * n + 1e-9, (p, m, n, S))
    R["tables"]["G2_max_M4_over_bound"] = worst
    R["timing"]["G"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- H. characters of order k
COS = {2: [Fr(1), Fr(-1)], 3: [Fr(1), Fr(-1, 2), Fr(-1, 2)], 4: [Fr(1), Fr(0), Fr(-1), Fr(0)],
       6: [Fr(1), Fr(1, 2), Fr(-1, 2), Fr(-1), Fr(-1, 2), Fr(1, 2)]}


def abs2_exact(Ncls, k):
    return sum(Ncls[i] * Ncls[j] * COS[k][(i - j) % k] for i in range(k) for j in range(k))


def theta_lo(w):
    """rational lower bound for theta = (1-u(max(w,1/4)))/2."""
    return (1 - u_hi(max(Fr(w), Fr(1, 4)))) / 2


def orderk_checks(p, k, A, B, dlog, key):
    """exact checks of section 5 for one pair (A,B) and the order-k character."""
    dk = (p - 1) // k; m = len(A); n = len(B)
    if not (m <= dk + 1 and m + dk - 1 <= p - 1):
        return
    Ncls = [0] * k; r = 0
    ebad = [[0] * n for _ in range(k)]
    for bi, b in enumerate(B):
        for a in A:
            x = (a + b) % p
            if x == 0:
                r += 1; continue
            c = dlog[x] % k
            Ncls[c] += 1
            for i in range(k):
                if i != c:
                    ebad[i][bi] += 1
    for i in range(k):
        Nb = sum(ebad[i])
        # Lemma 1.1 for the dilated set (per value i)
        for E in range(1, (m + 1) // 2 + 1):
            if Nb < E * n:
                check(key + "_jensen", (E * n - Nb) * (n * (2 * m + 3 - 3 * E) - Nb) <= 2 * n * E * (dk - E + 1 + r),
                      (p, k, A, B, i, E))
    wk = Fr(dk + r, m * n)
    if wk < 1:
        th = theta_lo(wk)
        M = m * n - r
        for i in range(k):
            check(key + "_Ni", Ncls[i] <= (1 - th) * M, (p, k, A, B, i))
        if k in COS:
            ck = COS[k][1]
            check(key + "_final", abs2_exact(Ncls, k) <= M * M * (1 - 2 * th * (1 - th) * (1 - ck)), (p, k, A, B))


def section_H():
    t = time.time()
    rng = random.Random(99)
    # H1/H2: exhaustive A containing 0 (small p), several B per A
    for (p, klist, exh) in [(13, [3, 4, 6], True), (19, [3, 6], True), (31, [3, 6], False), (37, [3, 4, 6], False),
                           (61, [3, 4, 6], False)]:
        dlog = dlog_table(p)
        for k in klist:
            dk = (p - 1) // k
            mmax = min(dk + 1, p - dk)
            if exh:
                sets = []
                for mm in range(1, mmax + 1):
                    for rest in itertools.combinations(range(1, p), mm - 1):
                        sets.append((0,) + rest)
                if len(sets) > 6000:
                    sets = rng.sample(sets, 6000)
            else:
                sets = [tuple(rng.sample(range(p), rng.randint(1, mmax))) for _ in range(400)]
            zeta = np.exp(2j * np.pi * np.arange(k) / k)
            for A in sets:
                Fb = np.zeros(p, dtype=complex)
                for a in A:
                    for b in range(p):
                        x = (a + b) % p
                        if x:
                            Fb[b] += zeta[dlog[x] % k]
                Bs = []
                for i in range(k):
                    score = (np.conj(zeta[i]) * Fb).real
                    order = np.argsort(-score, kind="stable")
                    for n in {max(1, (dk + 1) // len(A) + 1), 2 * dk // len(A) + 1, rng.randint(1, p)}:
                        if n <= p:
                            Bs.append([int(b) for b in order[:n]])
                Bs.append(rng.sample(range(p), rng.randint(1, p)))
                for B in Bs:
                    orderk_checks(p, k, list(A), B, dlog, "H")
    # H3: sharpness example A={0}, B = H_k plus x points of g H_k
    h3 = []
    for (p, k) in [(1009, 3), (1009, 4), (1009, 6), (10009, 3), (10009, 4)]:
        dk = (p - 1) // k; dlog = dlog_table(p)
        Hk = [x for x in range(1, p) if dlog[x] % k == 0]
        gH = [x for x in range(1, p) if dlog[x] % k == 1]
        check("H3_subgroup_size", len(Hk) == dk, (p, k))
        for lam in [Fr(0), Fr(1, 10), Fr(1, 2)]:
            x = int(lam * dk)
            B = Hk + gH[:x]
            Ncls = [0] * k
            for b in B:
                Ncls[dlog[b] % k] += 1
            mn = len(B)
            if k in COS:
                a2 = abs2_exact(Ncls, k)
                check("H3_value", a2 == dk * dk + x * x + 2 * dk * x * COS[k][1], (p, k, x))
                h3.append({"p": p, "k": k, "lambda": str(lam), "mn/d_k": mn / dk, "bias": math.sqrt(a2) / mn})
    R["tables"]["H3_sharpness"] = h3
    # H4: vertex lemma by brute force: max |sum_i zeta^i N_i|^2 over integer N >= 0, sum N = M, N_i <= (1-th) M
    for k in [3, 4, 6]:
        for M in [6, 10, 12]:
            for th in [Fr(1, 6), Fr(1, 4), Fr(1, 3), Fr(1, 2)]:
                cap = (1 - th) * M
                best = Fr(0)
                for comp in itertools.product(range(M + 1), repeat=k - 1):
                    s = sum(comp)
                    if s > M:
                        continue
                    Nv = list(comp) + [M - s]
                    if max(Nv) > cap:
                        continue
                    best = max(best, abs2_exact(Nv, k))
                check("H4_vertex", best <= M * M * (1 - 2 * th * (1 - th) * (1 - COS[k][1])), (k, M, th))
    R["timing"]["H"] = round(time.time() - t, 1)


# ---------------------------------------------------------------- main
def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, Fr):
        return str(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    return o


def main():
    cpu = {}
    for name, fn in [("A", section_A), ("B", section_B), ("C", section_C), ("D", section_D), ("E", section_E),
                     ("F", section_F), ("G", section_G), ("H", section_H)]:
        c0 = time.process_time()
        fn()
        cpu[name] = round(time.process_time() - c0, 1)
        print("section", name, "done; cpu", cpu[name], "s; failures so far", len(FAIL), flush=True)
    R["timing"]["cpu_seconds_by_section"] = cpu
    R["timing"]["cpu_seconds_total"] = round(time.process_time(), 1)
    R["timing"]["wall_seconds_total"] = round(time.time() - T0, 1)
    R["total_checks"] = sum(R["checks"].values())
    R["failures"] = len(FAIL)
    R["notes"]["u(w)"] = "u(w) = (sqrt(12w-3w^2)-w)/2; eta_star(kappa) = 1-u(1/(1+2kappa)) = 4k^2/3-40k^3/9+O(k^4)"
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(jsonable(R), f, indent=1)
    print("total checks", R["total_checks"], "failures", len(FAIL))
    print("cpu", R["timing"]["cpu_seconds_total"], "wall", R["timing"]["wall_seconds_total"])
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
