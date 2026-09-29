#!/usr/bin/env python3
"""Verifier for research/stepanov-hybrid-2026-09-29.md (Stepanov wave, worker `hybrid`).

Hybrid of the Hanson-Petridis/Stepanov count constraints (balanced note, Theorem 3.5) with
positive-semidefinite (spectral / Gram / moment) constraints.

Run:  /opt/miniconda3/bin/python3 experiments/stepanov_hybrid_2026_09_29.py
Writes results/stepanov_hybrid_2026_09_29.json.

Sections
  A  Paley matrix identities; validity and standalone value of the spectral Gram family (P1)
  B  level-set compressions of the Paley matrix on actual configurations (P3 identities)
  C  moment / Gram matrices of lifted sign vectors (P4): exact entries, Weil entries,
     symmetrisation identity, Krawtchouk-Weil inequalities on actual configurations
  D  the decoupling construction (Theorem 2.1): explicit T for arbitrary profiles
  E  hybrid LP certificates (balanced rows + Krawtchouk-Weil rows + P1 + P3), exact, at
     p = 1009, 10009, 40009, 100049, 1000033 (bicliques and cliques)
  F  independent numerical SDP solve of the level-set relaxation (cvxpy; Clarabel, SCS fallback)
  G  two-point (Delsarte-type) program P5: SDP, rational rounding, exact repair, exact LDL^T PSD checks
  G2 P5 pure-feasibility solves on two extreme one-point profiles (numerical)
  H  bipartite theta program, symmetrised (numerical)
Every exact check increments a counter; failures are recorded with witnesses.
Runtime about 6 minutes (5 of them in the P5 SDPs at p = 10009, 40009).
Requires numpy, scipy and cvxpy with the Clarabel and SCS solvers.
"""
import json
import math
import os
import sys
import time
import random
from fractions import Fraction as Fr
from math import comb, isqrt

import numpy as np

T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__))) if '__file__' in globals() else os.getcwd()
OUT = os.path.join(ROOT, 'results', 'stepanov_hybrid_2026_09_29.json')

CHECKS = {}
FAILS = []
RESULTS = {}


def check(name, cond, info=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not cond:
        FAILS.append({'check': name, 'info': repr(info)[:400]})
    return cond


# ------------------------------------------------------------------ basic arithmetic

def chi_table(p):
    t = [0] * p
    for x in range(1, p):
        t[x * x % p] = 1
    return [0] + [1 if t[x] else -1 for x in range(1, p)]


def is_prime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def krawtchouk_table(M):
    """K[t][j] = sum_i (-1)^i C(j,i) C(M-j,t-i) for 0<=t,j<=M (exact ints), via the
    three-term recurrence (t+1)K_{t+1}(j) = (M-2j)K_t(j) - (M-t+1)K_{t-1}(j)."""
    K = [[0] * (M + 1) for _ in range(M + 1)]
    for j in range(M + 1):
        K[0][j] = 1
        if M >= 1:
            K[1][j] = M - 2 * j
    for t in range(1, M):
        for j in range(M + 1):
            num = (M - 2 * j) * K[t][j] - (M - t + 1) * K[t - 1][j]
            q, r = divmod(num, t + 1)
            assert r == 0
            K[t + 1][j] = q
    return K


def sqrt_sign(X, Y, p):
    """sign of X + Y*sqrt(p) for rationals X, Y (exact)."""
    if X >= 0 and Y >= 0:
        return (X > 0 or Y > 0) and 1 or 0
    if X <= 0 and Y <= 0:
        return -1 if (X < 0 or Y < 0) else 0
    # opposite signs: compare X^2 with p Y^2
    a, b = X * X, p * Y * Y
    if X > 0:
        return 1 if a > b else (0 if a == b else -1)
    return -1 if a > b else (0 if a == b else 1)


def psd2_sqrt(a0, a1, b0, b1, c0, c1, p):
    """Is [[a0+a1 r, b0+b1 r],[b0+b1 r, c0+c1 r]] PSD, r = sqrt(p)?  (exact)"""
    if sqrt_sign(a0, a1, p) < 0 or sqrt_sign(c0, c1, p) < 0:
        return False
    # det = (a0+a1r)(c0+c1r) - (b0+b1r)^2 = X + Y r
    X = a0 * c0 + a1 * c1 * p - b0 * b0 - b1 * b1 * p
    Y = a0 * c1 + a1 * c0 - 2 * b0 * b1
    return sqrt_sign(X, Y, p) >= 0


# ------------------------------------------------------------------ the count LP (balanced note, Thm 3.5)
# Variables: x[j] = n_j (0<=j<=m): points x not in -A with e_x = j;
#            x[m+1+j] = n'_j (0<=j<=m-1): points of -A with e_x = j.
# f_A(x) = m - delta_x - 2 e_x.

def val(m, j, dl):
    return m - dl - 2 * j


def balanced_rows(p, m):
    """Re-implementation of the theorem rows of the balanced note (lp_rows, without target/spread).
    Returns list of (coeff dict, sense, rhs, name)."""
    d = (p - 1) // 2
    I = lambda j: j
    J = lambda j: m + 1 + j
    rows = []
    rows.append(({**{I(j): 1 for j in range(m + 1)}, **{J(j): 1 for j in range(m)}}, '=', p, 'total'))
    rows.append(({J(j): 1 for j in range(m)}, '=', m, 'negA'))
    rows.append(({**{I(j): val(m, j, 0) for j in range(m + 1)}, **{J(j): val(m, j, 1) for j in range(m)}}, '=', 0, 'moment1'))
    rows.append(({**{I(j): val(m, j, 0) ** 2 for j in range(m + 1)}, **{J(j): val(m, j, 1) ** 2 for j in range(m)}},
                 '=', m * (p - m), 'moment2'))
    for side in ('A', 'nuA'):
        ee0 = (lambda j: j) if side == 'A' else (lambda j: m - j)
        ee1 = (lambda j: j) if side == 'A' else (lambda j: m - 1 - j)
        for e in range(0, (m - 1) // 2 + 1):
            co = {}
            for j in range(m + 1):
                ee = ee0(j)
                if ee <= e:
                    co[I(j)] = co.get(I(j), 0) + Fr((e + 1 - ee) * (2 * m - (3 * e + ee)), 2)
            for j in range(m):
                ee = ee1(j)
                if ee <= e:
                    co[J(j)] = co.get(J(j), 0) + Fr((e + 1 - ee) * (2 * m - (3 * e + ee) - 2), 2)
            rows.append((co, '<=', (e + 1) * (d - e), 'star_%s_%d' % (side, e)))
        for t in range(1, m + 1):
            co = {}
            for j in range(m + 1):
                ee = ee0(j)
                if m - 1 - ee >= t - 1:
                    co[I(j)] = co.get(I(j), 0) + comb(m - 1 - ee, t - 1) * (m - ee)
            for j in range(m):
                ee = ee1(j)
                if m - 1 - ee >= t - 1 and m - ee - 1 > 0:
                    co[J(j)] = co.get(J(j), 0) + comb(m - 1 - ee, t - 1) * (m - ee - 1)
            rows.append((co, '<=', comb(m, t) * d, 'subsetHP_%s_%d' % (side, t)))
    rows.append(({I(0): 2 * m - 2, I(1): m - 2}, '<=', 2 * d - 2, 'pencil_A'))
    rows.append(({I(m): 2 * m - 2, I(m - 1): m - 2}, '<=', 2 * d - 2, 'pencil_nuA'))
    rows.append(({I(0): m - 1, I(m): m - 1, J(0): m - 2, J(m - 1): m - 2}, '<=', d - 1, 'derivHP'))
    sq = isqrt(p)
    for k in (2, 3, 4):
        df = 1
        for i in range(1, 2 * k, 2):
            df *= i
        rhs = df * m ** k * p + (2 * k - 1) * m ** (2 * k) * sq
        rows.append(({**{I(j): val(m, j, 0) ** (2 * k) for j in range(m + 1)},
                      **{J(j): val(m, j, 1) ** (2 * k) for j in range(m)}}, '<=', rhs, 'weil_moment_%d' % (2 * k)))
    for j in range(m + 1):
        co = {I(j): 1}
        if j < m:
            co[J(j)] = 1
        main = Fr(comb(m, j) * p, 2 ** m)
        err = Fr(comb(m, j) * m * (sq + 1), 2) + 2 * m
        rows.append((co, '<=', main + err, 'weil_count_up_%d' % j))
        rows.append(({k: -v for k, v in co.items()}, '<=', -(main - err), 'weil_count_lo_%d' % j))
    return rows


KCACHE = {}


def ktab(M):
    if M not in KCACHE:
        KCACHE[M] = krawtchouk_table(M)
    return KCACHE[M]


def krawtchouk_weil_rows(p, m):
    """P4 aggregated (Prop. 2.3): for 3 <= t <= m,
       | sum_{|U|=t} sum_x prod_{a in U} chi(x+a) | <= C(m,t) (t-1) sqrt(p),
    and the left side is  sum_j n_j K_t(j;m) + sum_j n'_j K_t(j;m-1)."""
    I = lambda j: j
    J = lambda j: m + 1 + j
    Km, Km1 = ktab(m), ktab(m - 1)
    sq = isqrt(p)
    rows = []
    for t in range(3, m + 1):
        co = {I(j): Km[t][j] for j in range(m + 1)}
        for j in range(m):
            co[J(j)] = Km1[t][j] if t <= m - 1 else 0
        co = {k: v for k, v in co.items() if v != 0}
        bound = comb(m, t) * (t - 1) * sq
        rows.append((co, '<=', bound, 'KW_up_%d' % t))
        rows.append(({k: -v for k, v in co.items()}, '<=', bound, 'KW_lo_%d' % t))
    return rows


def exact_eval_rows(rows, xr):
    """Exact verification with a common denominator (integer arithmetic).
    Returns (all_ok, list of (name, lhs, rhs, ok)) with lhs/rhs as Fractions."""
    D = 1
    for v in xr:
        D = D * v.denominator // math.gcd(D, v.denominator)
    X = [int(v * D) for v in xr]
    out = []
    ok_all = True
    for co, sense, rhs, name in rows:
        # coefficients may be Fractions (star rows have halves)
        cden = 1
        for c in co.values():
            if isinstance(c, Fr):
                cden = cden * c.denominator // math.gcd(cden, c.denominator)
        tot = 0
        for k, c in co.items():
            tot += int(c * cden) * X[k]
        lhs = Fr(tot, D * cden)
        r = Fr(rhs)
        ok = (lhs == r) if sense == '=' else (lhs <= r)
        ok_all &= ok
        out.append((name, lhs, r, ok))
    return ok_all, out


def solve_certificate(p, m, mode, target, band=0.4, extra_rows=True, Q=1000):
    """mode 'biclique': n_0 = target, n'_0 = 0 (r = 0).
       mode 'clique'  : n'_0 = m (all of -A' at level 0), n'_j = 0 (j>=1), n_0 = target (0).
    Rows: balanced theorem rows (+ Krawtchouk-Weil rows if extra_rows) + band restriction.
    Returns dict with exact certificate (list of Fractions) or None."""
    from scipy.optimize import linprog
    nv = 2 * m + 1
    thm = balanced_rows(p, m)
    if extra_rows:
        thm = thm + krawtchouk_weil_rows(p, m)
    fixed = {}
    if mode == 'biclique':
        fixed[0] = Fr(target)
        fixed[m + 1] = Fr(0)
    else:
        fixed[0] = Fr(target)
        fixed[m + 1] = Fr(m)
        for j in range(1, m):
            fixed[m + 1 + j] = Fr(0)
    lo, hi = (math.ceil(band * m), math.floor((1 - band) * m)) if band is not None else (1, m)
    for j in range(1, m + 1):
        if j < lo or j > hi:
            fixed[j] = Fr(0)
    if mode == 'biclique':
        for j in range(m):
            if j < lo or j > hi:
                fixed[m + 1 + j] = Fr(0)
    SLACKFREE = ('subsetHP_A_1', 'subsetHP_A_2', 'subsetHP_nuA_1', 'subsetHP_nuA_2')
    Aeq, beq, Aub, bub = [], [], [], []
    for co, sense, rhs, name in thm:
        v = np.zeros(nv + 1)
        sc = max([abs(c) for c in co.values()] + [abs(rhs), 1])
        const = all(k in fixed for k in co)
        for k, c in co.items():
            v[k] = ratio_float(c, sc)
        if sense == '=':
            Aeq.append(v); beq.append(ratio_float(rhs, sc))
        else:
            if name not in SLACKFREE and not const:
                v[nv] = 1.0
            Aub.append(v); bub.append(ratio_float(rhs, sc))
    for k, fv in fixed.items():
        v = np.zeros(nv + 1); v[k] = 1.0
        Aeq.append(v); beq.append(float(fv))
    c = np.zeros(nv + 1); c[nv] = -1
    res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * nv + [(0, 1)], method='highs')
    if res.status != 0:
        return {'status': 'infeasible_or_error', 'msg': res.message}
    x = res.x[:nv]
    xr = [Fr(int(round(v * Q)), Q) if v > 1e-7 else Fr(0) for v in x]
    for k, fv in fixed.items():
        xr[k] = fv
    # exact repair of equalities
    if mode == 'biclique':
        jp = max([j for j in range(m) if (m + 1 + j) not in fixed], key=lambda j: x[m + 1 + j])
        xr[m + 1 + jp] = Fr(m) - sum(xr[m + 1 + j] for j in range(m) if j != jp)
    free = [j for j in range(1, m + 1) if j not in fixed]
    piv = sorted(free, key=lambda j: -x[j])[:3]
    rest = [j for j in range(m + 1) if j not in piv]
    r1 = p - sum(xr[m + 1 + j] for j in range(m)) - sum(xr[j] for j in rest)
    r2 = 0 - sum(xr[m + 1 + j] * val(m, j, 1) for j in range(m)) - sum(xr[j] * val(m, j, 0) for j in rest)
    r3 = m * (p - m) - sum(xr[m + 1 + j] * val(m, j, 1) ** 2 for j in range(m)) - sum(xr[j] * val(m, j, 0) ** 2 for j in rest)
    V = [[Fr(1)] * 3, [Fr(val(m, j, 0)) for j in piv], [Fr(val(m, j, 0) ** 2) for j in piv]]
    sol = solve3(V, [r1, r2, r3])
    for k, j in enumerate(piv):
        xr[j] = sol[k]
    nonneg = all(v >= 0 for v in xr)
    rows_all = thm + [({k: 1}, '=', fv, 'fixed_%d' % k) for k, fv in fixed.items()]
    ok, out = exact_eval_rows(rows_all, xr)
    n_ineq = sum(1 for r in rows_all if r[1] == '<=')
    n_eq = len(rows_all) - n_ineq
    return {'status': 'ok', 'x': xr, 'nonneg': nonneg, 'exact_ok': ok and nonneg, 'rows': out,
            'n_rows': len(rows_all), 'n_ineq': n_ineq, 'n_eq': n_eq, 'lp_slack': float(res.x[nv]),
            'band': [lo, hi]}


def ratio_float(a, b):
    """correctly rounded float(a/b) for ints or Fractions (no overflow for huge ints)"""
    if isinstance(a, int) and isinstance(b, int):
        return a / b
    q = Fr(a) / Fr(b)
    return q.numerator / q.denominator


def solve3(V, r):
    """exact solve of a 3x3 rational system by Cramer's rule"""
    def det3(M):
        return (M[0][0] * (M[1][1] * M[2][2] - M[1][2] * M[2][1]) - M[0][1] * (M[1][0] * M[2][2] - M[1][2] * M[2][0])
                + M[0][2] * (M[1][0] * M[2][1] - M[1][1] * M[2][0]))
    D = det3(V)
    out = []
    for c in range(3):
        M = [row[:] for row in V]
        for i in range(3):
            M[i][c] = r[i]
        out.append(det3(M) / D)
    return out


# ------------------------------------------------------------------ exact PSD test (rational LDL^T)

def psd_exact(M):
    """Exact PSD test of a symmetric matrix of Fractions/ints by symmetric Gaussian elimination.
    Returns (is_psd, rank, min_pivot)."""
    n = len(M)
    A = [[Fr(M[i][j]) for j in range(n)] for i in range(n)]
    rank = 0
    minpiv = None
    for k in range(n):
        piv = A[k][k]
        if piv < 0:
            return False, rank, piv
        if piv == 0:
            if any(A[k][j] != 0 for j in range(k + 1, n)):
                return False, rank, piv
            continue
        rank += 1
        minpiv = piv if minpiv is None else min(minpiv, piv)
        for i in range(k + 1, n):
            if A[i][k] == 0:
                continue
            f = A[i][k] / piv
            for j in range(k + 1, n):
                A[i][j] -= f * A[k][j]
            A[i][k] = Fr(0)
        for j in range(k + 1, n):
            A[k][j] = Fr(0)
    return True, rank, minpiv


# ------------------------------------------------------------------ Section A: Paley matrices and P1

def paley_S(p, chi):
    return np.array([[chi[(x + y) % p] for y in range(p)] for x in range(p)], dtype=np.int64)


def complete_partners(A, p, chi):
    return [x for x in range(p) if all(chi[(x + a) % p] >= 0 for a in A)]


def p1_data(A, B, p, chi):
    Aset, Bset = set(A), set(B)
    m, n = len(A), len(B)
    sA = sum(chi[(a + a2) % p] for a in A for a2 in A)
    sB = sum(chi[(b + b2) % p] for b in B for b2 in B)
    tau = sum(chi[(a + b) % p] for a in A for b in B)
    s = len(Aset & Bset)
    return m, n, s, sA, sB, tau


def p1_psd(m, n, s, sA, sB, tau, p):
    """G_+- = [[m-m^2/p +- sA/sqrt p, s-mn/p +- tau/sqrt p],[., n-n^2/p +- sB/sqrt p]] >= 0 (exact).
    x/sqrt(p) = (x/p) sqrt(p)."""
    ok = True
    for sg in (1, -1):
        ok &= psd2_sqrt(Fr(m) - Fr(m * m, p), Fr(sg * sA, p), Fr(s) - Fr(m * n, p), Fr(sg * tau, p),
                        Fr(n) - Fr(n * n, p), Fr(sg * sB, p), p)
    return ok


def section_A():
    rng = random.Random(20260929)
    out = {'primes': [], 'p1_bicliques_checked': 0, 'p1_chung_bound_checked': 0, 'clique_theta_checked': 0}
    primes = [p for p in range(5, 140) if is_prime(p)]
    for p in primes:
        chi = chi_table(p)
        S = paley_S(p, chi)
        J = np.ones((p, p), dtype=np.int64)
        check('A.S1=0', not np.any(S.sum(axis=1)), p)
        check('A.S^2=pI-J', np.array_equal(S @ S, p * np.eye(p, dtype=np.int64) - J), p)
        check('A.trS=0', int(np.trace(S)) == 0, p)
        check('A.S_symmetric', np.array_equal(S, S.T), p)
        if p % 4 == 1:
            M = np.array([[chi[(x - y) % p] for y in range(p)] for x in range(p)], dtype=np.int64)
            check('A.M_symmetric', np.array_equal(M, M.T), p)
            check('A.M^2=pI-J', np.array_equal(M @ M, p * np.eye(p, dtype=np.int64) - J), p)
            # eigenvalues +-sqrt(p), each with multiplicity d (trace 0, S^2 = p on 1^perp)
        out['primes'].append(p)
        # P1 on actual bicliques: random A, B = complete partners (and random subsets of it)
        for trial in range(25 if p < 60 else 8):
            m = rng.randint(1, 4)
            A = rng.sample(range(p), m)
            Bfull = complete_partners(A, p, chi)
            if not Bfull:
                continue
            for B in (Bfull, rng.sample(Bfull, rng.randint(1, len(Bfull)))):
                mm, nn, s, sA, sB, tau = p1_data(A, B, p, chi)
                check('A.P1_psd_actual', p1_psd(mm, nn, s, sA, sB, tau, p), (p, A, B))
                out['p1_bicliques_checked'] += 1
                r = len(set(A) & set((-b) % p for b in B))
                check('A.biclique_tau=mn-r', tau == mm * nn - r, (p, A, B))
                # Chung form (Prop 1.2a): (mn - r)^2 p <= mn(p-m)(p-n)
                check('A.P1_chung_bound', (mm * nn - r) ** 2 * p <= mm * nn * (p - mm) * (p - nn), (p, A, B))
                out['p1_chung_bound_checked'] += 1
        if p % 4 == 1:
            # greedy random maximal cliques containing 0: P1 in difference form gives m <= sqrt(p)
            for trial in range(20):
                C = [0]
                cand = [x for x in range(1, p) if chi[x] == 1]
                rng.shuffle(cand)
                for x in cand:
                    if all(chi[(x - c) % p] == 1 for c in C):
                        C.append(x)
                m = len(C)
                ok = psd2_sqrt(Fr(m) - Fr(m * m, p), Fr(m * (m - 1), p), Fr(0), Fr(0), Fr(1), Fr(0), p) and \
                    psd2_sqrt(Fr(m) - Fr(m * m, p), Fr(-m * (m - 1), p), Fr(0), Fr(0), Fr(1), Fr(0), p)
                check('A.P1_clique_psd', ok, (p, C))
                check('A.clique_m^2<=p', m * m <= p, (p, C))
                out['clique_theta_checked'] += 1
    # standalone value of P1 for bicliques: the optimum tau* = sqrt(p m' n') is attained with
    # sigma_A = sigma_B = 0, s = 0 (exact check of feasibility at tau = floor(sqrt(p m'n')) and
    # infeasibility of every sigma pair at tau = ceil(...)+1 on a grid is not needed: the bound is proved).
    tab = []
    for p in (1009, 10009, 100049, 1000033):
        d = (p - 1) // 2
        # largest m (= n) with m^2 - m <= sqrt(p) m (1 - m/p)  [r = m worst case] and with r = 0
        def p1_ok(m, r):
            x = m * m - r
            return x >= 0 and x * x * p <= (m * (p - m)) ** 2
        mm0 = max(m for m in range(1, isqrt(p) + 3) if p1_ok(m, 0))
        mhp = max(m for m in range(1, isqrt(p) + 3) if m * m <= d + m)
        tab.append({'p': p, 'P1_max_balanced_m_r0': mm0, 'P1_m2_over_p': mm0 * mm0 / p,
                    'HP_max_balanced_m': mhp, 'HP_m2_over_p': mhp * mhp / p})
        # the supremum tau* = sqrt(p) m' is approached with sigma_A = sigma_B = -m^2/sqrt(p), s = 0:
        # exact check that (m, n) = (mm0 - 1, mm0 - 1), r = 0 is P1-feasible with integer sigma = round(-m^2/sqrt p)
        m = mm0 - 1
        sg = -round(m * m / math.sqrt(p))
        check('A.P1_value_attained', p1_psd(m, m, 0, sg, sg, m * m, p), (p, m, sg))
    out['P1_standalone_table'] = tab
    RESULTS['A'] = out


# ------------------------------------------------------------------ Section B: level-set compressions (P3) on real data

def level_of(x, A, p, chi):
    e = sum(1 for a in A if chi[(x + a) % p] == -1)
    dl = 1 if any((x + a) % p == 0 for a in A) else 0
    return (e, dl)


def section_B():
    rng = random.Random(7)
    out = {'configs': 0, 'min_eig_rel': 1.0}
    for p in [29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 13, 17, 19, 23, 31, 43, 47, 59, 67, 71, 79, 83]:
        chi = chi_table(p)
        S = paley_S(p, chi)
        for trial in range(6):
            m = rng.randint(2, 4)
            A = rng.sample(range(p), m)
            B = complete_partners(A, p, chi)
            if not B:
                continue
            n = len(B)
            labA = [level_of(x, A, p, chi) for x in range(p)]
            labB = [level_of(x, B, p, chi) for x in range(p)]
            joint = sorted(set(zip(labA, labB)))
            idx = {c: i for i, c in enumerate(joint)}
            L = len(joint)
            X = np.zeros((p, L), dtype=np.int64)
            for x in range(p):
                X[x, idx[(labA[x], labB[x])]] = 1
            T = X.T @ S @ X
            nu = X.sum(axis=0)
            fA = np.array([m - c[0][1] - 2 * c[0][0] for c in joint], dtype=np.int64)
            fB = np.array([n - c[1][1] - 2 * c[1][0] for c in joint], dtype=np.int64)
            aa = np.array([sum(1 for x in A if (labA[x], labB[x]) == c) for c in joint], dtype=np.int64)
            bb = np.array([sum(1 for x in B if (labA[x], labB[x]) == c) for c in joint], dtype=np.int64)
            check('B.T1_rowsum0', not np.any(T.sum(axis=1)), (p, A))
            check('B.T2_A', np.array_equal(T @ fA, p * aa - m * nu), (p, A))
            check('B.T2_B', np.array_equal(T @ fB, p * bb - n * nu), (p, A))
            check('B.moment2_A', int((nu * fA * fA).sum()) == m * (p - m), (p, A))
            check('B.joint_moment', int((nu * fA * fB).sum()) == p * len(set(A) & set(B)) - m * n, (p, A))
            # 1_A, 1_B in the span of the joint cells (B = B(A) maximal, A subset of B(B))
            check('B.A_subset_B(B)', all(chi[(a + b) % p] >= 0 for a in A for b in B), (p, A))
            # Gram under Pi_+- of the family {1_cells, 1_A, 1_B} is PSD (numerical)
            V = np.column_stack([X, np.isin(np.arange(p), A).astype(np.int64), np.isin(np.arange(p), B).astype(np.int64)]).astype(float)
            GI = V.T @ V - np.outer(V.sum(axis=0), V.sum(axis=0)) / p
            GS = V.T @ S.astype(float) @ V
            for sg in (1, -1):
                ev = np.linalg.eigvalsh((GI + sg * GS / math.sqrt(p)) / 2)
                rel = ev.min() / max(1.0, abs(ev).max())
                out['min_eig_rel'] = min(out['min_eig_rel'], rel)
                check('B.Gpm_psd_numeric', rel > -1e-9, (p, A, sg, rel))
            out['configs'] += 1
    RESULTS['B'] = out


# ------------------------------------------------------------------ Section C: P4 (moment matrices of lifted sign vectors)

def esym(vals, t):
    """elementary symmetric polynomial e_t of a list of ints"""
    e = [1] + [0] * t
    for v in vals:
        for k in range(t, 0, -1):
            e[k] += e[k - 1] * v
    return e[t]


def section_C():
    from itertools import combinations, permutations, product
    rng = random.Random(11)
    out = {'krawtchouk_tables_checked': 0, 'Z2_entries_checked': 0, 'symmetrisation_checked': 0,
           'KW_actual_checked': 0, 'KW_ratio_max': 0.0}
    # Krawtchouk table vs direct definition
    for M in range(1, 16):
        K = ktab(M)
        for t in range(M + 1):
            for j in range(M + 1):
                direct = sum((-1) ** i * comb(j, i) * comb(M - j, t - i) for i in range(0, t + 1))
                check('C.krawtchouk_table', K[t][j] == direct, (M, t, j))
        out['krawtchouk_tables_checked'] += 1
    for p in [29, 37, 41, 53, 61, 101, 109, 113, 23, 31, 43, 47, 59, 67]:
        chi = chi_table(p)
        for trial in range(4):
            m = rng.randint(3, 4)
            A = rng.sample(range(p), m)
            idxA = list(range(m))
            Ts = [()] + [(i,) for i in idxA] + list(combinations(idxA, 2))
            svec = [[chi[(x + a) % p] for a in A] for x in range(p)]
            def psi(sv):
                return [math.prod(sv[i] for i in T) for T in Ts]
            Z = np.zeros((len(Ts), len(Ts)), dtype=np.int64)
            for x in range(p):
                v = np.array(psi(svec[x]), dtype=np.int64)
                Z += np.outer(v, v)
            # entry formulas
            def S_U(U):
                return sum(math.prod(chi[(x + A[i]) % p] for i in U) for x in range(p))
            for a_, T in enumerate(Ts):
                for b_, T2 in enumerate(Ts):
                    U = tuple(sorted(set(T) ^ set(T2)))
                    I = set(T) & set(T2)
                    corr = sum(math.prod(chi[(A[i] - A[c]) % p] for i in U) for c in I)
                    val_ = S_U(U) - corr
                    check('C.Z2_entry_formula', Z[a_, b_] == val_, (p, A, T, T2))
                    if len(U) == 0:
                        check('C.S_empty=p', S_U(U) == p, p)
                    elif len(U) == 1:
                        check('C.S_1=0', S_U(U) == 0, p)
                    elif len(U) == 2:
                        check('C.S_2=-1', S_U(U) == -1, p)
                    else:
                        check('C.S_U_weil', S_U(U) ** 2 <= (len(U) - 1) ** 2 * p, (p, A, U))
                    out['Z2_entries_checked'] += 1
            # symmetrisation identity: average over S_m of P Z P^T = sum over profile of uniform moments
            prof = {}
            for x in range(p):
                sv = svec[x]
                key = (sum(1 for v in sv if v == -1), 1 if 0 in sv else 0)
                prof[key] = prof.get(key, 0) + 1
            Zbar = [[Fr(0)] * len(Ts) for _ in Ts]
            perms = list(permutations(idxA))
            for pi in perms:
                for x in range(p):
                    sv = [svec[x][pi[i]] for i in idxA]
                    v = psi(sv)
                    for a_ in range(len(Ts)):
                        if v[a_] == 0:
                            continue
                        for b_ in range(len(Ts)):
                            Zbar[a_][b_] += Fr(v[a_] * v[b_], len(perms))
            Zmod = [[Fr(0)] * len(Ts) for _ in Ts]
            for (j, dl), cnt in prof.items():
                words = []
                if dl == 0:
                    for pos in combinations(idxA, j):
                        words.append([(-1 if i in pos else 1) for i in idxA])
                else:
                    for z in idxA:
                        rest = [i for i in idxA if i != z]
                        for pos in combinations(rest, j):
                            words.append([0 if i == z else (-1 if i in pos else 1) for i in idxA])
                for w in words:
                    v = psi(w)
                    for a_ in range(len(Ts)):
                        for b_ in range(len(Ts)):
                            Zmod[a_][b_] += Fr(cnt * v[a_] * v[b_], len(words))
            check('C.symmetrisation_identity', Zbar == Zmod, (p, A))
            ok, rk, mp = psd_exact(Zmod)
            check('C.Zbar_psd_exact', ok, (p, A))
            out['symmetrisation_checked'] += 1
        # Krawtchouk-Weil on actual sets, larger m
        for trial in range(4):
            m = rng.randint(3, min(9, p - 1))
            A = rng.sample(range(p), m)
            n_ = [0] * (m + 1); n1 = [0] * m
            for x in range(p):
                e, dl = level_of(x, A, p, chi)
                if dl:
                    n1[e] += 1
                else:
                    n_[e] += 1
            Km, Km1 = ktab(m), ktab(m - 1)
            for t in range(1, m + 1):
                lhs = sum(n_[j] * Km[t][j] for j in range(m + 1)) + sum(n1[j] * (Km1[t][j] if t <= m - 1 else 0) for j in range(m))
                direct = sum(esym([chi[(x + a) % p] for a in A], t) for x in range(p))
                check('C.KW_identity', lhs == direct, (p, A, t))
                if t == 1:
                    check('C.KW_t1_exact', lhs == 0, (p, A))
                elif t == 2:
                    check('C.KW_t2_exact', lhs == -comb(m, 2), (p, A))
                else:
                    check('C.KW_bound', lhs * lhs <= (comb(m, t) * (t - 1)) ** 2 * p, (p, A, t))
                    out['KW_ratio_max'] = max(out['KW_ratio_max'], abs(lhs) / (comb(m, t) * (t - 1) * math.sqrt(p)))
                out['KW_actual_checked'] += 1
    RESULTS['C'] = out


# ------------------------------------------------------------------ Section D: the decoupling construction (Theorem 2.1)

def mat_inv(G):
    n = len(G)
    A = [[Fr(G[i][j]) for j in range(n)] + [Fr(int(i == j)) for j in range(n)] for i in range(n)]
    for c in range(n):
        piv = next(r for r in range(c, n) if A[r][c] != 0)
        A[c], A[piv] = A[piv], A[c]
        pv = A[c][c]
        A[c] = [v / pv for v in A[c]]
        for r in range(n):
            if r != c and A[r][c] != 0:
                f = A[r][c]
                A[r] = [A[r][k] - f * A[c][k] for k in range(2 * n)]
    return [row[n:] for row in A]


def decoupling_construction(p, nu, sets, fvals):
    """Theorem 2.1 construction in the level model L^2(Omega).
    nu: list of cell masses (Fractions); sets: dict name -> (placement list over cells, size, value list over cells)
    where value list gives the level function S~1_X on the cells; overlaps between sets are given
    in sets['_overlap'] (dict (X,Y) -> s).  Returns (Tcells, Tfamily, family names, report)."""
    L = len(nu)
    names = [k for k in sets if not k.startswith('_')]
    ov = sets.get('_overlap', {})
    # atoms: cells 0..L-1, then the sets
    natoms = L + len(names)
    Ga = [[Fr(0)] * natoms for _ in range(natoms)]
    for l in range(L):
        Ga[l][l] = Fr(nu[l])
    for i, X in enumerate(names):
        plc, size, _ = sets[X]
        for l in range(L):
            Ga[l][L + i] = Ga[L + i][l] = Fr(plc[l])
        Ga[L + i][L + i] = Fr(size)
        for k, Y in enumerate(names):
            if k != i:
                Ga[L + i][L + k] = Fr(ov.get((X, Y), ov.get((Y, X), 0)))
    one = [Fr(1)] * L + [Fr(0)] * len(names)
    gens, imgs = [], []
    for i, X in enumerate(names):
        plc, size, vals = sets[X]
        u = [-Fr(size, p)] * L + [Fr(0)] * len(names)
        u[L + i] = Fr(1)
        f = [Fr(v) for v in vals] + [Fr(0)] * len(names)
        gens += [u, f]
        imgs += [f, [p * c for c in u]]
    ip = lambda x, y: sum(x[i] * Ga[i][j] * y[j] for i in range(natoms) if x[i] for j in range(natoms) if y[j])
    Gg = [[ip(a, b) for b in gens] for a in gens]
    Gi = [[ip(a, b) for b in imgs] for a in imgs]
    Gx = [[ip(a, b) for b in imgs] for a in gens]    # <g_a, S g_b>
    rep = {}
    rep['isometry'] = all(Gi[a][b] == p * Gg[a][b] for a in range(len(gens)) for b in range(len(gens)))
    rep['symmetric'] = all(Gx[a][b] == Gx[b][a] for a in range(len(gens)) for b in range(len(gens)))
    rep['orth_to_1'] = all(ip(g, one) == 0 for g in gens)
    # basis of V0 (greedy by exact rank)
    basis = []
    for a in range(len(gens)):
        trial = basis + [a]
        M = [[Gg[i][j] for j in trial] for i in trial]
        ok, rk, _ = psd_exact(M)
        if rk == len(trial):
            basis.append(a)
    Gb = [[Gg[i][j] for j in basis] for i in basis]
    Gbi = mat_inv(Gb)
    Mb = [[Gx[i][j] for j in basis] for i in basis]
    rep['dimV0'] = len(basis)

    def coeff(w):
        rhs = [ip(gens[b], w) for b in basis]
        return [sum(Gbi[i][j] * rhs[j] for j in range(len(basis))) for i in range(len(basis))]
    atoms = [[Fr(int(i == k)) for i in range(natoms)] for k in range(natoms)]
    C = [coeff(w) for w in atoms]
    T = [[sum(C[x][i] * Mb[i][j] * C[y][j] for i in range(len(basis)) for j in range(len(basis)))
          for y in range(natoms)] for x in range(natoms)]
    return T, Ga, names, rep


def verify_construction(p, nu, sets, T, Ga, names, tag):
    L = len(nu)
    natoms = len(T)
    check('D.T_symmetric', all(T[i][j] == T[j][i] for i in range(natoms) for j in range(natoms)), tag)
    check('D.T1_rowsum0', all(sum(T[i][:L]) == 0 for i in range(natoms)), tag)
    for k, X in enumerate(names):
        plc, size, vals = sets[X]
        # T2: sum_l T_{i,l} f_X(l) = p <atom_i, 1_X> - size <atom_i, 1>
        for i in range(natoms):
            lhs = sum(T[i][l] * vals[l] for l in range(L))
            rhs = p * Ga[i][L + k] - size * (Ga[i][i] if i < L else sum(Ga[i][:L]))
            check('D.T2', lhs == rhs, (tag, X, i))
        # column of 1_X: <cell_l, S~ 1_X> = nu_l f_X(l)
        for l in range(L):
            check('D.S1_X=f_X', T[l][L + k] == nu[l] * vals[l], (tag, X, l))
    # Gram under Pi_+- (numerical): G = Ga - uu^T/p +- T/sqrt p, u_i = <atom_i, 1>
    u = np.array([float(sum(Ga[i][:L])) for i in range(natoms)])
    GA = np.array([[float(x) for x in row] for row in Ga])
    TT = np.array([[float(x) for x in row] for row in T])
    worst = 1.0
    for sg in (1, -1):
        G = (GA - np.outer(u, u) / p + sg * TT / math.sqrt(p)) / 2
        ev = np.linalg.eigvalsh(G)
        worst = min(worst, ev.min() / max(1.0, abs(ev).max()))
    check('D.Gpm_psd_numeric', worst > -1e-9, (tag, worst))
    return worst


def section_D_selftest():
    """Construction applied to the level profiles of real configurations (one- and two-sided)."""
    rng = random.Random(5)
    out = {'cases': 0}
    for p in [29, 37, 41, 53, 61, 101, 109, 113, 23, 31, 43]:
        chi = chi_table(p)
        for trial in range(3):
            m = rng.randint(2, 4)
            A = rng.sample(range(p), m)
            B = complete_partners(A, p, chi)
            if not B:
                continue
            n = len(B)
            labA = [level_of(x, A, p, chi) for x in range(p)]
            labB = [level_of(x, B, p, chi) for x in range(p)]
            joint = sorted(set(zip(labA, labB)))
            cnt = {c: 0 for c in joint}
            ca = {c: 0 for c in joint}
            cb = {c: 0 for c in joint}
            for x in range(p):
                c = (labA[x], labB[x])
                cnt[c] += 1
                ca[c] += x in A
                cb[c] += x in B
            nu = [Fr(cnt[c]) for c in joint]
            sets = {'A': ([ca[c] for c in joint], m, [m - c[0][1] - 2 * c[0][0] for c in joint]),
                    'B': ([cb[c] for c in joint], n, [n - c[1][1] - 2 * c[1][0] for c in joint]),
                    '_overlap': {('A', 'B'): len(set(A) & set(B))}}
            T, Ga, names, rep = decoupling_construction(p, nu, sets, None)
            check('D.isometry', rep['isometry'], (p, A))
            check('D.symmetric_on_V0', rep['symmetric'], (p, A))
            check('D.gens_orth_1', rep['orth_to_1'], (p, A))
            verify_construction(p, nu, sets, T, Ga, names, ('selftest', p, tuple(A)))
            out['cases'] += 1
    RESULTS['D_selftest'] = out


# ------------------------------------------------------------------ Section E: hybrid certificates

def lp_max_n0(p, m, with_kw=True, mode='biclique'):
    """max n_0 (real) subject to all theorem rows (+ KW rows), r = 0 (biclique) or n'_0 = m (clique)."""
    from scipy.optimize import linprog
    nv = 2 * m + 1
    rows = balanced_rows(p, m) + (krawtchouk_weil_rows(p, m) if with_kw else [])
    Aeq, beq, Aub, bub = [], [], [], []
    for co, sense, rhs, name in rows:
        v = np.zeros(nv)
        sc = max([abs(c) for c in co.values()] + [abs(rhs), 1])
        for k, c in co.items():
            v[k] = ratio_float(c, sc)
        (Aeq if sense == '=' else Aub).append(v)
        (beq if sense == '=' else bub).append(ratio_float(rhs, sc))
    v = np.zeros(nv); v[m + 1] = 1.0
    Aeq.append(v); beq.append(0.0 if mode == 'biclique' else float(m))
    c = np.zeros(nv); c[0] = -1.0
    res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * nv, method='highs')
    return float(-res.fun) if res.status == 0 else None


def profile_levels(x, m):
    """list of (label, mass, value f) for the nonzero entries of a profile"""
    out = []
    for j in range(m + 1):
        if x[j] > 0:
            out.append((('N', j), x[j], m - 2 * j))
    for j in range(m):
        if x[m + 1 + j] > 0:
            out.append((('A', j), x[m + 1 + j], m - 1 - 2 * j))
    return out


def nw_corner(ra, rb):
    """north-west corner coupling of two nonnegative rational vectors with equal sums"""
    ra, rb = list(ra), list(rb)
    pi = {}
    i = j = 0
    while i < len(ra) and j < len(rb):
        q = min(ra[i], rb[j])
        if q > 0:
            pi[(i, j)] = pi.get((i, j), 0) + q
        ra[i] -= q; rb[j] -= q
        if ra[i] == 0:
            i += 1
        else:
            j += 1
    return pi


def joint_coupling(p, m, n, xA, xB):
    """Two-sided joint profile for Theorem 2.1 (s = 0, r = 0). Returns (nu, sets, info) or None."""
    LA, LB = profile_levels(xA, m), profile_levels(xB, n)
    iA0 = next(i for i, l in enumerate(LA) if l[0] == ('N', 0))
    iB0 = next(i for i, l in enumerate(LB) if l[0] == ('N', 0))
    if LA[iA0][1] != n or LB[iB0][1] != m:
        return None
    candA = [i for i, l in enumerate(LA) if l[0][0] == 'N' and l[0][1] >= 1 and l[1] >= m]
    candB = [i for i, l in enumerate(LB) if l[0][0] == 'N' and l[0][1] >= 1 and l[1] >= n]
    sA = min(candA, key=lambda i: (abs(LA[i][2]), -LA[i][1]))
    sB = min(candB, key=lambda i: (abs(LB[i][2]), -LB[i][1]))
    RA = [l[1] for l in LA]; RB = [l[1] for l in LB]
    RA[iA0] -= n; RA[sA] -= m
    RB[iB0] -= m; RB[sB] -= n
    tot = sum(RA)
    assert tot == sum(RB) == p - m - n
    target = -m * n - m * n * (LA[sA][2] + LB[sB][2])
    prod = {(i, j): RA[i] * RB[j] / tot for i in range(len(LA)) for j in range(len(LB)) if RA[i] > 0 and RB[j] > 0}
    V = lambda pi: sum(q * LA[i][2] * LB[j][2] for (i, j), q in pi.items())
    oa = sorted(range(len(LA)), key=lambda i: LA[i][2])
    ob_desc = sorted(range(len(LB)), key=lambda j: -LB[j][2])
    ob_asc = sorted(range(len(LB)), key=lambda j: LB[j][2])
    anti = {(oa[i], ob_desc[j]): q for (i, j), q in nw_corner([RA[k] for k in oa], [RB[k] for k in ob_desc]).items()}
    mono = {(oa[i], ob_asc[j]): q for (i, j), q in nw_corner([RA[k] for k in oa], [RB[k] for k in ob_asc]).items()}
    v0, v1, v2 = V(prod), V(anti), V(mono)
    if v1 <= target <= v0:
        lam = (target - v0) / (v1 - v0) if v1 != v0 else Fr(0)
        other = anti
    elif v0 <= target <= v2:
        lam = (target - v0) / (v2 - v0) if v2 != v0 else Fr(0)
        other = mono
    else:
        return None
    pi = {}
    for k, q in prod.items():
        pi[k] = pi.get(k, 0) + (1 - lam) * q
    for k, q in other.items():
        pi[k] = pi.get(k, 0) + lam * q
    cells, nu, fa, fb, pa, pb = [], [], [], [], [], []
    cells.append(('B', iA0, sB)); nu.append(Fr(n)); fa.append(LA[iA0][2]); fb.append(LB[sB][2]); pa.append(0); pb.append(n)
    cells.append(('A', sA, iB0)); nu.append(Fr(m)); fa.append(LA[sA][2]); fb.append(LB[iB0][2]); pa.append(m); pb.append(0)
    for (i, j), q in sorted(pi.items()):
        if q > 0:
            cells.append(('R', i, j)); nu.append(Fr(q)); fa.append(LA[i][2]); fb.append(LB[j][2]); pa.append(0); pb.append(0)
    sets = {'A': (pa, m, fa), 'B': (pb, n, fb), '_overlap': {('A', 'B'): 0}}
    info = {'sigma_A': m * LA[sA][2], 'sigma_B': n * LB[sB][2], 'A_level': LA[sA][0], 'B_level': LB[sB][0],
            'lambda': float(lam), 'coupling': 'anti' if other is anti else 'mono', 'n_cells': len(cells)}
    # the hypotheses of Theorem 2.1, exactly
    check('E.joint_marginal_mass', sum(nu) == p, (p, m))
    check('E.joint_moment1', sum(v * a for v, a in zip(nu, fa)) == 0 and sum(v * b for v, b in zip(nu, fb)) == 0, (p, m))
    check('E.joint_moment2', sum(v * a * a for v, a in zip(nu, fa)) == m * (p - m) and
          sum(v * b * b for v, b in zip(nu, fb)) == n * (p - n), (p, m))
    check('E.joint_cross_moment', sum(v * a * b for v, a, b in zip(nu, fa, fb)) == -m * n, (p, m))
    check('E.joint_nonneg', all(v >= 0 for v in nu), (p, m))
    return nu, sets, info


def certificate_summary(r):
    return {('n_%d' % j if j <= (len(r['x']) - 1) // 2 else "n'_%d" % (j - (len(r['x']) - 1) // 2 - 1)): str(v)
            for j, v in enumerate(r['x']) if v != 0}


def kw_max_ratio(p, m, xr):
    Km, Km1 = ktab(m), ktab(m - 1)
    best = 0.0
    for t in range(3, m + 1):
        S_ = sum(xr[j] * Km[t][j] for j in range(m + 1) if xr[j]) + \
            sum(xr[m + 1 + j] * Km1[t][j] for j in range(m) if xr[m + 1 + j] and t <= m - 1)
        best = max(best, ratio_float(abs(S_), comb(m, t) * (t - 1) * isqrt(p)))
    return best


def section_E(instances):
    out = []
    certs = {}
    for (p, m, n, kind) in instances:
        t0 = time.time()
        d = (p - 1) // 2
        rec = {'p': p, 'm': m, 'kind': kind}
        if kind == 'biclique':
            rA = solve_certificate(p, m, 'biclique', n)
            check('E.cert_found', rA['status'] == 'ok', (p, m))
            check('E.cert_exact', rA['exact_ok'], (p, m, [o for o in rA['rows'] if not o[3]][:3]))
            if n != m:
                rB = solve_certificate(p, n, 'biclique', m)
                check('E.cert_found', rB['status'] == 'ok', (p, n))
                check('E.cert_exact', rB['exact_ok'], (p, n))
            else:
                rB = rA
            rec['KW_max_ratio'] = kw_max_ratio(p, m, rA['x'])
            rec.update({'n': n, 'mn_over_p': m * n / p, 'g': d - m * n, 'HP_real_n0': d / m,
                        'rows_A': rA['n_rows'], 'KW_rows': 2 * (m - 2), 'lp_common_slack': rA['lp_slack'],
                        'band': rA['band'], 'support_A': certificate_summary(rA),
                        'support_B': certificate_summary(rB) if n != m else 'same as A'})
            jc = joint_coupling(p, m, n, rA['x'], rB['x'])
            check('E.joint_coupling_found', jc is not None, (p, m))
            nu, sets, info = jc
            rec['joint'] = info
            ok = p1_psd(m, n, 0, info['sigma_A'], info['sigma_B'], m * n, p)
            check('E.P1_exact', ok, (p, m, info))
            T, Ga, names, rep = decoupling_construction(p, nu, sets, None)
            check('E.isometry', rep['isometry'], (p, m))
            check('E.symmetric_on_V0', rep['symmetric'], (p, m))
            check('E.gens_orth_1', rep['orth_to_1'], (p, m))
            rec['dimV0'] = rep['dimV0']
            rec['Gpm_min_rel_eig'] = verify_construction(p, nu, sets, T, Ga, names, ('E', p, m))
            certs[(p, m, kind)] = rA['x']
            if n != m:
                certs[(p, n, kind)] = rB['x']
        else:
            r = solve_certificate(p, m, 'clique', 0)
            check('E.cert_found', r['status'] == 'ok', (p, m))
            check('E.cert_exact', r['exact_ok'], (p, m, [o for o in r['rows'] if not o[3]][:3]))
            rec['KW_max_ratio'] = kw_max_ratio(p, m, r['x'])
            rec.update({'m(m-1)': m * (m - 1), 'd': d, 'm2_over_p': m * m / p, 'rows': r['n_rows'],
                        'KW_rows': 2 * (m - 2), 'lp_common_slack': r['lp_slack'], 'band': r['band'],
                        'support': certificate_summary(r)})
            # P1 in difference form (theta bound): m(m-1) <= sqrt(p)(m - m^2/p)
            ok = psd2_sqrt(Fr(m) - Fr(m * m, p), Fr(m * (m - 1), p), Fr(0), Fr(0), Fr(1), Fr(0), p) and \
                psd2_sqrt(Fr(m) - Fr(m * m, p), Fr(-m * (m - 1), p), Fr(0), Fr(0), Fr(1), Fr(0), p)
            check('E.P1_clique_exact', ok, (p, m))
            # one-sided construction with C = the cell of -A' (n'_0 = m), value m-1 there
            L = profile_levels(r['x'], m)
            nu = [l[1] for l in L]
            plc = [m if l[0] == ('A', 0) else 0 for l in L]
            sets = {'C': (plc, m, [l[2] for l in L])}
            T, Ga, names, rep = decoupling_construction(p, nu, sets, None)
            check('E.isometry', rep['isometry'], (p, m))
            check('E.symmetric_on_V0', rep['symmetric'], (p, m))
            rec['Gpm_min_rel_eig'] = verify_construction(p, nu, sets, T, Ga, names, ('Eclique', p, m))
            certs[(p, m, kind)] = r['x']
        rec['lp_max_n0_real_with_KW'] = lp_max_n0(p, m, True, 'biclique') if kind == 'biclique' else None
        rec['seconds'] = round(time.time() - t0, 2)
        out.append(rec)
    RESULTS['E'] = out
    return certs


# ------------------------------------------------------------------ Section G: the two-point (Delsarte-type) program P5

def Hset(j, k, M):
    if j < 0 or k < 0 or j > M or k > M:
        return []
    return [h for h in range(abs(j - k), min(j + k, 2 * M - j - k) + 1) if (h - j - k) % 2 == 0]


def tp_structure(p, m, xr):
    """Pair classes. levels: N_j (x not in -A, e_x = j) and A_j (x in -A).  Variables: ordered pairs
    (x in level a, y in level b), stored for a <= b for NN and AA and for (N, A) for NA, split by
    nu = #nonzero coordinates of s(x)s(y), h = #(-1) among them, and eps (value of x at y's zero),
    eps2 (value of y at x's zero; eps2 = chi(-1) eps)."""
    chim1 = 1 if p % 4 == 1 else -1
    levels = [('N', j, xr[j]) for j in range(m + 1) if xr[j] > 0] + \
             [('A', j, xr[m + 1 + j]) for j in range(m) if xr[m + 1 + j] > 0]
    L = len(levels)
    var = []
    for a in range(L):
        for b in range(L):
            ta, j, _ = levels[a]
            tb, k, _ = levels[b]
            if ta == 'N' and tb == 'N' and a <= b:
                for h in Hset(j, k, m):
                    var.append((a, b, m, h, 0, 0))
            elif ta == 'N' and tb == 'A':
                for eps in (1, -1):
                    for h in Hset(j - (eps == -1), k, m - 1):
                        var.append((a, b, m - 1, h, eps, 0))
            elif ta == 'A' and tb == 'A' and a <= b:
                for eps in (1, -1):
                    eps2 = chim1 * eps
                    for h in Hset(j - (eps == -1), k - (eps2 == -1), m - 2):
                        var.append((a, b, m - 2, h, eps, eps2))
    return levels, var


def tp_equalities(p, m, levels, var):
    """Linear equalities E P = e valid for every real configuration (exact coefficients).
    Returns list of (dict var->coef, rhs, name)."""
    L = len(levels)
    KT = {m: ktab(m), m - 1: ktab(m - 1), m - 2: ktab(m - 2)}
    eqs = []
    for a in range(L):
        for b in range(L):
            idx = [i for i, v in enumerate(var) if v[0] == a and v[1] == b]
            if not idx:
                continue
            na, nb = levels[a][2], levels[b][2]
            tot = na * nb - (na if (a == b and levels[a][0] == 'A') else 0)
            eqs.append(({i: 1 for i in idx}, tot, 'marg_%d_%d' % (a, b)))
    for a in range(L):
        ta, j, na = levels[a]
        if ta == 'N':
            idx = [i for i, v in enumerate(var) if v[0] == a and levels[v[1]][0] == 'A' and v[4] == -1]
        else:
            idx = [i for i, v in enumerate(var) if v[0] == a and levels[v[1]][0] == 'A' and v[4] == -1] + \
                  [i for i, v in enumerate(var) if v[1] == a and v[0] != a and levels[v[0]][0] == 'A' and v[5] == -1]
        if idx or j * na != 0:
            eqs.append(({i: 1 for i in idx}, j * na, 'epssplit_%d' % a))
    # t = 1 and t = 2 row identities: sum_b Gamma^(t)_{ab} = 0, resp. = -g^(2)_a
    for t in (1, 2):
        for a in range(L):
            co = {}
            const = Fr(0)
            for i, (x, y, nu, h, e1, e2) in enumerate(var):
                c = KT[nu][t][h] if t <= nu else 0
                if c == 0:
                    continue
                if x == a:
                    co[i] = co.get(i, 0) + c
                if y == a and x != a:
                    co[i] = co.get(i, 0) + c
            if levels[a][0] == 'A':
                const += levels[a][2] * comb(m - 1, t)
            ga = levels[a][2] * (KT[m][t][levels[a][1]] if levels[a][0] == 'N' else KT[m - 1][t][levels[a][1]])
            rhs = (0 if t == 1 else -ga) - const
            eqs.append((co, rhs, 'row_t%d_%d' % (t, a)))
    return eqs


def tp_gamma_exact(p, m, levels, var, P, t):
    L = len(levels)
    KT = {m: ktab(m), m - 1: ktab(m - 1), m - 2: ktab(m - 2)}
    G = [[Fr(0)] * L for _ in range(L)]
    for i, (a, b, nu, h, e1, e2) in enumerate(var):
        c = KT[nu][t][h] if t <= nu else 0
        if c and P[i]:
            G[a][b] += c * P[i]
    for a in range(L):
        for b in range(a):
            G[a][b] = G[b][a]
        if levels[a][0] == 'A':
            G[a][a] += levels[a][2] * comb(m - 1, t)
    g = [levels[a][2] * (KT[m][t][levels[a][1]] if levels[a][0] == 'N' else (KT[m - 1][t][levels[a][1]] if t <= m - 1 else 0))
         for a in range(L)]
    aug = [G[a] + [g[a]] for a in range(L)] + [g + [Fr(comb(m, t))]]
    return aug


def forced_nulls_num(levels, m, t):
    """null directions of aug^(t) forced by exact identities (used only to measure the SDP margin)"""
    L = len(levels)
    nulls = []
    if t == 1:
        nulls.append(np.concatenate([np.ones(L), [0.0]]))
    if t == 2:
        nulls.append(np.concatenate([np.ones(L), [1.0]]))
    for a in range(L):
        tp_, j, na = levels[a]
        if tp_ == 'N' and j == 0:
            v = np.zeros(L + 1); v[a] = 1.0; v[L] = -float(na); nulls.append(v)
        if tp_ == 'A' and float(na) == m:
            if j == 0:
                v = np.zeros(L + 1); v[a] = 1.0; v[L] = -float(m - t); nulls.append(v)
            if t == m - 1:
                v = np.zeros(L + 1); v[a] = 1.0; v[L] = -float((-1) ** j); nulls.append(v)
    return nulls


def tp_solve(p, m, xr, clique=False):
    """SDP (Clarabel) maximising the smallest eigenvalue of aug^(t) on the complement of the
    forced null directions; sparse formulation."""
    import cvxpy as cp
    import scipy.sparse as sp
    import warnings
    warnings.filterwarnings('ignore')
    levels, var = tp_structure(p, m, xr)
    L, nvar = len(levels), len(var)
    PS = float(p) ** 2
    Pn = cp.Variable(nvar, nonneg=True)
    cons = []
    eqs = tp_equalities(p, m, levels, var)
    rr, cc, vv, bb = [], [], [], []
    for k, (co, rhs, name) in enumerate(eqs):
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)) / PS, 1.0])
        for i_, c in co.items():
            rr.append(k); cc.append(i_); vv.append(float(c) / sc)
        bb.append(float(rhs) / PS / sc)
    Eq = sp.csr_matrix((vv, (rr, cc)), shape=(len(eqs), nvar))
    cons.append(Eq @ Pn == np.array(bb))
    for a in range(L):
        if levels[a][0] == 'N':
            idx = [i_ for i_, v in enumerate(var) if v[0] == a and v[1] == a and v[3] == 0]
            cons.append(cp.sum(Pn[idx]) >= float(levels[a][2]) / PS)
    KT = {m: ktab(m), m - 1: ktab(m - 1), m - 2: ktab(m - 2)}
    slack = cp.Variable()
    skip = set([m]) | (set([m - 1, m - 2]) if clique else set())
    D = L + 1
    for t in range(1, m + 1):
        sc = float(comb(m, t)) * PS
        rr, cc, vv = [], [], []
        for i_, (a, b, nu, h, e1, e2) in enumerate(var):
            c = KT[nu][t][h] if t <= nu else 0
            if c:
                w = float(c) * PS / sc
                rr.append(a * D + b); cc.append(i_); vv.append(w)
                if a != b:
                    rr.append(b * D + a); cc.append(i_); vv.append(w)
        Ct = sp.csr_matrix((vv, (rr, cc)), shape=(D * D, nvar))
        const = np.zeros((D, D))
        for a in range(L):
            if levels[a][0] == 'A':
                const[a, a] = float(levels[a][2]) * comb(m - 1, t) / sc
            ga = float(levels[a][2]) * (KT[m][t][levels[a][1]] if levels[a][0] == 'N' else (KT[m - 1][t][levels[a][1]] if t <= m - 1 else 0)) / sc
            const[a, L] = const[L, a] = ga
        const[L, L] = comb(m, t) / sc
        Mt = cp.reshape(Ct @ Pn, (D, D), order='C') + const
        Ms = (Mt + Mt.T) / 2
        if t >= 3:
            ones = np.zeros(D * D)
            for a in range(L):
                for b in range(L):
                    ones[a * D + b] = 1.0
            cons.append(ones @ (Ct @ Pn) + float(sum(const[a, a] for a in range(L))) <= 0.999 * comb(m, t) * (t - 1) ** 2 * p / sc)
        if t in skip:
            cons.append(Ms >> 0)
        else:
            nulls = forced_nulls_num(levels, m, t)
            if nulls:
                u_, s_, vt_ = np.linalg.svd(np.column_stack(nulls), full_matrices=True)
                rk = int((s_ > 1e-9).sum())
                Nb = u_[:, rk:]
            else:
                Nb = np.eye(D)
            cons.append(Nb.T @ Ms @ Nb >> slack * np.eye(Nb.shape[1]))
    prob = cp.Problem(cp.Maximize(slack), cons + [slack <= 1e-2])
    try:
        prob.solve(solver='CLARABEL')
    except Exception:
        prob.solve(solver='SCS', eps=1e-10, max_iters=200000)
    return levels, var, (None if Pn.value is None else Pn.value * PS), prob.status, float(slack.value) if slack.value is not None else None


def minnorm_repair(eqs, P0, nvar, support):
    """exact minimum-norm correction supported on `support` with E (P0 + delta) = e."""
    S = list(support)
    pos = {i: k for k, i in enumerate(S)}
    rows = []
    for co, rhs, name in eqs:
        r_ = Fr(rhs) - sum(Fr(c) * P0[i] for i, c in co.items())
        rows.append(({pos[i]: Fr(c) for i, c in co.items() if i in pos}, r_))
    # independent subset of rows (exact elimination on dense copies)
    dense = [[row.get(k, Fr(0)) for k in range(len(S))] for row, _ in rows]
    basis, red = [], []
    for ri, vec in enumerate(dense):
        v = list(vec)
        for (bi, bv, pc) in red:
            if v[pc] != 0:
                f = v[pc] / bv[pc]
                v = [x - f * y for x, y in zip(v, bv)]
        pc = next((k for k, x in enumerate(v) if x != 0), None)
        if pc is not None:
            red.append((ri, v, pc))
            basis.append(ri)
    EI = [dense[ri] for ri in basis]
    rI = [rows[ri][1] for ri in basis]
    G = [[sum(a * b for a, b in zip(EI[i], EI[j]) if a and b) for j in range(len(EI))] for i in range(len(EI))]
    Gi = mat_inv(G)
    lam = [sum(Gi[i][j] * rI[j] for j in range(len(EI))) for i in range(len(EI))]
    delta = [sum(EI[i][k] * lam[i] for i in range(len(EI))) for k in range(len(S))]
    P = list(P0)
    for k, i in enumerate(S):
        P[i] = P0[i] + delta[k]
    return P, len(basis)


def tp_certificate(p, m, xr, clique=False, Q=10 ** 6):
    t0 = time.time()
    levels, var, Pf, status, slack = tp_solve(p, m, xr, clique)
    out = {'p': p, 'm': m, 'kind': 'clique' if clique else 'biclique', 'n_levels': len(levels), 'n_pair_vars': len(var),
           'sdp_status': status, 'sdp_margin': slack}
    if Pf is None:
        out['exact_ok'] = False
        return out
    P0 = [Fr(int(round(v * Q)), Q) if v > 0 else Fr(0) for v in Pf]
    eqs = tp_equalities(p, m, levels, var)
    mx = max(Pf)
    support = [i for i in range(len(var)) if Pf[i] > 1e-6 * mx]
    P, rank = minnorm_repair(eqs, P0, len(var), support)
    consistent = True
    okE = all(sum(Fr(c) * P[i] for i, c in co.items()) == Fr(rhs) for co, rhs, name in eqs)
    check('G.equalities_exact', okE, (p, m))
    nonneg = all(v >= 0 for v in P)
    check('G.nonneg', nonneg, (p, m, min(P)))
    for a in range(len(levels)):
        if levels[a][0] == 'N':
            idx = [i for i, v in enumerate(var) if v[0] == a and v[1] == a and v[3] == 0]
            check('G.diagonal_pairs', sum(P[i] for i in idx) >= levels[a][2], (p, m, a))
    psd_all = True
    phi_ok = True
    ranks = []
    phimax = 0.0
    for t in range(1, m + 1):
        aug = tp_gamma_exact(p, m, levels, var, P, t)
        L = len(levels)
        if t >= 3:
            Phi = sum(aug[a][b] for a in range(L) for b in range(L))
            okphi = Phi <= comb(m, t) * (t - 1) ** 2 * p
            check('G.weil_Phi', okphi, (p, m, t))
            phi_ok &= okphi
            phimax = max(phimax, float(Phi / (comb(m, t) * (t - 1) ** 2 * p)))
        elif t == 1:
            check('G.Phi1=0', sum(aug[a][b] for a in range(L) for b in range(L)) == 0, (p, m))
        else:
            check('G.Phi2=C(m,2)', sum(aug[a][b] for a in range(L) for b in range(L)) == comb(m, 2), (p, m))
        ok, rk, mp = psd_exact(aug)
        check('G.aug_psd_exact', ok, (p, m, t))
        psd_all &= ok
        ranks.append(rk)
    out.update({'exact_ok': bool(okE and nonneg and psd_all and consistent and phi_ok), 'aug_ranks': ranks,
                'repair_rank': rank, 'max_Phi_over_weil_bound': phimax, 'seconds': round(time.time() - t0, 2)})
    return out


# ------------------------------------------------------------------ Section F: independent SDP solve of the level-set relaxation

def section_F(p=1009, m=22):
    """max n_0 over: count LP rows + KW rows, r = 0, placement alpha of A over the 2m+1 levels,
    a free symmetric compression T (levels x levels) with T1 (T 1 = 0), T2 (T f = p alpha - m nu),
    and T3: Gram of {1_levels, 1_A} under (I - J/p +- S/sqrt p) PSD (Schur form).  Numerical (Clarabel).
    Scaled variables: xt = nu/p, Tt = T/p^1.5, family w_l = 1_l/sqrt p, w_A = 1_A/sqrt m."""
    import cvxpy as cp
    import warnings
    warnings.filterwarnings('ignore')
    nv = 2 * m + 1
    rows = balanced_rows(p, m) + krawtchouk_weil_rows(p, m)
    xt = cp.Variable(nv, nonneg=True)
    al = cp.Variable(nv, nonneg=True)
    Tt = cp.Variable((nv, nv), symmetric=True)
    cons = [xt[m + 1] == 0, al <= p * xt, cp.sum(al) == m]
    for co, sense, rhs, name in rows:
        sc = max([abs(c) for c in co.values()] + [abs(rhs), 1])
        expr = sum(ratio_float(c, sc) * p * xt[k] for k, c in co.items())
        cons.append(expr == ratio_float(rhs, sc) if sense == '=' else expr <= ratio_float(rhs, sc))
    f = np.array([val(m, j, 0) for j in range(m + 1)] + [val(m, j, 1) for j in range(m)], dtype=float)
    cons.append(Tt @ np.ones(nv) == 0)
    cons.append(Tt @ f == (al - m * xt) / math.sqrt(p))
    sig = f @ al
    col = lambda v, k: cp.reshape(v, (k, 1), order='C')
    GI = cp.bmat([[cp.diag(xt), col(al / math.sqrt(p * m), nv)], [col(al / math.sqrt(p * m), nv).T, np.array([[1.0]])]])
    GSs = cp.bmat([[Tt, col(cp.multiply(f, xt) / math.sqrt(m), nv)],
                   [col(cp.multiply(f, xt) / math.sqrt(m), nv).T, cp.reshape(sig / (m * math.sqrt(p)), (1, 1), order='C')]])
    u = cp.hstack([xt, np.array([math.sqrt(m / p)])])
    for sg in (1, -1):
        M = cp.bmat([[GI + sg * GSs, col(u, nv + 1)], [col(u, nv + 1).T, np.array([[1.0]])]])
        cons.append((M + M.T) / 2 >> 0)
    prob = cp.Problem(cp.Maximize(xt[0]), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception:
        prob.solve(solver='SCS', eps=1e-9, max_iters=200000)
    lp = lp_max_n0(p, m, True)
    out = {'p': p, 'm': m, 'status': prob.status, 'sdp_max_n0': float(prob.value) * p, 'lp_max_n0': lp,
           'HP_d_over_m': ((p - 1) // 2) / m}
    check('F.sdp_equals_lp', abs(out['sdp_max_n0'] - lp) < 1e-4 * lp, out)
    RESULTS['F'] = out


def tp_feasibility(p, m, xr):
    """pure feasibility solve of P5 for a given profile (numerical, Clarabel status)"""
    import cvxpy as cp
    import scipy.sparse as sp
    import warnings
    warnings.filterwarnings('ignore')
    levels, var = tp_structure(p, m, xr)
    L, nvar = len(levels), len(var)
    PS = float(p) ** 2
    Pn = cp.Variable(nvar, nonneg=True)
    cons = []
    eqs = tp_equalities(p, m, levels, var)
    rr, cc, vv, bb = [], [], [], []
    for k, (co, rhs, name) in enumerate(eqs):
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)) / PS, 1.0])
        for i_, c in co.items():
            rr.append(k); cc.append(i_); vv.append(float(c) / sc)
        bb.append(float(rhs) / PS / sc)
    cons.append(sp.csr_matrix((vv, (rr, cc)), shape=(len(eqs), nvar)) @ Pn == np.array(bb))
    for a in range(L):
        if levels[a][0] == 'N':
            idx = [i_ for i_, v in enumerate(var) if v[0] == a and v[1] == a and v[3] == 0]
            cons.append(cp.sum(Pn[idx]) >= float(levels[a][2]) / PS)
    KT = {m: ktab(m), m - 1: ktab(m - 1), m - 2: ktab(m - 2)}
    D = L + 1
    for t in range(1, m + 1):
        sc = float(comb(m, t)) * PS
        rr, cc, vv = [], [], []
        for i_, (a, b, nu, h, e1, e2) in enumerate(var):
            c = KT[nu][t][h] if t <= nu else 0
            if c:
                w = float(c) * PS / sc
                rr.append(a * D + b); cc.append(i_); vv.append(w)
                if a != b:
                    rr.append(b * D + a); cc.append(i_); vv.append(w)
        Ct = sp.csr_matrix((vv, (rr, cc)), shape=(D * D, nvar))
        const = np.zeros((D, D))
        for a in range(L):
            if levels[a][0] == 'A':
                const[a, a] = float(levels[a][2]) * comb(m - 1, t) / sc
            ga = float(levels[a][2]) * (KT[m][t][levels[a][1]] if levels[a][0] == 'N' else (KT[m - 1][t][levels[a][1]] if t <= m - 1 else 0)) / sc
            const[a, L] = const[L, a] = ga
        const[L, L] = comb(m, t) / sc
        Mt = cp.reshape(Ct @ Pn, (D, D), order='C') + const
        if t >= 3:
            ones = np.zeros(D * D)
            for a in range(L):
                for b in range(L):
                    ones[a * D + b] = 1.0
            cons.append(ones @ (Ct @ Pn) + float(sum(const[a, a] for a in range(L))) <= comb(m, t) * (t - 1) ** 2 * p / sc)
        cons.append((Mt + Mt.T) / 2 >> 0)
    prob = cp.Problem(cp.Minimize(0), cons)
    try:
        prob.solve(solver='CLARABEL')
    except Exception:
        return 'solver_error'
    return prob.status


def section_G2(p=1009, m=22, n0=22):
    """P5 on extreme one-point-feasible profiles (theorem rows + KW, r = 0, n_0 = n0, no band)."""
    from scipy.optimize import linprog
    nv = 2 * m + 1
    rows = balanced_rows(p, m) + krawtchouk_weil_rows(p, m)
    Aeq, beq, Aub, bub = [], [], [], []
    for co, sense, rhs, name in rows:
        v = np.zeros(nv)
        sc = max([abs(c) for c in co.values()] + [abs(rhs), 1])
        for k, c in co.items():
            v[k] = ratio_float(c, sc)
        (Aeq if sense == '=' else Aub).append(v)
        (beq if sense == '=' else bub).append(ratio_float(rhs, sc))
    for k, fv in ((0, n0), (m + 1, 0)):
        v = np.zeros(nv); v[k] = 1.0
        Aeq.append(v); beq.append(float(fv))
    tests = {}
    c = np.zeros(nv)
    for j in range(m):
        c[m + 1 + j] = -j
    tests['max_sum_j_nprime_j'] = c
    c = np.zeros(nv)
    for j in range(m + 1):
        c[j] = -(m - 2 * j) ** 4
    tests['max_fourth_moment'] = c
    out = []
    for name, c in tests.items():
        res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                      bounds=[(0, None)] * nv, method='highs')
        x = res.x
        xr = [Fr(int(round(v * 1000)), 1000) if v > 1e-7 else Fr(0) for v in x]
        xr[0] = Fr(n0); xr[m + 1] = Fr(0)
        jp = max(range(1, m), key=lambda j: x[m + 1 + j])
        xr[m + 1 + jp] = Fr(m) - sum(xr[m + 1 + j] for j in range(m) if j != jp)
        piv = sorted(range(1, m + 1), key=lambda j: -x[j])[:3]
        rest = [j for j in range(m + 1) if j not in piv]
        r1 = p - sum(xr[m + 1 + j] for j in range(m)) - sum(xr[j] for j in rest)
        r2 = 0 - sum(xr[m + 1 + j] * val(m, j, 1) for j in range(m)) - sum(xr[j] * val(m, j, 0) for j in rest)
        r3 = m * (p - m) - sum(xr[m + 1 + j] * val(m, j, 1) ** 2 for j in range(m)) - sum(xr[j] * val(m, j, 0) ** 2 for j in rest)
        sol = solve3([[Fr(1)] * 3, [Fr(val(m, j, 0)) for j in piv], [Fr(val(m, j, 0) ** 2) for j in piv]], [r1, r2, r3])
        for k, j in enumerate(piv):
            xr[j] = sol[k]
        ok, _ = exact_eval_rows(rows, xr)
        ok = ok and all(v >= 0 for v in xr)
        check('G2.extreme_profile_exact', ok, name)
        st = tp_feasibility(p, m, xr) if ok else 'skipped'
        out.append({'objective': name, 'profile_exact_ok': ok, 'support': certificate_summary({'x': xr}),
                    'P5_feasibility_status': st})
    RESULTS['G2'] = out


# ------------------------------------------------------------------ Section H: bipartite theta (symmetrised), numerical

def section_H(primes=(101, 1009, 10009, 100049)):
    """theta-type bound for bicliques: max sum of the cross block of Y >= 0 (PSD), Tr Y_AA = Tr Y_BB = 1,
    Y_{x,y'} = 0 if chi(x+y) = -1 (optionally Y >= 0 entrywise).  Averaging over x -> qx+t, y -> qy-t
    (q in Q) reduces it to 7 variables and 2x2 PSD conditions at the three frequency classes 0, Q, N
    (p = 1 mod 4).  A biclique gives the feasible point v v^T, v = (1_A/sqrt m, 1_B/sqrt n), value sqrt(mn)."""
    import cvxpy as cp
    out = []
    for p in primes:
        if p % 4 != 1:
            continue
        d = (p - 1) // 2
        sq = math.sqrt(p)
        res = {}
        for nonneg in (False, True):
            a0, aQ, aN, b0, bQ, bN, g0, gQ = [cp.Variable() for _ in range(8)]
            def hat(z0, zQ, zN, eps):
                if eps == 0:
                    return z0 + d * zQ + d * zN
                return z0 + zQ * (-1 + eps * sq) / 2 + zN * (-1 - eps * sq) / 2
            # variables are p times the kernel values (a0 = b0 = 1 <-> trace 1); hats are p times the transforms,
            # divided by d for conditioning (PSD conditions are homogeneous)
            cons = [a0 == 1, b0 == 1]
            for eps in (0, 1, -1):
                A_ = hat(a0, aQ, aN, eps) / d; B_ = hat(b0, bQ, bN, eps) / d; G_ = hat(g0, gQ, 0, eps) / d
                cons.append(cp.bmat([[A_, G_], [G_, B_]]) >> 0)
            if nonneg:
                cons += [aQ >= 0, aN >= 0, bQ >= 0, bN >= 0, g0 >= 0, gQ >= 0]
            prob = cp.Problem(cp.Maximize(hat(g0, gQ, 0, 0) / d), cons)
            try:
                prob.solve(solver='CLARABEL')
            except Exception:
                prob.solve(solver='SCS', eps=1e-10, max_iters=100000)
            res['nonneg' if nonneg else 'plain'] = float(prob.value) * d
        out.append({'p': p, 'theta_bip': res['plain'], 'theta_bip_schrijver': res['nonneg'],
                    'theta2_over_p': res['plain'] ** 2 / p, 'theta2_schrijver_over_p': res['nonneg'] ** 2 / p})
    RESULTS['H'] = out


# ------------------------------------------------------------------ main

def jsonable(o):
    if isinstance(o, Fr):
        return str(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    return o


def main():
    timings = {}
    t = time.time(); section_A(); timings['A'] = round(time.time() - t, 2)
    t = time.time(); section_B(); timings['B'] = round(time.time() - t, 2)
    t = time.time(); section_C(); timings['C'] = round(time.time() - t, 2)
    t = time.time(); section_D_selftest(); timings['D'] = round(time.time() - t, 2)
    t = time.time()
    inst = []
    for (p, m, n) in [(1009, 22, 22), (10009, 71, 70), (40009, 141, 141), (100049, 224, 223), (1000033, 707, 707)]:
        inst.append((p, m, n, 'biclique'))
        inst.append((p, m, None, 'clique'))
    certs = section_E(inst)
    timings['E'] = round(time.time() - t, 2)
    t = time.time(); section_F(); timings['F'] = round(time.time() - t, 2)
    t = time.time()
    tp = []
    for (p, m, n) in TP_INSTANCES:
        for kind in ('biclique', 'clique'):
            xr = certs[(p, m, kind)]
            o = tp_certificate(p, m, xr, clique=(kind == 'clique'))
            check('G.certificate_exact', o['exact_ok'], (p, m, kind))
            tp.append(o)
    RESULTS['G'] = tp
    timings['G'] = round(time.time() - t, 2)
    t = time.time(); section_G2(); timings['G2'] = round(time.time() - t, 2)
    t = time.time(); section_H(); timings['H'] = round(time.time() - t, 2)
    RESULTS['timings'] = timings
    RESULTS['checks'] = CHECKS
    RESULTS['n_checks'] = sum(CHECKS.values())
    RESULTS['failures'] = FAILS
    RESULTS['n_failures'] = len(FAILS)
    RESULTS['elapsed'] = round(time.time() - T0, 1)
    with open(OUT, 'w') as fh:
        json.dump(jsonable(RESULTS), fh, indent=1)
    print('checks:', RESULTS['n_checks'], 'failures:', len(FAILS), 'elapsed:', RESULTS['elapsed'])
    for r in RESULTS['E']:
        print('E', r['p'], r['m'], r['kind'], r.get('mn_over_p', r.get('m2_over_p')), r.get('lp_max_n0_real_with_KW'), r['seconds'])
    for o in RESULTS['G']:
        print('G', o['p'], o['m'], o['kind'], o['exact_ok'], o.get('sdp_margin'), o.get('max_Phi_over_weil_bound'), o.get('seconds'))
    for o in RESULTS['G2']:
        print('G2', o['objective'], o['profile_exact_ok'], o['P5_feasibility_status'])
    if FAILS:
        print('FAILURES:', FAILS[:10])


TP_INSTANCES = [(1009, 22, 22), (10009, 71, 70), (40009, 141, 141), (100049, 224, 223)]

if __name__ == '__main__':
    main()
