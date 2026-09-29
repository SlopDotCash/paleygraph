#!/usr/bin/env python3
"""Verifier for research/stepanov-quartic-2026-09-29.md (worker `quartic`, Stepanov wave).

Run:  /opt/miniconda3/bin/python3 experiments/stepanov_quartic_2026_09_29.py
Writes results/stepanov_quartic_2026_09_29.json (check counts, failures, witnesses, certificates).

Sections
  A  identities and inequalities of section 1 of the note, exactly, at small primes:
     A1 Jacobi second-moment identities for all pairs of characters of order K in {4, 8};
     A2 |J(a,b)|^2 = p, K(a,-a) = -1;  A3 class rows (star)_k and subset HP_k, k in {2,4,8}
     (robust note Prop. 2.9), exhaustive at small p;  A4 pattern Hanson--Petridis (Prop. 1.6) with exact
     polynomial degrees and orders;  A5 W_k^{k/2} - 1 (Prop. 1.7);  A6 Weil fourth moments (Lemma 1.4);
     A7 Theorem 2.1(3) on random abstract column statistics;  A8 the exact Krawtchouk--Weil engine
     (Newton identities + reflection) against brute force over subsets, and the float engine against it.
  B  actual complete bicliques (p = 1 mod 4, p <= 1500) from results/stepanov_balanced_2026_09_29_data.json
     and results/sigma_biclique_2026_09_05_search.json: every row of the family evaluated exactly as a
     theorem; sign-matrix statistics; the exact scaling factor Theta of Theorem 2.1.
  C  Theorem 3.1: exact rational certificates -- one-sided refined profiles for A and B, an explicit sign
     matrix, a joint coupling -- every row re-checked exactly (the middle KW rows, 12 < t < M - 12, in
     double precision with an allowance, as stated in the note);  C2: one-sided clique certificates.
  D  the tight family A = {0,1}: column statistics of the sign matrix for every p = 1 mod 4 up to 3000.

Run with no argument to execute everything and write the results file (about 100 s); an argument such as
`C+` runs only the listed sections (and the first two certificates) without writing.

Only the standard library, numpy, scipy (LP solver, floating point; every certificate is re-checked
exactly with Fractions / Python integers) and sympy (primitive roots) are used.
"""
import cmath
import itertools
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr
from math import comb, isqrt

import numpy as np

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'results', 'stepanov_quartic_2026_09_29.json')

CHECKS = {}
FAILS = []
RESULTS = {}


def check(name, cond, info=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not cond:
        FAILS.append((name, repr(info)[:400]))
    return cond


def log(*a):
    print('[%7.1fs]' % (time.time() - T0), *a, flush=True)


# ============================================================ field, characters, Z[zeta_K]

_IND = {}


def ind_table(p):
    """generator g and discrete-log table ind[x] (ind[0] = -1)."""
    if p in _IND:
        return _IND[p]
    from sympy import primitive_root
    g = int(primitive_root(p))
    ind = np.full(p, -1, dtype=np.int64)
    x = 1
    for k in range(p - 1):
        ind[x] = k
        x = x * g % p
    _IND[p] = (g, ind)
    return g, ind


def zred(v, K):
    """element sum_j v_j zeta^j (0 <= j < K) of Z[zeta_K], K = 2^s >= 2, reduced with zeta^{K/2} = -1."""
    h = K // 2
    return tuple(v[j] - v[j + h] for j in range(h))


def zpow(j, K):
    v = [0] * K
    v[j % K] = 1
    return zred(v, K)


def zmul(u, v, K):
    w = [0] * K
    for i, a in enumerate(u):
        if a:
            for j, b in enumerate(v):
                if b:
                    w[i + j] += a * b
    return zred(w, K)


def zadd(u, v):
    return tuple(a + b for a, b in zip(u, v))


def zscale(u, s):
    return tuple(a * s for a in u)


def zconj(u, K):
    w = [0] * K
    for j, a in enumerate(u):
        w[(-j) % K] += a
    return zred(w, K)


def zcx(u, K):
    return sum(complex(float(a)) * cmath.exp(2j * math.pi * j / K) for j, a in enumerate(u))


def zint(n, K):
    return tuple([n] + [0] * (K // 2 - 1))


def zf(counts, j, K):
    """f_{phi^j} of a class-count vector: sum_c counts[c] zeta^{j c}."""
    v = [0] * K
    for c, nc in enumerate(counts):
        if nc:
            v[(j * c) % K] += nc
    return zred(v, K)


class Field:
    """p = 1 mod K, K in {2,4,8}; phi = character of order K with phi(g) = zeta_K."""

    def __init__(self, p, K):
        assert (p - 1) % K == 0
        self.p, self.K, self.d = p, K, (p - 1) // 2
        self.g, self.ind = ind_table(p)
        w = np.arange(2, p, dtype=np.int64)
        iw, i1w = self.ind[w], self.ind[(1 - w) % p]
        self.J = {}
        for a in range(1, K):
            for b in range(1, K):
                e = (a * iw + b * i1w) % K
                cnt = np.bincount(e, minlength=K)
                self.J[(a, b)] = zred([int(c) for c in cnt], K)
        im1 = (p - 1) // 2
        # K(a,b) := phi^a(-1) J(phi^a, phi^b); sum_x phi^a(x+u) phi^b(x+v) = K(a,b) phi^{a+b}(v-u), u != v
        self.Kc = {ab: zmul(zpow(ab[0] * im1, K), J, K) for ab, J in self.J.items()}

    def cls(self, x):
        """class of x != 0 modulo K (ind mod K); -1 for x = 0."""
        x %= self.p
        return -1 if x == 0 else int(self.ind[x] % self.K)


def class_counts(X, p, K, ind, xs=None):
    """for every x (or x in xs): delta_x = [x in -X] and counts n_c(x) = #{a in X: a+x != 0, class(a+x) = c}."""
    if xs is None:
        xs = np.arange(p, dtype=np.int64)
    Xa = np.array(sorted(X), dtype=np.int64)
    S = (xs[:, None] + Xa[None, :]) % p
    cl = np.where(S == 0, -1, ind[S] % K)
    delta = (S == 0).sum(axis=1)
    counts = np.stack([(cl == c).sum(axis=1) for c in range(K)], axis=1)
    return delta, counts


def star_row_lhs2(ebad, delta, m, e):
    """2 x LHS of (star) at level e: sum over x with ebad <= e of (e+1-ebad)(2m - 3e - ebad - 2 delta)."""
    sel = ebad <= e
    return int(((e + 1 - ebad[sel]) * (2 * m - 3 * e - ebad[sel] - 2 * delta[sel])).sum())


def dk(p, k):
    return (p - 1) // k


# ============================================================ Section A: identities

def section_A():
    log('Section A: identities')
    rng = random.Random(20260929)
    # ---- A1 + A2: Jacobi second moments, |J|^2 = p
    jac_primes = [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137]
    for p in jac_primes:
        for K in (4, 8):
            if (p - 1) % K:
                continue
            F = Field(p, K)
            for (a, b), J in F.J.items():
                if (a + b) % K:
                    check('A2.|J|^2=p', zmul(J, zconj(J, K), K) == zint(p, K), (p, K, a, b))
                else:
                    check('A2.K(a,-a)=-1', F.Kc[(a, b)] == zint(-1, K), (p, K, a, b))
            for trial in range(6):
                mX, mY = rng.randint(1, 7), rng.randint(1, 7)
                X = rng.sample(range(p), mX)
                Y = rng.sample(range(p), mY) if trial % 3 else X[:]  # includes X = Y
                if trial % 3 == 1:
                    Y = list(set(Y) | set(X[:2]))  # overlapping sets
                _, cX = class_counts(X, p, K, F.ind)
                _, cY = class_counts(Y, p, K, F.ind)
                G = cX.T.astype(object) @ cY.astype(object)  # G[c,c'] = sum_x n_c^X(x) n_c'^Y(x)
                # classes of differences v - u, u in X, v in Y, u != v
                dcl = [F.cls(v - u) for u in X for v in Y if u != v]
                inter = len(set(X) & set(Y))
                for a in range(1, K):
                    for b in range(1, K):
                        v = [0] * K
                        for c in range(K):
                            for c2 in range(K):
                                v[(a * c + b * c2) % K] += int(G[c, c2])
                        lhs = zred(v, K)
                        if (a + b) % K:
                            w = [0] * K
                            for c in dcl:
                                w[((a + b) * c) % K] += 1
                            rhs = zmul(F.Kc[(a, b)], zred(w, K), K)
                        else:
                            rhs = zint(-len(dcl) + (p - 1) * inter, K)
                        check('A1.jacobi_second_moment', lhs == rhs, (p, K, X, Y, a, b))

    # ---- A3: class rows (star)_k (robust note Prop 2.9), exhaustive at small p, random above
    def star_all(p, k, A, F):
        m = len(A)
        delta, cnt = class_counts(A, p, k, F.ind)
        d_k = dk(p, k)
        for c in range(k):
            ebad = m - delta - cnt[:, c]
            hist = {}
            for eb, dl in zip(ebad.tolist(), delta.tolist()):
                hist[(eb, dl)] = hist.get((eb, dl), 0) + 1
            hist_star_ok(hist, m, d_k, min((m - 1) // 2, d_k // 2), 'A3.star_k_class_row', (p, k, A, c))
            hist_subset_ok(hist, m, d_k, 'A3.subsetHP_k_class_row', (p, k, A, c))

    exh = [(13, 4, 6), (17, 4, 6), (17, 8, 6), (29, 4, 4), (37, 4, 4), (41, 8, 4)]
    nsets = 0
    for p, k, mmax in exh:
        F = Field(p, k)
        for m in range(1, mmax + 1):
            for rest in itertools.combinations(range(1, p), m - 1):
                A = (0,) + rest
                star_all(p, k, A, F)
                nsets += 1
    for p in (61, 73, 89, 97, 101, 113, 137, 193, 257, 401):
        for k in (4, 8):
            if (p - 1) % k:
                continue
            F = Field(p, k)
            for trial in range(8):
                m = rng.randint(2, min(14, (p - 1) // (2 * k) * 2 + 1))
                A = rng.sample(range(p), m)
                star_all(p, k, A, F)
                nsets += 1
    RESULTS['A3_sets_checked'] = nsets

    # ---- A4: pattern Hanson--Petridis (Prop 1.5) with exact degree and orders
    def lagrange_c(A, p):
        cs = []
        for i, a in enumerate(A):
            pr = 1
            for j, b in enumerate(A):
                if j != i:
                    pr = pr * (a - b) % p
            cs.append(pow(pr, p - 2, p))
        return cs

    for p in (13, 17, 29, 37, 41, 53, 61, 73, 89, 97):
        for k in (4, 8):
            if (p - 1) % k:
                continue
            F = Field(p, k)
            d_k = dk(p, k)
            for trial in range(10):
                m = rng.randint(2, 6)
                A = rng.sample(range(p), m)
                D = d_k + m - 1
                c = lagrange_c(A, p)
                pats = {}
                for x in range(p):
                    if any((x + a) % p == 0 for a in A):
                        continue
                    pats.setdefault(tuple(pow((x + a) % p, d_k, p) for a in A), []).append(x)
                for sig, xs in pats.items():
                    cp = [c[i] * pow(sig[i], p - 2, p) % p for i in range(m)]
                    # l(sigma) = min{l : sum_i cp_i a_i^l != 0}
                    l0 = None
                    for l in range(m):
                        if sum(cp[i] * pow(A[i], l, p) for i in range(m)) % p:
                            l0 = l
                            break
                    check('A4.pattern_l_exists', l0 is not None, (p, k, A, sig))
                    # exact degree of F_sigma = -1 + sum_i cp_i (x + a_i)^D: coefficient of x^{D-l}
                    coef = lambda l: comb(D, l) * sum(cp[i] * pow(A[i], l, p) for i in range(m)) % p
                    degF = D - l0
                    check('A4.pattern_degree', coef(l0) != 0 and all(coef(l) == 0 for l in range(l0)), (p, k, A))
                    for x0 in xs:
                        # orders: F(x0) = 0 and F^{(j)}(x0)/(D)_j = sum cp_i (x0+a_i)^{D-j} = 0, 1 <= j <= m-1
                        v0 = (-1 + sum(cp[i] * pow((x0 + A[i]) % p, D, p) for i in range(m))) % p
                        ok = v0 == 0 and all(sum(cp[i] * pow((x0 + A[i]) % p, D - j, p) for i in range(m)) % p == 0
                                             for j in range(1, m))
                        check('A4.pattern_order>=m', ok, (p, k, A, x0))
                    check('A4.pattern_HP', m * len(xs) <= degF, (p, k, A, sig, len(xs), degF))
                    if len(set(sig)) == 1:
                        check('A4.constant_pattern_l=m-1', l0 == m - 1, (p, k, A, sig))

    # ---- A5: Proposition 1.6: W = sum c_a (x+a)^{D_k}; W^{k/2} - 1 has degree exactly d and vanishes to order
    #      >= m - delta at every b whose column of M^(k) is constant; W(b)^{k/2} != 1 when it is not.
    def poly_W(A, p, Dk, c):
        # coefficients of W(x) = sum_i c_i (x + a_i)^Dk, degree <= Dk; returns list low->high
        co = [0] * (Dk + 1)
        for i, a in enumerate(A):
            apow = 1
            # coefficient of x^j is C(Dk, j) a^{Dk-j}
            pw = [1] * (Dk + 1)
            for t in range(1, Dk + 1):
                pw[t] = pw[t - 1] * a % p
            for j in range(Dk + 1):
                co[j] = (co[j] + c[i] * comb(Dk, j) % p * pw[Dk - j]) % p
        while len(co) > 1 and co[-1] == 0:
            co.pop()
        return co

    def pmul(a, b, p):
        r = [0] * (len(a) + len(b) - 1)
        for i, x in enumerate(a):
            if x:
                for j, y in enumerate(b):
                    r[i + j] = (r[i + j] + x * y) % p
        while len(r) > 1 and r[-1] == 0:
            r.pop()
        return r

    def ord_at(P, x0, p, cap):
        # multiplicity of root x0 (up to cap): repeated synthetic division
        o = 0
        Q = P[:]
        while o < cap:
            # evaluate
            v = 0
            for co in reversed(Q):
                v = (v * x0 + co) % p
            if v != 0 or len(Q) <= 1:
                break
            # divide by (x - x0)
            n_ = len(Q) - 1
            R = [0] * n_
            acc = 0
            for i in range(n_, 0, -1):
                acc = (acc * x0 + Q[i]) % p
                R[i - 1] = acc
            Q = R
            o += 1
        return o

    for p in (13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 113):
        for k in (4, 8):
            if (p - 1) % k:
                continue
            F = Field(p, k)
            d_k = dk(p, k)
            for trial in range(6):
                m = rng.randint(2, 5)
                A = rng.sample(range(p), m)
                Dk = d_k + m - 1
                c = lagrange_c(A, p)
                W = poly_W(A, p, Dk, c)
                check('A5.deg_W=d_k', len(W) - 1 == d_k, (p, k, A, len(W) - 1))
                Wp = [1]
                for _ in range(k // 2):
                    Wp = pmul(Wp, W, p)
                Wp[0] = (Wp[0] - 1) % p
                check('A5.deg_W^{k/2}-1=d', len(Wp) - 1 == (p - 1) // 2, (p, k, A))
                for b in range(p):
                    vals = [(b + a) % p for a in A]
                    if any(v != 0 and F.ind[v] % 2 for v in vals):
                        continue  # b not a complete point
                    dl = sum(1 for v in vals if v == 0)
                    cl = set(F.cls(v) for v in vals if v != 0)
                    const = len(cl) <= 1
                    o = ord_at(Wp, b, p, m + 1)
                    if const:
                        check('A5.constant_column=>order>=m-delta', o >= m - dl, (p, k, A, b, o))
                    if dl == 0:
                        check('A5.order>=m=>constant_column', (o >= m) == const, (p, k, A, b, o))

    # ---- A6: Weil fourth moment for characters of order >= 3 (Lemma 1.3), exact integer data
    for p in (29, 37, 41, 53, 61, 73, 89, 97, 101, 113, 137, 193, 257, 401, 1009):
        for K in (4, 8):
            if (p - 1) % K:
                continue
            F = Field(p, K)
            for trial in range(5):
                m = rng.randint(2, 12)
                A = rng.sample(range(p), m)
                _, cnt = class_counts(A, p, K, F.ind)
                for j in range(1, K):
                    if (2 * j) % K == 0:
                        continue  # order 2 (Legendre): not this lemma
                    tot = zint(0, K)
                    for row in cnt.tolist():
                        f = zf(row, j, K)
                        a2 = zmul(f, zconj(f, K), K)
                        tot = zadd(tot, zmul(a2, a2, K))
                    val = zcx(tot, K).real
                    bound = (2 * m * m - m) * p + 3 * (m ** 4 - 2 * m * m + m) * math.sqrt(p)
                    check('A6.weil_fourth_moment', val <= bound * (1 - 1e-12), (p, K, m, j, val, bound))
                    # main-term count: 2m^2 - m tuples with equal multisets
    # ---- A7: Theorem 2.1(iii) on random abstract column statistics (exact integer arithmetic)
    rs = np.random.RandomState(20260929)
    for trial in range(3000):
        k = int(rs.choice([4, 8]))
        p = int(rs.choice([1009, 10009, 40009, 100049]))
        d = (p - 1) // 2
        m = int(rs.randint(4, 61))
        n = int(rs.randint(1, d // m + 1))
        cols = rs.multinomial(m, [2.0 / k] * (k // 2), size=n)
        const = rs.random_sample(n) < rs.random_sample()  # a random fraction of constant columns
        cls_ = rs.randint(0, k // 2, size=n)
        cols[const] = 0
        cols[const, cls_[const]] = m
        d_k = dk(p, k)
        for ci in range(k // 2):
            tb = m - cols[:, ci]
            for e in range(0, (m - 1) // 2 + 1):
                sel = tb <= e
                W = int(m * sel.sum())
                if 2 * W * k <= 4 * (m * n) or W == 0:  # half-weight hypothesis (fraction 2/k), r = 0
                    lhs2 = int(((e + 1 - tb[sel]) * (2 * m - 3 * e - tb[sel])).sum())
                    check('A7.dichotomy_iii', lhs2 <= 2 * (e + 1) * (d_k - e), (k, p, m, n, e))
    # ---- A8: the exact KW engine (Newton identities + reflection) against brute force over subsets U
    for p, K, m in ((41, 8, 6), (37, 4, 7), (73, 8, 5), (29, 4, 6), (41, 8, 9), (37, 4, 10)):
        F = Field(p, K)
        A = rng.sample(range(p), m)
        if m >= 9:
            A = A[:-1] + [(-A[0]) % p] if (-A[0]) % p not in A else A  # include a point of -A structure
        dl, cnt = class_counts(A, p, K, F.ind)
        for j in kw_family(K):
            sums = kw_point_sums(dl, cnt, j, K, m, 3)
            for t, z in sums.items():
                w = [0] * K
                for U in itertools.combinations(A, t):
                    for x in range(p):
                        vals = [(x + a) % p for a in U]
                        if all(v != 0 for v in vals):
                            w[(j * sum(int(F.ind[v]) for v in vals)) % K] += 1
                check('A8.KW_engine_vs_bruteforce', z == zred(w, K), (p, K, j, t))
        # and the float normalized recurrence against the exact values
        prof = {}
        for dl_, row in zip(dl.tolist(), cnt.tolist()):
            prof[(dl_, tuple(row))] = prof.get((dl_, tuple(row)), 0) + 1
        fl = kw_float_rows(list(prof.items()), m, p, K, tmin=0)
        for j in kw_family(K):
            ex = kw_point_sums(dl, cnt, j, K, m, 3)
            for t, z in ex.items():
                if t >= 3:
                    check('A8.KW_float_vs_exact', abs(fl[(j, t)] * comb(m, t) - zcx(z, K)) <= 1e-6 * comb(m, t) * p,
                          (p, K, j, t))
    log('Section A done')




# ============================================================ Section C: exact certificates (Theorem 3.1)
#
# A *type* of a point x with respect to a set X (|X| = m) and the character phi of order K is
# (delta, counts) with delta = [x in -X] and counts[c] = #{a in X : a + x != 0, ind(a+x) = c mod K}.
# A *profile* is a nonnegative rational weight on types.  Every row below is a linear function of
# the profile.  Rotation c -> c + 2 (multiplication of X, x by a square g) acts on types; all
# profiles constructed here are rotation invariant (variables are rotation orbits).

def t_ebad(tp, k, c, K):
    """class-(k, c) bad count: #{a : a + x != 0, ind(a+x) != c mod k}."""
    dl, cnt = tp
    return sum(cnt) - sum(cnt[j] for j in range(K) if j % k == c)


def t_rot(tp, s, K):
    dl, cnt = tp
    return (dl, tuple(cnt[(j - s) % K] for j in range(K)))


def t_orbit(tp, K):
    out = []
    for s in range(0, K, 2):
        r = t_rot(tp, s, K)
        if r not in out:
            out.append(r)
    return tuple(sorted(out))


def t_fchi(tp, m, K):
    dl, cnt = tp
    return m - dl - 2 * sum(cnt[j] for j in range(1, K, 2))


def lcomb(n, k):
    if k < 0 or n < k:
        return None
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


QUAD_IDENT = ('subsetHP2_A_1', 'subsetHP2_A_2', 'subsetHP2_nuA_1', 'subsetHP2_nuA_2')


def quad_rows(p, m, n0, band=0.4):
    """theorem rows of the balanced note (Theorem 3.5 list) over level variables
    I(j) = #{x not in -X : e_x = j}, J(j) = #{x in -X : e_x = j}; plus fixings (B cell = n0, r = 0, band)."""
    d = (p - 1) // 2
    I = lambda j: j
    J = lambda j: m + 1 + j
    val = lambda j, dl: m - dl - 2 * j
    rows = []
    rows.append(({**{I(j): 1 for j in range(m + 1)}, **{J(j): 1 for j in range(m)}}, '=', p, 'total'))
    rows.append(({J(j): 1 for j in range(m)}, '=', m, 'negA'))
    rows.append(({**{I(j): val(j, 0) for j in range(m + 1)}, **{J(j): val(j, 1) for j in range(m)}}, '=', 0,
                 'moment1'))
    rows.append(({**{I(j): val(j, 0) ** 2 for j in range(m + 1)}, **{J(j): val(j, 1) ** 2 for j in range(m)}}, '=',
                 m * (p - m), 'moment2'))
    for side in ('A', 'nuA'):
        ee0 = (lambda j: j) if side == 'A' else (lambda j: m - j)
        ee1 = (lambda j: j) if side == 'A' else (lambda j: m - 1 - j)
        for e in range(0, (m - 1) // 2 + 1):
            co = {}
            for j in range(m + 1):
                ee = ee0(j)
                if ee <= e:
                    co[I(j)] = co.get(I(j), 0) + (e + 1 - ee) * (2 * m - (3 * e + ee)) / 2
            for j in range(m):
                ee = ee1(j)
                if ee <= e:
                    co[J(j)] = co.get(J(j), 0) + (e + 1 - ee) * (2 * m - (3 * e + ee) - 2) / 2
            rows.append((co, '<=', (e + 1) * (d - e), 'star2_%s_%d' % (side, e)))
        for t in range(1, m + 1):
            co = {}
            for j in range(m + 1):
                ee = ee0(j)
                if m - 1 - ee >= t - 1:
                    co[I(j)] = co.get(I(j), 0) + (m - ee) * math.exp(lcomb(m - 1 - ee, t - 1) - lcomb(m, t))
            for j in range(m):
                ee = ee1(j)
                if m - 1 - ee >= t - 1 and m - ee - 1 > 0:
                    co[J(j)] = co.get(J(j), 0) + (m - ee - 1) * math.exp(lcomb(m - 1 - ee, t - 1) - lcomb(m, t))
            rows.append((co, '<=', d, 'subsetHP2_%s_%d' % (side, t)))
    rows.append(({I(0): 2 * m - 2, I(1): m - 2}, '<=', 2 * d - 2, 'pencil_A'))
    rows.append(({I(m): 2 * m - 2, I(m - 1): m - 2}, '<=', 2 * d - 2, 'pencil_nuA'))
    rows.append(({I(0): m - 1, I(m): m - 1, J(0): m - 2, J(m - 1): m - 2}, '<=', d - 1, 'derivHP'))
    sq = isqrt(p)  # floor(sqrt p) <= sqrt p: every Weil row becomes stricter
    for k in (2, 3, 4):
        df = 1
        for i in range(1, 2 * k, 2):
            df *= i
        rhs = df * m ** k * p + (2 * k - 1) * m ** (2 * k) * sq
        rows.append(({**{I(j): float(val(j, 0) ** (2 * k)) for j in range(m + 1)},
                      **{J(j): float(val(j, 1) ** (2 * k)) for j in range(m)}}, '<=', float(rhs), 'weil2_moment_%d' % (2 * k)))
    for j in range(m + 1):
        co = {I(j): 1}
        if j < m:
            co[J(j)] = 1
        main = math.exp(lcomb(m, j) + math.log(p) - m * math.log(2))
        err = math.exp(lcomb(m, j) + math.log(m * (sq + 1) / 2)) + 2 * m
        if main - err > 0:
            rows.append(({k: -v for k, v in co.items()}, '<=', -(main - err), 'weil2_count_lo_%d' % j))
        if main + err < p:
            rows.append((co, '<=', main + err, 'weil2_count_up_%d' % j))
    fix = [({I(0): 1}, '=', n0, 'fix_B'), ({J(0): 1}, '=', 0, 'fix_r=0')]
    lo, hi = math.ceil(band * m), math.floor((1 - band) * m)
    for j in range(1, m + 1):
        if j < lo or j > hi:
            fix.append(({I(j): 1}, '=', 0, 'fix_band_%d' % j))
    return rows, fix, 2 * m + 1


def lp_solve_cuts(nv, rows, margin_pred, cutgen, rounds=40, Q=10 ** 6):
    """lp_solve with cutting planes: cutgen(x_float) returns new '<=' rows violated (or nearly) by x."""
    rows = list(rows)
    for it in range(rounds):
        x, mg = lp_solve(nv, rows, margin_pred, Q=Q)
        if x is None:
            return None, mg, rows
        new = cutgen([float(v) for v in x])
        if not new:
            return x, mg, rows
        rows.extend(new)
    return None, 'cutting planes did not converge', rows


def lp_solve(nv, rows, margin_pred, Q=10 ** 6):
    """maximise a common relative slack t (0 <= t <= 1) on the '<=' rows selected by margin_pred,
    subject to all rows; round to denominator Q and repair every '=' row exactly (pivots chosen
    greedily among the largest support variables).  Returns (exact solution, margin) or (None, msg)."""
    from scipy.optimize import linprog
    Aeq, beq, Aub, bub = [], [], [], []
    for co, sense, rhs, name in rows:
        v = np.zeros(nv + 1)
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)), 1.0])
        for k, c in co.items():
            v[k] += float(c) / sc
        if sense == '=':
            Aeq.append(v); beq.append(float(rhs) / sc)
        else:
            if margin_pred(name):
                v[nv] = 1.0
            Aub.append(v); bub.append(float(rhs) / sc)
    c = np.zeros(nv + 1); c[nv] = -1
    res = linprog(c, A_ub=np.array(Aub) if Aub else None, b_ub=np.array(bub) if Aub else None,
                  A_eq=np.array(Aeq), b_eq=np.array(beq), bounds=[(0, None)] * nv + [(0, 1)], method='highs')
    if res.status != 0:
        return None, res.message
    x = res.x[:nv]
    xr = [Fr(round(v * Q), Q) if v > 1e-9 else Fr(0) for v in x]
    eqs = [(co, rhs) for co, s, rhs, nm in rows if s == '=']
    Mrows = [[Fr(co.get(i, 0)) for i in range(nv)] + [Fr(rhs)] for co, rhs in eqs]
    R = [r[:] for r in Mrows]
    order = sorted([i for i in range(nv) if x[i] > 1e-7], key=lambda i: -x[i])
    piv, used = [], []
    for col in order:
        prow = next((ri for ri in range(len(R)) if ri not in used and R[ri][col] != 0), None)
        if prow is None:
            continue
        used.append(prow); piv.append(col)
        pv = R[prow][col]
        for ri in range(len(R)):
            if ri != prow and R[ri][col] != 0:
                f = R[ri][col] / pv
                R[ri] = [a - f * b for a, b in zip(R[ri], R[prow])]
        if len(used) == len(R):
            break
    free = [i for i in range(nv) if i not in piv]
    # exact solve of the pivot system (original rows `used`) by Fraction Gauss-Jordan
    Amat = [[Mrows[ri][ci] for ci in piv] + [Mrows[ri][-1] - sum(Mrows[ri][i] * xr[i] for i in free if xr[i])]
            for ri in used]
    nn = len(piv)
    for col in range(nn):
        pr = next(r for r in range(col, nn) if Amat[r][col] != 0)
        Amat[col], Amat[pr] = Amat[pr], Amat[col]
        pv = Amat[col][col]
        Amat[col] = [a / pv for a in Amat[col]]
        for r in range(nn):
            if r != col and Amat[r][col] != 0:
                f = Amat[r][col]
                Amat[r] = [a - f * b for a, b in zip(Amat[r], Amat[col])]
    for k, ci in enumerate(piv):
        xr[ci] = Amat[k][-1]
    return xr, float(res.x[-1])


def rows_ok(xr, rows):
    bad = []
    for co, sense, rhs, name in rows:
        lhs = sum(Fr(c) * xr[k] for k, c in co.items() if xr[k])
        if (sense == '=' and lhs != rhs) or (sense == '<=' and lhs > rhs):
            bad.append(name)
    return bad


# ------------------------------------------------ stage 2: refinement of levels into order-K class types

def orbit_features(orb, m, K):
    """rotation-orbit averages of |f_j|^2, |f_j|^4 (1 <= j < K/2) and f_i f_j (i + j = K/2)."""
    h = K // 2
    acc = {}
    for tp in orb:
        f = {j: zf(tp[1], j, K) for j in range(1, K)}
        for j in range(1, h):
            a2 = zmul(f[j], zconj(f[j], K), K)
            for key, val in ((('abs2', j), a2), (('abs4', j), zmul(a2, a2, K))):
                acc[key] = zadd(acc.get(key, zint(0, K)), val)
        for i in range(1, h):
            if i <= h - i:
                key = ('jac', i, h - i)
                acc[key] = zadd(acc.get(key, zint(0, K)), zmul(f[i], f[h - i], K))
    L = len(orb)
    return {k: tuple(Fr(x, L) for x in v) for k, v in acc.items()}


def class_coefs(tp, m, K, p):
    """coefficients of one type in the class rows (star_k, subsetHP_k), k in {4, 8} dividing K."""
    dl, cnt = tp
    out = {}
    for k in (4, 8):
        if K % k:
            continue
        d_k = (p - 1) // k
        emax = min((m - 1) // 2, d_k // 2)
        for c in range(k):
            eb = t_ebad(tp, k, c, K)
            for e in range(eb, emax + 1):
                out['star%d_c%d_e%d' % (k, c, e)] = (e + 1 - eb) * (2 * m - 3 * e - eb - 2 * dl) / 2
            if m - eb - dl > 0:
                for t in range(1, m - eb + 1):
                    out['sub%d_c%d_t%d' % (k, c, t)] = (m - eb - dl) * math.exp(lcomb(m - 1 - eb, t - 1) - lcomb(m, t))
    return out


def class_rhs(name, p):
    head = name.split('_')[0]
    k = int(head[4:]) if head.startswith('star') else int(head[3:])
    d_k = (p - 1) // k
    if head.startswith('star'):
        e = int(name.split('_')[2][1:])
        return (e + 1) * (d_k - e)
    return d_k


def refine_side(F, m, levels, fixed, S=40, seed=1):
    """levels: (delta, e) -> count (cells neg and rest); fixed: [(type, count)] for the complete cell.
    Candidate types: S seeded multinomial samples per level (bad partners spread uniformly over the
    odd classes, good ones over the even classes), as rotation orbits.  Returns [(cell, orbit, count)]."""
    p, K = F.p, F.K
    h = K // 2
    rs = np.random.RandomState(seed)
    cands = []
    for (dl, e), N in sorted(levels.items()):
        seen = set()
        for s in range(S):
            ev = rs.multinomial(m - dl - e, [1.0 / h] * h)
            od = rs.multinomial(e, [1.0 / h] * h)
            cnt = [0] * K
            for i in range(h):
                cnt[2 * i], cnt[2 * i + 1] = int(ev[i]), int(od[i])
            orb = t_orbit((dl, tuple(cnt)), K)
            if orb not in seen:
                seen.add(orb)
                cands.append(((dl, e), orb))
    nv = len(cands)
    sigma = sum(N * (m - dl - 2 * e) for (dl, e), N in levels.items() if dl == 1)
    const = {}
    for tp, cn in fixed:
        for key, v in orbit_features((tp,), m, K).items():
            const[key] = zadd(const.get(key, zint(0, K)), zscale(v, cn))
    ofs = [orbit_features(orb, m, K) for lv, orb in cands]
    rows = []
    for lv, N in sorted(levels.items()):
        rows.append(({i: 1 for i, (l2, o) in enumerate(cands) if l2 == lv}, '=', N, 'level_%d_%d' % lv))
    for j in range(1, h):
        key, tgt = ('abs2', j), zint(m * (p - m), K)
        for co in range(h):
            cf = {i: ofs[i][key][co] for i in range(nv) if ofs[i][key][co]}
            rhs = tgt[co] - const.get(key, zint(0, K))[co]
            if cf or rhs:
                rows.append((cf, '=', rhs, 'abs2_%d_c%d' % (j, co)))
    for i_ in range(1, h):
        if i_ > h - i_:
            continue
        key, tgt = ('jac', i_, h - i_), zscale(F.Kc[(i_, h - i_)], sigma)
        for co in range(h):
            cf = {i: ofs[i][key][co] for i in range(nv) if ofs[i][key][co]}
            rhs = tgt[co] - const.get(key, zint(0, K))[co]
            if cf or rhs:
                rows.append((cf, '=', rhs, 'jac_%d_%d_c%d' % (i_, h - i_, co)))
    bnd = (2 * m * m - m) * p + 3 * m ** 4 * isqrt(p)
    for j in range(1, h):
        key = ('abs4', j)
        rows.append(({i: zcx(ofs[i][key], K).real for i in range(nv)}, '<=',
                     bnd - zcx(const.get(key, zint(0, K)), K).real, 'weil4_%d' % j))
    crow, cconst = {}, {}
    for i, (lv, orb) in enumerate(cands):
        for tp in orb:
            for nm, v in class_coefs(tp, m, K, p).items():
                crow.setdefault(nm, {})
                crow[nm][i] = crow[nm].get(i, 0) + v / len(orb)
    for tp, cn in fixed:
        for nm, v in class_coefs(tp, m, K, p).items():
            cconst[nm] = cconst.get(nm, 0) + v * cn
    for nm in sorted(set(crow) | set(cconst)):
        rows.append((crow.get(nm, {}), '<=', class_rhs(nm, p) - cconst.get(nm, 0), nm))
    js = [j for j in kw_family(K) if 2 * j != K]
    E = {}
    for i, (lv, orb) in enumerate(cands):
        E[i] = {j: sum(ehat(tp[1], j, K) * kw_ratio(sum(tp[1]), m) for tp in orb) / len(orb) for j in js}
        E[i] = {j: np.concatenate([v, np.zeros(m + 1 - len(v))]) for j, v in E[i].items()}
    const = {j: np.zeros(m + 1, dtype=complex) for j in js}
    for tp, cn in fixed:
        for j in js:
            v = ehat(tp[1], j, K) * kw_ratio(sum(tp[1]), m)
            const[j][:len(v)] += float(cn) * v
    y, mg, rows = lp_solve_cuts(nv, rows, lambda nm: not nm.endswith('_t1'), kw_cutgen(E, const, m, p, K, js))
    if y is None or min(y) < 0 or rows_ok(y, [r for r in rows if r[1] == '=']):
        return None, mg
    byorb = {}
    for tp, cn in fixed:
        byorb.setdefault(t_orbit(tp, K), []).append((tp, cn))
    out = []
    for orb, lst in sorted(byorb.items()):
        cnts = dict(lst)
        assert all(cnts.get(tp, 0) == cnts[orb[0]] for tp in orb), 'complete cell not rotation invariant'
        out.append(('comp', orb, Fr(sum(cnts.values()))))
    for i, (lv, orb) in enumerate(cands):
        if y[i]:
            out.append(('neg' if lv[0] == 1 else 'rest', orb, y[i]))
    return out, mg


# ------------------------------------------------ the explicit sign matrix

def sign_matrix(K, m, n, seed):
    """M[a][b] = class index (even) of phi_K(a+b): base 2*((a+b) mod h), h = K/2, then seeded 2x2
    switches, which preserve every row and column class count."""
    h = K // 2
    M = [[2 * ((a + b) % h) for b in range(n)] for a in range(m)]
    rng = random.Random(seed)
    for _ in range(6 * m * n):
        a, a2 = rng.randrange(m), rng.randrange(m)
        b, b2 = rng.randrange(n), rng.randrange(n)
        u, v = M[a][b], M[a][b2]
        if a != a2 and b != b2 and u != v and M[a2][b2] == u and M[a2][b] == v:
            M[a][b], M[a][b2], M[a2][b], M[a2][b2] = v, u, u, v
    return M


def comp_orbits(M, K):
    """complete-cell profiles: columns (A side) and rows (B side) as types, aggregated by rotation orbit."""
    m, n = len(M), len(M[0])
    colt, rowt = {}, {}
    for b in range(n):
        cnt = [0] * K
        for a in range(m):
            cnt[M[a][b]] += 1
        tp = (0, tuple(cnt))
        colt[tp] = colt.get(tp, 0) + 1
    for a in range(m):
        cnt = [0] * K
        for b in range(n):
            cnt[M[a][b]] += 1
        tp = (0, tuple(cnt))
        rowt[tp] = rowt.get(tp, 0) + 1
    return colt, rowt


# ------------------------------------------------ the joint (two-sided) stage

def pair_orbits(oA, oB, K):
    tA = oA[0]
    res = set()
    for tB in oB:
        res.add(tuple(sorted(set((t_rot(tA, s, K), t_rot(tB, s, K)) for s in range(0, K, 2)))))
    return sorted(res)


def tgrid_for(m, n, full=30):
    if m <= full and n <= full:
        return [(a, b) for a in range(m + 1) for b in range(n + 1) if a + b >= 1]
    sa = sorted(set([0, 1, 2, 3, 4, 5, 6, 8, 10, m // 8, m // 4, m // 2, (3 * m) // 4, m]))
    sb = sorted(set([0, 1, 2, 3, 4, 5, 6, 8, 10, n // 8, n // 4, n // 2, (3 * n) // 4, n]))
    return [(a, b) for a in sa for b in sb if a + b >= 1]


def union_ident(nm):
    """union subset rows that are identities (t = 1 for every order; t = 2 for the quadratic class)."""
    if nm.endswith('_t1_0') or nm.endswith('_t0_1'):
        return True
    return nm.startswith('U_sub2_') and (nm.endswith('_t2_0') or nm.endswith('_t0_2') or nm.endswith('_t1_1'))


def union_float_coefs(mem, m, n, K, p, tg):
    from scipy.special import gammaln
    def lc(N, k):
        N = np.asarray(N, dtype=float); k = np.asarray(k, dtype=float)
        ok = (k >= 0) & (N >= k)
        v = gammaln(N + 1) - gammaln(np.where(ok, k, 0) + 1) - gammaln(np.where(ok, N - k, 0) + 1)
        return np.where(ok, v, -np.inf)
    out = {}
    M_ = m + n
    L = len(mem)
    ta = np.array([a for a, b in tg], dtype=float)
    tb = np.array([b for a, b in tg], dtype=float)
    base = lc(m, ta) + lc(n, tb)
    names = {}
    for tA, tB in mem:
        for k in (2, 4, 8):
            if K % k:
                continue
            d_k = (p - 1) // k
            for c in range(k):
                eA, eB = t_ebad(tA, k, c, K), t_ebad(tB, k, c, K)
                eU, dU = eA + eB, tA[0] + tB[0]
                for e in range(eU, min((M_ - 1) // 2, d_k // 2) + 1):
                    nm = 'U_star%d_c%d_e%d' % (k, c, e)
                    out[nm] = out.get(nm, 0) + (e + 1 - eU) * (2 * M_ - 3 * e - eU - 2 * dU) / 2 / L
                val = (ta + tb) * np.exp(lc(m - eA, ta) + lc(n - eB, tb) - base)
                if tA[0]:
                    val = val - np.where(ta >= 1, np.exp(lc(m - 1 - eA, ta - 1) + lc(n - eB, tb) - base), 0.0)
                if tB[0]:
                    val = val - np.where(tb >= 1, np.exp(lc(m - eA, ta) + lc(n - 1 - eB, tb - 1) - base), 0.0)
                for g_, (a_, b_) in enumerate(tg):
                    if val[g_]:
                        nm = 'U_sub%d_c%d_t%d_%d' % (k, c, a_, b_)
                        out[nm] = out.get(nm, 0) + val[g_] / L
    return out


def joint_stage(F, m, n, sideA, sideB):
    """joint cells: B (A-complete), A (B-complete), -A, -B, other; |A cap B| = 0, r = 0."""
    p, K = F.p, F.K
    h = K // 2
    A = {c: [(o, x) for (cc, o, x) in sideA if cc == c] for c in ('comp', 'neg', 'rest')}
    B = {c: [(o, x) for (cc, o, x) in sideB if cc == c] for c in ('comp', 'neg', 'rest')}
    cells = [('JB', 'comp', 'rest'), ('JA', 'rest', 'comp'), ('JnegA', 'neg', 'rest'), ('JnegB', 'rest', 'neg'),
             ('JO', 'rest', 'rest')]
    var = []
    for cell, ca, cb in cells:
        for ia, (oa, xa) in enumerate(A[ca]):
            for ib, (ob, xb) in enumerate(B[cb]):
                for mem in pair_orbits(oa, ob, K):
                    var.append((cell, (ca, ia), (cb, ib), mem))
    nv = len(var)

    def jfeat(mem):
        acc = {}
        for tA, tB in mem:
            fA = {j: zf(tA[1], j, K) for j in range(1, K)}
            fB = {j: zf(tB[1], j, K) for j in range(1, K)}
            cA, cB = t_fchi(tA, m, K), t_fchi(tB, n, K)
            vals = [(('J1',), zint(cA * cB, K)), (('chiA',), zint(cA, K)), (('chiB',), zint(cB, K))]
            for j in range(1, h):
                vals.append((('J2', j), zmul(fA[j], zconj(fB[j], K), K)))
                vals.append((('J3', j), zmul(fA[j], fB[h - j], K)))
            for key, v in vals:
                acc[key] = zadd(acc.get(key, zint(0, K)), v)
        return {k: tuple(Fr(x, len(mem)) for x in v) for k, v in acc.items()}

    jf = [jfeat(v[3]) for v in var]
    rows = []
    for ca in ('comp', 'neg', 'rest'):
        for ia, (oa, xa) in enumerate(A[ca]):
            rows.append(({i: 1 for i, v in enumerate(var) if v[1] == (ca, ia)}, '=', xa, 'margA_%s_%d' % (ca, ia)))
    for cb in ('comp', 'neg', 'rest'):
        for ib, (ob, xb) in enumerate(B[cb]):
            rows.append(({i: 1 for i, v in enumerate(var) if v[2] == (cb, ib)}, '=', xb, 'margB_%s_%d' % (cb, ib)))
    for cell, sz in (('JA', m), ('JnegB', n)):
        rows.append(({i: 1 for i, v in enumerate(var) if v[0] == cell}, '=', sz, 'size_%s' % cell))
    rows.append(({i: jf[i][('J1',)][0] for i in range(nv) if jf[i][('J1',)][0]}, '=', -m * n, 'J1'))
    for j in range(1, h):
        for co in range(h):
            rows.append(({i: jf[i][('J2', j)][co] for i in range(nv) if jf[i][('J2', j)][co]}, '=',
                         -m * n if co == 0 else 0, 'J2_%d_c%d' % (j, co)))
        KK = F.Kc[(j, h - j)]
        for co in range(h):
            cf = {}
            for i in range(nv):
                v = jf[i][('J3', j)][co]
                if var[i][0] == 'JnegA':
                    v -= KK[co] * jf[i][('chiB',)][0]
                if v:
                    cf[i] = v
            rows.append((cf, '=', 0, 'J3_%d_c%d' % (j, co)))
    cf = {}
    for i in range(nv):
        if var[i][0] == 'JnegA':
            cf[i] = cf.get(i, 0) + jf[i][('chiB',)][0]
        if var[i][0] == 'JnegB':
            cf[i] = cf.get(i, 0) - jf[i][('chiA',)][0]
    rows.append(({i: v for i, v in cf.items() if v}, '=', 0, 'Dchi_consistency'))
    tg = tgrid_for(m, n)
    urow = {}
    for i, v in enumerate(var):
        for nm, c in union_float_coefs(v[3], m, n, K, p, tg).items():
            urow.setdefault(nm, {})[i] = c
    for nm, cf in sorted(urow.items()):
        parts = nm.split('_')
        k = int(parts[1][4:]) if parts[1].startswith('star') else int(parts[1][3:])
        d_k = (p - 1) // k
        rhs = (int(parts[3][1:]) + 1) * (d_k - int(parts[3][1:])) if parts[1].startswith('star') else d_k
        rows.append((cf, '<=', rhs, nm))
    M_U = m + n
    js = kw_family(K)
    Ecache = {}

    def Evec(i):
        if i not in Ecache:
            mem = var[i][3]
            acc = {j: np.zeros(M_U + 1, dtype=complex) for j in js}
            for tA, tB in mem:
                cnt = tuple(a + b for a, b in zip(tA[1], tB[1]))
                for j in js:
                    v = ehat(cnt, j, K) * kw_ratio(sum(cnt), M_U)
                    acc[j][:len(v)] += v / len(mem)
            Ecache[i] = acc
        return Ecache[i]

    zeroc = {j: np.zeros(M_U + 1, dtype=complex) for j in js}

    def gen(x):
        tot = {j: np.zeros(M_U + 1, dtype=complex) for j in js}
        for i, xi in enumerate(x):
            if xi > 0:
                Ei = Evec(i)
                for j in js:
                    tot[j] += xi * Ei[j]
        viol = [(j, t) for j in js for t in range(3, M_U + 1)
                if kw_live(j, t, K) and abs(tot[j][t]) > 0.985 * (t - 1) * math.sqrt(p)]
        if not viol:
            return []
        Eall = {i: Evec(i) for i in range(nv)}
        out = []
        for j, t in viol:
            rot_ = cmath.exp(-1j * cmath.phase(tot[j][t]))
            co = {i: (rot_ * Eall[i][j][t]).real for i in range(nv) if abs(Eall[i][j][t]) > 0}
            out.append((co, '<=', 0.97 * (t - 1) * math.sqrt(p), 'U_KWcut_j%d_t%d' % (j, t)))
        return out

    z, mg, rows = lp_solve_cuts(nv, rows, lambda nm: nm.startswith('U_') and not union_ident(nm), gen)
    if z is None or min(z) < 0:
        return None, mg, tg
    eqbad = rows_ok(z, [r for r in rows if r[1] == '='])
    if eqbad:
        return None, 'equality repair failed: %s' % eqbad[:3], tg
    return [(var[i][0], var[i][3], z[i]) for i in range(nv) if z[i]], mg, tg


# ------------------------------------------------ exact verification of a certificate

def hist_star_ok(hist, M_, d_k, emax, tag, info):
    """hist: (ebad, delta) -> weight.  Checks (star) rows e = 0..emax exactly; returns min relative slack."""
    minsl = None
    items = sorted(hist.items())
    for e in range(0, emax + 1):
        lhs2 = sum(w * (e + 1 - eb) * (2 * M_ - 3 * e - eb - 2 * dl) for (eb, dl), w in items if eb <= e)
        rhs2 = 2 * (e + 1) * (d_k - e)
        check(tag, lhs2 <= rhs2, (info, e, float(lhs2), rhs2))
        sl = float((rhs2 - lhs2) / rhs2)
        minsl = sl if minsl is None else min(minsl, sl)
    return minsl


def hist_subset_ok(hist, M_, d_k, tag, info, tmax=None):
    minsl = None
    for t in range(1, (tmax or M_) + 1):
        lhs = sum(w * comb(M_ - 1 - eb, t - 1) * (M_ - eb - dl) for (eb, dl), w in hist.items()
                  if M_ - 1 - eb >= t - 1 and M_ - eb - dl > 0)
        rhs = comb(M_, t) * d_k
        check(tag, lhs <= rhs, (info, t))
        if t > 1:
            sl = float((rhs - lhs) / rhs)
            minsl = sl if minsl is None else min(minsl, sl)
    return minsl


def weil4_exact_ok(val, bound_int, K):
    """val in Z[zeta_K] real (Q-coefficients); K = 4: integer; K = 8: a + b*sqrt2 with sqrt2 = zeta - zeta^3.
    exact test val <= bound_int."""
    if K == 4:
        return val[1] == 0 and val[0] <= bound_int
    a, b, c2, b3 = val
    if not (c2 == 0 and b3 == -b):
        return False
    diff = bound_int - a  # need b*sqrt2 <= diff
    if b <= 0:
        return diff >= 0 or diff * diff <= 2 * b * b
    return diff >= 0 and diff * diff >= 2 * b * b


def et_values(cnt, j, K, T):
    """[z^t] prod_c (1 + zeta^{j c} z)^{cnt[c]} for t = 0..T, as elements of Z[zeta_K]:
    the elementary symmetric functions of the values phi^j(x + a), a in X (zeros contribute 1)."""
    poly = [zint(1, K)] + [zint(0, K)] * T
    for c, nc in enumerate(cnt):
        if not nc:
            continue
        u = zpow(j * c, K)
        fac = [zint(0, K)] * (T + 1)
        up = zint(1, K)
        for i in range(0, min(nc, T) + 1):
            fac[i] = zscale(up, comb(nc, i))
            up = zmul(up, u, K)
        new = [zint(0, K)] * (T + 1)
        for i1, a1 in enumerate(poly):
            if any(a1):
                for i2 in range(0, T + 1 - i1):
                    if any(fac[i2]):
                        new[i1 + i2] = zadd(new[i1 + i2], zmul(a1, fac[i2], K))
        poly = new
    return poly


def kw_rows_ok(prof, mX, p, K, T, tag, info):
    """Krawtchouk--Weil rows for every character phi^j of the family (1 <= j <= K/2) and 3 <= t <= T:
    |sum_x e_t(phi^j(x+a))_a|^2 <= (C(m,t)(t-1))^2 p  (Weil: sum_{|U|=t} sum_x phi^j(prod_{a in U}(x+a)))."""
    h = K // 2
    worst = 0.0
    for j in range(1, h + 1):
        tot = [zint(0, K)] * (T + 1)
        for tp, w in prof.items():
            ev = et_values(tp[1], j, K, T)
            for t in range(T + 1):
                if any(ev[t]):
                    tot[t] = zadd(tot[t], zscale(ev[t], w))
        for t in range(3, T + 1):
            z2 = zmul(tot[t], zconj(tot[t], K), K)
            R = (comb(mX, t) * (t - 1)) ** 2 * p
            check(tag, weil4_exact_ok(z2, R, K), (info, j, t))
            worst = max(worst, math.sqrt(max(zcx(z2, K).real, 0.0) / R))
    return worst


# ------------------------------------------------ Krawtchouk--Weil rows for every character of the family
#
# For a type (delta, cnt) the multiset of values phi^j(x + a), a in X, a + x != 0, is {zeta^{jc}: mult cnt[c]}.
# e_t = t-th elementary symmetric function; sum_x e_t = sum_{|U| = t} sum_x phi^j(prod_{a in U}(x + a)), so by
# Weil |sum_x e_t| <= C(m, t)(t - 1) sqrt(p) for t >= 2 (t = 2 is also an exact identity).  We use t >= 3.

_HCACHE = {}


def hyper_w(n1, n2):
    key = (n1, n2)
    if key not in _HCACHE:
        from scipy.special import gammaln as lg
        i = np.arange(n1 + 1)[:, None]
        k = np.arange(n2 + 1)[None, :]
        l1 = lg(n1 + 1) - lg(i + 1) - lg(n1 - i + 1)
        l2 = lg(n2 + 1) - lg(k + 1) - lg(n2 - k + 1)
        l3 = lg(n1 + n2 + 1) - lg(i + k + 1) - lg(n1 + n2 - i - k + 1)
        _HCACHE[key] = (np.exp(l1 + l2 - l3), (i + k).ravel())
    return _HCACHE[key]


def hyper_merge(a, b):
    """normalized convolution: if a_i = e_i(V)/C(|V|,i), b_k = e_k(W)/C(|W|,k) then the result is
    e_t(V u W)/C(|V|+|W|, t); every output is a convex combination of products (numerically stable)."""
    n1, n2 = len(a) - 1, len(b) - 1
    W, idx = hyper_w(n1, n2)
    prod = (a[:, None] * b[None, :] * W).ravel()
    re = np.bincount(idx, weights=prod.real, minlength=n1 + n2 + 1)
    im = np.bincount(idx, weights=prod.imag, minlength=n1 + n2 + 1)
    return re + 1j * im


_ECACHE = {}


def ehat(cnt, j, K):
    """normalized elementary symmetric functions e_t / C(M, t), t = 0..M, M = sum(cnt)."""
    key = (cnt, j, K)
    if key not in _ECACHE:
        e = np.ones(1, dtype=complex)
        for c, nc in enumerate(cnt):
            if nc:
                u = cmath.exp(2j * math.pi * ((j * c) % K) / K)
                e = hyper_merge(e, u ** np.arange(nc + 1))
        _ECACHE[key] = e
    return _ECACHE[key]


def kw_ratio(Mx, M):
    """C(Mx, t)/C(M, t) for t = 0..Mx (Mx <= M)."""
    if Mx == M:
        return 1.0
    from scipy.special import gammaln
    ts = np.arange(Mx + 1)
    return np.exp(gammaln(Mx + 1) - gammaln(Mx - ts + 1) - gammaln(M + 1) + gammaln(M - ts + 1))


def kw_cutgen(E, const, M, p, K, js, margin=0.97, trigger=0.985):
    """cutting planes for |sum_i y_i E_i[j][t] + const[j][t]| <= (t-1) sqrt p (normalized by C(M,t))."""
    def gen(x):
        out = []
        for j in js:
            tot = const[j].copy()
            for i, xi in enumerate(x):
                if xi > 0 and i in E:
                    tot += xi * E[i][j]
            for t in range(3, M + 1):
                if not kw_live(j, t, K):
                    continue
                R = (t - 1) * math.sqrt(p)
                z = tot[t]
                if abs(z) > trigger * R:
                    rot_ = cmath.exp(-1j * cmath.phase(z))
                    co = {i: (rot_ * E[i][j][t]).real for i in E if abs(E[i][j][t]) > 0}
                    rhs = margin * R - (rot_ * const[j][t]).real
                    out.append((co, '<=', rhs, 'KWcut_j%d_t%d_%d' % (j, t, len(out))))
        return out
    return gen


def kw_family(K):
    """characters phi^j of the family used in KW rows (one of each conjugate pair): j = 1..K/2."""
    return list(range(1, K // 2 + 1))


def kw_live(j, t, K):
    """rows killed by rotation symmetry (profile invariant under c -> c + 2) have e_t-sums multiplied by
    zeta^{2jt}; they vanish unless 2jt = 0 mod K."""
    return (2 * j * t) % K == 0


def kw_float_rows(prof_float_items, M, p, K, tmin=3):
    """returns {(j, t): complex value of sum_x e_t / C(M, t)} for a float/Fraction profile."""
    out = {}
    for j in kw_family(K):
        acc = np.zeros(M + 1, dtype=complex)
        for tp, w in prof_float_items:
            e = ehat(tp[1], j, K)
            Mx = len(e) - 1
            # e_t(x)/C(M,t) = ehat_t * C(Mx, t)/C(M, t)
            acc[:Mx + 1] += float(w) * e * kw_ratio(Mx, M)
        for t in range(tmin, M + 1):
            out[(j, t)] = acc[t]
    return out


def kw_exact_ends(cnt, j, K, T):
    """exact e_t in Z[zeta_K] for t <= T and, by e_{M-s}(v) = (prod v) conj(e_s(v)), for t >= M - T."""
    M = sum(cnt)
    low = et_values(cnt, j, K, min(T, M))
    pv = zpow(j * sum(c * nc for c, nc in enumerate(cnt)), K)
    out = {t: low[t] for t in range(len(low))}
    for s_ in range(0, min(T, M) + 1):
        out[M - s_] = zmul(pv, zconj(low[s_], K), K)
    return out


def kw_verify(prof, M, p, K, tag, info, T=12):
    """every KW row 3 <= t <= M, every character phi^j (1 <= j <= K/2): exact for t <= T and t >= M - T
    (all t if M <= 2T + 2); the middle rows are evaluated in floating point with an error allowance of
    p * 1e-9 on the normalized value (|e_t| / C(M,t) <= 1 for every point); returns the worst ratio."""
    worst_exact, worst_float = 0.0, 0.0
    Tn = T if M > 2 * T + 2 else M
    ends = {tp: {j: kw_exact_ends(tp[1], j, K, Tn) for j in kw_family(K)} for tp in prof}
    exact_ts = sorted(set(t for t in range(3, M + 1) if t <= Tn or t >= M - Tn))
    for j in kw_family(K):
        for t in exact_ts:
            z = zint(0, K)
            for tp, w in prof.items():
                v = ends[tp][j].get(t)
                if v is not None and any(v):
                    z = zadd(z, zscale(v, w))
            z2 = zmul(z, zconj(z, K), K)
            R = (comb(M, t) * (t - 1)) ** 2 * p
            check(tag + '_exact', weil4_exact_ok(z2, R, K), (info, j, t))
            worst_exact = max(worst_exact, math.sqrt(max(zcx(z2, K).real, 0.0) / R))
    mid = [t for t in range(3, M + 1) if t not in set(exact_ts)]
    if mid:
        vals = kw_float_rows(list(prof.items()), M, p, K)
        for j in kw_family(K):
            for t in mid:
                v = abs(vals[(j, t)]) + p * 1e-9
                ok = v <= (t - 1) * math.sqrt(p) * (1 - 1e-9)
                check(tag + '_float', ok, (info, j, t, v))
                worst_float = max(worst_float, v / ((t - 1) * math.sqrt(p)))
    return worst_exact, worst_float


def kraw_col(e, M, tmax):
    """K_t(e; M) = e_t of a +-1 vector of length M with e minus signs, t = 0..tmax (exact recurrence)."""
    row = [1, M - 2 * e]
    for t in range(1, tmax):
        nxt = (M - 2 * e) * row[t] - (M - t + 1) * row[t - 1]
        row.append(nxt // (t + 1))
    return row[:tmax + 1]


def kw2_cutgen(p, m):
    allrows = kw2_lp_rows(p, m)

    added = set()

    def gen(x):
        out = []
        for co, sense, rhs, nm in allrows:
            if nm in added:
                continue
            v = sum(c * x[k] for k, c in co.items() if x[k])
            if v > 0.995 * rhs:
                added.add(nm)
                out.append((co, sense, rhs * 0.99, nm))
        return out
    return gen


def kw2_lp_rows(p, m):
    """quadratic KW rows over the level variables of quad_rows (float coefficients, exact data)."""
    I = lambda j: j
    J = lambda j: m + 1 + j
    cols0 = [kraw_col(j, m, m) for j in range(m + 1)]
    cols1 = [kraw_col(j, m - 1, m - 1) for j in range(m)]
    rows = []
    sp = math.sqrt(p)
    for t in range(3, m + 1):
        C = comb(m, t)
        co = {I(j): cols0[j][t] / C for j in range(m + 1) if cols0[j][t]}
        if t <= m - 1:
            for j in range(m):
                if cols1[j][t]:
                    co[J(j)] = cols1[j][t] / C
        rows.append((co, '<=', (t - 1) * sp, 'KW2_up_%d' % t))
        rows.append(({k: -v for k, v in co.items()}, '<=', (t - 1) * sp, 'KW2_lo_%d' % t))
    return rows


def _zv_mul(u, v, K):
    h = K // 2
    out = np.zeros_like(u)
    for a in range(h):
        for b in range(h):
            e = a + b
            if e < h:
                out[:, e] += u[:, a] * v[:, b]
            else:
                out[:, e - h] -= u[:, a] * v[:, b]
    return out


def _zv_conj(u, K):
    h = K // 2
    out = np.zeros_like(u)
    for jj in range(h):
        e = (-jj) % K
        if e < h:
            out[:, e] += u[:, jj]
        else:
            out[:, e - h] -= u[:, jj]
    return out


def _zv_unit(E, K):
    """rows zeta^{E_x}"""
    h = K // 2
    out = np.zeros((len(E), h), dtype=object)
    for x, e in enumerate(E):
        e %= K
        if e < h:
            out[x, e] = 1
        else:
            out[x, e - h] = -1
    return out


def kw_point_sums(delta, cnt, j, K, m, T):
    """exact sum over points x of e_t(phi^j(x + a))_a for t <= T and for t >= m - T (reflection), from the
    per-point class counts; Newton's identities in Z[zeta_K] (exact integer arithmetic, dtype=object)."""
    h = K // 2
    N = cnt.shape[0]
    cnt = cnt.astype(object)
    P = []
    for i in range(1, T + 1):
        v = np.zeros((N, h), dtype=object)
        for c in range(K):
            e = (j * c * i) % K
            if e < h:
                v[:, e] += cnt[:, c]
            else:
                v[:, e - h] -= cnt[:, c]
        P.append(v)
    E = [np.zeros((N, h), dtype=object)]
    E[0][:, 0] = 1
    for t in range(1, T + 1):
        acc = np.zeros((N, h), dtype=object)
        for i in range(1, t + 1):
            term = _zv_mul(E[t - i], P[i - 1], K)
            acc = acc + term if i % 2 == 1 else acc - term
        assert all(x % t == 0 for x in acc.ravel())
        E.append(acc // t)
    out = {}
    for t in range(0, T + 1):
        out[t] = tuple(int(x) for x in E[t].sum(axis=0))
    # top end: point with M = m - delta nonzero coordinates: e_{M-s} = (prod v) conj(e_s)
    Etot = np.array([(j * sum(c * int(cn) for c, cn in enumerate(row))) % K for row in cnt.tolist()])
    U = _zv_unit(Etot, K)
    refl = [_zv_mul(U, _zv_conj(E[s_], K), K) for s_ in range(T + 1)]
    dlt = np.asarray(delta)
    for t in range(max(T + 1, m - T), m + 1):
        # every point contributes e_t = (prod v) conj(e_s) with s = (m - delta) - t, and s <= T here
        tot = zint(0, K)
        for dl in (0, 1):
            s_ = m - dl - t
            if s_ < 0:
                continue  # e_t = 0 for points with fewer than t nonzero coordinates
            sel = (dlt == dl)
            if sel.any():
                tot = zadd(tot, tuple(int(x) for x in refl[s_][sel].sum(axis=0)))
        out[t] = tot
    return out


def verify_side(F, prof, mX, nX, label, rep):
    """every one-sided row of the family, exactly, on a profile prof: type -> weight."""
    p, K, d = F.p, F.K, F.d
    h = K // 2
    info = (label, p, K, mX)
    tot = sum(prof.values())
    check('C.side.total=p', tot == p, info)
    check('C.side.delta_count=m', sum(w for tp, w in prof.items() if tp[0] == 1) == mX, info)
    comp = sum(w for tp, w in prof.items() if tp[0] == 0 and t_ebad(tp, 2, 0, K) == 0)
    check('C.side.complete_points=n', comp == nX, info)
    check('C.side.r=0', all(not (tp[0] == 1 and t_ebad(tp, 2, 0, K) == 0) for tp in prof), info)
    check('C.side.nonneg', all(w >= 0 for w in prof.values()), info)
    # quadratic rows of the balanced-note list not covered below ((star)_2 and subset HP_2 for A and
    # nu A are the classes k = 2, c = 0, 1 of the class loop; the moment identities are Jacobi pairs)
    lev = {}
    for tp, w in prof.items():
        key = (t_ebad(tp, 2, 0, K), tp[0])
        lev[key] = lev.get(key, 0) + w
    g = lambda e, dl: lev.get((e, dl), 0)
    check('C.side.pencil_A', (2 * mX - 2) * g(0, 0) + (mX - 2) * g(1, 0) <= 2 * d - 2, info)
    check('C.side.pencil_nuA', (2 * mX - 2) * g(mX, 0) + (mX - 2) * g(mX - 1, 0) <= 2 * d - 2, info)
    check('C.side.derivHP', (mX - 1) * (g(0, 0) + g(mX, 0)) + (mX - 2) * (g(0, 1) + g(mX - 1, 1)) <= d - 1, info)
    sq = isqrt(p)
    for kk in (2, 3, 4):
        df = 1
        for i_ in range(1, 2 * kk, 2):
            df *= i_
        lhs = sum(w * (mX - dl - 2 * e) ** (2 * kk) for (e, dl), w in lev.items())
        check('C.side.weil2_moment', lhs <= df * mX ** kk * p + (2 * kk - 1) * mX ** (2 * kk) * sq, (info, kk))
    for j in range(mX + 1):
        Nj = g(j, 0) + g(j, 1)
        main = Fr(comb(mX, j) * p, 2 ** mX)
        err = Fr(comb(mX, j) * mX * (sq + 1), 2) + 2 * mX
        check('C.side.weil2_count', main - err <= Nj <= main + err, (info, j))
    # identities: all pairs (i, j) of nontrivial characters phi^i, phi^j
    fv = {tp: {j: zf(tp[1], j, K) for j in range(1, K)} for tp in prof}
    for j in range(1, K):
        s = zint(0, K)
        for tp, w in prof.items():
            s = zadd(s, zscale(fv[tp][j], w))
        check('C.side.first_moment=0', s == zint(0, K), (info, j))
    for i in range(1, K):
        for j in range(i, K):
            lhs = zint(0, K)
            for tp, w in prof.items():
                lhs = zadd(lhs, zscale(zmul(fv[tp][i], fv[tp][j], K), w))
            if (i + j) % K == 0:
                rhs = zint(mX * (p - mX), K)
            else:
                l = (i + j) % K
                tneg = zint(0, K)
                for tp, w in prof.items():
                    if tp[0] == 1:
                        tneg = zadd(tneg, zscale(zf(tp[1], l, K) if l else zint(sum(tp[1]), K), w))
                rhs = zmul(F.Kc[(i, j)], tneg, K)
            check('C.side.jacobi_identity', lhs == rhs, (info, i, j))
    # class rows for every k | K, k >= 2, every class, and averaged subset HP_k
    minsl = {}
    for k in (2, 4, 8):
        if K % k:
            continue
        d_k = (p - 1) // k
        for c in range(k):
            hist = {}
            for tp, w in prof.items():
                key = (t_ebad(tp, k, c, K), tp[0])
                hist[key] = hist.get(key, 0) + w
            s1 = hist_star_ok(hist, mX, d_k, min((mX - 1) // 2, d_k // 2), 'C.side.star_k_row', (info, k, c))
            s2 = hist_subset_ok(hist, mX, d_k, 'C.side.subsetHP_k_row', (info, k, c))
            minsl['star%d' % k] = min(minsl.get('star%d' % k, 9), s1)
            if s2 is not None:
                minsl['subset%d' % k] = min(minsl.get('subset%d' % k, 9), s2)
    # Weil fourth moments for phi^j of order >= 3
    bnd = (2 * mX * mX - mX) * p + 3 * mX ** 4 * isqrt(p)
    for j in range(1, K):
        if (2 * j) % K == 0:
            continue
        s = zint(0, K)
        for tp, w in prof.items():
            a2 = zmul(fv[tp][j], zconj(fv[tp][j], K), K)
            s = zadd(s, zscale(zmul(a2, a2, K), w))
        check('C.side.weil4_row', weil4_exact_ok(s, bnd, K), (info, j))
    # Krawtchouk--Weil rows (t <= 12; all t when m <= 30)
    rep['max_KW_ratio_%s(exact,float)' % label] = kw_verify(prof, mX, p, K, 'C.side.KW_row', info)
    # mean squares off the complete cell (report)
    off = {tp: w for tp, w in prof.items() if not (tp[0] == 0 and t_ebad(tp, 2, 0, K) == 0)}
    wtot = sum(off.values())
    ms = {'chi': float(sum(w * t_fchi(tp, mX, K) ** 2 for tp, w in off.items()) / wtot)}
    for j in range(1, h):
        ms['phi^%d' % j] = float(sum(w * zcx(zmul(fv[tp][j], zconj(fv[tp][j], K), K), K).real
                                     for tp, w in off.items()) / wtot)
    rep['mean_square_off_complete_%s' % label] = ms
    rep['min_rel_slack_%s' % label] = minsl
    return True


def verify_joint(F, cert, m, n, tg, rep):
    p, K = F.p, F.K
    h = K // 2
    info = (p, K, m, n)
    # expand to member-level joint weights
    jw = {}
    for cell, mem, z in cert:
        for pr in mem:
            key = (cell, pr)
            jw[key] = jw.get(key, 0) + z / len(mem)
    check('C.joint.nonneg', all(w >= 0 for w in jw.values()), info)
    sizes = {}
    for (cell, pr), w in jw.items():
        sizes[cell] = sizes.get(cell, 0) + w
    check('C.joint.cell_sizes', sizes.get('JB') == n and sizes.get('JA') == m and sizes.get('JnegA') == m
          and sizes.get('JnegB') == n and sum(sizes.values()) == p, (info, {k: float(v) for k, v in sizes.items()}))
    # cell structure: JB points are A-complete, JA points B-complete; deltas
    for (cell, (tA, tB)), w in jw.items():
        okA = (t_ebad(tA, 2, 0, K) == 0 and tA[0] == 0) == (cell == 'JB')
        okB = (t_ebad(tB, 2, 0, K) == 0 and tB[0] == 0) == (cell == 'JA')
        okd = (tA[0] == 1) == (cell == 'JnegA') and (tB[0] == 1) == (cell == 'JnegB')
        check('C.joint.cell_structure', okA and okB and okd, (info, cell))
    # the joint identities for all pairs (i, j): sum_x f^A_i f^B_j = K(i,j) D_{i+j}  or  p|AcapB| - mn
    D = {}
    DB = {}
    for l in range(1, K):
        s, s2 = zint(0, K), zint(0, K)
        for (cell, (tA, tB)), w in jw.items():
            if cell == 'JnegA':
                s = zadd(s, zscale(zf(tB[1], l, K), w))
            if cell == 'JnegB':
                s2 = zadd(s2, zscale(zf(tA[1], l, K), w))
        D[l], DB[l] = s, s2
        # sum_{x in -B} f^A_l(x) = phi^l(-1) D_l
        check('C.joint.D_consistency', s2 == zmul(zpow(l * ((p - 1) // 2), K), s, K), (info, l))
    for i in range(1, K):
        for j in range(1, K):
            lhs = zint(0, K)
            for (cell, (tA, tB)), w in jw.items():
                lhs = zadd(lhs, zscale(zmul(zf(tA[1], i, K), zf(tB[1], j, K), K), w))
            rhs = zint(-m * n, K) if (i + j) % K == 0 else zmul(F.Kc[(i, j)], D[(i + j) % K], K)
            check('C.joint.jacobi_identity', lhs == rhs, (info, i, j))
    # union rows (A u B, size m + n, disjoint): star_k and joint subset HP, exactly
    M_ = m + n
    minsl = {}
    for k in (2, 4, 8):
        if K % k:
            continue
        d_k = (p - 1) // k
        for c in range(k):
            hist, jh = {}, {}
            for (cell, (tA, tB)), w in jw.items():
                eA, eB = t_ebad(tA, k, c, K), t_ebad(tB, k, c, K)
                key = (eA + eB, tA[0] + tB[0])
                hist[key] = hist.get(key, 0) + w
                key2 = (eA, eB, tA[0], tB[0])
                jh[key2] = jh.get(key2, 0) + w
            s1 = hist_star_ok(hist, M_, d_k, min((M_ - 1) // 2, d_k // 2), 'C.union.star_k_row', (info, k, c))
            minsl['U_star%d' % k] = min(minsl.get('U_star%d' % k, 9), s1)
            for (ta, tb) in tg:
                lhs = 0
                for (eA, eB, dA, dB), w in jh.items():
                    v = (ta + tb) * (comb(m - eA, ta) if m - eA >= ta else 0) * (comb(n - eB, tb) if n - eB >= tb else 0)
                    if dA and ta >= 1 and m - 1 - eA >= ta - 1 and n - eB >= tb:
                        v -= comb(m - 1 - eA, ta - 1) * comb(n - eB, tb)
                    if dB and tb >= 1 and m - eA >= ta and n - 1 - eB >= tb - 1:
                        v -= comb(m - eA, ta) * comb(n - 1 - eB, tb - 1)
                    lhs += w * v
                rhs = comb(m, ta) * comb(n, tb) * d_k
                check('C.union.subsetHP_k_row', lhs <= rhs, (info, k, c, ta, tb))
                if not union_ident('U_sub%d_c%d_t%d_%d' % (k, c, ta, tb)):
                    minsl['U_subset%d' % k] = min(minsl.get('U_subset%d' % k, 9), float((rhs - lhs) / rhs))
            if k == 2:
                N0 = hist.get((0, 0), 0)
                N1 = hist.get((1, 0), 0)
                check('C.union.pencil', (2 * M_ - 2) * N0 + (M_ - 2) * N1 <= 2 * F.d - 2, (info, c))
                if c == 0:
                    n0 = sum(w for (eb, dl), w in hist.items() if eb == 0)
                    r0 = sum(w for (eb, dl), w in hist.items() if eb == 0 and dl == 1)
                    nM = sum(w for (eb, dl), w in hist.items() if eb == M_ - dl)
                    rM = sum(w for (eb, dl), w in hist.items() if eb == M_ - 1 and dl == 1)
                    check('C.union.derivHP', (M_ - 1) * (n0 + nM) <= F.d - 1 + r0 + rM, info)
    rep['min_rel_slack_union'] = minsl
    uprof = {}
    for (cell, (tA, tB)), w in jw.items():
        tp = (tA[0] + tB[0], tuple(a + b for a, b in zip(tA[1], tB[1])))
        uprof[tp] = uprof.get(tp, 0) + w
    rep['max_KW_ratio_union(exact,float)'] = kw_verify(uprof, M_, p, K, 'C.union.KW_row', info)
    return jw


def verify_matrix(F, M, m, n, jw, rep):
    """the explicit sign matrix: entries even classes; its column / row types are exactly the complete
    cells; pattern HP_k and the monochromatic-rectangle HP_k bounds for every k | K, k >= 4."""
    p, K = F.p, F.K
    info = (p, K, m, n)
    check('C.matrix.entries_even', all(M[a][b] % 2 == 0 for a in range(m) for b in range(n)), info)
    colt, rowt = comp_orbits(M, K)
    profB, profA = {}, {}
    for (cell, (tA, tB)), w in jw.items():
        if cell == 'JB':
            profB[tA] = profB.get(tA, 0) + w
        if cell == 'JA':
            profA[tB] = profA.get(tB, 0) + w
    # rotation-invariant comparison: orbit totals
    def orbtot(d):
        o = {}
        for tp, w in d.items():
            key = t_orbit(tp, K)
            o[key] = o.get(key, 0) + w
        return o
    check('C.matrix.columns=complete_cell_A', orbtot(colt) == orbtot(profB) and
          {k: Fr(v) for k, v in colt.items()} == profB, info)
    check('C.matrix.rows=complete_cell_B', orbtot(rowt) == orbtot(profA) and
          {k: Fr(v) for k, v in rowt.items()} == profA, info)
    stats = {}
    for k in (4, 8):
        if K % k:
            continue
        d_k = (p - 1) // k
        colpat, rowpat = {}, {}
        for b in range(n):
            key = tuple(M[a][b] % k for a in range(m))
            colpat[key] = colpat.get(key, 0) + 1
        for a in range(m):
            key = tuple(M[a][b] % k for b in range(n))
            rowpat[key] = rowpat.get(key, 0) + 1
        check('C.matrix.pattern_HP_columns', m * max(colpat.values()) <= d_k + m - 1, (info, k))
        check('C.matrix.pattern_HP_rows', n * max(rowpat.values()) <= d_k + n - 1, (info, k))
        for c in range(0, k, 2):
            gc = max(sum(1 for a in range(m) if M[a][b] % k == c) for b in range(n))
            gr = max(sum(1 for b in range(n) if M[a][b] % k == c) for a in range(m))
            check('C.matrix.monochromatic_rectangle_HP', gc * gr <= d_k, (info, k, c, gc, gr))
        stats['k=%d' % k] = {'distinct_column_patterns': len(colpat), 'distinct_row_patterns': len(rowpat),
                             'max_column_pattern_multiplicity': max(colpat.values()),
                             'max_row_pattern_multiplicity': max(rowpat.values())}
    rep['sign_matrix'] = stats


CERT_PARAMS = [(1009, 8, 20, 25), (1013, 4, 24, 21), (10009, 8, 72, 69), (10037, 4, 70, 71),
               (40009, 8, 136, 147), (100049, 8, 224, 223), (100069, 4, 224, 223), (1000033, 8, 708, 706)]


def build_certificate(p, K, m, n, seed=7):
    F = Field(p, K)
    M = sign_matrix(K, m, n, seed)
    colt, rowt = comp_orbits(M, K)
    sides = []
    for (mX, nX, compt) in ((m, n, colt), (n, m, rowt)):
        rows, fix, nv = quad_rows(p, mX, nX)
        xr, mg1, _ = lp_solve_cuts(nv, rows + fix, lambda nm: not nm.startswith('fix') and nm not in QUAD_IDENT,
                                   kw2_cutgen(p, mX))
        assert xr is not None and min(xr) >= 0 and not rows_ok(xr, [r for r in rows + fix if r[1] == '=']), (p, mX, mg1)
        levels = {}
        for j in range(1, mX + 1):
            if xr[j]:
                levels[(0, j)] = xr[j]
        for j in range(mX):
            if xr[mX + 1 + j]:
                levels[(1, j)] = xr[mX + 1 + j]
        side, mg2 = refine_side(F, mX, levels, sorted(compt.items()), seed=seed)
        assert side is not None, (p, mX, mg2)
        sides.append((side, mg1, mg2, levels))
    cert, mg3, tg = joint_stage(F, m, n, sides[0][0], sides[1][0])
    assert cert is not None, (p, m, n, mg3)
    return F, M, sides, cert, mg3, tg


def section_C(params=CERT_PARAMS):
    log('Section C: certificates')
    out = []
    for (p, K, m, n) in params:
        t0 = time.time()
        F, M, sides, cert, mg3, tg = build_certificate(p, K, m, n)
        rep = {'p': p, 'K': K, 'm': m, 'n': n, 'mn': m * n, 'mn_over_p': m * n / p,
               'defect_g=d-mn': F.d - m * n, 'lp_margins': {'stage1_A': sides[0][1], 'stage1_B': sides[1][1],
                                                           'stage2_A': sides[0][2], 'stage2_B': sides[1][2],
                                                           'joint': mg3},
               'K(phi^i,phi^j)': {'%d,%d' % ij: list(v) for ij, v in F.Kc.items()}}
        for (a_, b_), J_ in F.J.items():
            if (a_ + b_) % K:
                check('C.|J|^2=p_at_certificate_prime', zmul(J_, zconj(J_, K), K) == zint(p, K), (p, K, a_, b_))
        rep['band_for_non_complete_points'] = [math.ceil(0.4 * m), math.floor(0.6 * m)]
        jw = verify_joint(F, cert, m, n, tg, rep)
        profA, profB = {}, {}
        for (cell, (tA, tB)), w in jw.items():
            profA[tA] = profA.get(tA, 0) + w
            profB[tB] = profB.get(tB, 0) + w
        verify_side(F, profA, m, n, 'A', rep)
        verify_side(F, profB, n, m, 'B', rep)
        verify_matrix(F, M, m, n, jw, rep)
        rep['levels_A'] = {'%s%d' % ('n' if dl == 0 else "n'", e): float(v) for (dl, e), v in sides[0][3].items()}
        rep['levels_B'] = {'%s%d' % ('n' if dl == 0 else "n'", e): float(v) for (dl, e), v in sides[1][3].items()}
        rep['support_sizes'] = {'A_types': len(profA), 'B_types': len(profB), 'joint_types': len(jw)}
        rep['certificate'] = [[cell, [[list(tA[1]) + [tA[0]], list(tB[1]) + [tB[0]]] for tA, tB in mem], str(z)]
                              for cell, mem, z in cert]
        rep['sign_matrix_rows'] = [''.join(str(x // 2) for x in row) for row in M] if m * n <= 3000 else None
        rep['seconds'] = time.time() - t0
        log('  certificate p=%d K=%d m=%d n=%d mn/p=%.4f margins %s  (%.1fs)' % (
            p, K, m, n, m * n / p, {k: round(v, 4) for k, v in rep['lp_margins'].items()}, rep['seconds']))
        out.append(rep)
    RESULTS['certificates'] = out


# ------------------------------------------------ C2: one-sided clique certificates (B = -A, r = m), K = 4

CLIQUE_PARAMS = [(1013, 21), (1009, 21), (10037, 71), (10009, 69), (100069, 223), (100049, 221), (1013, 23)]  # last: g = 0 (critical)


def clique_matrix(p, m):
    """circulant quartic sign matrix of an abstract m-clique on Z/m: g(t) = class of psi(c' - c), t = c' - c.
    psi(-1) = -1 (p = 5 mod 8): tournament g(-t) = g(t) + 2; psi(-1) = +1: symmetric, needs m = 1 mod 4.
    Every row has (m-1)/2 entries of each sign."""
    h = (m - 1) // 2
    if ((p - 1) // 2) % 4 == 2:
        return {t: (0 if t <= h else 2) for t in range(1, m)}
    assert (m - 1) % 4 == 0
    S = set(range(1, (m - 1) // 4 + 1))
    return {t: (0 if (t in S or (m - t) in S) else 2) for t in range(1, m)}


def clique_certificate(p, m):
    K = 4
    F = Field(p, K)
    band = 0.35 if m < 40 else 0.4
    rows, fix, nv = quad_rows(p, m, 0, band=band)
    fix = [r for r in fix if r[3] not in ('fix_B', 'fix_r=0')] + [({0: 1}, '=', 0, 'fix_n0=0'),
                                                                   ({m + 1: 1}, '=', m, 'fix_clique')]
    xr, mg1, _ = lp_solve_cuts(nv, rows + fix, lambda nm: not nm.startswith('fix') and nm not in QUAD_IDENT,
                               kw2_cutgen(p, m))
    assert xr is not None and min(xr) >= 0, (p, m, mg1)
    levels = {(0, j): xr[j] for j in range(1, m + 1) if xr[j]}
    assert all(xr[m + 1 + j] == 0 for j in range(1, m))
    g = clique_matrix(p, m)
    rowt = [0] * K
    for t, cl in g.items():
        rowt[cl] += 1
    comp = [((1, tuple(rowt)), m)]
    sigma = m * (m - 1)
    cands = []
    for (dl, e), N in sorted(levels.items()):
        E0, E1 = m - dl - e, e
        Rr = isqrt(8 * m) + 1
        for u in range(-min(E0, Rr), min(E0, Rr) + 1):
            if (u - E0) % 2:
                continue
            for v in range(-min(E1, Rr), min(E1, Rr) + 1):
                if (v - E1) % 2 or u * u + v * v > 8 * m:
                    continue
                tp = (dl, ((E0 + u) // 2, (E1 + v) // 2, (E0 - u) // 2, (E1 - v) // 2))
                orb = t_orbit(tp, K)
                if orb[0] == tp:
                    cands.append(((dl, e), orb))
    nvar = len(cands)
    ofs = [orbit_features(orb, m, K) for lv, orb in cands]
    const = {}
    for tp, cn in comp:
        for key, v in orbit_features((tp,), m, K).items():
            const[key] = zadd(const.get(key, zint(0, K)), zscale(v, cn))
    rows2 = []
    for lv, N in sorted(levels.items()):
        rows2.append(({i: 1 for i, (l2, o) in enumerate(cands) if l2 == lv}, '=', N, 'level'))
    key = ('abs2', 1)
    rows2.append(({i: ofs[i][key][0] for i in range(nvar) if ofs[i][key][0]}, '=',
                  m * (p - m) - const.get(key, zint(0, K))[0], 'abs2'))
    key, tgt = ('jac', 1, 1), zscale(F.Kc[(1, 1)], sigma)
    for co in range(2):
        rows2.append(({i: ofs[i][key][co] for i in range(nvar) if ofs[i][key][co]}, '=',
                      tgt[co] - const.get(key, zint(0, K))[co], 'jac_c%d' % co))
    bnd = (2 * m * m - m) * p + 3 * m ** 4 * isqrt(p)
    rows2.append(({i: zcx(ofs[i][('abs4', 1)], K).real for i in range(nvar)}, '<=',
                  bnd - zcx(const.get(('abs4', 1), zint(0, K)), K).real, 'weil4'))
    crow, cconst = {}, {}
    for i, (lv, orb) in enumerate(cands):
        for tp in orb:
            for nm, v in class_coefs(tp, m, K, p).items():
                crow.setdefault(nm, {})
                crow[nm][i] = crow[nm].get(i, 0) + v / len(orb)
    for tp, cn in comp:
        for nm, v in class_coefs(tp, m, K, p).items():
            cconst[nm] = cconst.get(nm, 0) + v * cn
    for nm in sorted(set(crow) | set(cconst)):
        rows2.append((crow.get(nm, {}), '<=', class_rhs(nm, p) - cconst.get(nm, 0), nm))
    E = {i: {1: np.concatenate([v, np.zeros(m + 1 - len(v))])} for i, v in
         ((i, sum(ehat(tp[1], 1, K) * kw_ratio(sum(tp[1]), m) for tp in orb) / len(orb))
          for i, (lv, orb) in enumerate(cands))}
    const = {1: np.zeros(m + 1, dtype=complex)}
    for tp, cn in comp:
        v = ehat(tp[1], 1, K) * kw_ratio(sum(tp[1]), m)
        const[1][:len(v)] += cn * v
    y, mg2, rows2 = lp_solve_cuts(nvar, rows2, lambda nm: not nm.endswith('_t1'), kw_cutgen(E, const, m, p, K, [1]))
    assert y is not None and min(y) >= 0, (p, m, mg2)
    prof = {tp: Fr(cn) for tp, cn in comp}
    for i, (lv, orb) in enumerate(cands):
        if y[i]:
            for tp in orb:
                prof[tp] = prof.get(tp, 0) + y[i] / len(orb)
    return F, g, prof, levels, mg1, mg2


def section_C2(params=CLIQUE_PARAMS):
    log('Section C2: clique certificates (one-sided)')
    out = []
    for (p, m) in params:
        t0 = time.time()
        F, g, prof, levels, mg1, mg2 = clique_certificate(p, m)
        K = 4
        check('C2.|J(psi,psi)|^2=p', zmul(F.J[(1, 1)], zconj(F.J[(1, 1)], K), K) == zint(p, K), p)
        rep = {'p': p, 'm': m, 'm(m-1)': m * (m - 1), 'd': F.d, 'm^2/p': m * m / p, 'psi(-1)': 1 if ((p - 1) // 2) % 4 == 0 else -1,
               'lp_margins': {'stage1': mg1, 'stage2': mg2}}
        # the clique cell: delta = 1 points are exactly the complete ones, all others have e >= 1
        check('C2.clique_cell', all((tp[0] == 1) == (t_ebad(tp, 2, 0, K) == 0) for tp in prof), (p, m))
        check('C2.clique_cell_size', sum(w for tp, w in prof.items() if tp[0] == 1) == m, (p, m))
        # all one-sided rows (verify_side's r = 0 test does not apply to cliques and is skipped)
        nf = len(FAILS)
        verify_side(F, prof, m, 0, 'clique', rep)
        extra = [f for f in FAILS[nf:] if f[0] != 'C.side.r=0']
        del FAILS[nf:]
        FAILS.extend(extra)
        CHECKS['C.side.r=0'] -= 1
        # the circulant sign matrix: (anti)symmetry matches psi(-1), rows are the clique-cell types,
        # pattern HP_4 for rows, and the monochromatic-rectangle HP_4 bound
        sgn = rep['psi(-1)']
        check('C2.matrix_symmetry', all(g[m - t] == ((g[t] + (0 if sgn == 1 else 2)) % 4) for t in range(1, m)), (p, m))
        rows = [tuple(-1 if c == c2 else g[(c2 - c) % m] for c2 in range(m)) for c in range(m)]
        check('C2.rows_distinct', len(set(rows)) == m, (p, m))
        gc = max(sum(1 for x in r if x == 0) for r in rows)
        check('C2.monochromatic_rectangle_HP', gc * gc <= (p - 1) // 4, (p, m, gc))
        # phase coherence forced by the Jacobi identity sum_x f_psi^2 = K(psi,psi) m(m-1)
        s2 = sum(w * zcx(zmul(zf(tp[1], 1, K), zf(tp[1], 1, K), K), K) for tp, w in prof.items())
        a2 = sum(w * zcx(zmul(zf(tp[1], 1, K), zconj(zf(tp[1], 1, K), K), K), K).real for tp, w in prof.items())
        rep['phase_coherence_|sum f_psi^2|/sum|f_psi|^2'] = abs(complex(s2)) / float(a2)
        rep['forced_coherence_sqrt(p)(m-1)/(p-m)'] = math.sqrt(p) * (m - 1) / (p - m)
        rep['levels'] = {'n%d' % e: float(v) for (dl, e), v in levels.items()}
        rep['support_types'] = len(prof)
        rep['certificate'] = [[list(tp[1]) + [tp[0]], str(w)] for tp, w in sorted(prof.items())]
        rep['seconds'] = time.time() - t0
        log('  clique p=%d m=%d m(m-1)/d=%.4f coherence %.4f (%.1fs)' % (
            p, m, m * (m - 1) / F.d, rep['phase_coherence_|sum f_psi^2|/sum|f_psi|^2'], rep['seconds']))
        out.append(rep)
    RESULTS['clique_certificates'] = out


# ============================================================ Section B: actual complete bicliques

def load_configs(pmax=1500):
    """complete bicliques from the two data files, p = 1 mod 4, p <= pmax (deduplicated)."""
    R = os.path.join(ROOT, 'results')
    sig = json.load(open(os.path.join(R, 'sigma_biclique_2026_09_05_search.json')))['table']
    bal = json.load(open(os.path.join(R, 'stepanov_balanced_2026_09_29_data.json')))['rows']
    confs = []
    for ps, row in sig.items():
        p = int(ps)
        if p % 4 != 1 or p > pmax:
            continue
        if row.get('clique'):
            C = row['clique']
            confs.append((p, 'clique', tuple(C), tuple((-c) % p for c in C)))
        if row.get('balanced') and row['balanced'].get('A'):
            confs.append((p, 'balanced', tuple(row['balanced']['A']), None))
        if row.get('sumclique_A'):
            confs.append((p, 'sumclique', tuple(row['sumclique_A']), tuple(row['sumclique_A'])))
        for w in row.get('profile_witnesses', []):
            if w['k'] >= 3 and w.get('A'):
                confs.append((p, 'profile_k=%d' % w['k'], tuple(w['A']), None))
        if row.get('greedy') and row['greedy'].get('A'):
            confs.append((p, 'greedy', tuple(row['greedy']['A']), tuple(row['greedy']['B'])))
    for ps, row in bal.items():
        p = int(ps)
        if p > pmax:
            continue
        for K_, ent in row.get('P', {}).items():
            if ent.get('A'):
                confs.append((p, 'P_K=%s' % K_, tuple(ent['A']), None))
    seen, out = set(), []
    for p, kind, A, B in confs:
        key = (p, tuple(sorted(A)), None if B is None else tuple(sorted(B)))
        if key not in seen:
            seen.add(key)
            out.append((p, kind, A, B))
    return out


def side_profile(F, X):
    """exact profile of the set X: type -> count (integers), and per-point types."""
    delta, cnt = class_counts(X, F.p, F.K, F.ind)
    prof = {}
    for dl, row in zip(delta.tolist(), cnt.tolist()):
        tp = (dl, tuple(row))
        prof[tp] = prof.get(tp, 0) + 1
    return prof, delta, cnt


def column_stats(F, A, B):
    """sign-matrix statistics: for every k | K, k >= 4, and every even class c of order k: the column
    bad counts t^c_b and weights; lopsidedness L^c = max_e W^c_e/(mn - r) over e <= (m-1)/2."""
    p, K = F.p, F.K
    m, n = len(A), len(B)
    r = sum(1 for a in A for b in B if (a + b) % p == 0)
    out = {}
    for k in (4, 8):
        if K % k:
            continue
        for c in range(0, k, 2):
            tb = []
            for b in B:
                t = sum(1 for a in A if (a + b) % p and F.cls(a + b) % k != c)
                dl = sum(1 for a in A if (a + b) % p == 0)
                tb.append((t, dl))
            L = 0.0
            for e in range(0, (m - 1) // 2 + 1):
                W = sum(m - dl for t, dl in tb if t <= e)
                L = max(L, W / (m * n - r))
            out[(k, c)] = (L, tb)
    return out, r


def check_B_restricted_dichotomy(F, A, B, stats, r, info):
    """Theorem 2.1 on an actual biclique: (a) half-weight => the B-restricted class rows follow from HP
    (checked as: B-restricted LHS <= (e+1)(d_k - e)); (b) the lopsided bound (a consequence of RHP_k)."""
    p, d = F.p, F.d
    m, n = len(A), len(B)
    for (k, c), (L, tb) in stats.items():
        d_k = (p - 1) // k
        for e in range(0, min((m - 1) // 2, d_k // 2) + 1):
            lhs2 = sum((e + 1 - t) * (2 * m - 3 * e - t - 2 * dl) for t, dl in tb if t <= e)
            check('B.B_restricted_class_row', lhs2 <= 2 * (e + 1) * (d_k - e), (info, k, c, e))
        for eta in (Fr(1, 16), Fr(1, 8)):
            e_ = math.floor(eta * m)
            W = sum(m - dl for t, dl in tb if t <= e_)
            if W == 0 or m < 8:
                continue
            phi = W / (m * n - r)
            bound = 2.0 / (k * phi) * (1 - math.sqrt(2 * eta)) ** -2 * d + m / phi
            check('B.lopsided_bound(Thm2.1b)', m * n - r <= bound * (1 + 1e-12), (info, k, c, float(eta)))


def section_B(pmax=1500):
    log('Section B: actual bicliques')
    confs = load_configs(pmax)
    table = []
    agg = {}
    for (p, kind, A, B) in confs:
        K = 8 if p % 8 == 1 else 4
        F = Field(p, K)
        A = tuple(a % p for a in A)
        if B is None:
            B = tuple(b for b in range(p) if all((a + b) % p == 0 or F.ind[(a + b) % p] % 2 == 0 for a in A))
        m, n = len(A), len(B)
        complete = all((a + b) % p == 0 or F.ind[(a + b) % p] % 2 == 0 for a in A for b in B)
        check('B.is_complete_biclique', complete, (p, kind))
        if not complete or m < 2 or n < 2:
            continue
        info = (p, kind, m, n)
        profA, dA, cA = side_profile(F, A)
        profB, dB, cB = side_profile(F, B)
        # every one-sided row (theorems) on the actual profiles: class rows for all k | K, identities
        for X, prof, mX in ((A, profA, m), (B, profB, n)):
            if mX > (p + 1) // 2:
                continue
            for k in (2, 4, 8):
                if K % k:
                    continue
                d_k = (p - 1) // k
                for c in range(k):
                    hist = {}
                    for tp, w in prof.items():
                        key = (t_ebad(tp, k, c, K), tp[0])
                        hist[key] = hist.get(key, 0) + w
                    hist_star_ok(hist, mX, d_k, min((mX - 1) // 2, d_k // 2), 'B.star_k_row', (info, k, c))
                    hist_subset_ok(hist, mX, d_k, 'B.subsetHP_k_row', (info, k, c), tmax=min(mX, 8))
            fv = {tp: {j: zf(tp[1], j, K) for j in range(1, K)} for tp in prof}
            for i in range(1, K):
                for j in range(i, K):
                    lhs = zint(0, K)
                    for tp, w in prof.items():
                        lhs = zadd(lhs, zscale(zmul(fv[tp][i], fv[tp][j], K), w))
                    if (i + j) % K == 0:
                        rhs = zint(mX * (p - mX), K)
                    else:
                        tneg = zint(0, K)
                        for tp, w in prof.items():
                            if tp[0] == 1:
                                tneg = zadd(tneg, zscale(zf(tp[1], (i + j) % K, K), w))
                        rhs = zmul(F.Kc[(i, j)], tneg, K)
                    check('B.jacobi_identity_one_sided', lhs == rhs, (info, i, j))
            if mX >= 3:
                dl_, cn_ = (dA, cA) if X is A else (dB, cB)
                for j in kw_family(K):
                    sums = kw_point_sums(dl_, cn_, j, K, mX, 5)
                    for t, zz in sums.items():
                        if t >= 3:
                            check('B.KW_row_exact', weil4_exact_ok(zmul(zz, zconj(zz, K), K),
                                                                    (comb(mX, t) * (t - 1)) ** 2 * p, K), (info, j, t))
            bnd = (2 * mX * mX - mX) * p + 3 * mX ** 4 * isqrt(p)
            for j in range(1, K):
                if (2 * j) % K:
                    s = zint(0, K)
                    for tp, w in prof.items():
                        a2 = zmul(fv[tp][j], zconj(fv[tp][j], K), K)
                        s = zadd(s, zscale(zmul(a2, a2, K), w))
                    check('B.weil4_row', weil4_exact_ok(s, bnd, K), (info, j))
        # two-sided Jacobi identities (all pairs)
        inter = len(set(A) & set(B))
        G = cA.T.astype(object) @ cB.astype(object)  # G[c, c'] = sum_x n^A_c(x) n^B_c'(x)
        for i in range(1, K):
            for j in range(1, K):
                v = [0] * K
                for c in range(K):
                    for c2 in range(K):
                        v[(i * c + j * c2) % K] += int(G[c, c2])
                lhs = zred(v, K)
                if (i + j) % K == 0:
                    rhs = zint(p * inter - m * n, K)
                else:
                    w = [0] * K
                    for a in A:
                        for b in B:
                            if a != b:
                                w[((i + j) * F.cls(b - a)) % K] += 1
                    rhs = zmul(F.Kc[(i, j)], zred(w, K), K)
                check('B.jacobi_identity_two_sided', lhs == rhs, (info, i, j))
        # sign matrix statistics
        stats, r = column_stats(F, A, B)
        statsT, _ = column_stats(F, B, A)
        check_B_restricted_dichotomy(F, A, B, stats, r, info)
        check_B_restricted_dichotomy(F, B, A, statsT, r, info + ('rows',))
        # column pattern HP (Prop 1.5)
        for k in (4, 8):
            if K % k:
                continue
            d_k = (p - 1) // k
            for X, Y, mX in ((A, B, m), (B, A, n)):
                pats = {}
                for y in Y:
                    if any((x + y) % p == 0 for x in X):
                        continue
                    key = tuple(F.cls(x + y) % k for x in X)
                    pats[key] = pats.get(key, 0) + 1
                if pats:
                    check('B.pattern_HP', mX * max(pats.values()) <= d_k + mX - 1, (info, k))
        M4 = [[0 if (a + b) % p == 0 else (1 if F.cls(a + b) % 4 == 0 else -1) for b in B] for a in A]
        S4 = sum(map(sum, M4))
        colsq = sum(sum(M4[i][j] for i in range(m)) ** 2 for j in range(n))
        rowsq = sum(sum(M4[i]) ** 2 for i in range(m))
        ent = {'p': p, 'kind': kind, 'm': m, 'n': n, 'r': r, 'mn_over_p': m * n / p,
               'S_psi_over_(mn-r)': S4 / (m * n - r),
               'colsq_over_nm^2': colsq / (n * m * m), 'rowsq_over_mn^2': rowsq / (m * n * n),
               'lopsided_cols': {'%d,%d' % kc: round(v[0], 4) for kc, v in stats.items()},
               'lopsided_rows': {'%d,%d' % kc: round(v[0], 4) for kc, v in statsT.items()}}
        # exact scaling factor Theta (Theorem 2.1): largest admissible (mn - r)/d for this column/row statistic
        theta = None
        for X_, st_, mX, nY in ((A, stats, m, n), (B, statsT, n, m)):
            for (k, c), (L, tb) in st_.items():
                d_k = (p - 1) // k
                for e in range(0, min((mX - 1) // 2, d_k // 2) + 1):
                    lhs2 = sum((e + 1 - t) * (2 * mX - 3 * e - t - 2 * dl) for t, dl in tb if t <= e)
                    if lhs2:
                        th = Fr(2 * (e + 1) * (d_k - e) * (m * n - r), F.d * lhs2)
                        if theta is None or th < theta[0]:
                            theta = (th, k, c, e, 'columns' if X_ is A else 'rows')
        ent['Theta'] = None if theta is None else float(theta[0])
        ent['Theta_at'] = None if theta is None else list(theta[1:])
        hw = all(v[0] <= 2.0 / kc[0] for kc, v in stats.items()) and all(v[0] <= 2.0 / kc[0] for kc, v in statsT.items())
        ent['half_weight_rows_and_columns'] = hw
        table.append(ent)
        kk = kind.split('=')[0]
        a = agg.setdefault(kk, {'configs': 0, 'half_weight_holds': 0, 'max_lopsided_quartic': 0.0,
                                'max_lopsided_octic': 0.0, 'max_|S_psi|/(mn-r)': 0.0, 'mean_colsq_ratio_x_m': 0.0})
        a['configs'] += 1
        a['half_weight_holds'] += int(hw)
        for (k, c), v in list(stats.items()) + list(statsT.items()):
            key = 'max_lopsided_quartic' if k == 4 else 'max_lopsided_octic'
            a[key] = max(a[key], v[0])
        a['max_|S_psi|/(mn-r)'] = max(a['max_|S_psi|/(mn-r)'], abs(S4) / (m * n - r))
        a.setdefault('Theta_values', []).append(float('inf') if ent['Theta'] is None else ent['Theta'])
        a['mean_colsq_ratio_x_m'] += colsq / (n * m * m) * m
    for a in agg.values():
        a['mean_colsq_ratio_x_m'] /= max(a['configs'], 1)
        tv = sorted(a.pop('Theta_values'))
        a['Theta_min'] = tv[0]
        a['Theta_median'] = tv[len(tv) // 2]
        a['Theta<1'] = sum(1 for x in tv if x < 1)
    RESULTS['B_actual_bicliques'] = {'n_configurations': len(table), 'aggregate_by_kind': agg,
                                     'balanced_and_cliques': [t for t in table if t['kind'] in
                                                              ('clique', 'balanced', 'sumclique', 'greedy')],
                                     'n_half_weight_fails': sum(1 for t in table if not t['half_weight_rows_and_columns']),
                                     'Theta<1': [t for t in table if t['Theta'] is not None and t['Theta'] < 1]}
    log('  %d configurations' % len(table))


# ============================================================ Section D: the tight family A = {0, 1}

def section_D(pmax=3000):
    log('Section D: A = {0,1}')
    from sympy import primerange
    rows = []
    worst = 0.0
    for p in primerange(5, pmax + 1):
        if p % 4 != 1:
            continue
        K = 8 if p % 8 == 1 else 4
        g, ind = ind_table(p)
        A = (0, 1)
        B = [b for b in range(p) if all((a + b) % p == 0 or ind[(a + b) % p] % 2 == 0 for a in A)]
        n = len(B)
        check('D.|B(0,1)|=(p+3)/4', n == (p + 3) // 4, p)
        r = 2
        cnt = {}
        for b in B:
            key = tuple(-1 if (a + b) % p == 0 else int(ind[(a + b) % p] % K) for a in A)
            cnt[key] = cnt.get(key, 0) + 1
        for k in (4, 8):
            if K % k:
                continue
            for c in range(0, k, 2):
                W0 = sum(v * (2 - sum(1 for x in key if x == -1)) for key, v in cnt.items()
                         if all(x == -1 or x % k == c for x in key))
                frac = W0 / (2 * n - r)
                worst = max(worst, frac * k / 2)
                check('D.half_weight', frac <= 2.0 / k, (p, k, c, frac))
                check('D.HP_k_class_row_e0', W0 <= (p - 1) // k, (p, k, c))
        if p in (13, 37, 41, 101, 1009, 1013, 2969) or p < 60:
            q = {}
            for key, v in cnt.items():
                if -1 not in key:
                    q[tuple((x % 4) // 2 for x in key)] = q.get(tuple((x % 4) // 2 for x in key), 0) + v
            rows.append({'p': p, 'n': n, 'quartic_column_patterns(0=+1,1=-1)': {str(k_): v for k_, v in sorted(q.items())}})
    RESULTS['D_A={0,1}'] = {'max_ratio_of_constant_column_weight_to_half_weight_threshold': worst, 'samples': rows}
    log('  max W_0/((2/k)(mn-r)) = %.4f' % worst)


def main(sections='AC'):
    if 'A' in sections:
        section_A()
    if 'B' in sections:
        section_B()
    if 'C' in sections:
        section_C(CERT_PARAMS if '+' not in sections else CERT_PARAMS[:2])
    if 'C' in sections:
        section_C2(CLIQUE_PARAMS if '+' not in sections else CLIQUE_PARAMS[:2])
    if 'D' in sections:
        section_D()
    RESULTS['checks'] = CHECKS
    RESULTS['n_checks'] = sum(CHECKS.values())
    RESULTS['failures'] = FAILS[:200]
    RESULTS['n_failures'] = len(FAILS)
    RESULTS['seconds'] = time.time() - T0
    return RESULTS


if __name__ == '__main__':
    sec = sys.argv[1] if len(sys.argv) > 1 else 'ABCD'
    main(sec)
    if len(sys.argv) <= 1:
        with open(OUT, 'w') as fh:
            json.dump(RESULTS, fh, indent=1, default=str)
    print(json.dumps({'n_checks': RESULTS['n_checks'], 'n_failures': len(FAILS), 'seconds': RESULTS['seconds']}))
    for f in FAILS[:20]:
        print(f)
    for kname, v in sorted(CHECKS.items()):
        print(kname, v)
