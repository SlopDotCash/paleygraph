#!/usr/bin/env python3
"""Verifier for research/stepanov-balanced-2026-09-29.md  (Stepanov wave, worker `balanced`).

Run:  /opt/miniconda3/bin/python3 experiments/stepanov_balanced_2026_09_29.py
Writes results/stepanov_balanced_2026_09_29.json (check counts, failures, witnesses, tables).

Sections (numbering follows the note):
  A  Theorem 1.1: F_A = lam * P1^m * P0^(m-1) * G, exact orders, deg G = g = d - mn + r,
     G(b) != 0 on B, same defect for F_B, converse (divisibility => biclique).
  B  Lemma 1.2 residue identities and the character constraint for even m.
  C  Proposition 1.3: derivative and Hankel-minor factorisations; the pencil degree
     deg R = n0(F_A); the Wronskian identity and "all roots have multiplicity <= m".
  D  Proposition 1.4 (derivative HP): (m-1)(N_0 + N_m) <= d - 1 + r_0 + r_m.
  E  Section 3 attempts: (b) truncation/abc bookkeeping, (c) (*)-slack for complete bicliques,
     (a) two-sided multiplicity bookkeeping, (d) rational roots of G.
  F  Task 2: factorisation type of G_A for extremal configurations; data table from
     results/stepanov_balanced_2026_09_29_data.json (produced by the companion driver
     experiments/stepanov_balanced_2026_09_29_data.py and the C helper), re-checked here.
"""
import json, os, sys, time, random, itertools, math
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
OUT = os.path.join(ROOT, 'results', 'stepanov_balanced_2026_09_29.json')
DATA = os.path.join(ROOT, 'results', 'stepanov_balanced_2026_09_29_data.json')
TABLE = os.path.join(ROOT, 'results', 'sigma_biclique_2026_09_05_search.json')

CHECKS = {}
FAILS = []
WIT = {}


def check(name, cond, info=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not cond:
        FAILS.append((name, info))
        if len(FAILS) < 20:
            print('FAIL', name, info)


# ---------------------------------------------------------------- polynomials mod p
# coefficient arrays, index = degree, dtype int64, entries in [0,p)

def trim(a):
    a = np.asarray(a, dtype=np.int64)
    nz = np.nonzero(a)[0]
    if len(nz) == 0:
        return np.zeros(1, dtype=np.int64)
    return a[:nz[-1] + 1].copy()


def deg(a):
    a = trim(a)
    return -1 if (len(a) == 1 and a[0] == 0) else len(a) - 1


def pmul(a, b, p):
    # exact for p < 2^15 and lengths < 2^20
    return trim(np.convolve(np.asarray(a, dtype=np.int64), np.asarray(b, dtype=np.int64)) % p)


def padd(a, b, p):
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] += b
    return trim(c % p)


def psub(a, b, p):
    n = max(len(a), len(b))
    c = np.zeros(n, dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] -= b
    return trim(c % p)


def pscale(a, s, p):
    return trim((np.asarray(a, dtype=np.int64) * (s % p)) % p)


def pdivmod(a, b, p):
    a = trim(a).copy()
    b = trim(b)
    db = len(b) - 1
    if db < 0 or (db == 0 and b[0] == 0):
        raise ZeroDivisionError
    inv = pow(int(b[-1]), p - 2, p)
    if len(a) - 1 < db:
        return np.zeros(1, dtype=np.int64), a
    q = np.zeros(len(a) - db, dtype=np.int64)
    for i in range(len(a) - 1, db - 1, -1):
        c = (int(a[i]) * inv) % p
        if c:
            q[i - db] = c
            a[i - db:i + 1] = (a[i - db:i + 1] - c * b) % p
    return trim(q), trim(a[:db] if db > 0 else np.zeros(1, dtype=np.int64))


def pgcd(a, b, p):
    a, b = trim(a), trim(b)
    while deg(b) >= 0:
        _, r = pdivmod(a, b, p)
        a, b = b, r
    if deg(a) < 0:
        return a
    return pscale(a, pow(int(a[-1]), p - 2, p), p)


def peval(a, x, p):
    r = 0
    for c in a[::-1]:
        r = (r * x + int(c)) % p
    return r


def pderiv(a, p):
    a = trim(a)
    if len(a) == 1:
        return np.zeros(1, dtype=np.int64)
    return trim((a[1:] * np.arange(1, len(a), dtype=np.int64)) % p)


def pfromroots(roots, p, mult=1):
    r = np.array([1], dtype=np.int64)
    for b in roots:
        r = pmul(r, np.array([(-b) % p, 1], dtype=np.int64), p)
    if mult != 1:
        r = ppow(r, mult, p)
    return r


def ppow(a, e, p):
    r = np.array([1], dtype=np.int64)
    base = trim(a)
    while e:
        if e & 1:
            r = pmul(r, base, p)
        e >>= 1
        if e:
            base = pmul(base, base, p)
    return r


def pmulmod(a, b, m, p):
    return pdivmod(pmul(a, b, p), m, p)[1]


def ppowmod(a, e, m, p):
    r = np.array([1], dtype=np.int64)
    base = pdivmod(a, m, p)[1]
    while e:
        if e & 1:
            r = pmulmod(r, base, m, p)
        e >>= 1
        if e:
            base = pmulmod(base, base, m, p)
    return r


# ---------------------------------------------------------------- basic arithmetic

def chi_table(p):
    c = [-1] * p
    c[0] = 0
    for x in range(1, p):
        c[x * x % p] = 1
    return c


_fact_cache = {}


def binom_mod(n, k, p):
    if k < 0 or k > n:
        return 0
    if p not in _fact_cache:
        f = [1] * p
        for i in range(1, p):
            f[i] = f[i - 1] * i % p
        _fact_cache[p] = f
    f = _fact_cache[p]
    assert n < p
    return f[n] * pow(f[k] * f[n - k] % p, p - 2, p) % p


def lagrange_c(A, p):
    cs = []
    for i, a in enumerate(A):
        pr = 1
        for j, b in enumerate(A):
            if j != i:
                pr = pr * (a - b) % p
        cs.append(pow(pr, p - 2, p))
    return cs


def h_series(A, N, p):
    """complete homogeneous symmetric polynomials h_0..h_N of A mod p."""
    h = np.zeros(N + 1, dtype=np.int64)
    h[0] = 1
    for a in A:
        a %= p
        # multiply by 1/(1 - a t): h_new[j] = h[j] + a h_new[j-1]
        for j in range(1, N + 1):
            h[j] = (h[j] + a * h[j - 1]) % p
    return h


def hp_poly(A, p):
    """F_A(x) = -1 + sum_k c_k (x + a_k)^D, D = d + m - 1, as coefficient array (degree d)."""
    m = len(A)
    d = (p - 1) // 2
    D = d + m - 1
    assert D <= p - 1
    h = h_series(A, d, p)
    F = np.zeros(d + 1, dtype=np.int64)
    for j in range(d + 1):
        F[d - j] = binom_mod(D, m - 1 + j, p) * int(h[j]) % p
    F[0] = (F[0] - 1) % p
    return trim(F)


def B_of_A(A, p, chi):
    return [b for b in range(p) if all((a + b) % p == 0 or chi[(a + b) % p] == 1 for a in A)]


# ---------------------------------------------------------------- Section A/B: Theorem 1.1, Lemma 1.2

def structure(A, B, p, chi, full=True):
    """Checks Theorem 1.1 and Lemma 1.2 for a complete biclique (A,B); returns G and data."""
    m, n = len(A), len(B)
    d = (p - 1) // 2
    D = d + m - 1
    negA = set((-a) % p for a in A)
    B0 = [b for b in B if b in negA]
    B1 = [b for b in B if b not in negA]
    r = len(B0)
    g = d - m * n + r
    check('A.hp_bound(mn-r<=d)', g >= 0, (p, A, B))
    F = hp_poly(A, p)
    check('A.degF=d', deg(F) == d, (p, A))
    lam = binom_mod(D, m - 1, p)
    check('A.leadF=binom(D,m-1)', int(F[-1]) == lam, (p, A))
    # direct evaluation of F at a few points against the definition
    cs = lagrange_c(A, p)
    for x in random.sample(range(p), min(p, 6)):
        v = (-1 + sum(c * pow((x + a) % p, D, p) for c, a in zip(cs, A))) % p
        check('A.F_eval_matches_definition', v == peval(F, x, p), (p, A, x))
    P1 = pfromroots(B1, p, m)
    P0 = pfromroots(B0, p, m - 1) if m >= 2 else np.array([1], dtype=np.int64)
    Pm = pmul(P1, P0, p)
    Q, R = pdivmod(F, Pm, p)
    check('A.divisibility', deg(R) < 0, (p, A, B))
    G = pscale(Q, pow(lam, p - 2, p), p)
    check('A.degG=g', deg(G) == g, (p, A, B, deg(G), g))
    check('A.G_monic', int(G[-1]) == 1, (p, A))
    for b in B:
        check('A.G(b)!=0 (exact orders m, m-1)', peval(G, b, p) != 0, (p, A, b))
    if not full:
        return G, g, r
    # Lemma 1.2 residue identities
    Pfull = pfromroots(B, p)
    dP = pderiv(Pfull, p)
    P0s = pfromroots(B0, p)
    P1s = pfromroots(B1, p)
    target = (pow(-1, m) * pow(2 * m, p - 2, p)) % p
    for b in B1:
        lhs = pow(peval(dP, b, p), m, p) * peval(G, b, p) % p
        for a in A:
            lhs = lhs * ((a + b) % p) % p
        lhs = lhs * pow(peval(P0s, b, p), p - 2, p) % p
        check('B.residue_identity_B1', lhs == target, (p, A, b))
        if m % 2 == 1:
            check('B.chi(G(b))=chi(-2m)chi(P\'(b)) (m odd)',
                  chi[peval(G, b, p)] == chi[(-2 * m) % p] * chi[peval(dP, b, p)], (p, A, b))
        if m % 2 == 0:
            check('B.chi(G(b))=chi(2m) (m even)', chi[peval(G, b, p)] == chi[(2 * m) % p], (p, A, b))
            check('B.chi(P0(b))=1', chi[peval(P0s, b, p)] == 1 or len(B0) == 0, (p, A, b))
    dP0 = pderiv(P0s, p)
    for b in B0:
        k0 = A.index((-b) % p)
        lhs = pow(peval(P1s, b, p), m, p) * pow(peval(dP0, b, p), m - 1, p) * peval(G, b, p) % p
        check('B.residue_identity_B0', lhs == (-cs[k0]) % p, (p, A, b))
    return G, g, r


def converse_check(A, p, chi, rng):
    """Theorem 1.1 converse: if P^m (resp. P^(m-1) at -A) divides F_A then b + A in Q u {0}."""
    m = len(A)
    F = hp_poly(A, p)
    negA = set((-a) % p for a in A)
    dF = [F]
    for s in range(m + 1):
        dF.append(pderiv(dF[-1], p))
    for b in range(p):
        need = m - 1 if b in negA else m
        vanish = all(peval(dF[s], b, p) == 0 for s in range(need))
        good = all((a + b) % p == 0 or chi[(a + b) % p] == 1 for a in A)
        check('A.converse(order>=m-delta iff complete)', vanish == good, (p, A, b))
        if good and b not in negA:
            check('A.exact_order_m', peval(dF[m], b, p) != 0, (p, A, b))
        if good and b in negA and m >= 2:
            check('A.exact_order_m-1', peval(dF[m - 1], b, p) != 0, (p, A, b))


# ---------------------------------------------------------------- Section C

def squarefree_info(G, p):
    """number of distinct roots of G in closure and the max multiplicity (deg G < p)."""
    if deg(G) <= 0:
        return 0, 0
    g1 = pgcd(G, pderiv(G, p), p)
    nd = deg(G) - deg(g1)
    # max multiplicity via repeated gcds (Yun-style)
    mult = 1
    cur = g1
    while deg(cur) > 0:
        mult += 1
        cur = pgcd(cur, pderiv(cur, p), p)
    return nd, mult


def section_C(A, B, p, chi, G, g, r, hankel=True):
    m, n = len(A), len(B)
    d = (p - 1) // 2
    D = d + m - 1
    negA = set((-a) % p for a in A)
    B0 = [b for b in B if b in negA]
    B1 = [b for b in B if b not in negA]
    F = hp_poly(A, p)
    # derivative factorisations: F^(s) = P1^(m-s) P0^(max(m-1-s,0)) Q_s, deg Q_s = g + s(n-1)
    Fs = F
    for s in range(0, m):
        if s > 0:
            Fs = pderiv(Fs, p)
        e1, e0 = m - s, max(m - 1 - s, 0)
        M = pmul(pfromroots(B1, p, e1), pfromroots(B0, p, e0) if e0 > 0 else np.array([1]), p)
        Qs, R = pdivmod(Fs, M, p)
        check('C.derivative_divisibility', deg(R) < 0, (p, A, s))
        expect = g + s * (n - 1) if s <= m - 2 else d - m - n + r + 1
        check('C.derivative_quotient_degree', deg(Qs) == expect, (p, A, s, deg(Qs), expect))
    # distinct roots and multiplicities
    nG, multG = squarefree_info(G, p)
    n0F = n + nG
    check('C.max_root_multiplicity_of_G<=m', multG <= m, (p, A, multG))
    # pencil map R = D F / F' - x : degree = n0(F)
    dF = pderiv(F, p)
    Vn = psub(pscale(F, D, p), pmul(np.array([0, 1]), dF, p), p)  # D F - x F'
    gg = pgcd(Vn, dF, p)
    num = pdivmod(Vn, gg, p)[0]
    den = pdivmod(dF, gg, p)[0]
    degR = max(deg(num), deg(den))
    check('C.pencil_degree=n0(F_A)', degR == n0F, (p, A, degR, n0F))
    check('C.pencil_degree<=n+g', degR <= n + g, (p, A))
    out = {'nG_distinct': nG, 'G_max_mult': multG, 'n0F': n0F, 'degR': degR}
    # Hankel minors e = 1, 2 (if allowed): H_{e+1} divisible by P1^{(e+1)(m-e)} P0^{(e+1)(m-1-e)}
    if hankel:
        u = [F]
        for s in range(1, 5):
            u.append(pderiv(u[-1], p))
        # normalise u_s = F^(s)/(D)_s
        us = []
        for s in range(5):
            ff = 1
            for i in range(s):
                ff = ff * (D - i) % p
            us.append(pscale(u[s], pow(ff, p - 2, p), p))
        for e in (1, 2):
            if 2 * e > m - 1:
                continue
            if e == 1:
                H = psub(pmul(us[0], us[2], p), pmul(us[1], us[1], p), p)
            else:
                # 3x3 Hankel det [u_{i+j}]
                a, b, c, dd, ee = us[0], us[1], us[2], us[3], us[4]
                t1 = pmul(a, psub(pmul(c, ee, p), pmul(dd, dd, p), p), p)
                t2 = pmul(b, psub(pmul(b, ee, p), pmul(c, dd, p), p), p)
                t3 = pmul(c, psub(pmul(b, dd, p), pmul(c, c, p), p), p)
                H = padd(psub(t1, t2, p), t3, p)
            check('C.hankel_degree=(e+1)(d-e)', deg(H) == (e + 1) * (d - e), (p, A, e))
            M = pmul(pfromroots(B1, p, (e + 1) * (m - e)),
                     pfromroots(B0, p, (e + 1) * (m - 1 - e)) if m - 1 - e > 0 else np.array([1]), p)
            K, R = pdivmod(H, M, p)
            check('C.hankel_divisibility', deg(R) < 0, (p, A, e))
            check('C.hankel_quotient_degree', deg(K) == (e + 1) * (g + e * (n - 1)), (p, A, e))
            # exact order at complete points not in -A: quotient nonzero there iff the
            # leading determinant det[(m)_{i+j}/(D)_{i+j}] is nonzero mod p
            ext = sum(1 for b in B1 if peval(K, b, p) == 0)
            out['hankel_extra_vanishing_e%d' % e] = ext
    return out


def wronskian_identity(A, p):
    """W(F, f_1..f_m) = kappa * prod (x+a_k)^(D-m), f_k = c_k (x+a_k)^D  (small m, exact)."""
    m = len(A)
    d = (p - 1) // 2
    D = d + m - 1
    cs = lagrange_c(A, p)
    F = hp_poly(A, p)
    fs = []
    for c, a in zip(cs, A):
        f = np.zeros(D + 1, dtype=np.int64)
        for l in range(D + 1):
            f[D - l] = binom_mod(D, l, p) * pow(a, l, p) % p * c % p
        fs.append(trim(f))
    cols = [F] + fs
    # matrix rows = derivative orders 0..m
    mat = []
    for col in cols:
        ds = [col]
        for _ in range(m):
            ds.append(pderiv(ds[-1], p))
        mat.append(ds)
    # determinant by Leibniz (m+1 <= 4)
    W = np.zeros(1, dtype=np.int64)
    idx = list(range(m + 1))
    for perm in itertools.permutations(idx):
        sgn = 1
        pl = list(perm)
        for i in range(len(pl)):
            for j in range(i + 1, len(pl)):
                if pl[i] > pl[j]:
                    sgn = -sgn
        term = np.array([1], dtype=np.int64)
        for colidx, row in enumerate(perm):
            term = pmul(term, mat[colidx][row], p)
        W = padd(W, term, p) if sgn == 1 else psub(W, term, p)
    target = np.array([1], dtype=np.int64)
    for a in A:
        target = pmul(target, ppow(np.array([a % p, 1]), D - m, p), p)
    q, rem = pdivmod(W, target, p)
    ok = deg(rem) < 0 and deg(q) == 0 and int(q[0]) != 0
    check('C.wronskian_identity', ok, (p, A))
    return ok


# ---------------------------------------------------------------- Section D: derivative HP

def section_D(A, p, chi):
    m = len(A)
    d = (p - 1) // 2
    negA = set((-a) % p for a in A)
    N0 = Nm = r0 = rm = 0
    for x in range(p):
        vals = [chi[(a + x) % p] for a in A]
        if all(v >= 0 for v in vals):
            N0 += 1
            r0 += x in negA
        if all(v <= 0 for v in vals):
            Nm += 1
            rm += x in negA
    lhs = (m - 1) * (N0 + Nm)
    rhs = d - 1 + r0 + rm
    check('D.derivative_HP', lhs <= rhs, (p, A, lhs, rhs))
    return lhs, rhs, N0, Nm


# ---------------------------------------------------------------- Section E

def e_profile(A, p, chi):
    """e_x (number of a with chi(a+x) = -1) and delta_x for all x."""
    negA = set((-a) % p for a in A)
    return [(sum(1 for a in A if chi[(a + x) % p] == -1), x in negA) for x in range(p)]


def star_slack(A, p, chi):
    """(*)_e of the robust note for all e <= (m-1)/2: returns list of (lhs, rhs)."""
    m = len(A)
    d = (p - 1) // 2
    prof = e_profile(A, p, chi)
    res = []
    for e in range(0, (m - 1) // 2 + 1):
        lhs = 0.0
        num = 0  # exact: 2*lhs integer
        for eb, dl in prof:
            if eb <= e:
                num += (e + 1 - eb) * (2 * m - (3 * e + eb) - 2 * dl)
        rhs2 = 2 * (e + 1) * (d - e)
        check('E.star_holds', num <= rhs2, (p, A, e))
        res.append((num, rhs2))
    return res


# ---------------------------------------------------------------- Section F: factorisation type

def roots_in_Fp(G, p):
    if deg(G) <= 0:
        return 0
    xp = ppowmod(np.array([0, 1]), p, G, p)
    h = psub(xp, np.array([0, 1]), p)
    return deg(pgcd(G, h, p))


def ddf_pattern(G, p, maxdeg=None):
    """distinct-degree factorisation of a squarefree part: returns {degree: number of factors}."""
    f = trim(G)
    if deg(f) <= 0:
        return {}
    f = pscale(f, pow(int(f[-1]), p - 2, p), p)
    # squarefree part
    g1 = pgcd(f, pderiv(f, p), p)
    f = pdivmod(f, g1, p)[0]
    pat = {}
    i = 0
    x = np.array([0, 1], dtype=np.int64)
    h = x
    # Frobenius via matrix: columns x^{p*j} mod f
    n = deg(f)
    xp = ppowmod(x, p, f, p)
    cols = [np.array([1], dtype=np.int64)]
    for j in range(1, n):
        cols.append(pmulmod(cols[-1], xp, f, p))
    Fm = np.zeros((n, n), dtype=np.int64)
    for j, c in enumerate(cols):
        Fm[:len(c), j] = c
    def frob(v):
        vv = np.zeros(n, dtype=np.int64)
        vv[:len(v)] = v
        return trim((Fm @ vv) % p)
    cur = f
    while deg(cur) >= 2 * (i + 1):
        i += 1
        if maxdeg and i > maxdeg:
            break
        h = frob(pdivmod(h, f, p)[1])
        gg = pgcd(cur, psub(h, x, p), p)
        if deg(gg) > 0:
            pat[i] = deg(gg) // i
            cur = pdivmod(cur, gg, p)[0]
    if deg(cur) > 0:
        pat[deg(cur)] = pat.get(deg(cur), 0) + 1
    return pat



# ---------------------------------------------------------------- Theorem 1.1(5): rational multiplicities

def rational_multiplicity_check(A, p, chi):
    m = len(A)
    F = hp_poly(A, p)
    negA = set((-a) % p for a in A)
    ders = [F]
    for s in range(m + 1):
        ders.append(pderiv(ders[-1], p))
    for x in range(p):
        if x in negA:
            continue
        ex = sum(1 for a in A if chi[(a + x) % p] == -1)
        mu = 0
        while mu <= m and peval(ders[mu], x, p) == 0:
            mu += 1
        if ex == 0:
            check('A.rational_order_complete=m', mu == m, (p, A, x))
        else:
            check('A.rational_order<=min(e-1,m-e)', mu <= min(ex - 1, m - ex), (p, A, x, mu, ex))


# ---------------------------------------------------------------- Remark 1.2': interpolation reach

def interpolation_check(A, B, p, chi, G, g):
    """G_A restricted to B equals the values prescribed by Lemma 1.2; if g <= n-2 the
    interpolant of the prescribed values has degree <= g (the n-1-g relations hold)."""
    m, n = len(A), len(B)
    negA = set((-a) % p for a in A)
    B0 = [b for b in B if b in negA]
    B1 = [b for b in B if b not in negA]
    P = pfromroots(B, p); dP = pderiv(P, p)
    P0 = pfromroots(B0, p); dP0 = pderiv(P0, p); P1 = pfromroots(B1, p)
    cs = lagrange_c(A, p)
    vals = {}
    for b in B1:
        den = pow(peval(dP, b, p), m, p)
        for a in A:
            den = den * ((a + b) % p) % p
        vals[b] = (pow(-1, m) * pow(2 * m, p - 2, p) * peval(P0, b, p) * pow(den, p - 2, p)) % p
    for b in B0:
        k0 = A.index((-b) % p)
        den = pow(peval(P1, b, p), m, p) * pow(peval(dP0, b, p), m - 1, p) % p
        vals[b] = (-cs[k0]) * pow(den, p - 2, p) % p
    # relations sum_b G(b) b^j / P'(b) = 0 for j <= n-2-g
    nrel = 0
    for j in range(0, n - 1 - g):
        tot = 0
        for b in B:
            tot = (tot + vals[b] * pow(b, j, p) * pow(peval(dP, b, p), p - 2, p)) % p
        check('B.interpolation_relations(g<=n-2)', tot == 0, (p, A, j))
        nrel += 1
    return nrel


# ---------------------------------------------------------------- Lemma 3.1: grid multiplicity SZ

def grid_sz_check(rng, trials=300):
    """Sum over A x B of mult(Phi) <= deg(Phi) * max(|A|,|B|) for products of linear forms."""
    for _ in range(trials):
        p = rng.choice([7, 11, 13, 17])
        t = rng.randint(1, 7)
        forms = []
        for _ in range(t):
            while True:
                al, be, ga = rng.randrange(p), rng.randrange(p), rng.randrange(p)
                if al or be:
                    break
            forms.append((al, be, ga))
        A = rng.sample(range(p), rng.randint(1, p))
        B = rng.sample(range(p), rng.randint(1, p))
        tot = 0
        for a in A:
            for b in B:
                tot += sum(1 for (al, be, ga) in forms if (al * a + be * b + ga) % p == 0)
        check('E.grid_multiplicity_SZ', tot <= t * max(len(A), len(B)), (p, forms, A, B))


# ---------------------------------------------------------------- Prop 1.3(2): Delta_e != 0

def sympy_rat(a, b):
    import sympy
    return sympy.Rational(a, b)


def delta_e(m, D, e, p):
    from fractions import Fraction
    import sympy
    M = sympy.Matrix(e + 1, e + 1, lambda i, j: sympy.Rational(math.perm(m, i + j), math.perm(D, i + j)))
    v = sympy.Rational(M.det())
    return v, (v.p % p != 0)


# ---------------------------------------------------------------- Theorem 3.5: LP certificate

def lp_rows(p, m, ntarget, band=0.4):
    from fractions import Fraction as Fr
    from math import comb, isqrt
    d = (p - 1) // 2
    I = lambda j: j
    J = lambda j: m + 1 + j
    val = lambda j, dl: m - dl - 2 * j
    rows = []
    rows.append(({**{I(j): 1 for j in range(m + 1)}, **{J(j): 1 for j in range(m)}}, '=', p, 'total'))
    rows.append(({J(j): 1 for j in range(m)}, '=', m, 'negA'))
    rows.append(({**{I(j): val(j, 0) for j in range(m + 1)}, **{J(j): val(j, 1) for j in range(m)}}, '=', 0, 'moment1'))
    rows.append(({**{I(j): val(j, 0) ** 2 for j in range(m + 1)}, **{J(j): val(j, 1) ** 2 for j in range(m)}}, '=', m * (p - m), 'moment2'))
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
    sq = isqrt(p)  # <= sqrt(p): makes every Weil row stricter
    for k in (2, 3, 4):
        df = 1
        for i in range(1, 2 * k, 2):
            df *= i
        rhs = df * m ** k * p + (2 * k - 1) * m ** (2 * k) * sq
        rows.append(({**{I(j): val(j, 0) ** (2 * k) for j in range(m + 1)},
                      **{J(j): val(j, 1) ** (2 * k) for j in range(m)}}, '<=', rhs, 'weil_moment_%d' % (2 * k)))
    # Weil counts (sharpen Lemma 3.1): |N_j - C(m,j)p/2^m| <= C(m,j)(m/2)(sqrt p + 1) + 2m, N_j = n_j + n'_j
    for j in range(m + 1):
        co = {I(j): 1}
        if j < m:
            co[J(j)] = 1
        main = Fr(comb(m, j) * p, 2 ** m)
        err = Fr(comb(m, j) * m * (sq + 1), 2) + 2 * m
        rows.append((co, '<=', main + err, 'weil_count_up_%d' % j))
        rows.append(({k: -v for k, v in co.items()}, '<=', -(main - err), 'weil_count_lo_%d' % j))
    rows.append(({I(0): -1}, '<=', -ntarget, 'target'))
    rows.append(({J(0): 1}, '<=', 0, 'r=0'))
    # 'spread' restriction (not a theorem; it makes the certificate stronger): every x outside
    # B has band*m <= e_x <= (1-band)*m, i.e. |f_A(x)| <= (1-2 band) m + 1 off B.
    lo, hi = math.ceil(band * m), math.floor((1 - band) * m)
    for j in range(1, m + 1):
        if j < lo or j > hi:
            rows.append(({I(j): 1}, '<=', 0, 'spread_%d' % j))
    for j in range(m):
        if j < lo or j > hi:
            rows.append(({J(j): 1}, '<=', 0, 'spreadneg_%d' % j))
    return rows, 2 * m + 1


def lp_certificate(p, m, ntarget):
    from fractions import Fraction as Fr
    from scipy.optimize import linprog
    rows, nv = lp_rows(p, m, ntarget)
    IDENT = ('subsetHP_A_1', 'subsetHP_A_2', 'subsetHP_nuA_1', 'subsetHP_nuA_2', 'target', 'r=0')
    Aeq, beq, Aub, bub = [], [], [], []
    for co, sense, rhs, name in rows:
        v = np.zeros(nv + 1)
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)), 1.0])
        for k, c in co.items():
            v[k] = float(c) / sc
        if sense == '=':
            Aeq.append(v); beq.append(float(rhs) / sc)
        else:
            if name not in IDENT and not name.startswith('spread'):
                v[nv] = 1.0
            Aub.append(v); bub.append(float(rhs) / sc)
    c = np.zeros(nv + 1); c[nv] = -1
    res = linprog(c, A_ub=np.array(Aub), b_ub=np.array(bub), A_eq=np.array(Aeq), b_eq=np.array(beq),
                  bounds=[(0, None)] * nv + [(0, 1)], method='highs-ipm')
    if res.status != 0:
        return None
    x = res.x[:nv]
    Q = 1000
    xr = [Fr(int(round(v * Q)), Q) if v > 1e-7 else Fr(0) for v in x]
    xr[0] = Fr(ntarget)
    xr[m + 1] = Fr(0)
    # exact repair of the four equalities: pivot n'_j (largest) and three largest n_j (j >= 1)
    jp = max(range(m), key=lambda j: x[m + 1 + j])
    xr[m + 1 + jp] = Fr(m) - sum(xr[m + 1 + j] for j in range(m) if j != jp)
    piv = sorted([j for j in range(1, m + 1)], key=lambda j: -x[j])[:3]
    val = lambda j, dl: m - dl - 2 * j
    rest = [j for j in range(m + 1) if j not in piv]
    r1 = p - sum(xr[m + 1 + j] for j in range(m)) - sum(xr[j] for j in rest)
    r2 = 0 - sum(xr[m + 1 + j] * val(j, 1) for j in range(m)) - sum(xr[j] * val(j, 0) for j in rest)
    r3 = m * (p - m) - sum(xr[m + 1 + j] * val(j, 1) ** 2 for j in range(m)) - sum(xr[j] * val(j, 0) ** 2 for j in rest)
    import sympy
    Mx = sympy.Matrix([[1, 1, 1], [val(piv[0], 0), val(piv[1], 0), val(piv[2], 0)],
                       [val(piv[0], 0) ** 2, val(piv[1], 0) ** 2, val(piv[2], 0) ** 2]])
    sol = Mx.LUsolve(sympy.Matrix([sympy.Rational(r1.numerator, r1.denominator),
                                   sympy.Rational(r2.numerator, r2.denominator),
                                   sympy.Rational(r3.numerator, r3.denominator)]))
    for k, j in enumerate(piv):
        xr[j] = Fr(int(sympy.fraction(sol[k])[0]), int(sympy.fraction(sol[k])[1]))
    ok = all(v >= 0 for v in xr)
    check('E.LP_certificate_nonneg', ok, (p, m))
    minslack = None
    for co, sense, rhs, name in rows:
        lhs = sum(Fr(cc) * xr[k] for k, cc in co.items())
        if sense == '=':
            check('E.LP_certificate_equality', lhs == rhs, (p, m, name))
        else:
            good = lhs <= rhs
            check('E.LP_certificate_inequality', good, (p, m, name, float(lhs), float(rhs)))
    support = {('n_%d' % j): float(xr[j]) for j in range(m + 1) if xr[j] != 0}
    support.update({('n\'_%d' % j): float(xr[m + 1 + j]) for j in range(m) if xr[m + 1 + j] != 0})
    d = (p - 1) // 2
    nspread = sum(1 for r_ in rows if r_[3].startswith('spread'))
    return {'p': p, 'm': m, 'n0': ntarget, 'm_n0_over_p': m * ntarget / p, 'g': d - m * ntarget,
            'n_rows': len(rows), 'n_theorem_rows': len(rows) - nspread, 'band': [math.ceil(0.4 * m), math.floor(0.6 * m)],
            'support': support, 'exact_ok': ok,
            'rest_mean_square': float((sum(xr[j] * (m - 2 * j) ** 2 for j in range(1, m + 1)) + sum(xr[m + 1 + j] * (m - 1 - 2 * j) ** 2 for j in range(m)))
                                       / (p - ntarget))}

# ---------------------------------------------------------------- main

def main():
    t0 = time.time()
    rng = random.Random(20260929)
    results = {'sections': {}}
    table = json.load(open(TABLE))['table']

    # ---- A/B/C/D exhaustive small primes
    exh = []
    NREL = [0]
    HX = [0]
    for p in (13, 17, 29, 37, 41):
        chi = chi_table(p)
        for m in range(2, 6):
            for rest in itertools.combinations(range(1, p), m - 1):
                A = [0] + list(rest)
                # normalisation: second element 1 or n0 is enough up to symmetry
                if rest[0] not in (1, next(x for x in range(2, p) if chi[x] == -1)):
                    continue
                if m <= 4 or p <= 29:
                    converse_check(A, p, chi, rng) if p <= 29 else None
                B = B_of_A(A, p, chi)
                if len(B) >= 1:
                    G, g, r = structure(A, B, p, chi)
                    if g <= len(B) - 2:
                        NREL[0] += interpolation_check(A, B, p, chi, G, g)
                    if p <= 29:
                        rational_multiplicity_check(A, p, chi)
                    if len(B) >= 2 and rng.random() < 0.5:
                        Bs = rng.sample(B, rng.randint(1, len(B)))
                        structure(A, Bs, p, chi, full=True)
                    info = section_C(A, B, p, chi, G, g, r, hankel=(p <= 29))
                    for kk in ('hankel_extra_vanishing_e1', 'hankel_extra_vanishing_e2'):
                        if info.get(kk) is not None:
                            HX[0] += 1
                            check('C.no_extra_hankel_vanishing_on_B1', info[kk] == 0, (p, A, kk))
                    # symmetric side: F_B has the same defect
                    if len(B) >= 2 and len(B) <= (p + 1) // 2:
                        GB, gB, rB = structure(B, A, p, chi, full=False)
                        check('A.same_defect_gA=gB', gB == g, (p, A, B))
                    exh.append((p, len(A), len(B), g, info.get('nG_distinct'), info.get('degR')))
                section_D(A, p, chi)
                if m <= 3 and p <= 17:
                    wronskian_identity(A, p)
    results['sections']['A-D_exhaustive'] = {'cases': len(exh), 'interpolation_relations_checked': NREL[0],
                                             'hankel_exact_order_cases': HX[0]}
    grid_sz_check(rng)
    # Delta_e != 0 mod p (exact order of Hankel minors on B_1)
    dz = []
    for p in (13, 29, 101, 509, 997, 1009, 10009):
        d = (p - 1) // 2
        for m in range(3, min(40, (p + 1) // 2) + 1):
            for e in range(1, min(5, (m - 1) // 2) + 1):
                v, nz = delta_e(m, d + m - 1, e, p)
                check('C.Delta_e_nonzero_mod_p', nz, (p, m, e))
                if e == 1:
                    D = d + m - 1
                    check('C.Delta_1_formula', v == sympy_rat(m * (m - D), D * D * (D - 1)), (p, m))
                if not nz:
                    dz.append((p, m, e))
    results['sections']['Delta_e_zero_cases'] = dz
    print('grid SZ, Delta_e done', time.time() - t0)
    print('A-D exhaustive done', time.time() - t0, len(exh))

    # ---- tight examples: m = 2 family and the three exceptional primes
    tight = []
    for p in [5, 13, 17, 29, 37, 41, 53, 61, 101, 197, 509, 997]:
        chi = chi_table(p)
        A = [0, 1]
        B = B_of_A(A, p, chi)
        G, g, r = structure(A, B, p, chi)
        check('A.m2_family_tight_g=0', g == 0 and len(B) == (p + 3) // 4, (p, len(B), g))
        if p % 4 == 1:
            lhs, rhs, N0, Nm = section_D(A, p, chi)
            check('D.derivative_HP_equality_A={0,1}', lhs == rhs and N0 == (p + 3) // 4 and Nm == (p - 1) // 4, (p, lhs, rhs))
        tight.append({'p': p, 'A': A, 'n': len(B), 'r': r, 'g': g})
    for p, A in ((13, [0, 1, 4]), (37, [0, 1, 11]), (41, [0, 1, 2, 10, 33])):
        chi = chi_table(p)
        B = B_of_A(A, p, chi)
        G, g, r = structure(A, B, p, chi)
        info = section_C(A, B, p, chi, G, g, r)
        tight.append({'p': p, 'A': A, 'B': B, 'n': len(B), 'r': r, 'g': g, 'C': info})
        check('A.exceptional_tight_g=0', g == 0, (p, A))
    results['tight_examples'] = tight

    # ---- table configurations: profile extremals, balanced witnesses, cliques, sum-cliques
    conf = []
    for ps in sorted(table, key=int):
        p = int(ps)
        if p > 1000 or p < 29:
            continue
        e = table[ps]
        chi = chi_table(p)
        items = []
        for w in e.get('profile_witnesses', []):
            if w['k'] >= 3:
                items.append(('profile', w['A']))
        for w in e.get('balanced_witnesses', [])[:1]:
            items.append(('balanced', w['A']))
        items.append(('clique', e['clique']))
        if e.get('sumclique_A'):
            items.append(('sumclique', e['sumclique_A']))
        for kind, A in items:
            A = [a % p for a in A]
            B = B_of_A(A, p, chi) if kind != 'clique' else [(-a) % p for a in A]
            if kind == 'clique':
                check('A.clique_is_biclique', all(chi[(a + b) % p] >= 0 for a in A for b in B), (p, A))
            if kind == 'sumclique':
                B = B_of_A(A, p, chi)
            m, n = len(A), len(B)
            if n < 1:
                continue
            G, g, r = structure(A, B, p, chi)
            info = section_C(A, B, p, chi, G, g, r, hankel=(p <= 200 and m >= 3))
            rows = {'p': p, 'kind': kind, 'm': m, 'n': n, 'r': r, 'g': g, 'g_over_d': g / ((p - 1) // 2),
                    'mn_over_p': m * n / p}
            rows.update(info)
            if n >= 2 and n <= (p + 1) // 2:
                GB, gB, rB = structure(B, A, p, chi, full=False)
                check('A.same_defect_gA=gB', gB == g, (p, kind))
                rows['nGB_distinct'] = squarefree_info(GB, p)[0]
            rows['G_roots_in_Fp'] = roots_in_Fp(G, p)
            if p <= 300 or (kind == 'profile' and p in (509, 997)):
                rows['G_pattern'] = ddf_pattern(G, p)
            conf.append(rows)
    results['configurations'] = conf
    print('table configurations done', time.time() - t0, len(conf))

    # ---- summary statistics of the defect polynomial over the configurations
    results['defect_summary'] = {
        'n_configs': len(conf),
        'G_squarefree_all': all(r.get('G_max_mult', 0) <= 1 for r in conf),
        'max_G_roots_in_Fp': max(r['G_roots_in_Fp'] for r in conf),
        'mean_G_roots_in_Fp': sum(r['G_roots_in_Fp'] for r in conf) / len(conf),
        'min_g_over_d_minK3_p>=100': min((r['g_over_d'], r['p'], r['kind']) for r in conf
                                          if min(r['m'], r['n']) >= 3 and r['p'] >= 100),
        'rational_root_histogram': sorted(__import__('collections').Counter(r['G_roots_in_Fp'] for r in conf).items()),
        'non_squarefree_G': [(r['p'], r['kind'], r['m'], r['n'], r['G_max_mult'], r['G_roots_in_Fp'])
                             for r in conf if r.get('G_max_mult', 0) > 1],
        'squarefree_by_kind': {k: (sum(1 for r in conf if r['kind'] == k and r.get('G_max_mult', 0) <= 1),
                                   sum(1 for r in conf if r['kind'] == k)) for k in ('profile', 'balanced', 'clique', 'sumclique')},
    }

    # ---- coefficient formula of section 3.4: coefficient of x^{d-j} in F_A/lam equals
    #      (-1)^j (1/2)^{(j)}/m^{(j)} h_j(A)  (rising factorials), for 1 <= j <= d
    for p in (13, 101, 509):
        d = (p - 1) // 2
        for A in ([0, 1, 4], [0, 1, 3, 7, 12], [0, 2, 5, 9]):
            m = len(A); D = d + m - 1
            F = hp_poly(A, p); lam = binom_mod(D, m - 1, p); h = h_series(A, d, p)
            il = pow(lam, p - 2, p); inv2 = pow(2, p - 2, p)
            for j in range(1, d + 1):
                num = 1; den = 1
                for i in range(j):
                    num = num * ((1 + 2 * i) * inv2) % p
                    den = den * (m + i) % p
                val = pow(-1, j, p) * num * pow(den, p - 2, p) % p * int(h[j]) % p
                coef = int(F[d - j]) * il % p if d - j > 0 else (int(F[0]) + 1) * il % p
                check('E.coefficient_formula', coef == val, (p, A, j))

    # ---- affine involution: A = {0,1,121,122} at p = 509 is fixed by x -> 122 - x;
    #      then F_A(-x-122) = F_A(x) (equivariance, section 3.4)
    for p, A, t in ((509, [0, 1, 121, 122], 122), (101, [0, 1, 3, 4], 4), (509, [0, 1, 2, 21, 490], 2)):
        F = hp_poly(A, p)
        check('E.affine_equivariance', all(peval(F, (-x - t) % p, p) == peval(F, x, p) for x in range(p)), (p, A))

    # ---- (c): points with few bad partners in the extremal configurations vs the threshold
    #      X(eta) >= (3/4) eta^2 p / (1 - 2 eta) needed for a constant gain from (*)_e
    xc = []
    for ps in ('509', '997'):
        p = int(ps); chi = chi_table(p)
        e = table[ps]
        for w in e['profile_witnesses'] + e.get('balanced_witnesses', [])[:1]:
            A = w['A']; m = len(A)
            if m < 5:
                continue
            prof = e_profile(A, p, chi)
            for eta in (0.2, 0.3, 0.4):
                ee = int(eta * m)
                X = sum(1 for (eb, dl) in prof if 1 <= eb <= ee and not dl)
                xc.append({'p': p, 'm': m, 'eta': eta, 'e': ee, 'X': X,
                           'needed': 0.75 * eta ** 2 * p / (1 - 2 * eta)})
            star_slack(A, p, chi)
    results['few_bad_partner_counts'] = xc

    # ---- Weil moments cannot see a balanced near-tight B: n m^{2k} vs (2k-1) m^{2k} sqrt(p)
    wr = []
    for p in (1009, 10009, 100049, 1000033):
        d = (p - 1) // 2; m = round(math.sqrt(p / 2)); n = d // m
        for k in (2, 3, 4):
            ratio = n * m ** (2 * k) / ((2 * k - 1) * m ** (2 * k) * math.sqrt(p))
            wr.append({'p': p, 'm': m, 'k': k, 'B_share_over_weil_error': ratio})
            check('E.weil_moment_blind', ratio < 1, (p, k))
    results['weil_moment_ratio'] = wr

    # ---- Theorem 3.5: exact rational LP certificates
    lps = []
    for p in (1009, 10009, 40009, 100049, 1000033):
        d = (p - 1) // 2; m = round(math.sqrt(p / 2)); nt = d // m
        cert = lp_certificate(p, m, nt)
        check('E.LP_certificate_found', cert is not None, p)
        if cert:
            lps.append(cert)
        print('LP', p, m, None if cert is None else cert['exact_ok'], time.time() - t0)
    results['lp_certificates'] = lps

    # ---- data file (task 2) re-check
    if os.path.exists(DATA):
        dat = json.load(open(DATA))['rows']
        drows = []
        for ps, row in sorted(dat.items(), key=lambda kv: int(kv[0])):
            p = int(ps)
            if 'error' in row:
                check('F.data_row_ok', False, (p, row['error']))
                continue
            chi = chi_table(p)
            for K, ent in row['P'].items():
                A = ent['A']; m = len(A)
                if p <= 3000:
                    B = B_of_A(A, p, chi)
                    check('F.witness_B_size', len(B) == ent['n'], (p, K, A))
                    check('F.witness_product', m * len(B) == ent['T_profile'], (p, K))
                    check('F.witness_min_side>=K', min(m, len(B)) >= int(K), (p, K))
                if ent.get('A_large'):
                    AL = ent['A_large']; BL = B_of_A(AL, p, chi)
                    check('F.large_witness', len(AL) * len(BL) >= ent['P'] and min(len(AL), len(BL)) >= int(K), (p, K))
                drows.append({'p': p, 'K': int(K), 'P': ent['P'], 'ratio': ent['P'] / p, 'exact': ent['exact'],
                              'argm': ent['argm'] if not ent.get('A_large') else len(ent['A_large'])})
            if 'mismatch' in row:
                check('F.profile_matches_table', False, (p, row['mismatch']))
            else:
                check('F.profile_matches_table', True)
        # independent pure-python M_3 for p <= 400
        for ps, row in dat.items():
            p = int(ps)
            if p > 400 or 3 not in [int(k) for k in row['M']]:
                continue
            chi = chi_table(p)
            n0 = next(x for x in range(2, p) if chi[x] == -1)
            best = 0
            for s_ in (1, n0):
                base = set(B_of_A([0, s_], p, chi))
                for c_ in range(1, p):
                    if c_ == s_:
                        continue
                    cnt = sum(1 for b in base if (c_ + b) % p == 0 or chi[(c_ + b) % p] == 1)
                    best = max(best, cnt)
            check('F.M3_pure_python', best == row['M'][str(3)] if str(3) in row['M'] else best == row['M'][3], (p, best))
        results['data_rows'] = drows

    results['checks'] = CHECKS
    results['n_checks'] = sum(CHECKS.values())
    results['failures'] = FAILS[:200]
    results['n_failures'] = len(FAILS)
    results['elapsed'] = time.time() - t0
    json.dump(results, open(OUT, 'w'), indent=1, default=str)
    print('checks', results['n_checks'], 'failures', len(FAILS), 'time', results['elapsed'])


if __name__ == '__main__':
    main()
