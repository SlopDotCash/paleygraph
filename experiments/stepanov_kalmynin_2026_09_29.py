#!/usr/bin/env python3
"""Verifier for research/stepanov-kalmynin-2026-09-29.md (worker `kalmynin`, Stepanov wave).

What is checked (exact arithmetic mod p, or exact rationals):
  A1  Taylor expansion of the Hanson-Petridis polynomial F_A at complete / anti-complete points
      (Theorem 2.1), all orders.
  A2  Relations X and Y of Kalmynin with the defect terms (Theorem 2.2, Corollary 2.3) on complete bicliques
      A + B in Q u {0}; the defect-free versions (Kalmynin's Lemma 10) are refuted for g > 0.
  A3  The Yip-Yoo global second-order identity with defect (Proposition 2.4): P*h | Z, deg W <= 2g-2.
  A4  Moment identities (Kalmynin Lemma 3 / Yip-Yoo (2.6)): true for exact decompositions,
      equivalent to uniform covering (Proposition 2.5), refuted under containment; A4b exhaustive
      list of critical r = 0 pairs (|A| <= 4) with their first non-vanishing moment (Prop 2.6).
  A5  Reach of coefficient comparison (Proposition 2.7): rank of the Hermite + top-coefficient map.
  A6  Moebius transfer for cliques (Kalmynin Lemma 7) and its profile map (Proposition 3.1).
  A7  Reciprocal transfer for bicliques (Rudnev-Tyrrell (1)) and its joint profile map
      (Proposition 3.2); twisted Hanson-Petridis (Proposition 3.3).
  A8  Hankel minors at complete points: next Taylor coefficient = universal constant * kappa (Prop 4.1),
      and the Hankel analogue of Relation X with defect K_1'/K_1 (e = 1).
  A9  Exhaustive search for Hanson-Petridis-tight pairs (g = 0) at small p, by type.
  B1  Clique LP + Moebius-transfer rows: exact rational certificates.
  B2  Joint biclique LP + reciprocal-transfer rows: exact rational certificates.

Run: /opt/miniconda3/bin/python3 experiments/stepanov_kalmynin_2026_09_29.py
Writes results/stepanov_kalmynin_2026_09_29.json
"""
import json
import math
import os
import random
import sys
import time
from fractions import Fraction as Fr

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'results', 'stepanov_kalmynin_2026_09_29.json')

CHECKS = {}
FAILS = []
WITNESS = {}


def check(name, cond, info=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not cond:
        FAILS.append((name, repr(info)[:400]))


def witness(name, item, cap=12):
    WITNESS.setdefault(name, [])
    if len(WITNESS[name]) < cap:
        WITNESS[name].append(item)


# ---------------------------------------------------------------- polynomials mod p (low -> high)

def trim(a):
    a = list(a)
    while a and a[-1] == 0:
        a.pop()
    return a


def deg(a):
    a = trim(a)
    return len(a) - 1


def padd(a, b, p):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def psub(a, b, p):
    n = max(len(a), len(b))
    return trim([((a[i] if i < len(a) else 0) - (b[i] if i < len(b) else 0)) % p for i in range(n)])


def pscale(a, s, p):
    return trim([(x * s) % p for x in a])


def pmul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return trim(r)


def pdivmod(a, b, p):
    a = trim(a)
    b = trim(b)
    if len(a) < len(b):
        return [], a
    inv = pow(b[-1], p - 2, p)
    q = [0] * (len(a) - len(b) + 1)
    r = list(a)
    for i in range(len(a) - len(b), -1, -1):
        c = (r[i + len(b) - 1] * inv) % p
        q[i] = c
        if c:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] - c * y) % p
    return trim(q), trim(r[:len(b) - 1])


def peval(a, x, p):
    v = 0
    for c in reversed(a):
        v = (v * x + c) % p
    return v


def pderiv(a, p):
    return trim([(i * a[i]) % p for i in range(1, len(a))])


def pfromroots(roots, p):
    r = [1]
    for z in roots:
        r = pmul(r, [(-z) % p, 1], p)
    return r


def ppow(a, e, p):
    r = [1]
    b = list(a)
    while e:
        if e & 1:
            r = pmul(r, b, p)
        b = pmul(b, b, p)
        e >>= 1
    return r


def taylor(a, b, p):
    """coefficients of a(b + Y) in Y (Horner shift)."""
    a = trim(a)
    c = list(a)
    n = len(c)
    for i in range(n):
        for j in range(n - 2, i - 1, -1):
            c[j] = (c[j] + b * c[j + 1]) % p
    return c


def inv(x, p):
    return pow(x % p, p - 2, p)


def series_inv(a, N, p):
    """power series inverse of a (a[0] != 0) to N terms."""
    a = list(a) + [0] * N
    r = [0] * N
    i0 = inv(a[0], p)
    r[0] = i0
    for k in range(1, N):
        s = 0
        for j in range(1, k + 1):
            s += a[j] * r[k - j]
        r[k] = (-s * i0) % p
    return r


def series_mul(a, b, N, p):
    r = [0] * N
    for i in range(min(N, len(a))):
        if a[i]:
            for j in range(min(N - i, len(b))):
                r[i + j] = (r[i + j] + a[i] * b[j]) % p
    return r


def binom_mod(n, k, p):
    if k < 0 or k > n:
        return 0
    return math.comb(n, k) % p


def chi_table(p):
    t = [0] * p
    for x in range(1, p):
        t[(x * x) % p] = 1
    return [0] + [1 if t[x] else -1 for x in range(1, p)]


def nonresidue(p, chi):
    for x in range(2, p):
        if chi[x] == -1:
            return x


def lagrange_c(A, p):
    c = []
    for a in A:
        prod = 1
        for a2 in A:
            if a2 != a:
                prod = prod * (a - a2) % p
        c.append(inv(prod, p))
    return c


def hp_poly(A, p, d=None):
    """F_A(x) = -1 + sum_k c_k (x + a_k)^D, D = d + m - 1 (d = (p-1)/2 by default)."""
    m = len(A)
    if d is None:
        d = (p - 1) // 2
    D = d + m - 1
    c = lagrange_c(A, p)
    F = [0] * (D + 1)
    for ck, a in zip(c, A):
        for i in range(D + 1):
            F[i] = (F[i] + ck * binom_mod(D, i, p) * pow(a, D - i, p)) % p
    F[0] = (F[0] - 1) % p
    return trim(F)


def complete_points(A, p, chi):
    """C1(A) = {x not in -A : chi(x+a) = 1 for all a};  C0(A) = {x in -A: x + A in Q u {0}}."""
    negA = set((-a) % p for a in A)
    C1, C0 = [], []
    for x in range(p):
        vals = [chi[(x + a) % p] for a in A]
        if x in negA:
            if all(v >= 0 for v in vals):
                C0.append(x)
        elif all(v == 1 for v in vals):
            C1.append(x)
    return C1, C0


def anti_complete_points(A, p, chi):
    negA = set((-a) % p for a in A)
    return [x for x in range(p) if x not in negA and all(chi[(x + a) % p] == -1 for a in A)]


def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            for j in range(i * i, n + 1, i):
                s[j] = False
    return [i for i in range(n + 1) if s[i]]


# ---------------------------------------------------------------- sampling helpers

def sample_sets(p, rng, sizes, per=6, exhaustive_below=0):
    """sets A containing 0 (translation normalised): all 2-sets {0,t}, plus random ones."""
    out = []
    for m in sizes:
        if m == 2:
            out += [[0, t] for t in range(1, p)]
            continue
        if math.comb(p - 1, m - 1) <= exhaustive_below:
            import itertools
            for T in itertools.combinations(range(1, p), m - 1):
                out.append([0] + list(T))
        else:
            seen = set()
            tries = 0
            while len(seen) < per and tries < 50 * per:
                tries += 1
                T = tuple(sorted(rng.sample(range(1, p), m - 1)))
                seen.add(T)
            out += [[0] + list(T) for T in seen]
    return out


# ---------------------------------------------------------------- A1: Taylor expansion at complete points

def sigma_s(A, x, s, p):
    """sigma_s(x) = sum_k c_k (x + a_k)^(-1-s)."""
    c = lagrange_c(A, p)
    return sum(ck * pow(inv(x + a, p), 1 + s, p) for ck, a in zip(c, A)) % p


def section_A1(rng):
    """Theorem 2.1: for x in C1(A): [Y^j]F_A(x+Y) = 0 (j<m), = C(D,m+s) sigma_s(x) (j = m+s),
    and sigma_s(x) = (-1)^(m-1) [Y^s] 1/h(x-Y),  h(X) = prod (X + a).
    For anti-complete x: [Y^0] = -2, [Y^j] = 0 (1<=j<m), [Y^(m+s)] = -C(D,m+s) sigma_s(x)."""
    stats = {'points': 0, 'coeffs': 0}
    for p in [13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 61]:
        chi = chi_table(p)
        d = (p - 1) // 2
        sets = sample_sets(p, rng, [2, 3, 4, 5], per=5)
        for A in sets:
            m = len(A)
            if m > (p + 1) // 2:
                continue
            D = d + m - 1
            F = hp_poly(A, p)
            check('A1.deg_F_eq_d', deg(F) == d, (p, A))
            h = pfromroots([(-a) % p for a in A], p)
            C1, C0 = complete_points(A, p, chi)
            AC = anti_complete_points(A, p, chi)
            for x, sgn in [(x, 1) for x in C1] + [(x, -1) for x in AC]:
                T = taylor(F, x, p) + [0] * (d + 2)
                # 1/h(x - Y) as a power series
                hx = taylor(h, x, p)
                hminus = [(c * (-1) ** i) % p for i, c in enumerate(hx)]
                ih = series_inv(hminus, d - m + 1, p)
                stats['points'] += 1
                if sgn == 1:
                    check('A1.low_order_zero', all(T[j] == 0 for j in range(m)), (p, A, x))
                else:
                    check('A1.anti_value', T[0] == p - 2 and all(T[j] == 0 for j in range(1, m)), (p, A, x))
                for s in range(0, d - m + 1):
                    sg = sigma_s(A, x, s, p)
                    check('A1.sigma_series', sg == ((-1) ** (m - 1) * ih[s]) % p, (p, A, x, s))
                    check('A1.taylor_coeff', T[m + s] == (sgn * binom_mod(D, m + s, p) * sg) % p, (p, A, x, s))
                    stats['coeffs'] += 1
            # points of -A with x + A in Q u {0}: order exactly m-1, explicit coefficients
            cA = lagrange_c(A, p)
            for x in C0:
                k0 = next(i for i, a in enumerate(A) if (x + a) % p == 0)
                T = taylor(F, x, p) + [0] * (d + 2)
                check('A1.C0_low_order_zero', all(T[j] == 0 for j in range(m - 1)), (p, A, x))
                check('A1.C0_leading', T[m - 1] == (-binom_mod(D, m - 1, p) * cA[k0]) % p, (p, A, x))
                for s in range(1, d - m + 2):
                    val = sum(cA[i] * pow(inv(x + a, p), s, p) for i, a in enumerate(A) if i != k0) % p
                    check('A1.C0_taylor_coeff', T[m - 1 + s] == binom_mod(D, m - 1 + s, p) * val % p, (p, A, x, s))
                    stats['coeffs'] += 1
    return stats


# ---------------------------------------------------------------- A2: Relations X, Y with defect

def biclique_structure(A, B, p, d=None):
    """F_A = lam * P1^m * P0^(m-1) * Gam (exact).  Returns dict or None if not divisible."""
    m = len(A)
    if d is None:
        d = (p - 1) // 2
    D = d + m - 1
    negA = set((-a) % p for a in A)
    B1 = [b for b in B if b not in negA]
    B0 = [b for b in B if b in negA]
    F = hp_poly(A, p, d)
    lam = binom_mod(D, m - 1, p)
    P1 = pfromroots(B1, p)
    P0 = pfromroots(B0, p)
    den = pscale(pmul(ppow(P1, m, p), ppow(P0, m - 1, p), p), lam, p)
    Gam, rem = pdivmod(F, den, p)
    if rem:
        return None
    return {'F': F, 'Gam': Gam, 'B1': B1, 'B0': B0, 'g': deg(Gam), 'r': len(B0), 'D': D, 'lam': lam}


def relation_quantities(A, B, b, st, p, d):
    """Taylor side and factor side of Relations X and Y at b in B1."""
    m = len(A)
    D = d + m - 1
    s0, s1, s2 = (sigma_s(A, b, s, p) for s in range(3))
    t0 = binom_mod(D, m, p) * s0 % p
    t1 = binom_mod(D, m + 1, p) * s1 % p
    t2 = binom_mod(D, m + 2, p) * s2 % p
    i0 = inv(t0, p)
    kX_taylor = t1 * i0 % p
    kY_taylor = (t2 * i0 - t1 * t1 * i0 * i0 * inv(2, p)) % p
    Gam = st['Gam']
    G0 = peval(Gam, b, p)
    G1 = peval(pderiv(Gam, p), b, p)
    G2 = peval(pderiv(pderiv(Gam, p), p), b, p)
    iG = inv(G0, p)
    S1 = sum(inv(b - bb, p) for bb in st['B1'] if bb != b) % p
    S2 = sum(inv(b - bb, p) ** 2 for bb in st['B1'] if bb != b) % p
    R1 = sum(inv(b - bb, p) for bb in st['B0']) % p
    R2 = sum(inv(b - bb, p) ** 2 for bb in st['B0']) % p
    defX = G1 * iG % p
    defY = (G2 * iG - G1 * G1 * iG * iG) * inv(2, p) % p
    kX_factor = (m * S1 + (m - 1) * R1 + defX) % p
    kY_factor = (-inv(2, p) * (m * S2 + (m - 1) * R2) + defY) % p
    hA1 = sum(inv(a + b, p) for a in A) % p          # h'(b)/h(b)
    hA2 = sum(inv(a + b, p) ** 2 for a in A) % p
    return dict(kX_taylor=kX_taylor, kY_taylor=kY_taylor, kX_factor=kX_factor, kY_factor=kY_factor,
                defX=defX, defY=defY, S1=S1, S2=S2, hA1=hA1, hA2=hA2, s0=s0, s1=s1)


def kalmynin_X_Y(A, B, b, p, d):
    """Kalmynin Lemma 10 (Relations X and Y), stated for A + B = mu_d, alpha = |A|:
       X: sum_a 1/(a+b) = alpha(alpha+1)/(d-1) * sum_{b' != b} 1/(b-b')
       Y: (sum 1/(a+b))^2 + sum 1/(a+b)^2 = alpha(alpha+1)(alpha+2)/((d-1)(d-2)) *
          (alpha (sum 1/(b-b'))^2 - sum 1/(b-b')^2)."""
    al = len(A)
    s1a = sum(inv(a + b, p) for a in A) % p
    s2a = sum(inv(a + b, p) ** 2 for a in A) % p
    s1b = sum(inv(b - bb, p) for bb in B if bb != b) % p
    s2b = sum(inv(b - bb, p) ** 2 for bb in B if bb != b) % p
    X = (s1a - al * (al + 1) * inv(d - 1, p) * s1b) % p == 0
    Y = (s1a * s1a + s2a - al * (al + 1) * (al + 2) * inv((d - 1) * (d - 2), p) * (al * s1b * s1b - s2b)) % p == 0
    return X, Y


def section_A2(rng):
    """Relations X_Gamma, Y_Gamma (Theorem 2.2) on complete bicliques; Kalmynin's defect-free
    Relations X, Y refuted for g > 0; the identity (d-1)/(m+1) sum_{a,b} 1/(a+b) = sum_b Gam'/Gam(b)."""
    stats = {'bicliques': 0, 'points': 0, 'X_defectfree_fail_g_pos': 0, 'X_defectfree_hold_g_pos': 0,
             'Y_defectfree_fail_g_pos': 0, 'by_g0': 0}
    for p in [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 19, 23, 31, 43, 47, 59, 67]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for A in sample_sets(p, rng, [2, 3, 4, 5, 6], per=6):
            m = len(A)
            C1, C0 = complete_points(A, p, chi)
            Bfull = C1 + C0
            if len(C1) < 2:
                continue
            choices = [Bfull, C1]
            if len(C1) >= 4:
                choices.append(rng.sample(C1, len(C1) // 2) + C0[:1])
            for B in choices:
                st = biclique_structure(A, B, p)
                check('A2.divisible', st is not None, (p, A, B))
                if st is None:
                    continue
                n, r, g = len(B), st['r'], st['g']
                check('A2.defect_formula', g == d - m * n + r, (p, A, B, g))
                stats['bicliques'] += 1
                # top coefficient (x^{d-1}):  (d/m) p1(A) = -m p1(B1) - (m-1) p1(B0) + [x^{g-1}] Gam
                gtop = st['Gam'][g - 1] if g >= 1 else 0
                lhs = d * inv(m, p) * sum(A) % p
                rhs = (-m * sum(st['B1']) - (m - 1) * sum(st['B0']) + gtop) % p
                check('A2.top_coefficient_M1', lhs == rhs, (p, A, B))
                if g == 0 and r == 0:
                    check('A2.critical_M1_zero', (n * sum(A) + m * sum(B)) % p == 0, (p, A, B))
                sumab = 0
                sumdef = 0
                for b in st['B1']:
                    q = relation_quantities(A, B, b, st, p, d)
                    stats['points'] += 1
                    check('A2.X_taylor_is_logderiv', q['kX_taylor'] == (d - 1) * inv(m + 1, p) * q['hA1'] % p, (p, A, b))
                    s2 = sigma_s(A, b, 2, p)
                    check('A2.sigma2_over_sigma0', s2 * inv(q['s0'], p) % p == (q['hA1'] ** 2 + q['hA2']) * inv(2, p) % p, (p, A, b))
                    check('A2.X_Gamma', q['kX_taylor'] == q['kX_factor'], (p, A, B, b))
                    check('A2.Y_Gamma', q['kY_taylor'] == q['kY_factor'], (p, A, B, b))
                    sumab += q['hA1']
                    sumdef += q['defX'] - (m - 1) * sum(inv(b - bb, p) for bb in st['B0'])
                    if r == 0:
                        X, Y = kalmynin_X_Y(A, B, b, p, d)
                        if g == 0:
                            check('A2.kalmynin_XY_at_g0', X and Y, (p, A, B, b))
                            stats['by_g0'] += 1
                        else:
                            if not X:
                                stats['X_defectfree_fail_g_pos'] += 1
                                witness('A2.relation_X_fails_under_containment', {'p': p, 'A': A, 'B': B, 'b': b, 'g': g})
                            else:
                                stats['X_defectfree_hold_g_pos'] += 1
                            if not Y:
                                stats['Y_defectfree_fail_g_pos'] += 1
                # summed relation: ((d-1)/(m+1)) sum_{a,b in B1} 1/(a+b)
                #   = sum_b [Gam'/Gam(b) + (m-1) sum_{B0} 1/(b-b0)] + m*0 - (m-1) ... (only r=0 recorded)
                if r == 0:
                    check('A2.summed_X', (d - 1) * inv(m + 1, p) * sumab % p == sumdef % p, (p, A, B))
                    if g == 0:
                        check('A2.critical_M_dminus1_zero', sumab % p == 0 and
                              sum(pow((a + b) % p, d - 1, p) for a in A for b in B) % p == 0, (p, A, B))
    return stats


def section_A2_mu4():
    """Kalmynin's Relations X, Y and the Yip-Yoo identity on the exact decompositions
    mu_4 = {0, -1-i} + {1, i} (p = 1 mod 4), with d = 4."""
    out = []
    for p in [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113]:
        i = next(x for x in range(2, p) if x * x % p == p - 1)
        A = [0, (-1 - i) % p]
        B = [1, i]
        H = sorted(((a + b) % p for a in A for b in B))
        mu4 = sorted([1, p - 1, i, p - i])
        check('A2.mu4_decomposition', H == mu4, p)
        d = 4
        st = biclique_structure(A, B, p, d)
        check('A2.mu4_critical_factorisation', st is not None and st['g'] == 0 and st['r'] == 0, p)
        for b in B:
            X, Y = kalmynin_X_Y(A, B, b, p, d)
            check('A2.mu4_relation_X', X, (p, b))
            check('A2.mu4_relation_Y', Y, (p, b))
        out.append({'p': p, 'A': A, 'B': B})
    return out


# ---------------------------------------------------------------- A3: Yip-Yoo second-order identity with defect

def yy_defect(A, B, p, d=None):
    """r = 0 only.  P = prod (X - b), h = prod (X + a), Gam_A = cofactor of F_A,
    GtB(X) = Gam_B(-X) with Gam_B the cofactor of F_B (HP polynomial of B, vanishing on A).
    Z = Q*Gam_A*GtB + 2(m+1) h P' Gam_A' GtB + 2(n+1) P h' Gam_A GtB',
    Q = m(m+1) h P'' + n(n+1) P h'' - 2(d-1) P' h'.   Claim: P*h | Z, deg(Z/(P h)) <= 2g - 2."""
    m, n = len(A), len(B)
    if d is None:
        d = (p - 1) // 2
    stA = biclique_structure(A, B, p, d)
    stB = biclique_structure(B, A, p, d)
    if stA is None or stB is None or stA['r'] != 0:
        return None
    g = stA['g']
    P = pfromroots(B, p)
    h = pfromroots([(-a) % p for a in A], p)
    GA = stA['Gam']
    GB = stB['Gam']
    GtB = [(c * (-1) ** i) % p for i, c in enumerate(GB)]
    dP, dh = pderiv(P, p), pderiv(h, p)
    ddP, ddh = pderiv(dP, p), pderiv(dh, p)
    Q = psub(padd(pscale(pmul(h, ddP, p), m * (m + 1), p), pscale(pmul(P, ddh, p), n * (n + 1), p), p),
             pscale(pmul(dP, dh, p), 2 * (d - 1), p), p)
    Z = pmul(pmul(Q, GA, p), GtB, p)
    Z = padd(Z, pscale(pmul(pmul(pmul(h, dP, p), pderiv(GA, p), p), GtB, p), 2 * (m + 1), p), p)
    Z = padd(Z, pscale(pmul(pmul(pmul(P, dh, p), GA, p), pderiv(GtB, p), p), 2 * (n + 1), p), p)
    W, rem = pdivmod(Z, pmul(P, h, p), p)
    return {'g': g, 'gB': stB['g'], 'rem_zero': not rem, 'degW': deg(W), 'Q_zero': not Q, 'W': W}


def section_A3(rng):
    stats = {'cases': 0, 'g0_cases': 0, 'yy29_holds_g0': 0, 'yy29_fails_gpos': 0, 'max_degW_minus_bound': -10 ** 9}
    for p in [13, 17, 29, 37, 41, 53, 61, 19, 23, 31, 43, 47]:
        chi = chi_table(p)
        for A in sample_sets(p, rng, [2, 3, 4, 5], per=6):
            C1, C0 = complete_points(A, p, chi)
            if len(C1) < 2:
                continue
            for B in [C1] + ([rng.sample(C1, max(2, len(C1) // 2))] if len(C1) >= 4 else []):
                res = yy_defect(A, B, p)
                if res is None:
                    continue
                stats['cases'] += 1
                check('A3.same_defect_both_sides', res['g'] == res['gB'], (p, A, B))
                check('A3.Ph_divides_Z', res['rem_zero'], (p, A, B))
                bound = 2 * res['g'] - 2
                check('A3.degW_bound', (not res['W']) or res['degW'] <= bound, (p, A, B, res['degW'], res['g']))
                stats['max_degW_minus_bound'] = max(stats['max_degW_minus_bound'], res['degW'] - bound if res['W'] else -99)
                if res['g'] == 0:
                    stats['g0_cases'] += 1
                    check('A3.yy29_at_g0', res['Q_zero'], (p, A, B))
                    stats['yy29_holds_g0'] += res['Q_zero']
                    witness('A3.yy29_holds_for_critical_containment_pairs', {'p': p, 'A': A, 'B': B})
                else:
                    if not res['Q_zero']:
                        stats['yy29_fails_gpos'] += 1
                        witness('A3.yy29_fails_under_containment', {'p': p, 'A': A, 'B': B, 'g': res['g']}, cap=6)
    return stats


# ---------------------------------------------------------------- A4: moment identities

def moments(A, B, p, J):
    return [sum(pow((a + b) % p, j, p) for a in A for b in B) % p for j in range(1, J + 1)]


def section_A4(rng):
    """(i) exact decompositions mu_4 = A + B satisfy M_j = sum (a+b)^j = 0, 1 <= j < d.
    (ii) Proposition 2.4: for a multiset supported on a subgroup H of order d (multiplicities < p),
    M_j = 0 for all 1 <= j <= d-1 iff the multiplicity function is constant.
    (iii) Hanson-Petridis-tight containment pairs (g = 0, r = 0) violate the moment identities and
    Kalmynin's alpha = beta conclusion (witnesses)."""
    stats = {'mu4': 0, 'uniform_equiv_tests': 0, 'critical_pairs_tested': 0, 'critical_pairs_moment_fail': 0,
             'critical_pairs_alpha_ne_beta': 0}
    for p in [13, 17, 29, 37, 41, 53]:
        i = next(x for x in range(2, p) if x * x % p == p - 1)
        A, B = [0, (-1 - i) % p], [1, i]
        check('A4.mu4_moments', all(v == 0 for v in moments(A, B, p, 3)), p)
        stats['mu4'] += 1
    # (ii) random multiplicity functions on H = Q
    for p in [13, 17, 29, 37]:
        chi = chi_table(p)
        H = [x for x in range(1, p) if chi[x] == 1]
        d = len(H)
        for trial in range(60):
            if trial % 3 == 0:
                c = rng.randint(1, 4)
                mu = {x: c for x in H}
            else:
                mu = {x: rng.randint(0, 3) for x in H}
            Ms = [sum(mu[x] * pow(x, j, p) for x in H) % p for j in range(1, d)]
            const = len(set(mu.values())) == 1
            check('A4.moments_iff_uniform', (all(v == 0 for v in Ms)) == const, (p, trial))
            stats['uniform_equiv_tests'] += 1
    # (iii) critical containment pairs
    for p in [13, 17, 19, 29, 37, 41, 53, 61, 73]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for A in sample_sets(p, rng, [2, 3, 4], per=30):
            C1, C0 = complete_points(A, p, chi)
            B = C1
            m, n = len(A), len(B)
            if n < 2 or m * n != d:
                continue
            stats['critical_pairs_tested'] += 1
            Ms = moments(A, B, p, d - 1)
            sums = sorted((a + b) % p for a in A for b in B)
            distinct = len(set(sums)) == len(sums)
            check('A4.critical_pair_not_decomposition', not distinct, (p, A, B))
            if any(Ms):
                stats['critical_pairs_moment_fail'] += 1
                j0 = 1 + next(j for j, v in enumerate(Ms) if v)
                witness('A4.critical_pair_moment_identity_fails', {'p': p, 'A': A, 'B': B, 'first_nonzero_j': j0, 'M_j': Ms[j0 - 1]})
            if m != n:
                stats['critical_pairs_alpha_ne_beta'] += 1
                witness('A4.critical_pair_alpha_ne_beta', {'p': p, 'A': A, 'B': B, 'm': m, 'n': n}, cap=6)
    return stats


def section_A4b():
    """exhaustive over A containing 0, |A| <= 4, p in a list: every critical pair (A, C1(A)) with r = 0
    (m |C1| = d, C0 = empty is automatic since m|C1| + (m-1)|C0| <= d); first non-vanishing moment."""
    import itertools
    from collections import Counter
    res = Counter()
    for p in [13, 17, 19, 29, 37, 41, 53, 61, 73]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for m in [2, 3, 4]:
            for T in itertools.combinations(range(1, p), m - 1):
                A = [0] + list(T)
                C1, C0 = complete_points(A, p, chi)
                if len(C1) < 2 or m * len(C1) != d:
                    continue
                check('A4b.r_zero', len(C0) == 0, (p, A))
                Ms = moments(A, C1, p, d - 1)
                check('A4b.M1_and_Mdm1_zero', Ms[0] == 0 and Ms[d - 2] == 0, (p, A))
                j0 = 1 + next(j for j, v in enumerate(Ms) if v)
                res['p=%d m=%d n=%d first_nonzero_j=%d' % (p, m, len(C1), j0)] += 1
    return dict(res)


# ---------------------------------------------------------------- A5: reach of coefficient comparison

def rank_mod_p(M, p):
    M = [list(r) for r in M]
    rk = 0
    cols = len(M[0]) if M else 0
    for c in range(cols):
        piv = None
        for i in range(rk, len(M)):
            if M[i][c] % p:
                piv = i
                break
        if piv is None:
            continue
        M[rk], M[piv] = M[piv], M[rk]
        iv = inv(M[rk][c], p)
        M[rk] = [(x * iv) % p for x in M[rk]]
        for i in range(len(M)):
            if i != rk and M[i][c] % p:
                f = M[i][c]
                M[i] = [(x - f * y) % p for x, y in zip(M[i], M[rk])]
        rk += 1
    return rk


def section_A5(rng):
    """Proposition 2.5: the linear map (c_0..c_{g-1}) -> (Taylor coefficients of order <= s at each
    b in B of sum c_i x^i ; c_{g-1},...,c_{g-J}) has rank min(g, (s+1)n + J) (J <= g)."""
    stats = {'cases': 0}
    for p in [101, 211, 409]:
        for trial in range(40):
            n = rng.randint(1, 6)
            s = rng.randint(0, 3)
            J = rng.randint(0, 4)
            g = rng.randint(J, (s + 1) * n + J + 6)
            B = rng.sample(range(p), n)
            rows = []
            for b in B:
                for t in range(s + 1):
                    # coefficient of Y^t in (b+Y)^i is C(i,t) b^(i-t)
                    rows.append([binom_mod(i, t, p) * pow(b, i - t, p) % p if i >= t else 0 for i in range(g)])
            for j in range(1, J + 1):
                rows.append([1 if i == g - j else 0 for i in range(g)])
            rk = rank_mod_p(rows, p) if g > 0 else 0
            check('A5.rank', rk == min(g, (s + 1) * n + J), (p, n, s, J, g, rk))
            stats['cases'] += 1
    return stats


# ---------------------------------------------------------------- A6: Moebius transfer for cliques

def find_cliques(p, chi, rng, count=20, target=None):
    """random greedy cliques of the Paley graph (p = 1 mod 4)."""
    out = []
    for _ in range(count * 5):
        order = list(range(p))
        rng.shuffle(order)
        C = []
        for x in order:
            if all(chi[(x - y) % p] == 1 for y in C):
                C.append(x)
        if len(C) >= 3:
            out.append(sorted(C))
        if len(out) >= count:
            break
    return out


def clique_profile(A, p, chi):
    """difference form: e_x = #{a in A : chi(x - a) = -1} for x not in A."""
    Aset = set(A)
    m = len(A)
    prof = [0] * (m + 2)
    for x in range(p):
        if x in Aset:
            continue
        prof[sum(1 for a in A if chi[(x - a) % p] == -1)] += 1
    return prof


def moebius_T(prof, m):
    """averaged transferred profile: T(n)_j = [(m-j) n_j + (m+1-j) n_{m+1-j}] / m."""
    n = list(prof) + [0] * 3
    return [Fr((m - j) * n[j] + ((m + 1 - j) * n[m + 1 - j] if 1 <= m + 1 - j <= m else 0), m) for j in range(m + 1)]


def section_A6(rng):
    stats = {'cliques': 0, 'transfers': 0, 'points': 0}
    for p in [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157]:
        chi = chi_table(p)
        for A in find_cliques(p, chi, rng, count=8):
            m = len(A)
            stats['cliques'] += 1
            tot = [Fr(0)] * (m + 1)
            for a0 in A:
                Aa = sorted([0] + [inv(a0 - a, p) for a in A if a != a0])
                check('A6.transfer_is_clique', all(chi[(x - y) % p] == 1 for x in Aa for y in Aa if x != y), (p, A, a0))
                check('A6.transfer_size', len(set(Aa)) == m, (p, A, a0))
                stats['transfers'] += 1
                Aas = set(Aa)
                for x in range(p):
                    if x in A:
                        continue
                    y = inv(a0 - x, p)
                    check('A6.point_map_off_clique', y not in Aas, (p, A, a0, x))
                    ex = sum(1 for a in A if chi[(x - a) % p] == -1)
                    ey = sum(1 for a in Aa if chi[(y - a) % p] == -1)
                    pred = ex if chi[(x - a0) % p] == 1 else m + 1 - ex
                    check('A6.profile_rule', ey == pred, (p, A, a0, x))
                    stats['points'] += 1
                pa = clique_profile(Aa, p, chi)
                for j in range(m + 1):
                    tot[j] += pa[j]
            T = moebius_T(clique_profile(A, p, chi), m)
            check('A6.averaged_profile', all(tot[j] / m == T[j] for j in range(m + 1)), (p, A))
    return stats


# ---------------------------------------------------------------- A7: reciprocal transfer for bicliques

def joint_counts(A, B, p, chi):
    """for every x: (u_x, v_x, delta_x, eps_x); u = #{a: chi(x+a) = -1}, v = #{b: chi(x-b) = -1}."""
    negA = set((-a) % p for a in A)
    Bs = set(B)
    out = []
    for x in range(p):
        u = sum(1 for a in A if chi[(x + a) % p] == -1)
        v = sum(1 for b in B if chi[(x - b) % p] == -1)
        out.append((u, v, int(x in negA), int(x in Bs)))
    return out


def sum_profile(A, p, chi):
    """sum form: n[j] (x not in -A), n1[j] (x in -A), j = #{a : chi(x+a) = -1}."""
    m = len(A)
    negA = set((-a) % p for a in A)
    n0 = [0] * (m + 1)
    n1 = [0] * (m + 1)
    for x in range(p):
        j = sum(1 for a in A if chi[(x + a) % p] == -1)
        (n1 if x in negA else n0)[j] += 1
    return n0, n1


def transferred_A_profile(JC, m, n):
    """averaged profile of S_b = 1/(A+b) over b in B (Proposition 3.2), from the joint counts."""
    t0 = [Fr(0)] * (m + 1)
    t1 = [Fr(0)] * (m + 1)
    t0[0] += 1
    for (u, v, dl, ep) in JC:
        adj = n - ep - v          # #{b != x : chi(x-b) = 1}
        tgt = t1 if dl else t0
        tgt[u] += Fr(adj, n)
        tgt[m - dl - u] += Fr(v, n)
    return t0, t1


def transferred_B_profile(JC, m, n):
    """averaged sum-form profile of T_a = 1/(B+a) over a in A (p = 1 mod 4), from the joint counts
    at x (the point y' = 1/(y-a) with y = -x); the extra point y' = 0 is spread as 1/m per a."""
    t0 = [Fr(0)] * (n + 1)
    t1 = [Fr(0)] * (n + 1)
    t0[0] += 1
    for (u, v, dl, ep) in JC:
        adj = m - dl - u
        tgt = t1 if ep else t0
        tgt[v] += Fr(adj, m)
        tgt[n - ep - v] += Fr(u, m)
    return t0, t1


def section_A7(rng):
    """(a) chi(x' + s) = chi(x+a) chi(x-b) for x' = 1/(x-b), s = 1/(a+b), b in C1(A);
    (b) the averaged transferred profile formula (Prop 3.2), for p = 1 mod 4 and r = 0;
    (c) Chung: sum_x f_A(x) g_B(x) = p r - m n;
    (d) twisted HP (Rudnev-Tyrrell Lemma 3.4, containment form): deg = d-1, orders >= m on T, >= m-1 at 0,
        and T_max = {1/(x-b) : x in C1(A), x != b}."""
    stats = {'bicliques': 0, 'point_checks': 0, 'twisted': 0}
    for p in [13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 19, 23, 31, 43]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for A in sample_sets(p, rng, [2, 3, 4, 5], per=4):
            m = len(A)
            C1, C0 = complete_points(A, p, chi)
            if len(C1) < 2:
                continue
            for B in [C1, C1 + C0]:
                n = len(B)
                stats['bicliques'] += 1
                negA = set((-a) % p for a in A)
                r = sum(1 for b in B if b in negA)
                JC = joint_counts(A, B, p, chi)
                fg = sum((m - dl - 2 * u) * (n - ep - 2 * v) for (u, v, dl, ep) in JC)
                check('A7.chung_joint', fg == p * r - m * n, (p, A, B))
                for b in [bb for bb in B if bb not in negA][:3]:
                    S = [inv(a + b, p) for a in A]
                    check('A7.S_in_Q', all(chi[s] == 1 for s in S), (p, A, b))
                    for x in range(p):
                        if x == b:
                            continue
                        xp = inv(x - b, p)
                        for a, s in zip(A, S):
                            check('A7.sign_rule', chi[(xp + s) % p] == chi[(x + a) % p] * chi[(x - b) % p], (p, A, b, x))
                            stats['point_checks'] += 1
                if p % 4 == 1 and r == 0:
                    t0, t1 = transferred_A_profile(JC, m, n)
                    a0 = [Fr(0)] * (m + 1)
                    a1 = [Fr(0)] * (m + 1)
                    for b in B:
                        S = [inv(a + b, p) for a in A]
                        s0, s1 = sum_profile(S, p, chi)
                        for j in range(m + 1):
                            a0[j] += Fr(s0[j], n)
                            a1[j] += Fr(s1[j], n)
                    check('A7.averaged_transfer', a0 == t0 and a1 == t1, (p, A, B))
                    u0, u1 = transferred_B_profile(JC, m, n)
                    c0 = [Fr(0)] * (n + 1)
                    c1 = [Fr(0)] * (n + 1)
                    for a in A:
                        T = [inv(b + a, p) for b in B]
                        s0, s1 = sum_profile(T, p, chi)
                        for k in range(n + 1):
                            c0[k] += Fr(s0[k], m)
                            c1[k] += Fr(s1[k], m)
                    check('A7.averaged_transfer_B', c0 == u0 and c1 == u1, (p, A, B))
                # twisted HP
                b = C1[0]
                S = [inv(a + b, p) for a in A]
                Tmax = [t for t in range(1, p) if all(chi[(s + t) % p] == chi[t] for s in S)]
                pred = sorted(inv(x - b, p) for x in C1 if x != b)
                check('A7.Tmax_is_complete_points', sorted(Tmax) == pred, (p, A, b))
                D = d + m - 1
                cS = lagrange_c(S, p)
                k = [c * inv(s, p) % p for c, s in zip(cS, S)]
                CS = sum(k) % p
                Ftw = [0] * (D + 1)
                for kk, s in zip(k, S):
                    for i in range(D + 1):
                        Ftw[i] = (Ftw[i] + kk * binom_mod(D, i, p) * pow(s, D - i, p)) % p
                Ftw[D] = (Ftw[D] - CS) % p
                Ftw = trim(Ftw)
                check('A7.twisted_degree', deg(Ftw) == d - 1, (p, A, b))
                T0 = taylor(Ftw, 0, p) + [0] * m
                check('A7.twisted_order_at_0', all(T0[j] == 0 for j in range(m - 1)), (p, A, b))
                for t in Tmax:
                    Tt = taylor(Ftw, t, p) + [0] * m
                    check('A7.twisted_order_on_T', all(Tt[j] == 0 for j in range(m)), (p, A, b, t))
                check('A7.twisted_count_equals_HP', m * len(Tmax) + m - 1 <= d - 1, (p, A, b))
                stats['twisted'] += 1
    return stats


# ---------------------------------------------------------------- A8: Hankel minors at complete points

def ff(x, k):
    r = 1
    for i in range(k):
        r *= (x - i)
    return r


def frac_det(M):
    import sympy
    return sympy.Matrix(M).det()


def hankel_constants(m, D, e):
    import sympy
    N0 = [[sympy.Rational(ff(m, i + j), ff(D, i + j)) for j in range(e + 1)] for i in range(e + 1)]
    N1 = [[sympy.Rational(ff(m + 1, i + j), ff(D, i + j)) for j in range(e + 1)] for i in range(e + 1)]
    Delta = sympy.Matrix(N0).det()
    E = 0
    for i in range(e + 1):
        M = [row[:] for row in N0]
        M[i] = N1[i][:]
        E += sympy.Matrix(M).det()
    return sympy.Rational(Delta), sympy.Rational(E)


def rat_mod(q, p):
    import sympy
    q = sympy.Rational(q)
    return int(q.p) % p * inv(int(q.q), p) % p


def series_det(M, N, p):
    """determinant of a small matrix of truncated power series (lists of length N) mod p."""
    k = len(M)
    import itertools
    tot = [0] * N
    for perm in itertools.permutations(range(k)):
        sgn = 1
        pl = list(perm)
        for i in range(k):
            for j in range(i + 1, k):
                if pl[i] > pl[j]:
                    sgn = -sgn
        term = [1] + [0] * (N - 1)
        for i in range(k):
            term = series_mul(term, M[i][perm[i]], N, p)
        for t in range(N):
            tot[t] = (tot[t] + sgn * term[t]) % p
    return tot


def section_A8(rng):
    """Proposition 4.1: at b in C1(A), H_{e+1}(b+Y) = phi^{e+1} Y^{(e+1)(m-e)} [Delta_e + kappa E_e Y + O(Y^2)]
    with phi = [Y^m]F(b+Y), kappa = [Y^{m+1}]/[Y^m] = ((d-1)/(m+1)) h'(b)/h(b)."""
    stats = {'points': 0}
    for p in [29, 37, 41, 53, 61, 101, 109]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for A in sample_sets(p, rng, [3, 4, 5, 6], per=3):
            m = len(A)
            D = d + m - 1
            F = hp_poly(A, p)
            C1, C0 = complete_points(A, p, chi)
            for b in C1[:2]:
                T = taylor(F, b, p) + [0] * (3 * m)
                for e in (1, 2):
                    if 2 * e > m - 1:
                        continue
                    order = (e + 1) * (m - e)
                    N = order + 2
                    # u_s(b+Y) = sum_j T_j (j)_s/(D)_s Y^{j-s}
                    u = []
                    for s in range(2 * e + 1):
                        cs = inv(ff(D, s) % p, p)
                        u.append([T[j + s] * (ff(j + s, s) % p) * cs % p for j in range(N)])
                    M = [[u[i + j] for j in range(e + 1)] for i in range(e + 1)]
                    H = series_det(M, N, p)
                    Dl, El = hankel_constants(m, D, e)
                    phi = T[m]
                    kap = T[m + 1] * inv(T[m], p) % p
                    check('A8.low_orders_zero', all(H[t] == 0 for t in range(order)), (p, A, b, e))
                    check('A8.leading', H[order] == pow(phi, e + 1, p) * rat_mod(Dl, p) % p, (p, A, b, e))
                    check('A8.next', H[order + 1] == pow(phi, e + 1, p) * kap % p * rat_mod(El, p) % p, (p, A, b, e))
                    check('A8.kappa_is_relation_X', kap == (d - 1) * inv(m + 1, p) * sum(inv(a + b, p) for a in A) % p, (p, A, b))
                    stats['points'] += 1
    # Hankel analogue of X_Gamma: H_2 = P1^{2(m-1)} K_1 with B = C1(A); at b in B:
    #   (E_1/Delta_1) kappa(b) = 2(m-1) sum_{b' != b} 1/(b-b') + K_1'(b)/K_1(b)
    stats['hankel_X'] = 0
    for p in [29, 37, 41, 53, 61]:
        chi = chi_table(p)
        d = (p - 1) // 2
        for A in sample_sets(p, rng, [3, 4, 5], per=3):
            m = len(A)
            D = d + m - 1
            C1, C0 = complete_points(A, p, chi)
            if len(C1) < 2:
                continue
            F = hp_poly(A, p)
            u = [F]
            for s_ in range(1, 3):
                u.append(pderiv(u[-1], p))
            u = [pscale(u[k], inv(ff(D, k) % p, p), p) for k in range(3)]
            H2 = psub(pmul(u[0], u[2], p), pmul(u[1], u[1], p), p)
            P1 = pfromroots(C1, p)
            K1, rem = pdivmod(H2, ppow(P1, 2 * (m - 1), p), p)
            check('A8.H2_divisible', not rem, (p, A))
            Dl, El = hankel_constants(m, D, 1)
            ratio = rat_mod(El, p) * inv(rat_mod(Dl, p), p) % p
            dK = pderiv(K1, p)
            for b in C1:
                kap = (d - 1) * inv(m + 1, p) * sum(inv(a + b, p) for a in A) % p
                Kb = peval(K1, b, p)
                check('A8.K1_nonzero_on_B', Kb != 0, (p, A, b))
                rhs = (2 * (m - 1) * sum(inv(b - bb, p) for bb in C1 if bb != b) + peval(dK, b, p) * inv(Kb, p)) % p
                check('A8.hankel_X_Gamma', ratio * kap % p == rhs, (p, A, b))
                stats['hankel_X'] += 1
    return stats


# ---------------------------------------------------------------- A9: exhaustive search for HP-tight pairs

def section_A9(pmax=61, mmax=5, budget=400000):
    """All A containing 0 with 2 <= |A| <= mmax (as budget allows) for primes 11 <= p <= pmax:
    Hanson-Petridis-tight pairs (A, C(A)), i.e. m|C1| + (m-1)|C0| = d, recorded by type."""
    import itertools
    found = []
    summary = {}
    for p in primes_upto(pmax):
        if p < 11:
            continue
        chi = chi_table(p)
        d = (p - 1) // 2
        good = np.array([chi[x] == 1 for x in range(p)])
        good0 = np.array([chi[x] >= 0 for x in range(p)])
        xs = np.arange(p)
        for m in range(2, mmax + 1):
            if math.comb(p - 1, m - 1) > budget:
                continue
            cnt_r0 = cnt_rpos = 0
            for T in itertools.combinations(range(1, p), m - 1):
                A = (0,) + T
                maskC = np.ones(p, dtype=bool)
                for a in A:
                    maskC &= good0[(xs + a) % p]
                if not maskC.any():
                    continue
                negA = np.zeros(p, dtype=bool)
                negA[[(-a) % p for a in A]] = True
                c0 = int((maskC & negA).sum())
                c1 = int((maskC & ~negA).sum())
                n = c0 + c1
                if n >= 2 and m * c1 + (m - 1) * c0 == d:
                    if c0 == 0:
                        cnt_r0 += 1
                    else:
                        cnt_rpos += 1
                    if min(m, n) >= 3 and len(found) < 400:
                        B = [int(x) for x in np.nonzero(maskC)[0]]
                        found.append({'p': p, 'A': list(A), 'B': B, 'm': m, 'n': n, 'r': c0})
            summary['%d_%d' % (p, m)] = {'tight_r0': cnt_r0, 'tight_rpos': cnt_rpos}
    return summary, found


# ---------------------------------------------------------------- B: LP machinery

def U_rows(p, m, P0, P1, tag):
    """The 'universal' count constraints of research/stepanov-balanced-2026-09-29.md, Theorem 3.5,
    valid for every set of size m in F_p (sum form): P0[j] = #{x not in -A : e_x = j} (j = 0..m),
    P1[j] = #{x in -A : e_x = j} (j = 0..m-1).  P0, P1 are lists of variable indices.
    Returns rows (coef dict, sense, rhs, name, is_theorem_inequality)."""
    from math import comb, isqrt
    d = (p - 1) // 2
    val = lambda j, dl: m - dl - 2 * j
    rows = []

    def row(co0, co1, sense, rhs, name, thm=True):
        co = {}
        for j, c in co0.items():
            if c:
                co[P0[j]] = co.get(P0[j], 0) + c
        for j, c in co1.items():
            if c:
                co[P1[j]] = co.get(P1[j], 0) + c
        rows.append((co, sense, rhs, tag + ':' + name, thm and sense == '<='))

    row({j: 1 for j in range(m + 1)}, {j: 1 for j in range(m)}, '=', p, 'total')
    row({}, {j: 1 for j in range(m)}, '=', m, 'negA')
    row({j: val(j, 0) for j in range(m + 1)}, {j: val(j, 1) for j in range(m)}, '=', 0, 'moment1')
    row({j: val(j, 0) ** 2 for j in range(m + 1)}, {j: val(j, 1) ** 2 for j in range(m)}, '=', m * (p - m), 'moment2')
    for side in ('A', 'nuA'):
        ee0 = (lambda j: j) if side == 'A' else (lambda j: m - j)
        ee1 = (lambda j: j) if side == 'A' else (lambda j: m - 1 - j)
        for e in range(0, (m - 1) // 2 + 1):
            c0, c1 = {}, {}
            for j in range(m + 1):
                ee = ee0(j)
                if ee <= e:
                    c0[j] = Fr((e + 1 - ee) * (2 * m - (3 * e + ee)), 2)
            for j in range(m):
                ee = ee1(j)
                if ee <= e:
                    c1[j] = Fr((e + 1 - ee) * (2 * m - (3 * e + ee) - 2), 2)
            row(c0, c1, '<=', (e + 1) * (d - e), 'star_%s_%d' % (side, e))
        for t in range(1, m + 1):
            c0, c1 = {}, {}
            for j in range(m + 1):
                ee = ee0(j)
                if m - 1 - ee >= t - 1:
                    c0[j] = comb(m - 1 - ee, t - 1) * (m - ee)
            for j in range(m):
                ee = ee1(j)
                if m - 1 - ee >= t - 1 and m - ee - 1 > 0:
                    c1[j] = comb(m - 1 - ee, t - 1) * (m - ee - 1)
            row(c0, c1, '<=', comb(m, t) * d, 'subsetHP_%s_%d' % (side, t), thm=(t > 2))
    row({0: 2 * m - 2, 1: m - 2}, {}, '<=', 2 * d - 2, 'pencil_A')
    row({m: 2 * m - 2, m - 1: m - 2}, {}, '<=', 2 * d - 2, 'pencil_nuA')
    row({0: m - 1, m: m - 1}, {0: m - 2, m - 1: m - 2}, '<=', d - 1, 'derivHP')
    sq = isqrt(p)
    for k in (2, 3, 4):
        df = 1
        for i in range(1, 2 * k, 2):
            df *= i
        rhs = df * m ** k * p + (2 * k - 1) * m ** (2 * k) * sq
        row({j: val(j, 0) ** (2 * k) for j in range(m + 1)}, {j: val(j, 1) ** (2 * k) for j in range(m)}, '<=', rhs,
            'weil_moment_%d' % (2 * k))
    for j in range(m + 1):
        main = Fr(comb(m, j) * p, 2 ** m)
        err = Fr(comb(m, j) * m * (sq + 1), 2) + 2 * m
        c1 = {j: 1} if j < m else {}
        row({j: 1}, c1, '<=', main + err, 'weil_count_up_%d' % j)
        row({j: -1}, {k: -v for k, v in c1.items()}, '<=', -(main - err), 'weil_count_lo_%d' % j)
    return rows


def solve_lp(nv, rows, extra_bounds=None, method='highs-ipm'):
    """maximise the common normalised slack s of theorem rows; returns (x, s) or (None, None)."""
    from scipy.optimize import linprog
    from scipy.sparse import lil_matrix
    eq = [r for r in rows if r[1] == '=']
    ub = [r for r in rows if r[1] == '<=']
    Aeq = lil_matrix((len(eq), nv + 1))
    beq = np.zeros(len(eq))
    Aub = lil_matrix((len(ub), nv + 1))
    bub = np.zeros(len(ub))
    for i, (co, s, rhs, name, thm) in enumerate(eq):
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)), 1.0])
        for k, c in co.items():
            Aeq[i, k] = float(c) / sc
        beq[i] = float(rhs) / sc
    for i, (co, s, rhs, name, thm) in enumerate(ub):
        sc = max([abs(float(c)) for c in co.values()] + [abs(float(rhs)), 1.0])
        for k, c in co.items():
            Aub[i, k] = float(c) / sc
        if thm:
            Aub[i, nv] = 1.0
        bub[i] = float(rhs) / sc
    c = np.zeros(nv + 1)
    c[nv] = -1
    bounds = [(0, None)] * nv + [(0, 1)]
    if extra_bounds:
        for k, (lo, hi) in extra_bounds.items():
            bounds[k] = (lo, hi)
    res = linprog(c, A_ub=Aub.tocsr(), b_ub=bub, A_eq=Aeq.tocsr() if len(eq) else None,
                  b_eq=beq if len(eq) else None, bounds=bounds, method=method)
    if res.status != 0:
        return None, None, res.status
    return res.x[:nv], res.x[nv], 0


def exact_check_rows(rows, xr, label):
    """exact rational check of every row; returns (all_ok, min theorem slack, n_rows)."""
    ok = True
    worst = None
    for co, sense, rhs, name, thm in rows:
        lhs = sum(Fr(c) * xr[k] for k, c in co.items())
        if sense == '=':
            good = lhs == rhs
        else:
            good = lhs <= rhs
            if thm:
                sl = float(Fr(rhs) - lhs) / max(1.0, abs(float(rhs)))
                worst = sl if worst is None else min(worst, sl)
        check(label, good, (name, float(lhs), float(rhs)))
        ok = ok and good
    return ok, worst, len(rows)


def clique_lp(p, m, K=3, band=0.4, Q=1000):
    """B1: profile of an m-clique A of the Paley graph (difference form = sum form of -A):
    n_j (x not in A), n'_j (x in A; forced n'_0 = m).  Rows: U(m) on (n, n') and on
    (T^k n, n') for k = 1..K (Moebius transfer, Proposition 3.1); band restriction on n.
    Returns an exact rational certificate record or the LP status."""
    import sympy
    d = (p - 1) // 2
    nI = lambda j: j
    nJ = lambda j: m + 1 + j
    base = 2 * m + 1
    Tk = lambda k, j: base + (k - 1) * (m + 1) + j
    nv = base + K * (m + 1)
    rows = []
    rows += U_rows(p, m, [nI(j) for j in range(m + 1)], [nJ(j) for j in range(m)], 'A')
    for k in range(1, K + 1):
        prev = (lambda j: nI(j)) if k == 1 else (lambda j, kk=k - 1: Tk(kk, j))
        for j in range(m + 1):
            co = {Tk(k, j): m}
            co[prev(j)] = co.get(prev(j), 0) - (m - j)
            if 1 <= m + 1 - j <= m:
                co[prev(m + 1 - j)] = co.get(prev(m + 1 - j), 0) - (m + 1 - j)
            rows.append((co, '=', 0, 'link_T%d_%d' % (k, j), False))
        rows += U_rows(p, m, [Tk(k, j) for j in range(m + 1)], [nJ(j) for j in range(m)], 'T%d' % k)
    lo = math.ceil(band * m)
    hi = m + 1 - lo
    bounds = {nJ(0): (m, m)}
    for j in range(1, m):
        bounds[nJ(j)] = (0, 0)
    for j in range(1, m + 1):
        if j < lo or j > hi:
            bounds[nI(j)] = (0, 0)
    x, s, st = solve_lp(nv, rows, bounds)
    rec = {'p': p, 'm': m, 'K': K, 'm(m-1)/p': m * (m - 1) / p, 'defect_g': d - m * (m - 1), 'band': [lo, hi],
           'lp_status': st}
    if x is None:
        return rec
    rec['float_slack'] = float(s)
    xr = [Fr(0)] * nv
    for j in range(m + 1):
        xr[nI(j)] = Fr(int(round(x[nI(j)] * Q)), Q) if x[nI(j)] > 1e-9 else Fr(0)
    xr[nJ(0)] = Fr(m)
    # exact repair of total, moment1, moment2 using the three largest band entries
    val = lambda j: m - 2 * j
    piv = sorted([j for j in range(lo, hi + 1)], key=lambda j: -x[nI(j)])[:3]
    rest = [j for j in range(m + 1) if j not in piv]
    r1 = p - m - sum(xr[nI(j)] for j in rest)
    r2 = -m * (m - 1) - sum(xr[nI(j)] * val(j) for j in rest)
    r3 = m * (p - m) - m * (m - 1) ** 2 - sum(xr[nI(j)] * val(j) ** 2 for j in rest)
    Mx = sympy.Matrix([[1, 1, 1], [val(j) for j in piv], [val(j) ** 2 for j in piv]])
    sol = Mx.LUsolve(sympy.Matrix([sympy.Rational(r.numerator, r.denominator) for r in (r1, r2, r3)]))
    for kk, j in enumerate(piv):
        q = sympy.Rational(sol[kk])
        xr[nI(j)] = Fr(int(q.p), int(q.q))
    prev = [xr[nI(j)] for j in range(m + 1)]
    for k in range(1, K + 1):
        cur = []
        for j in range(m + 1):
            v = (m - j) * prev[j] + ((m + 1 - j) * prev[m + 1 - j] if 1 <= m + 1 - j <= m else 0)
            cur.append(Fr(v, m))
        for j in range(m + 1):
            xr[Tk(k, j)] = cur[j]
        prev = cur
    nonneg = all(v >= 0 for v in xr)
    check('B1.nonneg', nonneg, (p, m))
    ok, worst, nrows = exact_check_rows(rows, xr, 'B1.row')
    # band check (restriction, not theorem)
    bandok = all(xr[nI(j)] == 0 for j in range(1, m + 1) if j < lo or j > hi)
    check('B1.band', bandok, (p, m))
    rec.update({'exact_ok': ok and nonneg and bandok, 'n_rows': nrows, 'min_theorem_slack': worst,
                'support': {('n_%d' % j): float(xr[nI(j)]) for j in range(m + 1) if xr[nI(j)] != 0},
                'T1_support': {('T1_%d' % j): float(xr[Tk(1, j)]) for j in range(m + 1) if xr[Tk(1, j)] != 0} if K else {}})
    return rec


def joint_lp(p, m, n, band=0.4, Q=1000, with_transfer=True):
    """B2: joint profile of a complete biclique A + B in Q u {0} (p = 1 mod 4, r = 0):
    X[j,k] = #{x not in B u -A : u_x = j, v_x = k}, Y[k] = #{b in B : v_b = k}, Z[j] = #{x in -A : u_x = j}
    (u_x = #{a : chi(x+a) = -1}, v_x = #{b : chi(x-b) = -1}); all in the band [band*size, (1-band)*size].
    Rows: U(m) on the A-profile, U(n) on the (-B)-profile and, if with_transfer, U(m) on the averaged
    reciprocal transfer of A (Prop 3.2) and U(n) on that of B; the Chung joint identity is TA:moment1."""
    import sympy
    d = (p - 1) // 2
    loA, hiA = math.ceil(band * m), math.floor((1 - band) * m)
    loB, hiB = math.ceil(band * n), math.floor((1 - band) * n)
    bA = list(range(loA, hiA + 1))
    bB = list(range(loB, hiB + 1))
    idx = {}
    for j in bA:
        for k in bB:
            idx[('X', j, k)] = len(idx)
    for k in bB:
        idx[('Y', k)] = len(idx)
    for j in bA:
        idx[('Z', j)] = len(idx)
    nprim = len(idx)
    for name, size in (('PA0', m + 1), ('PA1', m), ('PB0', n + 1), ('PB1', n), ('TA0', m + 1), ('TA1', m), ('TB0', n + 1), ('TB1', n)):
        for j in range(size):
            idx[(name, j)] = len(idx)
    nv = len(idx)
    V = lambda *key: idx[key]
    links = {}      # aux var -> dict(prim var -> coefficient)

    def add(aux, prim, c):
        links.setdefault(aux, {})
        links[aux][prim] = links[aux].get(prim, 0) + c

    for j in bA:
        for k in bB:
            x = V('X', j, k)
            add(V('PA0', j), x, 1)
            add(V('PB0', k), x, 1)
            add(V('TA0', j), x, Fr(n - k, n))
            add(V('TA0', m - j), x, Fr(k, n))
            add(V('TB0', k), x, Fr(m - j, m))
            add(V('TB0', n - k), x, Fr(j, m))
    for k in bB:
        y = V('Y', k)
        add(V('PA0', 0), y, 1)
        add(V('PB1', k), y, 1)
        add(V('TA0', 0), y, Fr(n - k, n))
        add(V('TA0', m), y, Fr(k, n))
        add(V('TB1', k), y, 1)
    for j in bA:
        z = V('Z', j)
        add(V('PA1', j), z, 1)
        add(V('PB0', 0), z, 1)
        add(V('TA1', j), z, 1)
        add(V('TB0', 0), z, Fr(m - j, m))
        add(V('TB0', n), z, Fr(j, m))
    rows = []
    for key, i in idx.items():
        if key[0] in ('X', 'Y', 'Z'):
            continue
        co = {i: 1}
        for k2, c in links.get(i, {}).items():
            co[k2] = co.get(k2, 0) - c
        rows.append((co, '=', 0, 'link_%s_%d' % key, False))
    rows += U_rows(p, m, [V('PA0', j) for j in range(m + 1)], [V('PA1', j) for j in range(m)], 'PA')
    rows += U_rows(p, n, [V('PB0', k) for k in range(n + 1)], [V('PB1', k) for k in range(n)], 'PB')
    if with_transfer:
        rows += U_rows(p, m, [V('TA0', j) for j in range(m + 1)], [V('TA1', j) for j in range(m)], 'TA')
        rows += U_rows(p, n, [V('TB0', k) for k in range(n + 1)], [V('TB1', k) for k in range(n)], 'TB')
    chung = {}
    for j in bA:
        for k in bB:
            chung[V('X', j, k)] = (m - 2 * j) * (n - 2 * k)
    for k in bB:
        chung[V('Y', k)] = m * (n - 1 - 2 * k)
    for j in bA:
        chung[V('Z', j)] = n * (m - 1 - 2 * j)
    rows.append((chung, '=', -m * n, 'chung', False))
    rows.append(({V('Y', k): 1 for k in bB}, '=', n, 'target_B', False))
    x, s, st = solve_lp(nv, rows)
    rec = {'p': p, 'm': m, 'n': n, 'mn/p': m * n / p, 'defect_g': d - m * n, 'bandA': [loA, hiA], 'bandB': [loB, hiB],
           'with_transfer': with_transfer, 'n_primary_vars': nprim, 'lp_status': st}
    if x is None:
        return rec
    rec['float_slack'] = float(s)
    xr = [Fr(0)] * nv
    for key, i in idx.items():
        if key[0] in ('X', 'Y', 'Z'):
            xr[i] = Fr(int(round(x[i] * Q)), Q) if x[i] > 1e-9 else Fr(0)
    # integer-exact Y and Z sums: fix the largest entries
    for tag, tot, keys in (('Y', n, [('Y', k) for k in bB]), ('Z', m, [('Z', j) for j in bA])):
        big = max(keys, key=lambda kk: x[idx[kk]])
        xr[idx[big]] += tot - sum(xr[idx[kk]] for kk in keys)
    # exact repair: total, PA-m1, PA-m2, PB-m1, PB-m2, chung via six X pivots
    Xkeys = [kk for kk in idx if kk[0] == 'X']
    piv = sorted(Xkeys, key=lambda kk: -x[idx[kk]])[:40]

    def forms(kk):
        j, k = kk[1], kk[2]
        return [1, m - 2 * j, (m - 2 * j) ** 2, n - 2 * k, (n - 2 * k) ** 2, (m - 2 * j) * (n - 2 * k)]

    # choose six pivots with an invertible system
    chosen = []
    for kk in piv:
        trial = chosen + [kk]
        Mx = sympy.Matrix([forms(c) for c in trial]).T
        if Mx.rank() == len(trial):
            chosen = trial
        if len(chosen) == 6:
            break
    rec['repair_pivots'] = len(chosen)
    if len(chosen) < 6:
        rec['exact_ok'] = False
        return rec
    Ys = sum(xr[idx[('Y', k)]] * 0 for k in bB)
    # targets for sum over X of each form
    tgt = [p - n - m,
           -(m * n) - sum(xr[idx[('Z', j)]] * (m - 1 - 2 * j) for j in bA),
           m * (p - m) - n * m * m - sum(xr[idx[('Z', j)]] * (m - 1 - 2 * j) ** 2 for j in bA),
           -(m * n) - sum(xr[idx[('Y', k)]] * (n - 1 - 2 * k) for k in bB),
           n * (p - n) - m * n * n - sum(xr[idx[('Y', k)]] * (n - 1 - 2 * k) ** 2 for k in bB),
           -m * n - sum(xr[idx[('Y', k)]] * m * (n - 1 - 2 * k) for k in bB) - sum(xr[idx[('Z', j)]] * n * (m - 1 - 2 * j) for j in bA)]
    rest = [kk for kk in Xkeys if kk not in chosen]
    rhs = []
    for t in range(6):
        rhs.append(tgt[t] - sum(xr[idx[kk]] * forms(kk)[t] for kk in rest))
    Mx = sympy.Matrix([forms(c) for c in chosen]).T
    sol = Mx.LUsolve(sympy.Matrix([sympy.Rational(r.numerator, r.denominator) for r in rhs]))
    for t, kk in enumerate(chosen):
        q = sympy.Rational(sol[t])
        xr[idx[kk]] = Fr(int(q.p), int(q.q))
    # aux variables exactly from links
    for aux, co in links.items():
        xr[aux] = sum(Fr(c) * xr[k2] for k2, c in co.items())
    nonneg = all(v >= 0 for v in xr)
    check('B2.nonneg', nonneg, (p, m, n))
    ok, worst, nrows = exact_check_rows(rows, xr, 'B2.row')
    rec.update({'exact_ok': ok and nonneg, 'n_rows': nrows, 'min_theorem_slack': worst,
                'Y': {k: float(xr[idx[('Y', k)]]) for k in bB if xr[idx[('Y', k)]]},
                'Z': {j: float(xr[idx[('Z', j)]]) for j in bA if xr[idx[('Z', j)]]},
                'X_support_size': sum(1 for kk in Xkeys if xr[idx[kk]] != 0),
                'X_support': {'%d,%d' % (kk[1], kk[2]): str(xr[idx[kk]]) for kk in Xkeys if xr[idx[kk]] != 0},
                'TA_complete_anticomplete': [float(xr[idx[('TA0', 0)]]), float(xr[idx[('TA0', m)]])]})
    return rec


# ---------------------------------------------------------------- A10: Kalmynin's clique relations (1), (2)

def clique_relations(A, p):
    al = len(A)
    out = []
    for a in A:
        s1 = sum(inv(a - b, p) for b in A if b != a) % p
        s2 = sum(inv(a - b, p) ** 2 for b in A if b != a) % p
        s3 = sum(inv(a - b, p) ** 3 for b in A if b != a) % p
        r1 = (s2 - inv(al, p) * s1 * s1) % p == 0
        r2 = (s3 - inv(al * al, p) * s1 ** 3) % p == 0
        out.append((r1, r2))
    return out


def section_A10(rng):
    """Kalmynin Corollary 2 relations (1), (2) hold for the critical clique at p = 41 and fail for
    near-critical cliques (g > 0); Kalmynin Theorem 5's operator identity D F_a = 0 and the
    factorisation of D(alpha, l) (sympy)."""
    import sympy
    stats = {}
    A41 = [0, 1, 9, 32, 40]
    chi = chi_table(41)
    check('A10.p41_clique', all(chi[(x - y) % 41] == 1 for x in A41 for y in A41 if x != y))
    check('A10.p41_critical', 5 * 4 == 20)
    rel = clique_relations(A41, 41)
    check('A10.p41_relations_hold', all(r1 and r2 for r1, r2 in rel))
    fails = 0
    tested = 0
    for p in [53, 61, 73, 89, 97, 101, 109, 113, 137, 149, 157, 173, 181, 193, 197]:
        chi = chi_table(p)
        d = (p - 1) // 2
        cl = find_cliques(p, chi, rng, count=10)
        cl.sort(key=len, reverse=True)
        for A in cl[:3]:
            g = d - len(A) * (len(A) - 1)
            rel = clique_relations(A, p)
            tested += 1
            if not all(r1 for r1, r2 in rel):
                fails += 1
                witness('A10.relation1_fails_near_critical_clique', {'p': p, 'A': A, 'g': g}, cap=6)
    stats['near_critical_cliques_tested'] = tested
    stats['relation1_fails'] = fails
    # Theorem 5: D F_a = 0 for F_a = (x-a)(1+s(x-a))^alpha, alpha symbolic via T = s(x-a)
    x, a, s, al, l = sympy.symbols('x a s alpha l')
    for alv in range(2, 9):
        Fa = (x - a) * (1 + s * (x - a)) ** alv
        Dg = lambda g: 4 * alv * (alv - 2) * sympy.diff(g, x) * sympy.diff(g, x, 3) - 3 * (alv - 1) * (alv - 2) * sympy.diff(g, x, 2) ** 2 - alv * (alv + 1) * g * sympy.diff(g, x, 4)
        check('A10.thm5_DFa_zero', sympy.expand(Dg(Fa)) == 0, alv)
    ff_ = lambda y, k: sympy.prod([y - i for i in range(k)]) if k > 0 else 1
    Dal = 4 * al * (al - 2) * (al * ff_(l, 3) + l * ff_(al, 3)) - 6 * (al - 1) * (al - 2) * ff_(al, 2) * ff_(l, 2) - al * (al + 1) * (ff_(al, 4) + ff_(l, 4))
    fac = -al * (al - l + 1) * (al - l) * (al - l - 1) * (al ** 2 - (l + 5) * al - l + 6)
    check('A10.thm5_D_alpha_l_factorisation', sympy.expand(Dal - fac) == 0)
    return stats


# ---------------------------------------------------------------- main

def critical_r0_followup(found_all):
    """for Hanson-Petridis-tight pairs with r = 0 (from A9 and the m = 2 family): Yip-Yoo (2.9) holds,
    the sums are not distinct, the moment identity fails; record."""
    out = []
    for rec in found_all:
        p, A, B = rec['p'], rec['A'], rec['B']
        if rec['r'] != 0:
            continue
        res = yy_defect(A, B, p)
        check('A9.r0_tight_yy29_holds', res is not None and res['g'] == 0 and res['Q_zero'], (p, A, B))
        d = (p - 1) // 2
        Ms = moments(A, B, p, d - 1)
        sums = [(a + b) % p for a in A for b in B]
        out.append({'p': p, 'A': A, 'B': B, 'distinct_sums': len(set(sums)) == len(sums),
                    'vanishing_moments_j': [j + 1 for j, v in enumerate(Ms) if v == 0]})
        check('A9.r0_tight_not_decomposition', len(set(sums)) < len(sums), (p, A, B))
    return out


def main():
    t0 = time.time()
    rng = random.Random(20260929)
    R = {'meta': {'script': 'experiments/stepanov_kalmynin_2026_09_29.py', 'note': 'research/stepanov-kalmynin-2026-09-29.md',
                  'python': sys.version.split()[0]}}
    timings = {}

    def run(name, fn, *a):
        t = time.time()
        R[name] = fn(*a)
        timings[name] = round(time.time() - t, 2)
        print(name, timings[name], 's', flush=True)

    run('A1_taylor', section_A1, rng)
    run('A2_relations_with_defect', section_A2, rng)
    run('A2_mu4_exact_decompositions', section_A2_mu4)
    run('A3_yip_yoo_with_defect', section_A3, rng)
    run('A4_moments', section_A4, rng)
    run('A4b_critical_pairs_exhaustive', section_A4b)
    run('A5_reach_rank', section_A5, rng)
    run('A6_moebius_cliques', section_A6, rng)
    run('A7_reciprocal_transfer', section_A7, rng)
    run('A8_hankel_next_coefficient', section_A8, rng)
    t = time.time()
    summ, found = section_A9(pmax=53, mmax=5, budget=300000)
    R['A9_tight_pairs_summary'] = summ
    R['A9_tight_pairs_min3'] = found[:60]
    R['A9_r0_tight_followup'] = critical_r0_followup(found)
    timings['A9'] = round(time.time() - t, 2)
    print('A9', timings['A9'], flush=True)
    run('A10_clique_relations', section_A10, rng)
    # B1: clique LP with Moebius-transfer rows
    B1 = []
    t = time.time()
    for p in [1009, 10009, 40009, 100049, 1000033]:
        d = (p - 1) // 2
        m = max(mm for mm in range(2, 5000) if mm * (mm - 1) <= d)
        rec = clique_lp(p, m, K=3)
        check('B1.certificate_found', rec.get('exact_ok', False), (p, m))
        if p > 20000:
            rec.pop('support', None)
            rec.pop('T1_support', None)
        B1.append(rec)
        print('B1', p, rec.get('exact_ok'), round(time.time() - t, 1), flush=True)
    # sanity: one more than the Hanson-Petridis clique bound is infeasible
    for p in [1009, 10009]:
        d = (p - 1) // 2
        m = max(mm for mm in range(2, 5000) if mm * (mm - 1) <= d) + 1
        rec = clique_lp(p, m, K=1)
        check('B1.sanity_infeasible_above_HP', rec['lp_status'] == 2, (p, m, rec['lp_status']))
        B1.append({'p': p, 'm': m, 'sanity_above_HP_lp_status': rec['lp_status']})
    R['B1_clique_lp'] = B1
    timings['B1'] = round(time.time() - t, 2)
    # B2: joint biclique LP with reciprocal-transfer rows
    B2 = []
    t = time.time()
    for p, m in [(1009, 22), (10009, 71), (40009, 141), (100049, 224), (1000033, 707)]:
        d = (p - 1) // 2
        n = d // m
        rec = joint_lp(p, m, n)
        check('B2.certificate_found', rec.get('exact_ok', False), (p, m, n))
        B2.append(rec)
        print('B2', p, rec.get('exact_ok'), round(time.time() - t, 1), flush=True)
    R['B2_joint_lp'] = B2
    timings['B2'] = round(time.time() - t, 2)
    R['timings_s'] = timings
    R['total_checks'] = sum(CHECKS.values())
    R['checks'] = CHECKS
    R['failures'] = FAILS[:200]
    R['n_failures'] = len(FAILS)
    R['witnesses'] = WITNESS
    R['elapsed_s'] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f:
        json.dump(R, f, indent=1, default=str)
    print('total checks', R['total_checks'], 'failures', len(FAILS), 'elapsed', R['elapsed_s'])


if __name__ == '__main__':
    main()
