#!/usr/bin/env python3
"""Sigma pass, direction `subgroup` (2026-09-05).

Audit of the pass-7 dyadic subgroup Gauss-period bound
    M(p,n) = max_{b!=0} |sum_{h in mu_n} e(bh/p)| <= (17+log N)^{1/9} N^{8/9}
on the pass-5 sixth-energy prime class, together with exact energy scans over
ALL primes p = 1 (mod n) in stated windows (not only class primes), the
classical moment identities, the Karatsuba--Konyagin inequality
    M <= (p E_k E_l)^{1/(2kl)} n^{1-1/k-1/l}      (all integers k,l >= 1),
an elementary level-descent bound, the cyclotomic norm estimate (N') at
orders 4 and 6, an exact-rational exponent ledger, and the ratio
M / sqrt(n log(p/n)) of the primary target.

Arithmetic conventions.  All additive energies E_k(H) = #{x_1+..+x_k =
y_1+..+y_k} are exact integers (numpy int64 sort/segment sums; every count is
bounded by n^{2k-1} <= 2^{42}).  Gauss periods eta are evaluated in IEEE
double precision from a primitive-root table; their absolute error is far
below 1e-9 at every prime used here, and every inequality involving a period
is tested with an explicit tolerance TOL = 1e-9 relative to its right-hand
side and recorded with its margin.  Exact rational arithmetic is used for the
norm-estimate certificates (Bareiss determinants over Z) and for the ledger.

Standard library + numpy only.  Writes results/sigma_subgroup_2026_09_05.json.
Set SIGMA_QUICK=1 for a short smoke run (same code path, fewer primes).
"""
import json
import math
import os
import random
import sys
import time
from fractions import Fraction
from hashlib import sha256
from multiprocessing import Pool
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / 'results' / 'sigma_subgroup_2026_09_05.json'
WORKERS = max(1, min(12, os.cpu_count() or 1))
TOL = 1e-9
T0 = time.time()
QUICK = os.environ.get('SIGMA_QUICK') == '1'


# ----------------------------------------------------------------------------
# elementary number theory helpers
# ----------------------------------------------------------------------------
def sieve(limit):
    s = np.ones(limit + 1, dtype=bool)
    s[:2] = False
    for i in range(2, int(limit ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.flatnonzero(s)


def distinct_prime_factors(x):
    out = []
    d = 2
    while d * d <= x:
        if x % d == 0:
            out.append(d)
            while x % d == 0:
                x //= d
        d += 1
    if x > 1:
        out.append(x)
    return out


def element_of_order(p, n):
    """An element of exact order n (n a power of two dividing p-1)."""
    e = (p - 1) // n
    for x in range(2, p):
        h = pow(x, e, p)
        if n == 1 or pow(h, n // 2, p) != 1:
            return h
    raise ValueError


def primitive_root(p):
    fs = distinct_prime_factors(p - 1)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs):
            return g
    raise ValueError


def subgroup(p, n):
    h = element_of_order(p, n)
    H = sorted(pow(h, j, p) for j in range(n))
    assert len(set(H)) == n
    return H


def T6(t):
    """Complex intrinsic six-term zero-sum count for dyadic t (pass 5)."""
    return 15 * t ** 3 - 45 * t ** 2 + 40 * t


def T4(t):
    return 3 * t ** 2 - 3 * t


# ----------------------------------------------------------------------------
# exact energies
# ----------------------------------------------------------------------------
def _segments(flat, weights):
    order = np.argsort(flat, kind='stable')
    fs = flat[order]
    ws = weights[order]
    starts = np.flatnonzero(np.r_[True, fs[1:] != fs[:-1]])
    return fs[starts], np.add.reduceat(ws, starts)


def energies23(p, H):
    """Exact E_2, E_3 and the supports of r_2, r_3 for H subset F_p."""
    Ha = np.array(H, dtype=np.int64)
    n = len(H)
    s2 = (Ha[:, None] + Ha[None, :]) % p
    v2, c2 = np.unique(s2.ravel(), return_counts=True)
    c2 = c2.astype(np.int64)
    E2 = int((c2 * c2).sum())
    s3 = (v2[:, None] + Ha[None, :]) % p
    v3, r3 = _segments(s3.ravel(), np.repeat(c2, n))
    E3 = int((r3 * r3).sum())
    return E2, E3, (v2, c2), (v3, r3)


def energy4(p, v2, c2):
    """Exact E_4 from r_4 = r_2 * r_2 on the support of r_2."""
    s4 = (v2[:, None] + v2[None, :]) % p
    w = c2[:, None] * c2[None, :]
    v4, r4 = _segments(s4.ravel(), w.ravel())
    return int((r4 * r4).sum())


# ----------------------------------------------------------------------------
# Gauss periods (double precision) via a primitive-root table
# ----------------------------------------------------------------------------
def power_table(g, p, length):
    """g^j mod p for 0 <= j < length, as int64 (block product, O(length))."""
    B = 1 << 12
    nb = (length + B - 1) // B
    base = np.empty(B, dtype=np.int64)
    x = 1
    for i in range(B):
        base[i] = x
        x = x * g % p
    gB = pow(g, B, p)
    big = np.empty(nb, dtype=np.int64)
    y = 1
    for i in range(nb):
        big[i] = y
        y = y * gB % p
    return ((big[:, None] * base[None, :]) % p).ravel()[:length]


def periods(p, n, g):
    """eta on the m = (p-1)/n cosets g^j H, j = 0..m-1, with H = <g^m>."""
    m = (p - 1) // n
    base_m = power_table(g, p, m)
    gm = pow(g, m, p)
    acc_c = np.zeros(m)
    acc_s = np.zeros(m)
    mult = 1
    scale = 2 * math.pi / p
    for _ in range(n):
        row = (mult * base_m) % p
        ang = row * scale
        acc_c += np.cos(ang)
        acc_s += np.sin(ang)
        mult = mult * gm % p
    return acc_c + 1j * acc_s, base_m


def coset_indices(vals, p, g, n, m, base_m):
    """Coset index j (x in g^j H) for each nonzero x in vals (vectorised)."""
    order = np.argsort(base_m)
    sorted_base = base_m[order]
    vals = np.asarray(vals, dtype=np.int64)
    res = np.full(len(vals), -1, dtype=np.int64)
    gm_inv = pow(pow(g, m, p), p - 2, p)
    mult = 1
    for _ in range(n):
        w = (vals * mult) % p
        idx = np.searchsorted(sorted_base, w)
        idx[idx == len(sorted_base)] = 0
        hit = (sorted_base[idx] == w) & (res < 0)
        res[hit] = order[idx[hit]]
        mult = mult * gm_inv % p
    assert (res >= 0).all()
    return res


# ----------------------------------------------------------------------------
# per-prime worker
# ----------------------------------------------------------------------------
KL_PAIRS = [(1, 1), (1, 2), (2, 2), (1, 3), (2, 3), (3, 3), (1, 4), (2, 4), (3, 4), (4, 4)]
RS_PAIRS = [(2, 1), (2, 2), (2, 3), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2)]


def scan_prime(args):
    p, N, want_periods, want_E4 = args
    rec = {'p': p, 'N': N}
    # energies at every dyadic level t | N (class condition needs all levels)
    levels = {}
    t = 2
    while t <= N:
        H_t = subgroup(p, t)
        E2, E3, (v2, c2), (v3, r3) = energies23(p, H_t)
        levels[t] = (E2, E3)
        if t == N:
            HN, v2N, c2N, v3N, r3N = H_t, v2, c2, v3, r3
        t *= 2
    E2, E3 = levels[N]
    rec['E2'], rec['E3'] = E2, E3
    rec['U_neg'] = [t for t, (a, b) in levels.items() if b < T6(t) or a < T4(t)]
    rec['class_sum'] = sum(Fraction(levels[t][1] - T6(t), t ** 3) for t in levels)
    rec['in_class'] = rec['class_sum'] <= math.log(N)
    rec['E3_bound_ok'] = E3 <= (15 + math.log(N)) * N ** 3
    rec['levels_E3'] = {t: levels[t][1] for t in levels}
    # elementary level-descent bound  min_j 2^{j/3} (E3(H_t)/t^3)^{1/9} N^{8/9}
    ld_terms = {int(math.log2(N // t)): 2 ** (math.log2(N // t) / 3) * (levels[t][1] / t ** 3) ** (1 / 9) * N ** (8 / 9)
                for t in levels}
    ld = min(ld_terms.values())
    rec['level_descent_bound'] = ld
    rec['level_descent_argmin_j'] = min(ld_terms, key=ld_terms.get)
    E4 = None
    if want_E4:
        E4 = energy4(p, v2N, c2N)
        rec['E4'] = E4
    if not want_periods:
        return rec
    n = N
    m = (p - 1) // n
    g = primitive_root(p)
    eta, base_m = periods(p, n, g)
    absmag = np.abs(eta)
    jstar = int(np.argmax(absmag))
    M = float(absmag[jstar])
    rec['M'] = M
    rec['ratio'] = M / math.sqrt(n * math.log(p / n))
    rec['g'] = g
    # moment identities: (n^{2k} + n sum_c |eta_c|^{2k}) / p  vs  E_k (exact)
    Ek = {1: n, 2: E2, 3: E3}
    if E4 is not None:
        Ek[4] = E4
    rec['moment_relerr'] = {}
    rec['moment_bound_ok'] = {}
    for k, E in Ek.items():
        spec = (n ** (2 * k) + n * float((absmag ** (2 * k)).sum())) / p
        rec['moment_relerr'][k] = abs(spec - E) / E
        # M^{2k} <= (p E_k - n^{2k}) / n
        rhs = (p * E - n ** (2 * k)) / n
        rec['moment_bound_ok'][k] = M ** (2 * k) <= rhs * (1 + TOL)
    # constancy on cosets: three random nonzero b evaluated directly
    rng = random.Random(p)
    bs = [rng.randrange(1, p) for _ in range(3)]
    js = coset_indices(bs, p, g, n, m, base_m)
    Hn = np.array(HN, dtype=np.int64)
    err = 0.0
    for b, j in zip(bs, js):
        direct = np.exp(2j * math.pi * ((b * Hn) % p) / p).sum()
        err = max(err, abs(direct - eta[j]))
    rec['coset_constancy_err'] = err
    # Karatsuba--Konyagin:  M <= (p E_k E_l)^{1/(2kl)} n^{1-1/k-1/l}
    kon = {}
    for k, l in KL_PAIRS:
        if k in Ek and l in Ek:
            rhs = (p * Ek[k] * Ek[l]) ** (1 / (2 * k * l)) * n ** (1 - 1 / k - 1 / l)
            kon[f'{k},{l}'] = (rhs - M) / rhs
    rec['konyagin_margin'] = kon
    # pass-7 gate (7):  Delta^{rs} <= 2/n + sqrt(p E_r E_s)/n^{r+s}
    gate = {}
    Delta = M / n
    for r, s in RS_PAIRS:
        if r in Ek and s in Ek:
            rhs = 2 / n + math.sqrt(p * Ek[r] * Ek[s]) / n ** (r + s)
            gate[f'{r},{s}'] = (rhs - Delta ** (r * s)) / rhs
    rec['gate_margin'] = gate
    # pass-7 (13):  M <= (B+2)^{1/9} n^{8/9}, B = E3/n^3, valid when p <= n^4
    B = E3 / n ** 3
    rec['pass7_conditional_ok'] = (p > n ** 4) or M <= (B + 2) ** (1 / 9) * n ** (8 / 9) * (1 + TOL)
    rec['pass7_class_bound_ok'] = M <= (17 + math.log(N)) ** (1 / 9) * N ** (8 / 9) * (1 + TOL)
    rec['konyagin33_ok'] = (p > n ** 4) or M <= B ** (1 / 9) * n ** (8 / 9) * (1 + TOL)
    rec['level_descent_ok'] = M <= ld * (1 + TOL)
    # Jensen step at the maximising coset: T = sum_y mu_3(y)|theta(ay)|^3 >= Delta^9
    u3 = float(r3N[v3N == 0].sum()) / n ** 3 if (v3N == 0).any() else 0.0
    a3 = (1 + n * float((absmag ** 3).sum() / n ** 3)) / p
    rec['a3_le_1_over_n'] = a3 <= 1 / n + TOL
    rec['u3_le_1_over_n'] = u3 <= 1 / n + TOL
    mask = v3N != 0
    cj = coset_indices(v3N[mask], p, g, n, m, base_m)
    theta_abs = absmag / n
    T = u3 + float((r3N[mask] * theta_abs[(jstar + cj) % m] ** 3).sum()) / n ** 3
    rec['jensen_margin'] = (T - Delta ** 9) / T
    # centred (SG) statistic at depth r = ceil(log m):  K_min = (Q_r/m)^{1/r} / (r n)
    r = max(1, math.ceil(math.log(m)))
    Qr = n * float((absmag ** (2 * r)).sum())
    rec['SG_depth'] = r
    rec['SG_Kmin'] = (Qr / m) ** (1 / r) / (r * n)
    # smallest depth at which the centred moment bound  M <= (Q_k/1)^{1/(2k)}  beats n
    kmin = None
    for k in range(1, 40):
        Qk = n * float((absmag ** (2 * k)).sum())
        if Qk ** (1 / (2 * k)) < n:
            kmin = k
            break
    rec['moment_depth_needed'] = kmin
    return rec


# ----------------------------------------------------------------------------
# cyclotomic norm estimate (N') at orders 4 and 6, exact
# ----------------------------------------------------------------------------
def bareiss_det(M):
    """Exact integer determinant (fraction-free elimination)."""
    A = [row[:] for row in M]
    n = len(A)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if A[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if A[i][k] != 0), None)
            if sw is None:
                return 0
            A[k], A[sw] = A[sw], A[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                A[i][j] = (A[i][j] * A[k][k] - A[i][k] * A[k][j]) // prev
        prev = A[k][k]
    return sign * A[n - 1][n - 1]


def norm_estimate_check(n, k):
    """Exact check of (N') for mu_n, k-term zero sums, over ALL primes p = 1 (n).

    Returns the list of primes with extra relations, the sum of U(p) log p and
    the (N') right-hand side, and the intrinsic count T_k(n)."""
    N = n // 2
    primes_seen = set()
    intrinsic = 0
    total = 0
    max_det = 1
    for idx in range(n ** (k - 1)):
        a = []
        x = idx
        for _ in range(k - 1):
            a.append(x % n)
            x //= n
        vec = [0] * N
        vec[0] += 1
        for e in a:
            q, rr = divmod(e, N)
            vec[rr] += (-1) ** q
        total += 1
        if not any(vec):
            intrinsic += 1
            continue
        # multiplication-by-f matrix on Z[X]/(X^N+1): column j = X^j * f
        cols = []
        for j in range(N):
            col = [0] * N
            for i, c in enumerate(vec):
                q, rr = divmod(i + j, N)
                col[rr] += c * (-1) ** q
            cols.append(col)
        Mx = [[cols[j][i] for j in range(N)] for i in range(N)]
        D = abs(bareiss_det(Mx))
        assert D != 0
        max_det = max(max_det, D)
        for q in distinct_prime_factors(D):
            if q % n == 1:
                primes_seen.add(q)
    Tk = intrinsic * n
    assert Tk == (T4(n) if k == 4 else T6(n))
    lhs = 0.0
    per_prime = {}
    for p in sorted(primes_seen):
        H = subgroup(p, n)
        E2, E3, _, _ = energies23(p, H)
        Z = E2 if k == 4 else E3
        U = Z - Tk
        assert U >= 0
        per_prime[p] = U
        lhs += U * math.log(p)
    # every prime with an extra relation divides some D_a (so it was seen):
    # check against a direct scan of all primes p = 1 mod n up to max_det
    scan_extra = []
    for p in sieve(int(max_det)):
        p = int(p)
        if p % n == 1:
            H = subgroup(p, n)
            E2, E3, _, _ = energies23(p, H)
            Z = E2 if k == 4 else E3
            if Z > Tk:
                scan_extra.append(p)
    assert scan_extra == sorted(q for q in per_prime if per_prime[q] > 0)
    y = n ** k - Tk
    rhs = y / 2 * math.log(k * n ** k / y)
    return {'n': n, 'k': k, 'T_k': Tk, 'normalized_tuples': total, 'intrinsic_normalized': intrinsic,
            'primes_with_extra_relations': {str(q): per_prime[q] for q in sorted(per_prime) if per_prime[q] > 0},
            'sum_U_log_p': lhs, 'norm_rhs': rhs, 'crude_rhs_half_nk_log_k': n ** k / 2 * math.log(k),
            'holds': lhs <= rhs <= n ** k / 2 * math.log(k)}


# ----------------------------------------------------------------------------
# exact-rational exponent ledger (Tasks 2 and 3 of the sigma brief)
# ----------------------------------------------------------------------------
def kon_exponent(k, l, ek, el, pe=Fraction(4)):
    """log_n of the Karatsuba--Konyagin bound when p = n^pe and E_j = n^{e_j}."""
    k, l, ek, el = Fraction(k), Fraction(l), Fraction(ek), Fraction(el)
    return (pe + ek + el) / (2 * k * l) + 1 - 1 / k - 1 / l


def ledger():
    F = Fraction
    out = {}
    # (a) moment bound at depth k with p <= n^4 and E_k <= n^{2k-3-2kc} gives M <= n^{1-c}
    for k in range(1, 13):
        for c in [F(1, 2), F(1, 8), F(1, 9), F(1, 100)]:
            assert (F(4) + (2 * k - 3 - 2 * k * c) - 1) / (2 * k) == 1 - c
    # what that requires against the diagonal k! n^k and the principal term n^{2k}/p = n^{2k-4}
    need_vs_diag = {k: F(k - 3, 2 * k) for k in range(1, 13)}          # c <= (k-3)/(2k)
    need_vs_principal = {k: F(1, 2 * k) for k in range(1, 13)}        # c <= 1/(2k)
    best_moment = max(min(need_vs_diag[k], need_vs_principal[k]) for k in range(1, 13))
    assert best_moment == F(1, 8) and min(need_vs_diag[4], need_vs_principal[4]) == F(1, 8)
    out['moment_route'] = {'E_k_needed_for_saving_c': 'E_k <= n^{2k-3-2kc} (p <= n^4)',
                           'c_max_vs_diagonal_k_factorial_n^k': {k: str(v) for k, v in need_vs_diag.items()},
                           'c_max_vs_principal_n^{2k-4}': {k: str(v) for k, v in need_vs_principal.items()},
                           'best_fixed_depth_saving': str(best_moment), 'attained_at_k': 4}
    # (b) Heath-Brown--Konyagin E_2 << n^{5/2} at p = n^4 is trivial in the moment bound and in (2,2)
    assert kon_exponent(2, 1, F(5, 2), 1) == F(11, 8) > 1
    assert kon_exponent(2, 2, F(5, 2), F(5, 2)) == F(9, 8) > 1
    # nontrivial range of the k=2 moment bound with E_2 <= n^{5/2}:  p < n^{5/2}
    assert (F(5, 2) + F(5, 2) - 1) / 4 == 1                                  # p = n^{5/2} boundary
    # (c) Konyagin k=l=3 with E_3 <= B n^3 and p <= n^4:  exponent 8/9 (equals pass 7 (13))
    assert kon_exponent(3, 3, 3, 3) == F(8, 9)
    # with E_3 <= n^4 (log):  boundary exactly at p = n^4  (Konyagin 2002 threshold n > p^{1/4})
    assert kon_exponent(3, 3, 4, 4) == F(1)
    assert kon_exponent(3, 3, 4, 4, pe=F(4) - F(1, 100)) < 1
    # (d) pass-7 (16) and Konyagin (4,4) / (1,4) / (2,4) with E_4 <= D n^4:  exponent 7/8
    assert kon_exponent(4, 4, 4, 4) == F(7, 8) == kon_exponent(1, 4, 1, 4) == kon_exponent(2, 4, 2, 4)
    # (e) optimum of Konyagin's inequality over all k,l with ideal energies e_j = max(j, 2j-4)
    def ideal(j):
        return max(F(j), F(2 * j - 4))
    best = F(0)
    argbest = []
    for k in range(1, 41):
        for l in range(1, 41):
            s = 1 - kon_exponent(k, l, ideal(k), ideal(l))
            if s > best:
                best, argbest = s, [(k, l)]
            elif s == best:
                argbest.append((k, l))
    assert best == F(1, 8)
    # closed forms covering all k,l (proof in the note): k,l>=4: 2/(kl); k<=3<l... etc.
    for k in range(4, 41):
        for l in range(4, 41):
            assert 1 - kon_exponent(k, l, ideal(k), ideal(l)) == F(2, k * l)
    for k in range(1, 4):
        for l in range(4, 41):
            assert 1 - kon_exponent(k, l, ideal(k), ideal(l)) == F(1, 2 * l)
    for k in range(1, 4):
        for l in range(1, 4):
            assert 1 - kon_exponent(k, l, ideal(k), ideal(l)) == F(k + l - 4, 2 * k * l)
    out['konyagin_ideal_optimum'] = {'best_saving': str(best), 'argmax_pairs': argbest,
                                     'closed_forms': 'k,l>=4: 2/(kl); min(k,l)<=3<max: 1/(2 max); k,l<=3: (k+l-4)/(2kl)'}
    # (f) Di Benedetto et al. Theorem 3.1 at p <= n^4:  2689/2880 + 4/72 = 2849/2880, saving 31/2880
    assert F(2689, 2880) + F(4, 72) == F(2849, 2880)
    assert 1 - F(2849, 2880) == F(31, 2880)
    # the paper's own 'in particular' form H^{1-31/2880} for H > p^{1/4} is the same number
    assert F(1) - F(31, 2880) == F(2849, 2880)
    # (g) level-descent exponent bookkeeping: 2^j (2^{4j} t^{10})^{1/18} t^{1/3} = 2^{j/3} N^{8/9}, t = N/2^j
    for j in range(0, 7):
        assert F(j) + F(4 * j, 18) - F(8 * j, 9) == F(j, 3)
    # (h) trivial energy chain E_3 <= n^2 E_2 <= n^2 * n^{5/2} in (3,3): exponent > 1
    assert kon_exponent(3, 3, F(9, 2), F(9, 2)) == (F(4) + 9) / 18 + F(1, 3) > 1
    # (i) recurrence-only descent: E_2(H_{2k}) <= 16 E_2(H_k) => E_2(H_N) <= 6 N^4 / 16 ... trivial (>= n^3)
    #     (with B <= E_2(K), T <= E_2(K): 2+6+8 = 16)
    assert 2 + 6 + 8 == 16
    out['di_benedetto_at_quartic'] = {'exponent': '2849/2880', 'saving': '31/2880', 'decimal': float(F(31, 2880))}
    # (j) cited unconditional energy inputs at n ~ p^{1/4} (log factors dropped) are all trivial in (2.1) at p = n^4:
    #     HBK e_2 = 5/2; MRSS Thm 3 e_2 = 49/20; MRSS Cor 7 e_3 = 4; Konyagin 2002 / Shkredov Thm 6 e_d = 2d-2+2^{1-d}
    kon_d = lambda d: F(2 * d - 2) + F(1, 2 ** (d - 1))
    assert kon_d(3) == F(17, 4) and kon_d(4) == F(49, 8)
    table = {'(2,1) HBK': kon_exponent(2, 1, F(5, 2), 1), '(2,2) HBK': kon_exponent(2, 2, F(5, 2), F(5, 2)),
             '(2,2) MRSS': kon_exponent(2, 2, F(49, 20), F(49, 20)), '(2,3) MRSS': kon_exponent(2, 3, F(49, 20), 4),
             '(3,3) MRSS': kon_exponent(3, 3, 4, 4), '(3,4) MRSS+Sh14': kon_exponent(3, 4, 4, kon_d(4)),
             '(4,4) Sh14': kon_exponent(4, 4, kon_d(4), kon_d(4)), '(1,3) MRSS': kon_exponent(1, 3, 1, 4)}
    assert table['(2,1) HBK'] == F(11, 8) and table['(2,2) MRSS'] == F(89, 80) and table['(2,3) MRSS'] == F(83, 80)
    assert table['(3,3) MRSS'] == 1 and table['(3,4) MRSS+Sh14'] == F(193, 192) and table['(4,4) Sh14'] == F(129, 128)
    assert all(v >= 1 for v in table.values())
    # exhaustive: with e_2 = 49/20, e_3 = 4, e_d = kon_d(d) (d >= 4), e_1 = 1, every pair k,l <= 40 is trivial
    def e_cited(j):
        return {1: F(1), 2: F(49, 20), 3: F(4)}.get(j, kon_d(j))
    assert all(kon_exponent(k, l, e_cited(k), e_cited(l)) >= 1 for k in range(1, 41) for l in range(1, 41))
    out['konyagin_with_cited_energies'] = {k: str(v) for k, v in table.items()}
    # (k) Shkredov arXiv:1705.09703 Cor. 16 at |Gamma| = p^{1/4}: k = ceil(2 log p/log|Gamma| + 4) = 12,
    #     proof line saving 1/2^{k+2} = 1/2^14; stated form p^{-delta/2^{7+2/delta}} at delta = 1/4 gives |Gamma|^{-1/2^15}
    k16 = math.ceil(2 * 4 + 4)
    assert k16 == 12
    assert F(1, 2 ** (k16 + 2)) == F(1, 16384)
    assert F(4) * F(1, 4) / F(2 ** (7 + 8)) == F(1, 32768)
    assert F(1, 16384) < F(31, 2880) and F(1, 32768) < F(31, 2880)
    assert F(175, 9437184) < F(1, 16384) < F(31, 2880)          # BG09 < Sh19 (proof line, delta = 1/4) < DGGGST
    out['shkredov_cor16_at_quarter'] = {'k': k16, 'saving_proof_line': '1/16384', 'saving_stated_form': '1/32768',
                                        'bourgain_garaev_2009': '175/9437184', 'di_benedetto_2020': '31/2880'}
    # (l) recurrence-only descent E_2(H_{2k}) <= 16 E_2(H_k):  E_2(H_N) <= 6*16^{log2 N - 1} = (3/8) N^4 >= N^3 for N >= 4
    for m in range(2, 8):
        Nn = 2 ** m
        assert F(6) * F(16) ** (m - 1) == F(3, 8) * Nn ** 4 and F(3, 8) * Nn ** 4 >= Nn ** 3
    return out


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def prime_list(lo, hi, n):
    ps = sieve(int(hi))
    ps = ps[(ps >= lo) & (ps % n == 1)]
    return [int(p) for p in ps]


def plan():
    """Windows and sampling.  (lo, hi) for energies; periods/E4 sampling rules."""
    jobs = []
    meta = {}
    for N in [2, 4, 8, 16, 32, 64]:
        lo = math.ceil(N ** 3.5)
        hi = math.floor(N ** 4.5)
        if N == 64:
            lo, hi = N ** 4 // 4, N ** 4          # quartic window only
        if N == 32:
            hi = min(hi, 5_931_641)
        ps = prime_list(lo, hi, N)
        ql, qh = N ** 4 // 4, N ** 4
        rule = {2: (1, 1), 4: (1, 1), 8: (1, 1), 16: (1, 1), 32: (1, 6), 64: (12, 12)}[N]
        if QUICK:
            ps = ps[:60] if N >= 32 else ps[:400]
            rule = (1, 1) if N < 64 else (5, 5)
        pin, pout = rule
        must = {6700417} if N == 64 else set()
        sel = []
        for i, p in enumerate(ps):
            inq = ql <= p <= qh
            step = pin if inq else pout
            want_periods = (i % step == 0) or p in must
            if N == 64:
                want_E4 = want_periods and (i % (12 * 8) == 0 or p in must)
            elif N == 32:
                want_E4 = want_periods and (i % 4 == 0)
            else:
                want_E4 = want_periods
            sel.append((p, N, want_periods, want_E4))
        jobs.extend(sel)
        meta[N] = {'energy_window': [lo, hi], 'quartic_window': [ql, qh], 'primes_scanned': len(ps),
                   'periods_computed': sum(1 for s in sel if s[2]), 'E4_computed': sum(1 for s in sel if s[3]),
                   'sampling_rule': f'energies: all primes p=1 mod {N} in [{lo},{hi}]; periods: every {pin}th prime in '
                                    f'the quartic window and every {pout}th outside' + (
                                        '; p=6700417 forced' if N == 64 else '')}
    return jobs, meta


def gaussian_expected_max_abs(m, h=0.002, xmax=12.0):
    """E[max_{i<=m}|Z_i|] for iid standard real Gaussians: int_0^inf (1 - erf(x/sqrt2)^m) dx (trapezoid)."""
    tot = 0.0
    x = 0.0
    prev = 1.0
    while x < xmax:
        x += h
        cur = 1.0 - math.erf(x / math.sqrt(2)) ** m
        tot += 0.5 * (prev + cur) * h
        prev = cur
    return tot


def summarize(N, recs, meta):
    log_N = math.log(N)
    inwin = [r for r in recs if meta['quartic_window'][0] <= r['p'] <= meta['quartic_window'][1]]
    with_M = [r for r in recs if 'M' in r]
    out = dict(meta)
    out['n'] = N
    out['E3_ge_T6_all_levels_violations'] = [r['p'] for r in recs if r['U_neg']]
    out['class_failures'] = [r['p'] for r in inwin if not r['in_class']]
    out['E3_endpoint_bound_failures'] = [r['p'] for r in inwin if not r['E3_bound_ok']]
    out['markov_bound_on_class_failures'] = (8 * N ** 3 - 8) * math.log(6) / (14 * log_N * math.log(N ** 4 / 4))
    out['quartic_primes'] = len(inwin)
    out['max_E2_quartic'] = max(r['E2'] for r in inwin) if inwin else None
    out['max_E3_quartic'] = max(r['E3'] for r in inwin) if inwin else None
    out['max_E2_over_n2_quartic'] = max(r['E2'] for r in inwin) / N ** 2 if inwin else None
    out['max_E3_over_n3_quartic'] = max(r['E3'] for r in inwin) / N ** 3 if inwin else None
    e4 = [r['E4'] for r in inwin if 'E4' in r]
    out['max_E4_over_n4_quartic'] = max(e4) / N ** 4 if e4 else None
    out['E4_count_quartic'] = len(e4)
    out['diag_k_factorial'] = {k: math.factorial(k) for k in (2, 3, 4)}
    out['principal_n^{2k}/p_at_window_ends'] = {k: [N ** (2 * k) / meta['quartic_window'][1], N ** (2 * k) / meta['quartic_window'][0]]
                                               for k in (2, 3, 4)}
    out['max_class_sum_quartic'] = float(max(r['class_sum'] for r in inwin)) if inwin else None
    out['log_N'] = log_N
    out['E2_intrinsic_count_quartic'] = sum(1 for r in inwin if r['E2'] == T4(N))
    out['E3_intrinsic_count_quartic'] = sum(1 for r in inwin if r['E3'] == T6(N))
    # over the whole energy window (all primes scanned), not only the quartic window
    out['max_E3_over_n3_all'] = max(r['E3'] for r in recs) / N ** 3 if recs else None
    out['max_E2_over_n2_all'] = max(r['E2'] for r in recs) / N ** 2 if recs else None
    fails = {}
    for key in ['pass7_conditional_ok', 'pass7_class_bound_ok', 'konyagin33_ok', 'level_descent_ok',
                'a3_le_1_over_n', 'u3_le_1_over_n']:
        fails[key] = [r['p'] for r in with_M if not r[key]]
    fails['moment_bound'] = [(r['p'], k) for r in with_M for k, ok in r['moment_bound_ok'].items() if not ok]
    fails['konyagin_any'] = [(r['p'], kl) for r in with_M for kl, mg in r['konyagin_margin'].items() if mg < -TOL]
    fails['gate_any'] = [(r['p'], rs) for r in with_M for rs, mg in r['gate_margin'].items() if mg < -TOL]
    fails['jensen'] = [r['p'] for r in with_M if r['jensen_margin'] < -TOL]
    out['inequality_failures'] = fails
    # the class bound (17+log N)^{1/9} N^{8/9} is only a claim for class primes in the quartic window;
    # record separately whether it happens to hold at NON-class quartic primes (data, not a theorem)
    out['pass7_class_bound_at_nonclass_quartic'] = [(r['p'], r['M'], r['pass7_class_bound_ok']) for r in with_M
                                                    if not r['in_class'] and meta['quartic_window'][0] <= r['p'] <= meta['quartic_window'][1]]
    out['pass7_class_bound_value'] = (17 + log_N) ** (1 / 9) * N ** (8 / 9)
    if with_M:
        out['max_moment_identity_relerr'] = max(max(r['moment_relerr'].values()) for r in with_M)
        out['max_coset_constancy_err'] = max(r['coset_constancy_err'] for r in with_M)
        rat = sorted(r['ratio'] for r in with_M)
        best = max(with_M, key=lambda r: r['ratio'])
        out['ratio_stats'] = {'count': len(rat), 'max': rat[-1], 'argmax_p': best['p'], 'M_at_argmax': best['M'],
                              'E2_at_argmax': best['E2'], 'E3_at_argmax': best['E3'],
                              'median': rat[len(rat) // 2], 'mean': sum(rat) / len(rat),
                              'p90': rat[int(0.9 * (len(rat) - 1))], 'min': rat[0]}
        ratq = sorted(r['ratio'] for r in with_M if meta['quartic_window'][0] <= r['p'] <= meta['quartic_window'][1])
        out['ratio_max_quartic'] = ratq[-1] if ratq else None
        out['ratio_count_quartic'] = len(ratq)
        out['max_M_over_n'] = max(r['M'] for r in with_M) / N
        out['max_M'] = max(r['M'] for r in with_M)
        out['min_konyagin33_margin'] = min(r['konyagin_margin']['3,3'] for r in with_M)
        out['min_gate33_margin'] = min(r['gate_margin']['3,3'] for r in with_M)
        out['min_jensen_margin'] = min(r['jensen_margin'] for r in with_M)
        out['max_SG_Kmin'] = max(r['SG_Kmin'] for r in with_M)
        out['SG_depth_range'] = [min(r['SG_depth'] for r in with_M), max(r['SG_depth'] for r in with_M)]
        depths = [r['moment_depth_needed'] for r in with_M if r['moment_depth_needed']]
        out['moment_depth_needed_max'] = max(depths) if depths else None
        out['moment_depth_needed_min'] = min(depths) if depths else None
        out['level_descent_exponent_max'] = max(math.log(r['level_descent_bound']) / log_N for r in with_M)
        out['level_descent_exponent_min'] = min(math.log(r['level_descent_bound']) / log_N for r in with_M)
        out['actual_exponent_max'] = max(math.log(r['M']) / log_N for r in with_M)
        # best fixed-(k,l) Konyagin exponent actually achieved by exact energies (max over primes of the bound)
        kb = {}
        for kl in KL_PAIRS:
            key = f'{kl[0]},{kl[1]}'
            vals = [math.log(r['M'] / (1 - r['konyagin_margin'][key])) / log_N for r in with_M if key in r['konyagin_margin']]
            if vals:
                kb[key] = {'max_exponent_of_bound': max(vals), 'count': len(vals)}
        out['konyagin_bound_exponents'] = kb
        top = sorted(with_M, key=lambda r: -r['ratio'])[:10]
        out['top_ratio_table'] = [{'p': r['p'], 'M': r['M'], 'ratio': r['ratio'], 'E2': r['E2'], 'E3': r['E3'],
                                   'E4': r.get('E4'), 'in_class': r['in_class'], 'class_sum': float(r['class_sum'])}
                                  for r in top]
        top = sorted(with_M, key=lambda r: -r['E3'])[:5]
        out['top_E3_table'] = [{'p': r['p'], 'M': r['M'], 'ratio': r['ratio'], 'E2': r['E2'], 'E3': r['E3'],
                                'in_class': r['in_class'], 'class_sum': float(r['class_sum']),
                                'level_descent_bound': r['level_descent_bound']} for r in top]
    topE3 = sorted(recs, key=lambda r: -r['E3'])[:5]
    out['top_E3_all_scanned'] = [{'p': r['p'], 'E2': r['E2'], 'E3': r['E3'], 'in_class': r['in_class'],
                                  'class_sum': float(r['class_sum'])} for r in topE3]
    out['class_failure_details'] = [{'p': r['p'], 'E2': r['E2'], 'E3': r['E3'], 'E3_over_n3': r['E3'] / N ** 3,
                                     'class_sum': float(r['class_sum']), 'M': r.get('M'), 'ratio': r.get('ratio'),
                                     'levels_E3_over_t3': {t: v / t ** 3 for t, v in r['levels_E3'].items()}}
                                    for r in inwin if not r['in_class']]
    out['nonclass_outside_quartic'] = [{'p': r['p'], 'E2': r['E2'], 'E3': r['E3'], 'class_sum': float(r['class_sum'])}
                                       for r in recs if not r['in_class'] and r not in inwin]
    if with_M:
        # HEURISTIC comparison: if the m nonprincipal coset periods were iid real N(0, n), the mean of
        # M/sqrt(n log(p/n)) would be E[max_m |Z|]/sqrt(log(p/n)); evaluated at the median-p prime.
        ps = sorted(r['p'] for r in with_M)
        pm = ps[len(ps) // 2]
        mm = (pm - 1) // N
        gm = gaussian_expected_max_abs(mm)
        obs = [r['ratio'] for r in with_M]
        out['gaussian_model'] = {'median_p': pm, 'm': mm, 'E_max_abs_Z': gm,
                                 'predicted_mean_ratio': gm / math.sqrt(math.log(pm / N)),
                                 'observed_mean_ratio': sum(obs) / len(obs),
                                 'observed_median_ratio': sorted(obs)[len(obs) // 2],
                                 'sqrt2_limit': math.sqrt(2)}
        # HEURISTIC: the max over S sampled primes of the max over m cosets is the max over m*S iid |Z|
        S = len(with_M)
        gS = gaussian_expected_max_abs(mm * S)
        out['gaussian_model']['sample_count'] = S
        out['gaussian_model']['predicted_max_ratio_over_sample'] = gS / math.sqrt(math.log(pm / N))
        out['gaussian_model']['observed_max_ratio'] = max(obs)
    out['level_descent_argmin_j_counts'] = {}
    for r in recs:
        j = r['level_descent_argmin_j']
        out['level_descent_argmin_j_counts'][j] = out['level_descent_argmin_j_counts'].get(j, 0) + 1
    return out


def main():
    jobs, meta = plan()
    print(f'{len(jobs)} prime jobs on {WORKERS} workers (QUICK={QUICK})', flush=True)
    recs = {N: [] for N in meta}
    with Pool(WORKERS) as pool:
        for i, rec in enumerate(pool.imap_unordered(scan_prime, jobs, chunksize=8)):
            recs[rec['N']].append(rec)
            if i % 5000 == 0:
                print(f'  {i} done, {time.time() - T0:.0f}s', flush=True)
    print(f'scan finished at {time.time() - T0:.0f}s', flush=True)
    summary = {N: summarize(N, recs[N], meta[N]) for N in meta}
    for N in meta:
        s = summary[N]
        print(f'n={N}: {s["primes_scanned"]} primes, class failures {len(s["class_failures"])} '
              f'(Markov {s["markov_bound_on_class_failures"]:.1f}), E3 endpoint failures {len(s["E3_endpoint_bound_failures"])}, '
              f'max E3/n^3 {s["max_E3_over_n3_all"]:.3f}, ratio max {s.get("ratio_stats", {}).get("max")}',
              flush=True)
    norm = [norm_estimate_check(4, 4), norm_estimate_check(8, 4), norm_estimate_check(16, 4),
            norm_estimate_check(4, 6), norm_estimate_check(8, 6)]
    assert all(x['holds'] for x in norm)
    print(f'norm estimates done at {time.time() - T0:.0f}s', flush=True)
    led = ledger()
    compact = {}
    for N in meta:
        compact[N] = [[r['p'], r['E2'], r['E3'], r.get('E4'), r.get('M'), r.get('ratio'), int(r['in_class']),
                       float(r['class_sum'])] for r in sorted(recs[N], key=lambda r: r['p'])]
    total_checks = sum(len(recs[N]) for N in meta)
    total_ineq = sum(len(r['konyagin_margin']) + len(r['gate_margin']) + len(r['moment_bound_ok']) + 7
                     for N in meta for r in recs[N] if 'M' in r)
    all_fail = {N: {k: v for k, v in summary[N]['inequality_failures'].items() if v} for N in meta}
    witnesses = {N: summary[N]['E3_endpoint_bound_failures'] for N in meta}
    inputs = ['research/parallel7-multiplier-2026-09-04.md', 'research/parallel7-subgroup-2026-09-04.md',
              'research/parallel5-subgroup-2026-09-04.md', 'research/parallel4-subgroup-2026-09-04.md',
              'research/cyclotomic-prime-average.md', 'research/subgroup-target.md',
              'experiments/sigma_subgroup_2026_09_05.py']
    result = {
        'status': ('All exact and float-tolerance checks passed; no used inequality failed at any scanned prime.'
                   if not any(all_fail[N] for N in meta) else 'FAILURES RECORDED'),
        'quick_mode': QUICK,
        'seconds': time.time() - T0,
        'workers': WORKERS,
        'tolerance': TOL,
        'conventions': __doc__,
        'per_order': summary,
        'norm_estimate_checks': norm,
        'exponent_ledger': led,
        'inequality_failures': all_fail,
        'E3_endpoint_bound_failures_by_order': witnesses,
        'class_failures_by_order': {N: summary[N]['class_failures'] for N in meta},
        'counts': {'primes_scanned': total_checks,
                   'periods_computed': sum(1 for N in meta for r in recs[N] if 'M' in r),
                   'E4_computed': sum(1 for N in meta for r in recs[N] if 'E4' in r),
                   'pointwise_inequality_checks': total_ineq,
                   'norm_estimate_certificates': len(norm)},
        'compact_records': {'columns': ['p', 'E2', 'E3', 'E4_or_null', 'M_or_null', 'ratio_or_null', 'in_class', 'class_sum'],
                            'by_order': compact},
        'input_sha256': {name: sha256((ROOT / name).read_bytes()).hexdigest() for name in inputs},
    }
    OUT.write_text(json.dumps(result, indent=1, default=float) + '\n')
    print(json.dumps({'status': result['status'], 'seconds': result['seconds'], 'counts': result['counts'],
                      'failures': all_fail,
                      'E3_endpoint_failures': witnesses,
                      'class_failures': result['class_failures_by_order']}, indent=1, default=float), flush=True)


if __name__ == '__main__':
    main()
