#!/usr/bin/env python3
"""Stepanov wave, direction `algebra` (2026-09-26): exact verifier.

Every identity and inequality stated in research/stepanov-algebra-2026-09-26.md
is checked here by exact modular / integer arithmetic at small primes.
Writes results/stepanov_algebra_2026_09_26.json.

Sections
  S1a  exact lifting of S(A,B) from one element of F_p (and the 2p-lift)
  S1b  binom(d,j) = (-1/4)^j C(2j,j) mod p, the kernel formula for S,
       (1+u)^d = (1+u)^(-1/2) mod u^p, chi(1+u) = T_d(u)
  S1c  gcd form of HP, N_e as a point-on-many-polynomials count
  S2   symmetric-function form: F = h_d(A+x) - 1, derivative = shift,
       HP lemma, syndrome lemma, order bounds, Hankel (e=1,2) polynomials,
       unweighted Phi_A has no multiplicity, multiset witness,
       sign-pattern HP  mu(u) <= (d+m-1)/m
  S3   list-decoding view: N_e versus HP, second moment, Johnson variants,
       distinct-word Johnson, Hankel bound
  S4   F_{p^2}: the second-moment bound is attained by the subfield
Standard library + numpy only.
"""
import json
import math
import os
import random
import sys
import time

import numpy as np

T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_algebra_2026_09_26.json")
RNG = random.Random(20260926)

RESULTS = {"checks": {}, "witnesses": {}, "tables": {}, "failures": []}


def bump(key, n=1):
    RESULTS["checks"][key] = RESULTS["checks"].get(key, 0) + n


def fail(key, info):
    RESULTS["failures"].append({"check": key, "info": info})


def require(cond, key, info=None):
    bump(key)
    if not cond:
        fail(key, info)
    return cond


# ---------------------------------------------------------------- basics

def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]


def legendre_table(p):
    d = (p - 1) // 2
    chi = [0] * p
    for x in range(1, p):
        chi[x] = 1 if pow(x, d, p) == 1 else -1
    return chi


def lar(x, p):
    """least absolute residue in (-p/2, p/2)"""
    x %= p
    return x - p if x > p // 2 else x


def inv(x, p):
    return pow(x % p, p - 2, p)


def fact_tables(p):
    f = [1] * p
    for i in range(1, p):
        f[i] = f[i - 1] * i % p
    fi = [1] * p
    fi[p - 1] = inv(f[p - 1], p)
    for i in range(p - 1, 0, -1):
        fi[i - 1] = fi[i] * i % p
    return f, fi


def binom_mod(n, k, f, fi, p):
    if k < 0 or k > n:
        return 0
    return f[n] * fi[k] % p * fi[n - k] % p


# polynomials: list of coefficients mod p, index = degree
def ptrim(a):
    while a and a[-1] == 0:
        a.pop()
    return a


def padd(a, b, p):
    n = max(len(a), len(b))
    return ptrim([((a[i] if i < len(a) else 0) + (b[i] if i < len(b) else 0)) % p for i in range(n)])


def pscale(a, c, p):
    return ptrim([x * c % p for x in a])


def pmul(a, b, p):
    if not a or not b:
        return []
    r = [0] * (len(a) + len(b) - 1)
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                r[i + j] = (r[i + j] + x * y) % p
    return ptrim(r)


def pdivmod(a, b, p):
    a = a[:]
    q = [0] * max(len(a) - len(b) + 1, 1)
    ib = inv(b[-1], p)
    while len(a) >= len(b) and a:
        c = a[-1] * ib % p
        s = len(a) - len(b)
        q[s] = c
        for i, y in enumerate(b):
            a[s + i] = (a[s + i] - c * y) % p
        ptrim(a)
    return ptrim(q), a


def pgcd(a, b, p):
    a, b = ptrim(a[:]), ptrim(b[:])
    while b:
        _, r = pdivmod(a, b, p)
        a, b = b, r
    if a:
        a = pscale(a, inv(a[-1], p), p)
    return a


def pderiv(a, p):
    return ptrim([(i * a[i]) % p for i in range(1, len(a))])


def peval(a, x, p):
    r = 0
    for c in reversed(a):
        r = (r * x + c) % p
    return r


def taylor_shift(a, b, p):
    """coefficients of a(x+b)"""
    r = []
    for c in reversed(a):
        # r <- r*(x+b) + c
        nr = [0] * (len(r) + 1)
        for i, y in enumerate(r):
            nr[i + 1] = (nr[i + 1] + y) % p
            nr[i] = (nr[i] + y * b) % p
        nr[0] = (nr[0] + c) % p
        r = nr
    return ptrim(r)


def order_at(a, b, p):
    """multiplicity of the root b of the nonzero polynomial a (a has degree < p)"""
    s = taylor_shift(a, b, p)
    for i, c in enumerate(s):
        if c:
            return i
    return None


def translate_power(a, N, f, fi, p):
    """(x+a)^N, N < p"""
    return ptrim([binom_mod(N, i, f, fi, p) * pow(a, N - i, p) % p for i in range(N + 1)])


def hp_weights(A, p):
    c = []
    for k, ak in enumerate(A):
        den = 1
        for l, al in enumerate(A):
            if l != k:
                den = den * (ak - al) % p
        c.append(inv(den, p))
    return c


def hp_poly(A, p, f, fi, extra=0):
    """F(x) = -1 + sum_k c_k (x+a_k)^(d+M-1+extra)"""
    d = (p - 1) // 2
    M = len(A)
    D = d + M - 1 + extra
    c = hp_weights(A, p)
    F = [p - 1]
    for ck, ak in zip(c, A):
        F = padd(F, pscale(translate_power(ak, D, f, fi, p), ck, p), p)
    return F, c, D


def complete_h(ys, imax, p):
    """h_0..h_imax of ys mod p via e_j recurrence"""
    M = len(ys)
    e = [1]
    for y in ys:  # prod (1 + y t)
        ne = e + [0]
        for j in range(len(e), 0, -1):
            ne[j] = (ne[j] + e[j - 1] * y) % p
        e = ne
    h = [1]
    for i in range(1, imax + 1):
        s = 0
        for j in range(1, min(i, M) + 1):
            term = e[j] * h[i - j]
            s = s + term if j % 2 == 1 else s - term
        h.append(s % p)
    return h


def bad_counts(A, p, chi):
    """e_b and zero indicator for all b"""
    eb = [0] * p
    z = [0] * p
    for b in range(p):
        s = 0
        for a in A:
            v = chi[(a + b) % p]
            if v == -1:
                s += 1
            elif v == 0:
                z[b] = 1
        eb[b] = s
    return eb, z


def save():
    RESULTS["elapsed_seconds"] = round(time.time() - T0, 2)
    RESULTS["total_checks"] = sum(RESULTS["checks"].values())
    RESULTS["n_failures"] = len(RESULTS["failures"])
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(RESULTS, fh, indent=1, sort_keys=True)


def rand_set(p, m, exclude=()):
    pool = [x for x in range(p) if x not in exclude]
    return sorted(RNG.sample(pool, m))


# ---------------------------------------------------------------- S1a
def section_1a():
    """S(A,B) = lar(sum (a+b)^d mod p) when mn < p/2;
    S = unique integer in [-mn,mn] with S = R (mod p), S = mn - r (mod 2) when mn < p."""
    key1, key2 = "S1a_lift_half", "S1a_lift_2p"
    wit_fail_half = None
    for p in primes_upto(400)[3:]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        for _ in range(12):
            m = RNG.randint(1, max(1, int(math.isqrt(p)) + 2))
            n = RNG.randint(1, max(1, (p - 1) // m))  # mn < p
            A = rand_set(p, m)
            B = rand_set(p, n)
            S = sum(chi[(a + b) % p] for a in A for b in B)
            R = sum(pow((a + b) % p, d, p) for a in A for b in B) % p
            r = len(set(B) & {(-a) % p for a in A})
            require(R == S % p, "S1a_euler_congruence", (p, A, B))
            if m * n < p / 2:
                require(lar(R, p) == S, key1, (p, A, B))
            # 2p-lift with parity, valid for mn < p
            cands = [x for x in range(-m * n, m * n + 1)
                     if (x - R) % p == 0 and (x - (m * n - r)) % 2 == 0]
            require(cands == [S], key2, (p, A, B, cands, S))
            if m * n > p / 2 and lar(R, p) != S and wit_fail_half is None:
                wit_fail_half = {"p": p, "A": A, "B": B, "S": S, "lar": lar(R, p), "mn": m * n}
    RESULTS["witnesses"]["S1a_lar_differs_from_S_when_mn_gt_p_over_2"] = wit_fail_half
    # structured check: complete bicliques A={0,t}, B = {b: chi(b)=chi(b+t)=1}, m n > p/2
    for p in [101, 103, 107, 109, 113]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        A = [0, 1]
        B = [b for b in range(p) if chi[b] == 1 and chi[(b + 1) % p] == 1]
        S = sum(chi[(a + b) % p] for a in A for b in B)
        R = sum(pow((a + b) % p, d, p) for a in A for b in B) % p
        require(S == 2 * len(B), "S1a_biclique_S_equals_mn")
        if 2 * len(B) > p / 2:
            RESULTS["witnesses"].setdefault("S1a_biclique_lar", []).append(
                {"p": p, "mn": 2 * len(B), "S": S, "lar": lar(R, p)})
    save()


def section_1a_witness():
    """explicit witness that the least absolute residue is not S once S > p/2
    (p/2 < mn < p); the 2p-lift with parity still recovers S."""
    out = []
    for p in [101, 103, 197, 401]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        A = [0, 1, 2]
        score = sorted(range(p), key=lambda b: -sum(chi[(a + b) % p] for a in A))
        n = (p - 1) // 3
        B = sorted(score[:n])
        S = sum(chi[(a + b) % p] for a in A for b in B)
        R = sum(pow((a + b) % p, d, p) for a in A for b in B) % p
        r = len(set(B) & {(-a) % p for a in A})
        cands = [x for x in range(-3 * n, 3 * n + 1)
                 if (x - R) % p == 0 and (x - (3 * n - r)) % 2 == 0]
        require(cands == [S], "S1a_lift_2p_witness")
        require(lar(R, p) != S, "S1a_lar_differs_witness")
        out.append({"p": p, "A": A, "n": n, "mn": 3 * n, "S": S, "lar": lar(R, p), "r": r})
    RESULTS["witnesses"]["S1a_lar_differs_from_S_when_mn_gt_p_over_2"] = out
    save()


# ---------------------------------------------------------------- S1b
def section_1b():
    # (i) binom(d,j) = (-1/4)^j C(2j,j) mod p for 0<=j<=d; C(2j,j)=0 mod p for d<j<p;
    #     symmetry beta_j = beta_{d-j}
    for p in primes_upto(3000):
        d = (p - 1) // 2
        i4 = inv(-4, p)
        bd, cc, pw = 1, 1, 1  # binom(d,j), C(2j,j), (-1/4)^j
        betas = []
        ok = True
        for j in range(0, p):
            if j <= d:
                ok &= (bd == cc * pw % p)
                betas.append(bd)
            else:
                ok &= (cc == 0)
            # advance
            bd = bd * (d - j) % p * inv(j + 1, p) % p if j + 1 < p else 0
            cc = cc * 2 * (2 * j + 1) % p * inv(j + 1, p) % p if j + 1 < p else 0
            pw = pw * i4 % p
        require(ok, "S1b_binom_central", p)
        require(all(betas[j] == betas[d - j] for j in range(d + 1)), "S1b_beta_symmetry", p)
    # (ii) power series: T(u)^2 (1+u) = 1 mod u^p with T = sum_{j<p} (-1/4)^j C(2j,j) u^j,
    #      and (1+u)^d == T mod u^p
    for p in primes_upto(160):
        d = (p - 1) // 2
        f, fi = fact_tables(p)
        T = [binom_mod(2 * j, j, f, fi, p) * pow(inv(-4, p), j, p) % p if 2 * j < p else 0
             for j in range(p)]
        # C(2j,j) for 2j >= p is 0 mod p (j<p); computed above as 0 when 2j>=p
        onepd = [binom_mod(d, j, f, fi, p) for j in range(p)]
        require(ptrim(T[:]) == ptrim(onepd[:]), "S1b_series_equals_truncation", p)
        sq = pmul(pmul(T, T, p), [1, 1], p)
        require(ptrim(sq[:p]) == [1], "S1b_series_square_inverse", p)
    # (iii) chi(1+u) = T_d(u) for all u in F_p; chi(x0+u) = chi(x0) T_d(u/x0)
    for p in primes_upto(700):
        d = (p - 1) // 2
        chi = legendre_table(p)
        f, fi = fact_tables(p)
        Td = [binom_mod(2 * j, j, f, fi, p) * pow(inv(-4, p), j, p) % p for j in range(d + 1)]
        ok = all((peval(Td, u, p) - chi[(1 + u) % p]) % p == 0 for u in range(p))
        require(ok, "S1b_chi_is_truncated_binomial_series", p)
        for _ in range(5):
            x0 = RNG.randrange(1, p)
            u = RNG.randrange(p)
            require((chi[x0] * peval(Td, u * inv(x0, p) % p, p) - chi[(x0 + u) % p]) % p == 0,
                    "S1b_local_expansion")
        # derivative identity (x^d)^{(j)}(x0) = (d)_j x0^{d-j} = (-1/2)_j chi(x0) x0^{-j}
        for _ in range(5):
            x0 = RNG.randrange(1, p)
            j = RNG.randrange(0, d + 1)
            lhs = 1
            for i in range(j):
                lhs = lhs * (d - i) % p
            lhs = lhs * pow(x0, d - j, p) % p
            rhs = 1
            for i in range(j):
                rhs = rhs * (-inv(2, p) - i) % p
            rhs = rhs * chi[x0] * pow(inv(x0, p), j, p) % p
            require(lhs == rhs, "S1b_derivative_minus_half")
    # (iv) kernel formula S = sum_j beta_j P_j(A) P_{d-j}(B) mod p, 0^0 = 1
    for p in primes_upto(500)[5:]:
        d = (p - 1) // 2
        chi = legendre_table(p)
        f, fi = fact_tables(p)
        beta = [binom_mod(2 * j, j, f, fi, p) * pow(inv(-4, p), j, p) % p for j in range(d + 1)]
        for _ in range(3):
            m = RNG.randint(1, 12)
            n = RNG.randint(1, 12)
            A = rand_set(p, m)
            B = rand_set(p, n)
            PA = [sum(pow(a, j, p) if (a or j == 0) else 0 for a in A) % p for j in range(d + 1)]
            PB = [sum(pow(b, j, p) if (b or j == 0) else 0 for b in B) % p for j in range(d + 1)]
            val = sum(beta[j] * PA[j] % p * PB[d - j] for j in range(d + 1)) % p
            S = sum(chi[(a + b) % p] for a in A for b in B)
            require(val == S % p, "S1b_kernel_formula", (p, A, B))
    save()


# ---------------------------------------------------------------- S1c
def section_1c():
    """deg gcd_{a in T}((x+a)^d - 1) = #common roots <= d/|T|;
    with zeros allowed: gcd of (x+a)^{d+1}-(x+a) has degree n' with |T| n' <= d + |B' cap -T|;
    N_e(A) = #{b : b is a root of >= m-e of the h_a = (x+a)^{d+1}-(x+a)} = #{b: e_b <= e}."""
    worst = []
    for p in [p for p in primes_upto(160) if p >= 7]:
        d = (p - 1) // 2
        f, fi = fact_tables(p)
        chi = legendre_table(p)
        for size in range(1, 6):
            for _ in range(3):
                T = rand_set(p, size)
                g = None
                g0 = None
                for a in T:
                    fa = padd(translate_power(a, d, f, fi, p), [p - 1], p)
                    ha = padd(translate_power(a, d + 1, f, fi, p), [(-a) % p, p - 1], p)
                    g = fa if g is None else pgcd(g, fa, p)
                    g0 = ha if g0 is None else pgcd(g0, ha, p)
                common = [b for b in range(p) if all(chi[(a + b) % p] == 1 for a in T)]
                common0 = [b for b in range(p) if all(chi[(a + b) % p] >= 0 for a in T)]
                r0 = len(set(common0) & {(-a) % p for a in T})
                require(len(g) - 1 == len(common), "S1c_gcd_degree_is_common_roots", (p, T))
                require(size * (len(g) - 1) <= d, "S1c_gcd_HP_bound", (p, T))
                require(len(g0) - 1 == len(common0), "S1c_gcd0_degree", (p, T))
                require(size * (len(g0) - 1) <= d + r0, "S1c_gcd0_HP_bound", (p, T))
                worst.append((size * (len(g) - 1) / d, p, T))
        # N_e as a count of points on many polynomials
        for _ in range(4):
            m = RNG.randint(2, 7)
            A = rand_set(p, m)
            eb, z = bad_counts(A, p, chi)
            on = [sum(1 for a in A if pow((b + a) % p, d + 1, p) == (b + a) % p) for b in range(p)]
            for e in range(0, m):
                require(sum(1 for b in range(p) if eb[b] <= e) == sum(1 for b in range(p) if on[b] >= m - e),
                        "S1c_Ne_equals_points_on_many")
    worst.sort(reverse=True)
    RESULTS["tables"]["S1c_max_ratio_size_times_deg_over_d"] = [
        {"ratio": round(w[0], 4), "p": w[1], "T": w[2]} for w in worst[:5]]
    save()


# ---------------------------------------------------------------- S2
def order_synth(a, b, p, cap=None):
    """multiplicity of root b of nonzero a, by repeated synthetic division"""
    k = 0
    a = a[:]
    while True:
        # divide by (x-b)
        n = len(a)
        if n == 0:
            return None
        q = [0] * (n - 1)
        acc = 0
        for i in range(n - 1, -1, -1):
            acc = (acc * b + a[i]) % p
            if i > 0:
                q[i - 1] = acc
        if acc != 0:
            return k
        k += 1
        a = q
        if cap is not None and k >= cap:
            return k


def section_2():
    stats = {"ord_at_bad": {}, "phi_order_at_biclique": {}, "W1_min_orders": [], "W2_min_orders": []}
    prime_list = [p for p in primes_upto(140) if p >= 29]
    for p in prime_list:
        d = (p - 1) // 2
        f, fi = fact_tables(p)
        chi = legendre_table(p)
        for trial in range(3):
            M = RNG.randint(3, min(8, int(math.isqrt(p)) + 1))
            if trial == 0:
                # a near-biclique: greedy set with many common square partners
                A = [0]
                while len(A) < M:
                    best = max((x for x in range(p) if x not in A),
                               key=lambda x: sum(1 for b in range(p)
                                                 if all(chi[(a + b) % p] == 1 for a in A + [x])))
                    A.append(best)
                A.sort()
            else:
                A = rand_set(p, M)
            F, c, D = hp_poly(A, p, f, fi)
            require(len(F) - 1 == d, "S2_degF_equals_d", (p, A))
            # (a) sum_k c_k y_k^N = h_{N-M+1}(y)
            x0 = RNG.randrange(p)
            ys = [(a + x0) % p for a in A]
            h = complete_h(ys, D + 1, p)
            for N in range(0, D + 1):
                lhs = sum(ck * pow(y, N, p) if (y or N == 0) else 0 for ck, y in zip(c, ys)) % p
                rhs = h[N - M + 1] if N - M + 1 >= 0 else 0
                require(lhs == rhs, "S2_divided_difference_is_h", (p, A, N))
            # (b) F = h_d(A+x)-1 and F^{(j)} = (D)_j h_{d-j}(A+x), pointwise at every x
            derivs = [F]
            for j in range(1, M + 2):
                derivs.append(pderiv(derivs[-1], p))
            falling = [1]
            for j in range(1, M + 2):
                falling.append(falling[-1] * (D - j + 1) % p)
            for x in range(p):
                hx = complete_h([(a + x) % p for a in A], d, p)
                require((peval(F, x, p) - (hx[d] - 1)) % p == 0, "S2_F_is_hd_minus_1")
                for j in range(1, M + 2):
                    require((peval(derivs[j], x, p) - falling[j] * hx[d - j]) % p == 0,
                            "S2_derivative_is_shift")
            eb, z = bad_counts(A, p, chi)
            negA = {(-a) % p for a in A}
            ords = {}
            for b in range(p):
                ys = [(a + b) % p for a in A]
                o = order_synth(F, b, p, cap=M + 3)
                ords[b] = o
                if b not in negA:
                    E = [k for k in range(M) if chi[ys[k]] == -1]
                    # (d) syndrome lemma
                    for j in range(M):
                        lhs = peval(derivs[j], b, p) * inv(falling[j], p) % p
                        rhs = -2 * sum(c[k] * pow(ys[k], M - 1 - j, p) for k in E) % p
                        require(lhs == rhs, "S2_syndrome_lemma", (p, A, b, j))
                    # (c) twisted HP lemma h_{d-j}(y) = sum c_k chi(y_k) y_k^{M-1-j}
                    hy = complete_h(ys, d, p)
                    for j in range(M):
                        rhs = sum(c[k] * chi[ys[k]] * pow(ys[k], M - 1 - j, p) for k in range(M)) % p
                        require(hy[d - j] == rhs, "S2_twisted_HP_lemma")
                    # (g) periodicity defect: h_{i+d} - h_i = -2 sum_{E} c_k y_k^{i+M-1}, 0<=i<=d-? 
                    hy2 = complete_h(ys, d + M + 3, p)
                    for i in range(0, M + 3):
                        lhs = (hy2[i + d] - hy2[i]) % p
                        rhs = -2 * sum(c[k] * pow(ys[k], i + M - 1, p) for k in E) % p
                        require(lhs == rhs, "S2_periodicity_defect")
                    # (e) orders
                    if eb[b] == 0:
                        require(o >= M, "S2_order_at_complete_point", (p, A, b))
                    else:
                        require(o <= eb[b] - 1, "S2_order_at_bad_point_le_e_minus_1", (p, A, b, o, eb[b]))
                        key = f"e={min(eb[b], 4)}"
                        stats["ord_at_bad"].setdefault(key, {})
                        stats["ord_at_bad"][key][o] = stats["ord_at_bad"][key].get(o, 0) + 1
                else:
                    if eb[b] == 0:
                        require(o >= M - 1, "S2_order_at_zero_point", (p, A, b))
            # (f) unweighted Phi_A = sum (x+a)^d : order of Phi_A - m at complete points
            Phi = [0]
            for a in A:
                Phi = padd(Phi, translate_power(a, d, f, fi, p), p)
            Phi_m = padd(Phi, [(-M) % p], p)
            for b in range(p):
                if b not in negA and eb[b] == 0:
                    o = order_synth(Phi_m, b, p, cap=M + 2)
                    stats["phi_order_at_biclique"][o] = stats["phi_order_at_biclique"].get(o, 0) + 1
                    # jet identity Phi^{(1)}(b) = d * sum_a (a+b)^{-1}
                    require(peval(pderiv(Phi, p), b, p) == d * sum(inv(a + b, p) for a in A) % p,
                            "S2_phi_jet_negative_power_sum")
            # (h) Hankel polynomials W1 (e=1) and W2 (e=2)
            z = [pscale(derivs[j], inv(falling[j], p), p) for j in range(5)]
            W1 = padd(pmul(z[0], z[2], p), pscale(pmul(z[1], z[1], p), p - 1, p), p)
            require(len(W1) - 1 == 2 * d - 2, "S2_W1_degree_exact", (p, A))
            lead = (binom_mod(D, M - 1, f, fi, p) * binom_mod(D - 2, M - 1, f, fi, p)
                    - binom_mod(D - 1, M - 1, f, fi, p) ** 2) % p
            require(lead == W1[-1] and lead != 0, "S2_W1_leading_coefficient")
            n0 = sum(1 for b in range(p) if b not in negA and eb[b] == 0)
            n1 = sum(1 for b in range(p) if b not in negA and eb[b] == 1)
            mo = {0: 10 ** 9, 1: 10 ** 9}
            for b in range(p):
                if b not in negA and eb[b] <= 1:
                    o = order_synth(W1, b, p, cap=2 * M + 2)
                    mo[eb[b]] = min(mo[eb[b]], o)
                    require(o >= (2 * M - 2 if eb[b] == 0 else M - 2), "S2_W1_orders", (p, A, b, o))
            require((M - 2) * n1 + (2 * M - 2) * n0 <= 2 * d - 2, "S2_W1_count_inequality")
            stats["W1_min_orders"].append({"p": p, "M": M, "min_ord_e0": mo[0], "min_ord_e1": mo[1],
                                           "n0": n0, "n1": n1})
            if M >= 5:
                # W2 = det [[z0,z1,z2],[z1,z2,z3],[z2,z3,z4]]
                def m3(a, b_, c_):
                    return pmul(pmul(a, b_, p), c_, p)
                t = [0]
                t = padd(t, m3(z[0], z[2], z[4]), p)
                t = padd(t, pscale(m3(z[0], z[3], z[3]), p - 1, p), p)
                t = padd(t, pscale(m3(z[1], z[1], z[4]), p - 1, p), p)
                t = padd(t, pscale(m3(z[1], z[2], z[3]), 2, p), p)
                t = padd(t, pscale(m3(z[2], z[2], z[2]), p - 1, p), p)
                W2 = t
                degW2 = len(W2) - 1
                require(degW2 <= 3 * d - 6, "S2_W2_degree_bound")
                mo2 = {}
                n2 = 0
                for b in range(p):
                    if b not in negA and eb[b] <= 2:
                        n2 += 1
                        o = order_synth(W2, b, p, cap=3 * M + 2)
                        mo2[eb[b]] = min(mo2.get(eb[b], 10 ** 9), o)
                        require(o >= M - 4, "S2_W2_order_ge_M_minus_4", (p, A, b, o))
                stats["W2_min_orders"].append({"p": p, "M": M, "deg": degW2, "min_orders": mo2, "N2": n2,
                                               "bound_3(d-2)/(M-4)": round(3 * (d - 2) / (M - 4), 2)})
    # (i) multiset witness: A = {0} with multiplicity m, B = Q -> "S" = m d
    p = 101
    chi = legendre_table(p)
    d = 50
    Bq = [b for b in range(1, p) if chi[b] == 1]
    m = 10
    Smult = m * sum(chi[b] for b in Bq)
    require(Smult == m * d and m * len(Bq) > p, "S2_multiset_witness")
    RESULTS["witnesses"]["S2_multiset_breaks_HP"] = {"p": p, "A": "{0} with multiplicity 10", "B": "Q",
                                                     "mn": m * len(Bq), "S": Smult}
    stats["ord_at_bad"] = {k: {str(kk): vv for kk, vv in v.items()} for k, v in stats["ord_at_bad"].items()}
    stats["phi_order_at_biclique"] = {str(k): v for k, v in stats["phi_order_at_biclique"].items()}
    RESULTS["tables"]["S2"] = stats
    save()


# ---------------------------------------------------------------- S2j sign-pattern HP
def section_2j():
    """For every u in {+-1}^A: mu(u) = #{b not in -A : chi(a+b) = u_a for all a} <= (d+m-1)/m,
    and <= d/m for constant u.  Exhaustive for small p, m=3,4 (a_1 = 0 by translation); random beyond."""
    import itertools
    table = []
    for p in [p for p in primes_upto(47) if p >= 11]:
        d = (p - 1) // 2
        chi = np.array(legendre_table(p), dtype=np.int64)
        for m in (3, 4):
            if m > (p + 1) // 2:
                continue
            worst_mixed, worst_const = 0.0, 0.0
            wm = None
            for rest in itertools.combinations(range(1, p), m - 1):
                A = (0,) + rest
                X = chi[(np.array(A)[:, None] + np.arange(p)[None, :]) % p]  # m x p
                ok = np.all(X != 0, axis=0)
                cols = X[:, ok]
                keys = ((cols + 1) // 2 * (1 << np.arange(m))[:, None]).sum(axis=0)
                counts = np.bincount(keys, minlength=1 << m)
                const_max = max(counts[0], counts[(1 << m) - 1])
                mixed_max = counts[1:(1 << m) - 1].max()
                bump("S2j_sign_pattern_sets", 1)
                require(m * mixed_max <= d + m - 1, "S2j_mixed_bound", (p, A))
                require(m * const_max <= d, "S2j_constant_bound", (p, A))
                if m * mixed_max / d > worst_mixed:
                    worst_mixed = m * mixed_max / d
                    wm = A
                worst_const = max(worst_const, m * const_max / d)
            table.append({"p": p, "m": m, "max m*mu(mixed)/d": round(worst_mixed, 4),
                          "bound (d+m-1)/d": round((d + m - 1) / d, 4), "argmax": list(wm),
                          "max m*mu(const)/d": round(worst_const, 4)})
    # random larger
    for p in [101, 211, 401, 1009]:
        d = (p - 1) // 2
        chi = np.array(legendre_table(p), dtype=np.int64)
        for _ in range(30):
            m = RNG.randint(2, 12)
            A = rand_set(p, m)
            X = chi[(np.array(A)[:, None] + np.arange(p)[None, :]) % p]
            ok = np.all(X != 0, axis=0)
            keys = ((X[:, ok] + 1) // 2 * (1 << np.arange(m))[:, None]).sum(axis=0)
            counts = np.bincount(keys, minlength=1 << m)
            require(m * counts.max() <= d + m - 1, "S2j_mixed_bound_random", (p, A))
    # the auxiliary polynomial F_u vanishes to order m at every b in the fibre (spot check)
    for p in [53, 61, 97]:
        d = (p - 1) // 2
        f, fi = fact_tables(p)
        chil = legendre_table(p)
        for _ in range(4):
            m = RNG.randint(2, 5)
            A = rand_set(p, m)
            b0 = RNG.choice([b for b in range(p) if all(chil[(a + b) % p] != 0 for a in A)])
            u = [chil[(a + b0) % p] for a in A]
            c = hp_weights(A, p)
            D = d + m - 1
            Fu = [p - 1]
            for ck, uk, ak in zip(c, u, A):
                Fu = padd(Fu, pscale(translate_power(ak, D, f, fi, p), ck * uk % p, p), p)
            require(len(Fu) > 0, "S2j_Fu_nonzero")
            for b in range(p):
                if all(chil[(a + b) % p] == uk for a, uk in zip(A, u)):
                    require(order_synth(Fu, b, p, cap=m + 1) >= m, "S2j_Fu_order", (p, A, b))
    RESULTS["tables"]["S2j_sign_pattern"] = table
    save()


# ---------------------------------------------------------------- S3 list decoding
def greedy_set(chi_mat, p, m, e, rng, start=None):
    """greedily add elements maximising #{b not in -A : e_b <= e} (ties random)"""
    neg = (chi_mat == -1).astype(np.int32)   # neg[a, b] = [chi(a+b) = -1]
    zero = (chi_mat == 0)
    A = [start if start is not None else rng.randrange(p)]
    eb = neg[A[0]].copy()
    zb = zero[A[0]].copy()
    while len(A) < m:
        cand = np.array([x for x in range(p) if x not in A])
        new_e = eb[None, :] + neg[cand]
        new_z = zb[None, :] | zero[cand]
        score = ((new_e <= e) & (~new_z)).sum(axis=1)
        best = np.flatnonzero(score == score.max())
        x = int(cand[rng.choice(list(best))])
        A.append(x)
        eb = eb + neg[x]
        zb = zb | zero[x]
    return sorted(A)


def list_stats(A, p, chi_arr, emax=4):
    d = (p - 1) // 2
    m = len(A)
    X = chi_arr[(np.array(A)[:, None] + np.arange(p)[None, :]) % p]   # m x p
    off = np.all(X != 0, axis=0)          # b not in -A
    eb = (X == -1).sum(axis=0)
    W = X[:, off].T.astype(np.int64)      # words, one per b not in -A
    out = {"p": p, "m": m, "A": A}
    # Gram of the code, distinct words
    words, inv_idx, mult = np.unique(W, axis=0, return_inverse=True, return_counts=True)
    G = words @ words.T
    np.fill_diagonal(G, -10 ** 9)
    lam = G.max() / m if len(words) > 1 else -1.0
    out["distinct_words"] = int(len(words))
    out["code_min_distance"] = int((m - G.max()) // 2) if len(words) > 1 else None
    out["lambda_max_corr_distinct"] = round(float(lam), 4)
    # spectral Johnson: largest eigenvalue of W^T W
    lam_spec = float(np.linalg.eigvalsh((W.T @ W).astype(float)).max())
    out["spectral_WtW"] = round(lam_spec, 3)
    require(lam_spec <= p + 1e-6, "S3_spectral_WtW_le_p", (p, A))
    rows = []
    for e in range(0, min(emax, (m - 1) // 2) + 1):
        in_ball = eb[off] <= e
        Ne = int(in_ball.sum())
        full = int((eb <= e).sum())
        on_negA = int(((eb <= e) & ~off).sum())
        rhp_lhs = m * full - on_negA          # m|B| - |B cap -A| with B = {b : e_b <= e}
        tau = 1 - 2 * e / m
        sm_b = m * (p - m) / (m - 2 * e) ** 2
        john_list = p * m / ((m - 2 * e) ** 2 + m)
        john_spec = lam_spec * m / (m - 2 * e) ** 2
        # distinct words inside the ball and their multiplicities
        ball_words = np.flatnonzero(((words == -1).sum(axis=1)) <= e)
        s_act = int(len(ball_words))
        mu_max = int(mult[ball_words].max()) if s_act else 0
        if tau ** 2 > lam:
            s_j = (1 - lam) / (tau ** 2 - lam)
        else:
            s_j = float("inf")
        row = {"e": e, "N_e_x": Ne, "rhp_lhs_over_d": round(rhp_lhs / d, 4),
               "mN_over_d": round(m * Ne / d, 4),
               "second_moment_over_d": round(m * sm_b / d, 4),
               "johnson_listvar_over_d": round(m * john_list / d, 4),
               "johnson_spectral_over_d": round(m * john_spec / d, 4),
               "distinct_in_ball": s_act, "max_mult_in_ball": mu_max,
               "distinct_johnson_s_bound": (round(s_j, 3) if s_j != float("inf") else "inf"),
               "distinct_johnson_times_fibre_over_d": (round(s_j * (d + m - 1) / d, 3)
                                                       if s_j != float("inf") else "inf")}
        # rigorous inequalities that must hold
        require(Ne <= sm_b + 1e-9, "S3_second_moment_b", (p, A, e))
        require(Ne <= john_list + 1e-9, "S3_johnson_listvar", (p, A, e))
        require(Ne <= john_spec + 1e-9, "S3_johnson_spectral", (p, A, e))
        if e == 0:
            require(rhp_lhs <= d, "S3_HP", (p, A))
        if e == 1 and m >= 3:
            require(Ne <= (2 * d - 2) / (m - 2) + 1e-9, "S3_hankel_e1", (p, A))
            row["hankel_e1_over_d"] = round(m * (2 * d - 2) / (m - 2) / d, 4)
        require(s_act * 1.0 <= s_j + 1e-9, "S3_distinct_johnson", (p, A, e))
        require(mu_max * m <= d + m - 1, "S3_fibre_bound", (p, A, e))
        rows.append(row)
    out["rows"] = rows
    return out


def section_3():
    rng = random.Random(7)
    tables = []
    for p in [101, 197, 401, 809]:
        chi_l = legendre_table(p)
        chi_arr = np.array(chi_l, dtype=np.int64)
        chi_mat = chi_arr[(np.arange(p)[:, None] + np.arange(p)[None, :]) % p]
        mmax = int(math.isqrt(p))
        for m in sorted({4, 6, mmax // 2 + 1, mmax}):
            fams = {"random": rand_set(p, m), "interval": list(range(1, m + 1))}
            for e in (0, 1, 2):
                if 2 * e < m:
                    fams[f"greedy_e{e}"] = greedy_set(chi_mat, p, m, e, rng)
            for name, A in fams.items():
                st = list_stats(A, p, chi_arr)
                st["family"] = name
                tables.append(st)
    RESULTS["tables"]["S3_list_decoding"] = tables
    # summary: best (smallest) Johnson-type constant versus the truth
    summ = []
    for st in tables:
        for row in st["rows"]:
            if row["e"] >= 1:
                summ.append({"p": st["p"], "m": st["m"], "family": st["family"], "e": row["e"],
                             "truth mN/d": row["mN_over_d"],
                             "best Johnson-type /d": min(row["second_moment_over_d"],
                                                         row["johnson_listvar_over_d"],
                                                         row["johnson_spectral_over_d"])})
    RESULTS["tables"]["S3_summary"] = summ
    save()


# ---------------------------------------------------------------- S4 F_{p^2}
def section_4():
    """Over F_q, q=p^2, A=F_p: the second-moment inequality is attained up to a factor (p-1)/p,
    and m|B| - r = p^2 - p = (2 - o(1)) (q-1)/2 with A+B in Q_q u {0}.  So f(0) >= 2 for any
    statement valid over all finite fields of odd order, and the second moment is optimal there."""
    out = []
    for p in [7, 11, 13, 17, 19, 23]:
        nu = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
        q = p * p
        dq = (q - 1) // 2

        def mul(u, v):
            return ((u[0] * v[0] + nu * u[1] * v[1]) % p, (u[0] * v[1] + u[1] * v[0]) % p)

        def pw(u, e):
            r = (1, 0)
            while e:
                if e & 1:
                    r = mul(r, u)
                u = mul(u, u)
                e >>= 1
            return r

        elems = [(x, y) for x in range(p) for y in range(p)]
        chiq = {}
        for z in elems:
            if z == (0, 0):
                chiq[z] = 0
            else:
                t = pw(z, dq)
                chiq[z] = 1 if t == (1, 0) else -1
        require(all(chiq[(x, 0)] == 1 for x in range(1, p)), "S4_subfield_all_squares", p)
        A = [(x, 0) for x in range(p)]
        FA = {}
        for b in elems:
            FA[b] = sum(chiq[((a[0] + b[0]) % p, (a[1] + b[1]) % p)] for a in A)
        m = p
        second = sum(v * v for v in FA.values())
        require(second == m * (q - m), "S4_second_moment_identity_Fq", p)
        sub = sum(FA[(x, 0)] ** 2 for x in range(p))
        good = [b for b in elems if all(chiq[((a[0] + b[0]) % p, (a[1] + b[1]) % p)] >= 0 for a in A)]
        r = sum(1 for b in good if b[1] == 0)   # b in -A = F_p
        lhs = m * len(good) - r
        require(lhs >= p * p - p, "S4_subfield_biclique_size", p)
        require(lhs > dq, "S4_HP_analogue_fails_over_Fp2", p)
        out.append({"p": p, "q": q, "m": m, "n": len(good), "r": r, "mn_minus_r": lhs,
                    "(q-1)/2": dq, "ratio_to_(q-1)/2": round(lhs / dq, 4),
                    "second_moment_share_of_subfield": f"{sub}/{second}"})
    RESULTS["tables"]["S4_Fp2"] = out
    save()


# ---------------------------------------------------------------- S3b coding-parameter facts
def section_3b():
    """(i) the only polynomials of degree <= d with values in {+-1} on F_p minus one point and 0 there
    are +-(x+a)^d; with no zero, only +-1 (exhaustive p = 7, 11);
    (ii) each translate agrees with the all-ones word in exactly d < d+1 = dim positions;
    (iii) closed form of the W_1 leading coefficient."""
    import itertools
    for p in (7, 11):
        d = (p - 1) // 2
        xs = np.arange(p)
        V = np.stack([xs ** j % p for j in range(d + 1)], axis=1)  # p x (d+1)
        found_one_zero, found_no_zero = [], []
        # enumerate coefficient vectors in chunks
        total = p ** (d + 1)
        idx = np.arange(total, dtype=np.int64)
        coeffs = np.stack([(idx // p ** j) % p for j in range(d + 1)], axis=1)
        vals = (coeffs @ V.T) % p                      # total x p
        is_pm = (vals == 1) | (vals == p - 1)
        nz = (vals == 0).sum(axis=1)
        ok1 = (nz == 1) & ((is_pm | (vals == 0)).all(axis=1))
        ok0 = (nz == 0) & is_pm.all(axis=1)
        found_one_zero = coeffs[ok1]
        found_no_zero = coeffs[ok0]
        require(len(found_no_zero) == 2, "S3b_pm1_valued_are_constants", p)
        require(len(found_one_zero) == 2 * p, "S3b_one_zero_are_translates", p)
        f, fi = fact_tables(p)
        expected = set()
        for a in range(p):
            t = translate_power(a, d, f, fi, p)
            t = t + [0] * (d + 1 - len(t))
            expected.add(tuple(t))
            expected.add(tuple((-x) % p for x in t))
        require({tuple(int(x) for x in c) for c in found_one_zero} == expected,
                "S3b_one_zero_are_translates_exactly", p)
        bump("S3b_polynomials_enumerated", total)
    for p in primes_upto(300)[2:]:
        chi = legendre_table(p)
        d = (p - 1) // 2
        a = RNG.randrange(p)
        require(sum(1 for y in range(p) if chi[(a + y) % p] == 1) == d, "S3b_agreement_is_d")
        f, fi = fact_tables(p)
        for M in range(2, min(12, d)):
            D = d + M - 1
            lhs = (binom_mod(D, M - 1, f, fi, p) * binom_mod(D - 2, M - 1, f, fi, p)
                   - binom_mod(D - 1, M - 1, f, fi, p) ** 2) % p
            rhs = (-binom_mod(D - 1, M - 1, f, fi, p) ** 2 * (M - 1) * inv(d * (D - 1), p)) % p
            require(lhs == rhs and lhs != 0, "S3b_W1_leading_closed_form", (p, M))
    save()


def section_1b_fq():
    """Lemma 1.4 is field-agnostic: for q = p^k, binom(d_q, j) = (-1/4)^j C(2j,j) mod p for j <= d_q
    and C(2j,j)(-1/4)^j = 0 mod p for d_q < j < q (integer computation)."""
    for p, k in [(3, 2), (3, 3), (5, 2), (7, 2), (11, 2), (13, 2), (5, 3)]:
        q = p ** k
        dq = (q - 1) // 2
        i4 = inv(-4, p)
        ok = True
        for j in range(q):
            beta = math.comb(2 * j, j) % p * pow(i4, j, p) % p
            if j <= dq:
                ok &= (math.comb(dq, j) % p == beta)
            else:
                ok &= (beta == 0)
        require(ok, "S1b_Fq_kernel_field_agnostic", (p, k))
    save()


def main():
    section_1a()
    section_1a_witness()
    section_1b()
    section_1b_fq()
    section_1c()
    section_2()
    section_2j()
    section_3()
    section_3b()
    section_4()
    RESULTS["status"] = "all checks passed" if not RESULTS["failures"] else "FAILURES PRESENT"
    save()
    print(json.dumps({"total_checks": RESULTS["total_checks"], "n_failures": RESULTS["n_failures"],
                      "elapsed_seconds": RESULTS["elapsed_seconds"]}))
    return 0 if not RESULTS["failures"] else 1


if __name__ == "__main__":
    sys.exit(main())
