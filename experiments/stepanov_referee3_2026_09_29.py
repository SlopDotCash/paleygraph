#!/usr/bin/env python3
"""Independent referee verifier (`referee3`) for research/stepanov-sharpen-2026-09-27.md
and for Proposition 2.9 of research/stepanov-robust-2026-09-26.md.

Written from scratch by the referee; the worker's verifier was not imported.
Run:  /opt/miniconda3/bin/python3 experiments/stepanov_referee3_2026_09_29.py [sections]
Writes results/stepanov_referee3_2026_09_29.json.

Sections
  A  Prop. 2.9 of [R] ((star)_k): exhaustive over A containing 0 for small (p,k), random and
     structured A for larger p, and a direct computation of the Hankel determinant H_{e+1}
     over F_p (degree, leading coefficient, vanishing orders) for a sample.
  B  Lemma 1.1, Corollary 1.2, Theorem 1.4 and Theorem 1.5 of [Sh]: exhaustive over ALL
     pairs (A,B) for p in {7,11,13,17,19,23} (A normalised by an affine map to contain 0,1;
     for every (m,n,r) the extremal B is computed exactly).  Exact integer arithmetic.
  C  Theorem 1.3 of [Sh]: analytic constants, an independent rational-interval check of the
     finite range, extended exact checks, float scan.
  D  Proposition 2.1 of [Sh] (the LP): exact check of the feasible profiles (ii), float LPs.
  E  Corollary 3.3 of [Sh]: beta_m in exact rational arithmetic.
  F  Section 4 of [Sh]: Lemma 3.1, Lemma 4.2, Theorem 4.3, Prop. 4.1 at primes where the
     hypotheses hold (numpy, exact integers).
  G  Section 5 of [Sh]: Theorem 5.1 / Cor. 5.2 for characters of order k (exhaustive over all
     (A,B) at p = 7, 13; random at larger p), vertex lemma, sharpness example.
  H  Numbers quoted in [Sh] and in Section 9 of paper/robust-hanson-petridis.md.
"""
import itertools, json, math, os, random, sys, time
from fractions import Fraction as Fr
import numpy as np

T0 = time.time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_referee3_2026_09_29.json")
RES = {"checks": {}, "failures": [], "witnesses": {}, "data": {}, "timing": {}}


def chk(key, cond, info=None, n=1):
    RES["checks"][key] = RES["checks"].get(key, 0) + n
    if not cond:
        RES["failures"].append({"key": key, "info": str(info)[:800]})
    return cond


def note(key, val):
    RES["data"][key] = val


# ----------------------------------------------------------------------------- basics
def is_prime(n):
    if n < 2:
        return False
    for q in range(2, int(n ** 0.5) + 1):
        if n % q == 0:
            return False
    return True


def next_prime(n):
    while not is_prime(n):
        n += 1
    return n


def chi_table(p):
    d = (p - 1) // 2
    t = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        t[x] = 1 if pow(x, d, p) == 1 else -1
    return t


def prim_root(p):
    fs = [q for q in range(2, p) if (p - 1) % q == 0 and is_prime(q)]
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs):
            return g


def dlog_table(p):
    g = prim_root(p)
    L = np.full(p, -1, dtype=np.int64)
    x = 1
    for i in range(p - 1):
        L[x] = i
        x = x * g % p
    return L, g


def isqrt_frac_lo(X, bits=200):
    """rational lower bound for sqrt(X), X >= 0 rational."""
    X = Fr(X)
    if X <= 0:
        return Fr(0)
    S = 1 << bits
    return Fr(math.isqrt(X.numerator * S * S // X.denominator), S)


def isqrt_frac_hi(X, bits=200):
    X = Fr(X)
    if X <= 0:
        return Fr(0)
    S = 1 << bits
    q = -(-X.numerator * S * S // X.denominator)
    r = math.isqrt(q)
    if r * r < q:
        r += 1
    return Fr(r, S)


def u_float(w):
    return (math.sqrt(12 * w - 3 * w * w) - w) / 2


def le_u(x, w):
    """exact: x <= u(w) = (sqrt(12w-3w^2)-w)/2 for rationals x, w in (0,1]."""
    x, w = Fr(x), Fr(w)
    L = 2 * x + w
    return L <= 0 or L * L <= 12 * w - 3 * w * w


def u_lo(w, bits=200):
    w = Fr(w)
    return (isqrt_frac_lo(12 * w - 3 * w * w, bits) - w) / 2


def u_hi(w, bits=200):
    w = Fr(w)
    return (isqrt_frac_hi(12 * w - 3 * w * w, bits) - w) / 2


# ============================================================================ Section A
def bad_matrix_k(p, k):
    """BAD[a][b] = 1 iff a+b != 0 and (a+b)^{(p-1)/k} != 1."""
    dk = (p - 1) // k
    inH = np.zeros(p, dtype=bool)
    for x in range(1, p):
        inH[x] = pow(x, dk, p) == 1
    Bm = np.zeros((p, p), dtype=np.int16)
    for a in range(p):
        for b in range(p):
            y = (a + b) % p
            Bm[a, b] = 1 if (y != 0 and not inH[y]) else 0
    return Bm


def all_sets_with0(p, rows):
    """doubling construction over subsets of {1..p-1}; returns (E, DEL, m) for A = {0} u T,
    E[i,b] = sum_{a in A} rows[a][b], DEL[i,b] = [b in -A]."""
    E = rows[0][None, :].astype(np.int16).copy()
    DEL = np.zeros((1, p), dtype=bool)
    DEL[0, 0] = True
    M = np.array([1], dtype=np.int16)
    for a in range(1, p):
        E = np.concatenate([E, E + rows[a][None, :]], axis=0)
        D2 = DEL.copy()
        D2[:, (-a) % p] = True
        DEL = np.concatenate([DEL, D2], axis=0)
        M = np.concatenate([M, M + 1])
    return E, DEL, M


def star_k_check(E, DEL, M, p, dk, key, extra=None):
    """(star)_k: sum_{e_b<=e} (e+1-e_b)(2m-3e-e_b-2delta_b) <= 2(e+1)(dk-e)  (doubled)."""
    n_ok = 0
    best = (-1.0, None)
    eq_count = {}
    Mi = M.astype(np.int64)
    Ei = E.astype(np.int64)
    Di = DEL.astype(np.int64)
    emax = int(min((int(M.max()) - 1) // 2, dk // 2))
    for e in range(0, emax + 1):
        adm = (2 * e <= np.minimum(Mi - 1, dk)) & (Mi + dk - 1 <= p - 1)
        if not adm.any():
            continue
        Ea, Da, Ma = Ei[adm], Di[adm], Mi[adm]
        mask = Ea <= e
        term = (e + 1 - Ea) * (2 * Ma[:, None] - 3 * e - Ea - 2 * Da)
        term = np.where(mask, term, 0)
        neg = (term < 0).sum()
        lhs = term.sum(axis=1)
        rhs = 2 * (e + 1) * (dk - e)
        viol = np.nonzero(lhs > rhs)[0]
        chk(key, neg == 0 and len(viol) == 0,
            {"e": e, "neg_terms": int(neg), "viol_rows": viol[:3].tolist()}, n=int(adm.sum()))
        n_ok += int(adm.sum())
        if rhs > 0:
            rat = lhs / rhs
            i = int(np.argmax(rat))
            if rat[i] > best[0]:
                best = (float(rat[i]), {"e": e, "m": int(Ma[i]), "ratio": float(rat[i])})
            eq_count[e] = int((lhs == rhs).sum())
    return n_ok, best, eq_count


def section_A():
    t = time.time()
    out = {}
    cases = [(7, 3), (7, 6), (11, 5), (11, 10), (13, 3), (13, 4), (13, 6), (13, 12),
             (17, 4), (17, 8), (17, 16), (19, 3), (19, 6), (19, 9), (19, 18),
             (7, 2), (11, 2), (13, 2)]  # k = 2 is (star) itself (control)
    for p, k in cases:
        dk = (p - 1) // k
        Bm = bad_matrix_k(p, k)
        E, DEL, M = all_sets_with0(p, Bm)
        n, best, eq = star_k_check(E, DEL, M, p, dk, "A1_star_k_exhaustive")
        out[f"p{p}_k{k}"] = {"instances(A,e)": n, "max_ratio": best[1], "equality_counts_by_e": eq}
        del E, DEL, M
    note("A1_exhaustive", out)

    # A2: random and structured A at larger p
    rng = random.Random(20260929)
    out2 = {}
    for p in [29, 31, 37, 41, 43, 61, 73, 97]:
        for k in [k for k in range(2, p) if (p - 1) % k == 0 and k <= 12]:
            dk = (p - 1) // k
            Bm = bad_matrix_k(p, k).astype(np.int64)
            L, g = dlog_table(p)
            H = [x for x in range(1, p) if pow(x, dk, p) == 1]
            sets = []
            mmax = p - dk
            for _ in range(150):
                m = rng.randint(1, mmax)
                sets.append(rng.sample(range(p), m))
            for m in range(1, mmax + 1, max(1, mmax // 8)):
                sets.append(list(range(m)))                                   # interval
                sets.append([(3 * i) % p for i in range(m)])                  # AP
                sets.append([pow(g, i, p) for i in range(m)])                 # geometric
            sets.append(H)                                                    # subgroup
            sets.append(sorted(set(H) | set((g * h) % p for h in H)))         # two cosets
            sets.append([0] + H[: max(0, min(len(H), mmax - 1))])
            sets = [s for s in sets if 1 <= len(s) <= mmax]
            n_inst = 0
            best = 0.0
            for A in sets:
                m = len(A)
                e_b = Bm[A, :].sum(axis=0)
                dl = np.zeros(p, dtype=np.int64)
                for a in A:
                    dl[(-a) % p] = 1
                for e in range(0, min((m - 1) // 2, dk // 2) + 1):
                    msk = e_b <= e
                    term = np.where(msk, (e + 1 - e_b) * (2 * m - 3 * e - e_b - 2 * dl), 0)
                    lhs = int(term.sum())
                    rhs = 2 * (e + 1) * (dk - e)
                    chk("A2_star_k_random_structured", (term >= 0).all() and lhs <= rhs,
                        {"p": p, "k": k, "A": A, "e": e, "lhs2": lhs, "rhs2": rhs})
                    n_inst += 1
                    if rhs:
                        best = max(best, lhs / rhs)
            out2[f"p{p}_k{k}"] = {"instances": n_inst, "max_ratio": round(best, 4)}
    note("A2_random_structured", out2)

    # A3: the Hankel determinant itself (degree, leading coefficient, vanishing orders)
    out3 = hankel_direct(rng)
    note("A3_hankel_direct", out3)
    RES["timing"]["A"] = round(time.time() - t, 1)


def pmul(f, g, p):
    r = np.zeros(len(f) + len(g) - 1, dtype=np.int64)
    for i, c in enumerate(f):
        if c:
            r[i:i + len(g)] = (r[i:i + len(g)] + c * g) % p
    return r


def padd(f, g, p):
    n = max(len(f), len(g))
    r = np.zeros(n, dtype=np.int64)
    r[:len(f)] += f
    r[:len(g)] += g
    return r % p


def ptrim(f):
    f = np.array(f, dtype=np.int64)
    nz = np.nonzero(f)[0]
    return f[: nz[-1] + 1] if len(nz) else np.zeros(1, dtype=np.int64)


def binom_mod(n, k, p):
    return math.comb(n, k) % p


def ord_at(f, b, p):
    """multiplicity of the root b of f in F_p[x] (f != 0)."""
    f = [int(c) for c in ptrim(f)]
    cnt = 0
    while True:
        if len(f) == 1:
            return cnt
        # synthetic division by (x - b); coefficients low->high
        n = len(f) - 1
        q = [0] * n
        acc = f[n]
        q[n - 1] = acc
        for i in range(n - 1, 0, -1):
            acc = (f[i] + b * acc) % p
            q[i - 1] = acc
        rem = (f[0] + b * acc) % p
        if rem != 0:
            return cnt
        cnt += 1
        f = q


def hankel_direct(rng):
    """For (p,k,A,e): H_{e+1} = det[u_{i+j}] with u_s = F^{(s)}/(D)_s,
    F = -1 + sum_l c_l (x+a_l)^D, D = d_k+m-1.  Check deg = (e+1)(d_k-e), the leading coefficient
    equals Lambda = det[C(D-i-j, m-1)] mod p, and ord_b H >= T_b for all b with e_b <= e."""
    out = {"cases": 0, "b_checks": 0, "strict_excess_orders": 0}
    plist = [(13, 3), (13, 4), (13, 2), (19, 3), (19, 6), (31, 3), (31, 5), (37, 4), (37, 6), (41, 5)]
    for p, k in plist:
        dk = (p - 1) // k
        Bm = bad_matrix_k(p, k).astype(np.int64)
        H = [x for x in range(1, p) if pow(x, dk, p) == 1]
        for trial in range(16):
            m = rng.randint(2, min(p - dk, 2 * dk + 1, 9))
            A = rng.sample(range(p), m)
            if trial % 4 == 0:
                # configurations with many small e_b: A inside {0} u H, or inside H
                A = ([0] + H)[: min(m, len(H) + 1)]
            elif trial % 4 == 1:
                A = rng.sample(H, min(m, len(H)))
            elif trial % 4 == 2:
                A = rng.sample([0] + H, min(m, len(H) + 1))
            m = len(A)
            if m < 1:
                continue
            D = dk + m - 1
            if D > p - 1:
                continue
            c = []
            for l in range(m):
                pr = 1
                for l2 in range(m):
                    if l2 != l:
                        pr = pr * (A[l] - A[l2]) % p
                c.append(pow(pr, p - 2, p))
            emax = min((m - 1) // 2, dk // 2)
            U = {}
            for s in range(0, 2 * emax + 1):
                N = D - s
                poly = np.zeros(N + 1, dtype=np.int64)
                for l in range(m):
                    a = A[l] % p
                    for i in range(N + 1):
                        poly[i] = (poly[i] + c[l] * binom_mod(N, i, p) * pow(a, N - i, p)) % p
                if s == 0:
                    poly[0] = (poly[0] - 1) % p
                U[s] = poly
            e_b = Bm[A, :].sum(axis=0)
            dl = np.zeros(p, dtype=np.int64)
            for a in A:
                dl[(-a) % p] = 1
            for e in range(0, emax + 1):
                n = e + 1
                Hp = np.zeros(1, dtype=np.int64)
                for perm in itertools.permutations(range(n)):
                    sgn = 1
                    for i in range(n):
                        for j in range(i + 1, n):
                            if perm[i] > perm[j]:
                                sgn = -sgn
                    term = np.ones(1, dtype=np.int64)
                    for i in range(n):
                        term = pmul(term, U[i + perm[i]], p)
                    Hp = padd(Hp, (sgn * term) % p, p)
                Hp = ptrim(Hp)
                degH = len(Hp) - 1
                # Lambda mod p
                Lm = [[math.comb(D - i - j, m - 1) for j in range(n)] for i in range(n)]
                Lam = 0
                for perm in itertools.permutations(range(n)):
                    sgn = 1
                    for i in range(n):
                        for j in range(i + 1, n):
                            if perm[i] > perm[j]:
                                sgn = -sgn
                    pr = sgn
                    for i in range(n):
                        pr *= Lm[i][perm[i]]
                    Lam += pr
                ok_deg = degH == (e + 1) * (dk - e) and int(Hp[-1]) == Lam % p and Lam % p != 0
                chk("A3_hankel_degree_leading", ok_deg,
                    {"p": p, "k": k, "A": A, "e": e, "deg": degH, "Lam_mod_p": Lam % p})
                out["cases"] += 1
                for b in range(p):
                    if e_b[b] <= e:
                        T2 = (e + 1 - e_b[b]) * (2 * m - 3 * e - e_b[b] - 2 * dl[b])  # 2*T_b
                        o = ord_at(Hp, b, p)
                        chk("A3_hankel_orders", 2 * o >= T2,
                            {"p": p, "k": k, "A": A, "e": e, "b": b, "ord": o, "2T": int(T2)})
                        out["b_checks"] += 1
                        if 2 * o > T2 + 1:
                            out["strict_excess_orders"] += 1
    return out


# ============================================================================ Section B
def all_F_01(p):
    """All A containing {0,1}: F[i,b] = sum_{a in A} chi(a+b), INS[i,b] = [b in -A], m."""
    chi = chi_table(p)
    rows = [np.roll(chi, -a) for a in range(p)]  # rows[a][b] = chi(a+b)
    F = (rows[0] + rows[1])[None, :].astype(np.int8)
    INS = np.zeros((1, p), dtype=bool)
    INS[0, 0] = True
    INS[0, p - 1] = True
    M = np.array([2], dtype=np.int16)
    for a in range(2, p):
        F = np.concatenate([F, F + rows[a][None, :].astype(np.int8)], axis=0)
        I2 = INS.copy()
        I2[:, (-a) % p] = True
        INS = np.concatenate([INS, I2], axis=0)
        M = np.concatenate([M, M + 1])
    return F, INS, M, chi


def extremal_S(p):
    """maxS[(m,n,r)], minS[(m,n,r)] over all A (up to affine maps) and all B with |B|=n,
    |B cap -A| = r; exhaustive.  Also per-(m,E,r) maxima of sum_B W_E (doubled)."""
    F, INS, M, chi = all_F_01(p)
    maxS, minS = {}, {}
    W = {}
    d = (p - 1) // 2
    # m = 1 : A = {0}
    vals_in = np.array([0])
    vals_out = np.array([int(chi[b]) for b in range(1, p)])
    groups = [(1, np.zeros((1, 1), dtype=np.int64), vals_out[None, :].astype(np.int64))]
    for m in range(2, p + 1):
        idx = np.nonzero(M == m)[0]
        if len(idx) == 0:
            continue
        Fm = F[idx].astype(np.int64)
        Im = INS[idx]
        ins = Fm[Im].reshape(len(idx), m)
        outs = Fm[~Im].reshape(len(idx), p - m)
        groups.append((m, ins, outs))
    npairs = 0
    bic = {}
    for m, ins, outs in groups:
        nA = ins.shape[0]
        # biclique data: N_0 (zeros allowed at b in -A) and N_0 without b in -A, for A and nu*A
        n0 = np.maximum(((ins == m - 1).sum(axis=1) + (outs == m).sum(axis=1)),
                        ((ins == -(m - 1)).sum(axis=1) + (outs == -m).sum(axis=1)))
        n0s = np.maximum((outs == m).sum(axis=1), (outs == -m).sum(axis=1))
        bic[m] = (int(n0.max()), int(n0s.max()))
        it = -np.sort(-ins, axis=1)
        ot = -np.sort(-outs, axis=1)
        ib = np.sort(ins, axis=1)
        ob = np.sort(outs, axis=1)
        z = np.zeros((nA, 1), dtype=np.int64)
        CIt = np.concatenate([z, np.cumsum(it, axis=1)], axis=1)
        COt = np.concatenate([z, np.cumsum(ot, axis=1)], axis=1)
        CIb = np.concatenate([z, np.cumsum(ib, axis=1)], axis=1)
        COb = np.concatenate([z, np.cumsum(ob, axis=1)], axis=1)
        for r in range(0, m + 1):
            for no in range(0, p - m + 1):
                n = r + no
                if n == 0:
                    continue
                maxS[(m, n, r)] = int((CIt[:, r] + COt[:, no]).max())
                minS[(m, n, r)] = int((CIb[:, r] + COb[:, no]).min())
                npairs += nA * math.comb(m, r) * math.comb(p - m, no)
        # Lemma 1.1 data (only m <= (p+1)/2): e-values e = (m - delta - F)/2, both signs of F
        if m <= (p + 1) // 2:
            for sgn in (1, -1):
                ein = (m - 1 - sgn * ins) // 2
                eout = (m - sgn * outs) // 2
                for E in range(1, (m + 1) // 2 + 1):
                    Win = np.where(ein < E, (E - ein) * (2 * m + 3 - 3 * E - ein), 0)
                    Wout = np.where(eout < E, (E - eout) * (2 * m + 3 - 3 * E - eout), 0)
                    so = Wout.sum(axis=1)
                    Wi = -np.sort(-Win, axis=1)
                    CW = np.concatenate([z, np.cumsum(Wi, axis=1)], axis=1)
                    for r in range(0, m + 1):
                        key = (m, E, r)
                        W[key] = max(W.get(key, -1), int((so + CW[:, r]).max()))
    return maxS, minS, W, npairs, bic


def thm14_ok(S, m, n, r, d):
    """S <= u(max(w',1/4)) m n with w' = (d+r)/(mn) < 1 (exact integers)."""
    Q = m * n
    P = d + r
    if 4 * P >= Q:
        L = 2 * S + P  # (2S + P)/Q vs sqrt(12PQ-3P^2)/Q
        return L <= 0 or L * L <= 12 * P * Q - 3 * P * P
    L = 8 * S + Q
    return L <= 0 or L * L <= 45 * Q * Q


def cor12_ok(S, m, n, r, d):
    """Corollary 1.2 for the pair data: with N = (mn - r - S)/2 (so 1-2tbar = (mn-2N)/mn),
    1-2tbar <= 1 - 2(1-w')/m  and  1-2tbar <= h_m(E,w') for every E <= (m+1)/2."""
    N2 = m * n - r - S  # = 2N
    ok = True
    # HP line: (mn - 2N)/(mn) <= 1 - 2(1 - (d+r)/(mn))/m   <=>  -2N m <= -2(mn - d - r)
    ok &= (-N2 * m) <= -2 * (m * n - d - r)
    for E in range(1, (m + 1) // 2 + 1):
        L = 2 * m * n - N2 + n * (3 - 2 * E)
        if L > 0:
            ok &= L * L <= n * ((2 * m - 4 * E + 3) ** 2 * n + 8 * E * (d + r))
    return ok


def w1_of(p, kap):
    m1 = math.isqrt(int(math.floor((Fr(1, 2) + kap) * p)))
    while Fr(m1 * m1) < (Fr(1, 2) + kap) * p:
        m1 += 1
    while m1 > 0 and Fr((m1 - 1) ** 2) >= (Fr(1, 2) + kap) * p:
        m1 -= 1
    return m1, Fr(p - 1 + 2 * m1) / ((1 + 2 * kap) * p)


def section_B():
    t = time.time()
    out = {}
    for p in [7, 11, 13, 17, 19, 23]:
        d = (p - 1) // 2
        maxS, minS, W, npairs, bic = extremal_S(p)
        st = {"pairs_(A,B)_up_to_affine": npairs, "mnr_triples": len(maxS)}
        # [Sh] section 4 G1 data: M_m(p) = max_{|A|=m} N_0(A); HP m*M_m <= d + m
        st["M_m(p)_zeros_allowed"] = [bic[m][0] for m in sorted(bic)]
        st["M_m(p)_no_zeros"] = [bic[m][1] for m in sorted(bic)]
        for m in bic:
            if m <= (p + 1) // 2:
                chk("B0_HP_biclique", m * bic[m][0] <= d + m, {"p": p, "m": m})
        best_bic = max(((m * bic[m][0], m, bic[m][0]) for m in bic if m >= 3 and bic[m][0] >= 3), default=(0, 0, 0))
        best_bic2 = max(((m * bic[m][1], m, bic[m][1]) for m in bic if m >= 3 and bic[m][1] >= 3), default=(0, 0, 0))
        st["largest_biclique_min>=3_(mn/p,m,n)_zeros_allowed"] = [best_bic[0] / p, best_bic[1], best_bic[2]]
        st["largest_biclique_min>=3_(mn/p,m,n)_no_zeros"] = [best_bic2[0] / p, best_bic2[1], best_bic2[2]]
        # validation of the affine normalisation by brute force at p = 7, 11
        if p in (7, 11):
            chi = chi_table(p)
            allB = np.array([[(bm >> b) & 1 for b in range(p)] for bm in range(1 << p)], dtype=np.int64)
            nB = allB.sum(axis=1)
            brute = {}
            for am in range(1, 1 << p):
                A = [a for a in range(p) if (am >> a) & 1]
                Fv = np.zeros(p, dtype=np.int64)
                dl = np.zeros(p, dtype=bool)
                for a in A:
                    Fv += np.roll(chi, -a)
                    dl[(-a) % p] = True
                Sv = allB @ Fv
                rv = allB[:, dl].sum(axis=1)
                key = nB * 100 + rv
                for kk, v in zip(key.tolist(), np.abs(Sv).tolist()):
                    t3 = (len(A), kk // 100, kk % 100)
                    if brute.get(t3, -1) < v:
                        brute[t3] = v
            agree = all(brute[t3] == max(maxS[t3], -minS[t3]) for t3 in brute if t3[1] > 0)
            chk("B0_normalisation_bruteforce", agree and all(t3 in brute for t3 in maxS), {"p": p})
        # Lemma 1.1
        nl = 0
        for (m, E, r), v in W.items():
            chk("B1_lemma11_mean_form", v <= 2 * E * (d - E + 1 + r), {"p": p, "m": m, "E": E, "r": r, "lhs2": v})
            nl += 1
        st["lemma11_(m,E,r)"] = nl
        worst14 = 0.0
        n14 = 0
        for (m, n, r) in maxS:
            if m > (p + 1) // 2 or r > min(m, n):
                continue
            for S in (maxS[(m, n, r)], -minS[(m, n, r)]):
                chk("B2_cor12_all_E_and_HP", cor12_ok(S, m, n, r, d),
                    {"p": p, "m": m, "n": n, "r": r, "S": S})
                if (d + r) < m * n:
                    ok = thm14_ok(S, m, n, r, d)
                    chk("B3_thm14_pair_form", ok, {"p": p, "m": m, "n": n, "r": r, "S": S})
                    n14 += 1
                    w = max(Fr(d + r, m * n), Fr(1, 4))
                    worst14 = max(worst14, S / (m * n) / u_float(float(w)))
        st["thm14_checks"] = n14
        st["thm14_max_S_over_u_mn"] = round(worst14, 4)
        # Theorem 1.5
        maxabs = {}
        for (m, n, r), v in maxS.items():
            maxabs[(m, n)] = max(maxabs.get((m, n), -10 ** 9), v, -minS[(m, n, r)])
        n15 = 0
        worst15 = (0.0, None)
        nontriv = 0
        for (m, n), S in maxabs.items():
            kmax = min(Fr(3, 2), Fr(m * n, p) - Fr(1, 2))
            if kmax <= 0:
                continue
            cands = {kmax}
            j = 1
            while True:
                kj = Fr(j * j, p) - Fr(1, 2)
                if kj > kmax:
                    break
                if kj > 0:
                    cands.add(kj)
                j += 1
            for kap in cands:
                m1, w1 = w1_of(p, kap)
                wb = min(Fr(1), w1)
                ok = le_u(Fr(S, m * n), wb)
                chk("B4_thm15_exhaustive", ok, {"p": p, "m": m, "n": n, "kappa": str(kap), "S": S})
                n15 += 1
                if wb < 1:
                    nontriv += 1
                    rr = S / (m * n) / u_float(float(wb))
                    if rr > worst15[0]:
                        worst15 = (rr, {"m": m, "n": n, "kappa": float(kap), "S": S, "u(w1)": u_float(float(wb))})
                # second inequality of Thm 1.5 (float, margin reported)
                kf = float(kap)
                rhs = u_float(1 / (1 + 2 * kf)) + 2.14 * (math.sqrt(2 * p) + 1) / p
                chk("B5_thm15_second_ineq_float", u_float(float(wb)) <= rhs + 1e-12,
                    {"p": p, "kappa": kf})
        # [Sh] section 3 exhaustive data: max bias over |A||B| >= (1/2+kappa)p
        mb = {}
        for kap in [Fr(0), Fr(1, 20), Fr(1, 10), Fr(3, 20), Fr(1, 5), Fr(1, 4), Fr(3, 10), Fr(7, 20), Fr(1, 2), Fr(1)]:
            best = (Fr(-1), None)
            for (m, n), S in maxabs.items():
                if m * n >= (Fr(1, 2) + kap) * p:
                    v = Fr(S, m * n)
                    if v > best[0]:
                        best = (v, (m, n))
            mb[str(kap)] = [float(best[0]), best[1]]
        st["max_bias_by_kappa_(value,(m,n))"] = mb
        st["thm15_checks"] = n15
        st["thm15_nontrivial_(w1<1)"] = nontriv
        st["thm15_worst_ratio"] = worst15
        out[f"p{p}"] = st
        del maxS, minS, W
    note("B_exhaustive", out)
    out6 = section_B6()
    note("B6_structured_larger_p", out6)
    RES["timing"]["B"] = round(time.time() - t, 1)


def greedy_clique_sum(p, chi, rng, size):
    """greedy A with chi(a+a') = 1 for a != a' in A ... used with B = A (sum form biclique)."""
    A = []
    cand = list(range(p))
    rng.shuffle(cand)
    for x in cand:
        if all(chi[(x + a) % p] == 1 for a in A) and chi[(2 * x) % p] == 1:
            A.append(x)
        if len(A) >= size:
            break
    return A


def section_B6():
    rng = random.Random(29)
    out = {}
    for p in [29, 37, 41, 53, 61, 101, 197]:
        chi = chi_table(p)
        d = (p - 1) // 2
        sets = []
        for m in sorted(set([2, 3, 4, 5, 6, 8, 10, int(math.sqrt(p)), int(math.sqrt(2 * p)), d // 2, d, d + 1])):
            if m < 1 or m > p:
                continue
            sets.append(("interval", list(range(m))))
            sets.append(("ap7", [(7 * i) % p for i in range(m)]))
            sets.append(("random", rng.sample(range(p), m)))
            qr = [x for x in range(1, p) if chi[x] == 1]
            sets.append(("QRsub", qr[:m] if m <= len(qr) else qr))
        for size in [3, 4, 5, 6, 8]:
            A = greedy_clique_sum(p, chi, rng, size)
            if len(A) >= 2:
                sets.append(("sumclique", A))
        worst14, n14, n15, worst15 = 0.0, 0, 0, 0.0
        for kind, A in sets:
            A = sorted(set(A))
            m = len(A)
            Fv = np.zeros(p, dtype=np.int64)
            dl = np.zeros(p, dtype=bool)
            for a in A:
                Fv += np.roll(chi, -a)
                dl[(-a) % p] = True
            ins, outs = Fv[dl], Fv[~dl]
            for sgn in (1, -1):
                it = np.sort(sgn * ins)[::-1]
                ot = np.sort(sgn * outs)[::-1]
                CI = np.concatenate([[0], np.cumsum(it)])
                CO = np.concatenate([[0], np.cumsum(ot)])
                for r in range(0, m + 1):
                    for no in range(0, p - m + 1):
                        n = r + no
                        if n == 0:
                            continue
                        S = int(CI[r] + CO[no])
                        if m <= (p + 1) // 2:
                            chk("B6_cor12_structured", cor12_ok(S, m, n, r, d), {"p": p, "A": kind, "m": m, "n": n, "r": r})
                            if d + r < m * n:
                                chk("B6_thm14_structured", thm14_ok(S, m, n, r, d), {"p": p, "A": kind, "m": m, "n": n, "r": r, "S": S})
                                n14 += 1
                                w = max((d + r) / (m * n), 0.25)
                                worst14 = max(worst14, S / (m * n) / u_float(w))
                        kmax = min(Fr(3, 2), Fr(m * n, p) - Fr(1, 2))
                        if kmax > 0 and p >= 11:
                            m1, w1 = w1_of(p, kmax)
                            chk("B6_thm15_structured_at_kappa_max", le_u(Fr(S, m * n), min(Fr(1), w1)),
                                {"p": p, "A": kind, "m": m, "n": n})
                            n15 += 1
        out[f"p{p}"] = {"sets": len(sets), "thm14_checks": n14, "thm14_max_ratio": round(worst14, 4), "thm15_checks": n15}
    return out


# ============================================================================ Section C
def h_opt(m, E, w):
    a = 1 - (2 * E - 1.5) / m
    beta = 1.5 / m
    return -a - beta + 2 * math.sqrt(a * a - a * w + (1 + beta) * w)


def h_opt_hi(m, E, w, bits=120):
    """rational upper bound for h_m(E,w)."""
    a = 1 - Fr(4 * E - 3, 2 * m)
    beta = Fr(3, 2 * m)
    return -a - beta + 2 * isqrt_frac_hi(a * a - a * w + (1 + beta) * w, bits)


def best_opt_hi(m, w, bits=120):
    w = Fr(w)
    best = 1 - 2 * (1 - w) / m
    for E in range(1, (m + 1) // 2 + 1):
        best = min(best, h_opt_hi(m, E, w, bits))
    return best


def section_C():
    t = time.time()
    out = {}
    # C1 analytic constants, exact
    eps = Fr(9, 20)
    c = (1 - eps) * (3 + eps) / 4
    val = (1 + Fr(2, 3) / isqrt_frac_lo(c)) * (3 - 2 * eps) / (2 * (1 - eps))
    chk("C1_eps_MJ_at_9/20_le_3.757", val <= Fr(3757, 1000), float(val))
    out["eps_MJ(9/20)_upper"] = float(val)
    eps = Fr(3, 4)
    c = (1 - eps) * (3 + eps) / 4
    first = 1 + Fr(2, 3) / isqrt_frac_lo(c)
    sec = max((3 - 2 * e) / (2 * e * (1 - e)) for e in (Fr(9, 20), Fr(3, 4)))
    chk("C1_MJ_le_10.085", first * sec <= Fr(10085, 1000), float(first * sec))
    out["MJ_upper_on_[1/4,11/20]"] = float(first * sec)
    # v >= 0.40 at w = 1/4 (v increasing in w)
    v14 = (Fr(3, 4) + isqrt_frac_lo(Fr(3) - Fr(3, 16))) / 6
    chk("C1_v_ge_0.40", v14 >= Fr(2, 5), float(v14))
    # symbolic facts
    import sympy as sp
    W, EPS, A_ = sp.symbols("w epsilon a", positive=True)
    u = (sp.sqrt(12 * W - 3 * W ** 2) - W) / 2
    chk("C1_u_second_derivative", sp.simplify(sp.diff(u, W, 2) + 2 * sp.sqrt(3) / (W * (4 - W)) ** sp.Rational(3, 2)) == 0)
    chk("C1_u(1)=1", sp.simplify(u.subs(W, 1) - 1) == 0)
    ser = sp.series(1 - u.subs(W, 1 / (1 + 2 * EPS)), EPS, 0, 4).removeO()
    chk("C1_eta_series", sp.simplify(ser - (sp.Rational(4, 3) * EPS ** 2 - sp.Rational(40, 9) * EPS ** 3)) == 0, str(ser))
    LAM = sp.symbols("lambda", positive=True)
    th = sp.series((1 - u.subs(W, 1 / (1 + LAM))) / 2, LAM, 0, 3).removeO()
    chk("C1_theta_lambda_series", sp.simplify(th - LAM ** 2 / 6) == 0, str(th))
    # (F1): 3v^2-3vw+w^2-w = 0 and h0(v) = u
    v = (3 * W + sp.sqrt(12 * W - 3 * W ** 2)) / 6
    chk("C1_F1", sp.simplify(sp.expand(3 * v ** 2 - 3 * v * W + W ** 2 - W)) == 0)
    Pv = v ** 2 - v * W + W
    chk("C1_F1b", sp.simplify(sp.expand(Pv - (2 * v - W) ** 2)) == 0)
    # (F5) identity v - w = 2w(1-w)/(sqrt(12w-3w^2)+3w)
    chk("C1_F5", sp.simplify(v - W - 2 * W * (1 - W) / (sp.sqrt(12 * W - 3 * W ** 2) + 3 * W)) == 0)
    # derivative claims: eps*MJ increasing in eps on (0, 9/20]; (3-2e)/(e(1-e)) critical point
    g2 = (3 - 2 * EPS) / (EPS * (1 - EPS))
    crit = sp.solve(sp.diff(g2, EPS), EPS)
    out["crit_points_(3-2e)/(e(1-e))"] = [str(x) for x in crit]
    chk("C1_crit_point", any(sp.simplify(x - (6 - sp.sqrt(12)) / 4) == 0 for x in crit))
    # (F6): 1-u <= eps^2/2 for eps <= sqrt3-1 : check exactly on a rational grid + symbolic
    for i in range(1, 201):
        e = Fr(i, 200) * Fr(732, 1000)
        # real test: 1 - u(w) <= e^2/2  <=>  u(w) >= 1 - e^2/2
        wv = 1 - e
        x = 1 - e * e / 2
        # u(w) >= x  <=>  2x + w <= sqrt(12w-3w^2)
        L = 2 * x + wv
        chk("C1_F6_exact", L <= 0 or L * L <= 12 * wv - 3 * wv * wv, float(e))
    # C2: independent rational-interval check of the finite range m <= 10, w in [1/4, 11/20]
    boxes = 0
    min_margin = Fr(10)
    for m in range(1, 11):
        stack = [(Fr(1, 4) + Fr(i, 40) * Fr(3, 10), Fr(1, 4) + Fr(i + 1, 40) * Fr(3, 10)) for i in range(40)]
        while stack:
            w1, w2 = stack.pop()
            hi = best_opt_hi(m, w2)
            lo = u_lo(w1, 120)
            if hi <= lo:
                boxes += 1
                min_margin = min(min_margin, lo - hi)
            else:
                if w2 - w1 < Fr(1, 10 ** 6):
                    chk("C2_interval_m_le_10", False, {"m": m, "w1": float(w1), "w2": float(w2)})
                    continue
                mid = (w1 + w2) / 2
                stack += [(w1, mid), (mid, w2)]
        chk("C2_interval_m_le_10", True)
    out["C2_boxes"] = boxes
    out["C2_min_margin"] = float(min_margin)
    # C2b: the same interval method, as an extra cross-check, for 11 <= m <= 60 on [1/4, 0.9]
    boxes2 = 0
    for m in range(11, 61):
        stack = [(Fr(1, 4) + Fr(i, 60) * Fr(13, 20), Fr(1, 4) + Fr(i + 1, 60) * Fr(13, 20)) for i in range(60)]
        while stack:
            w1, w2 = stack.pop()
            if best_opt_hi(m, w2) <= u_lo(w1, 120):
                boxes2 += 1
            else:
                if w2 - w1 < Fr(1, 10 ** 7):
                    chk("C2b_interval_extended", False, {"m": m, "w1": float(w1)})
                    continue
                mid = (w1 + w2) / 2
                stack += [(w1, mid), (mid, w2)]
        chk("C2b_interval_extended", True)
    out["C2b_boxes_m11_60_w_to_0.9"] = boxes2
    # C3: exact point checks of Theorem 1.3 for large m and w close to 1 (rigorous enclosures)
    rng = random.Random(7)
    npt = 0
    for _ in range(400):
        m = rng.choice([rng.randint(1, 50), rng.randint(50, 2000), rng.randint(2000, 20000)])
        w = Fr(rng.randint(250000, 999999), 10 ** 6)
        # evaluate only the best few E near the optimum to keep it fast
        a_target = u_float(float(w))
        vv = (3 * float(w) + math.sqrt(12 * float(w) - 3 * float(w) ** 2)) / 6
        Ec = int(round((1 - vv) * m / 2 + 0.75))
        best = 1 - 2 * (1 - w) / m
        for E in range(max(1, Ec - 3), min((m + 1) // 2, Ec + 3) + 1):
            best = min(best, h_opt_hi(m, E, w, 160))
        ok = best <= u_lo(w, 160)
        chk("C3_thm13_point_exact", ok, {"m": m, "w": str(w)})
        npt += 1
    out["C3_points"] = npt
    # C4 float scan and the size of the gap
    worst = -1.0
    for m in list(range(1, 400)) + list(range(400, 3001, 50)):
        for i in range(0, 151):
            w = 0.25 + i * (0.999 - 0.25) / 150
            b = 1 - 2 * (1 - w) / m
            for E in range(1, (m + 1) // 2 + 1):
                b = min(b, h_opt(m, E, w))
            worst = max(worst, b - u_float(w))
    chk("C4_float_scan", worst < 0, worst)
    out["C4_max(min_option - u)"] = worst
    note("C_thm13", out)
    RES["timing"]["C"] = round(time.time() - t, 1)


# ============================================================================ Section D
def star_constraints_hold(n, m, J, d):
    """all-e_b = J profile with n points, r = 0: (star)_e for all e <= (m-1)/2 (doubled)."""
    for e in range(J, (m - 1) // 2 + 1):
        if n * (e + 1 - J) * (2 * m - 3 * e - J) > 2 * (e + 1) * (d - e):
            return False, e
    return True, None


def section_D():
    t = time.time()
    out = {"profiles": []}
    # D1 profiles of Prop 2.1(ii),(iii) : exact
    for (w, m, d) in [(Fr(9, 10), 4000, 10 ** 10), (Fr(2, 3), 2000, 10 ** 10), (Fr(1, 2), 1000, 10 ** 9),
                      (Fr(1, 4), 400, 10 ** 8), (Fr(4, 5), 800, 10 ** 9), (Fr(19, 20), 3000, 10 ** 11),
                      (Fr(1, 3), 200, 10 ** 7), (Fr(3, 5), 60, 10 ** 6)]:
        ustar_lo = u_lo(w)
        tstar_hi = (1 - ustar_lo) / 2          # upper bound for t*
        delta = 5 * (Fr(3, 2 * m) + Fr(m, 2 * d))
        J = math.ceil((tstar_hi + delta) * m)
        n = int(d / (m * w))
        ok, bad = star_constraints_hold(n, m, J, d)
        # hypothesis of (ii): (t - t*)(1/2 - 1/(2m) - t) >= 3/(2m) + m/(2d), with t* <= tstar_hi
        tt = Fr(J, m)
        X = (tt - tstar_hi) * (Fr(1, 2) - Fr(1, 2 * m) - tt)
        hyp = X >= Fr(3, 2 * m) + Fr(m, 2 * d)
        chk("D1_profile_feasible_exact", ok, {"w": str(w), "m": m, "d": d, "bad_e": bad})
        chk("D1_profile_hypothesis_ii", hyp or m < 40, {"w": str(w), "m": m})
        # also: maximal J' <= J still feasible (the profile family is not an artefact of slack)
        Jmin = J
        while Jmin > 0 and star_constraints_hold(n, m, Jmin - 1, d)[0]:
            Jmin -= 1
        out["profiles"].append({"w": float(w), "m": m, "d": d, "n": n, "J": J,
                                "profile_bias": float(1 - Fr(2 * J, m)), "u(w)": u_float(float(w)),
                                "smallest_feasible_J": Jmin, "bias_at_Jmin": float(1 - Fr(2 * Jmin, m)),
                                "d/(mn)": float(Fr(d, m * n)), "hyp_ii": bool(hyp)})
    # D2 float LPs (d -> infinity, r = 0 normalisation): value vs u(w)
    try:
        from scipy.optimize import linprog
        lp = {}
        for w in [0.5, 0.8, 0.9, 0.95]:
            for m in [20, 50, 100, 200, 400]:
                J = m + 1
                c = -np.array([1 - 2 * j / m for j in range(J)])
                rows, rhs = [], []
                for e in range(0, (m - 1) // 2 + 1):
                    row = [((e + 1 - j) * (m - (3 * e + j) / 2) if j <= e else 0.0) for j in range(J)]
                    rows.append(row)
                    rhs.append((e + 1) * m * w)
                res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs), A_eq=np.ones((1, J)), b_eq=[1.0],
                              bounds=[(0, None)] * J, method="highs")
                V = -res.fun
                chk("D2_lp_le_u_float", V <= u_float(w) + 1e-9, {"w": w, "m": m, "V": V})
                lp[f"w{w}_m{m}"] = {"V": V, "u": u_float(w), "m*(u-V)": m * (u_float(w) - V)}
        out["D2_lp"] = lp
        # D3 finite LP with r variables at p in {101, 1009}: V/(mn) <= u(max(w'',1/4))
        cnt = 0
        worst = 0.0
        for p in [101, 1009]:
            d = (p - 1) // 2
            for kap in [Fr(1, 20), Fr(1, 4), Fr(1, 2), Fr(1)]:
                for m in range(2, 25):
                    n = math.ceil((Fr(1, 2) + kap) * p / m)
                    if m > (p + 1) // 2:
                        continue
                    J = m + 1
                    c = -np.array([1 - 2 * j / m for j in range(J)] + [(m - 1 - 2 * j) / m for j in range(J)])
                    rows, rhs = [], []
                    for e in range(0, (m - 1) // 2 + 1):
                        row = [((e + 1 - j) * (m - (3 * e + j) / 2) if j <= e else 0.0) for j in range(J)]
                        row += [((e + 1 - j) * (m - (3 * e + j) / 2 - 1) if j <= e else 0.0) for j in range(J)]
                        rows.append(row)
                        rhs.append((e + 1) * (d - e))
                    rows.append([0.0] * J + [1.0] * J)
                    rhs.append(min(m, n))
                    res = linprog(c, A_ub=np.array(rows), b_ub=np.array(rhs),
                                  A_eq=np.ones((1, 2 * J)), b_eq=[float(n)],
                                  bounds=[(0, None)] * (2 * J), method="highs")
                    V = -res.fun / n
                    w2 = (d + min(m, n)) / (m * n)
                    if w2 < 1:
                        ub = u_float(max(w2, 0.25))
                        chk("D3_lp_finite_le_u", V <= ub + 1e-9, {"p": p, "m": m, "n": n, "V": V})
                        worst = max(worst, V / ub)
                        cnt += 1
        out["D3_count"] = cnt
        out["D3_max_V_over_u"] = worst
    except ImportError:
        out["scipy"] = "missing"
    note("D_lp", out)
    RES["timing"]["D"] = round(time.time() - t, 1)


# ============================================================================ Section E
def beta_exact(m, q):
    """mean of 1-2j/m over the lowest mass q of Binomial(m,1/2) (q rational in (0,1])."""
    rem = Fr(q)
    tot = Fr(0)
    for j in range(m + 1):
        pj = Fr(math.comb(m, j), 2 ** m)
        take = min(pj, rem)
        tot += take * (1 - Fr(2 * j, m))
        rem -= take
        if rem == 0:
            break
    return tot / Fr(q)


def section_E():
    t = time.time()
    out = {}
    # beta_2 = w for w >= 1/3 ; beta_5 = 3/5 + w/8 for w in [8/15, 1]
    for i in range(0, 101):
        w = Fr(1, 3) + Fr(i, 100) * Fr(2, 3)
        chk("E1_beta2_eq_w", beta_exact(2, 1 / (4 * w)) == w)
        w5 = Fr(8, 15) + Fr(i, 100) * Fr(7, 15)
        chk("E1_beta5_formula", beta_exact(5, 1 / (10 * w5)) == Fr(3, 5) + w5 / 8)
        w1 = Fr(1, 2) + Fr(i, 200)
        chk("E1_beta1", beta_exact(1, 1 / (2 * w1)) == 2 * w1 - 1)
    chk("E1_crossover_24/35", Fr(3, 5) + Fr(24, 35) / 8 == Fr(24, 35) and (1 / Fr(24, 35) - 1) / 2 == Fr(11, 48))
    # beta_m <= w for 3 <= m < 200, w in [24/35, 1]: piecewise linear in w -> check breakpoints
    nbp = 0
    worst = Fr(-1)
    for m in range(3, 200):
        cum = Fr(0)
        bps = {Fr(24, 35), Fr(1)}
        for j in range(m + 1):
            cum += Fr(math.comb(m, j), 2 ** m)
            if cum > 0:
                wb = 1 / (2 * m * cum)
                if Fr(24, 35) <= wb <= 1:
                    bps.add(wb)
        for wb in bps:
            q = 1 / (2 * m * wb)
            if q > 1:
                continue
            val = beta_exact(m, q) - wb
            ok = val <= 0 if not (m == 5 and wb == Fr(24, 35)) else val == 0
            chk("E2_beta_m_le_w_breakpoints", val <= 0, {"m": m, "w": str(wb), "diff": float(val)})
            worst = max(worst, val)
            nbp += 1
    out["E2_breakpoints"] = nbp
    out["E2_max(beta_m - w)"] = float(worst)
    # Hoeffding tail for m >= 200:   1/2 + 2m e^{-m/8} < 24/35
    chk("E2_hoeffding_m200", 0.5 + 2 * 200 * math.exp(-25) < 24 / 35)
    chk("E2_hoeffding_monotone", all((m + 1) * math.exp(-(m + 1) / 8) <= m * math.exp(-m / 8) for m in range(200, 2000)))
    # E3: sup_m beta_m at the table's kappas (m up to 400 exact; tail by Hoeffding)
    tab = {}
    for kap in [Fr(1, 100), Fr(1, 20), Fr(1, 10), Fr(1, 4), Fr(1, 2), Fr(1), Fr(3, 2)]:
        w = 1 / (1 + 2 * kap)
        best, arg = Fr(-1), None
        for m in range(1, 401):
            q = 1 / (2 * m * w)
            if q > 1:
                continue
            b = beta_exact(m, q)
            if b > best:
                best, arg = b, m
        # Hoeffding: for m > 400, beta_m <= x0 + (2mw) exp(-m x0^2/2) with x0 = 0.3
        tail = 0.3 + 2 * 400 * math.exp(-400 * 0.09 / 2)
        chk("E3_tail_below_sup", tail < float(best), {"kappa": str(kap), "tail": tail, "best": float(best)})
        tab[str(kap)] = {"sup_beta": float(best), "argmax_m": arg, "saving_upper": float(1 - best),
                         "2k/(1+2k)": float(2 * kap / (1 + 2 * kap)), "eta_star": 1 - u_float(float(w)),
                         "exact_saving_upper": str(1 - best)}
    out["E3_table"] = tab
    chk("E3_kappa_half_29/80", Fr(tab["1/2"]["exact_saving_upper"]) == Fr(29, 80))
    # E4: realisation at a large prime (Weil counts): random A, B = smallest e_b
    p = 100003
    chi = chi_table(p)
    rng = random.Random(3)
    real = {}
    for m, kap in [(2, Fr(1, 10)), (5, Fr(1, 2)), (5, Fr(1, 4)), (7, Fr(1)), (6, Fr(3, 2))]:
        w = 1 / (1 + 2 * kap)
        n = math.ceil((Fr(1, 2) + kap) * p / m)
        A = rng.sample(range(p), m)
        Fv = np.zeros(p, dtype=np.int64)
        for a in A:
            Fv += np.roll(chi, -a)
        top = np.sort(Fv)[::-1][:n]
        S = int(top.sum())
        real[f"m{m}_k{kap}"] = {"S/(mn)": S / (m * n), "beta_m": float(beta_exact(m, 1 / (2 * m * w)))}
        chk("E4_realisation_close", abs(S / (m * n) - float(beta_exact(m, Fr(n, p)))) < 0.02, real[f"m{m}_k{kap}"])
    out["E4_realised"] = real
    note("E_upper_bounds", out)
    RES["timing"]["E"] = round(time.time() - t, 1)


# ============================================================================ Section F
def section_F():
    t = time.time()
    out = {}
    rng = random.Random(11)
    # F1 Lemma 3.1 exact counts at large p
    worst = 0.0
    for p in [10007, 100003]:
        chi = chi_table(p)
        for m in range(1, 11):
            for trial in range(3):
                if trial == 0:
                    A = list(range(m))
                elif trial == 1:
                    A = [(i * i) % p for i in range(1, m + 1)]
                else:
                    A = rng.sample(range(p), m)
                A = sorted(set(A))
                mm = len(A)
                e = np.zeros(p, dtype=np.int64)
                for a in A:
                    e += (np.roll(chi, -a) == -1)
                cnt = np.bincount(e, minlength=mm + 1)
                for j in range(mm + 1):
                    main = math.comb(mm, j) * p / 2 ** mm
                    bound = math.comb(mm, j) * (mm / 2) * (math.sqrt(p) + 1) + 2 * mm
                    dev = abs(cnt[j] - main)
                    # exact version: |2^m N_j - C p| <= 2^m bound  (bound irrational: compare squares)
                    lhs = abs(2 ** mm * int(cnt[j]) - math.comb(mm, j) * p)
                    # 2^m*bound = C*m*2^(m-1)*(sqrt p + 1) + 2^(m+1)*m ; check lhs - rational part <= C m 2^(m-1) sqrt p
                    R = lhs - math.comb(mm, j) * mm * 2 ** (mm - 1) - 2 ** (mm + 1) * mm
                    K = math.comb(mm, j) * mm * 2 ** (mm - 1)
                    ok = R <= 0 or R * R <= K * K * p
                    chk("F1_lemma31_exact", ok, {"p": p, "A": A, "j": j})
                    worst = max(worst, dev / bound)
    out["F1_max_dev_over_bound"] = worst
    # F2 Lemma 4.2 fourth moment, exact integers
    worst4 = 0.0
    for p in [1009, 10007]:
        chi = chi_table(p)
        for m in [3, 5, 10, 20, 40, 80]:
            for kind in ["interval", "squares", "random", "ap"]:
                if kind == "interval":
                    A = list(range(m))
                elif kind == "squares":
                    A = sorted(set((i * i) % p for i in range(1, m + 1)))
                elif kind == "ap":
                    A = [(7 * i) % p for i in range(m)]
                else:
                    A = rng.sample(range(p), m)
                mm = len(A)
                Fv = np.zeros(p, dtype=np.int64)
                for a in A:
                    Fv += np.roll(chi, -a)
                M4 = int((Fv ** 4).sum())
                # M4 <= 3 m^2 p + 3 m^4 sqrt p  (exact: M4 - 3m^2p <= 3m^4 sqrt p)
                R = M4 - 3 * mm * mm * p
                ok = R <= 0 or R * R <= 9 * mm ** 8 * p
                chk("F2_fourth_moment_exact", ok, {"p": p, "A": kind, "m": mm})
                worst4 = max(worst4, M4 / (3 * mm * mm * p + 3 * mm ** 4 * math.sqrt(p)))
    out["F2_max_ratio"] = worst4
    # F3 Theorem 4.3 at primes where the hypotheses hold
    res43 = []
    for k, kap in [(2, Fr(1)), (3, Fr(1)), (4, Fr(1)), (3, Fr(1, 2))]:
        tau = Fr(k, 2 ** k) + kap
        M = math.ceil(12 / tau)
        pmin = max(25, math.ceil((2 * M * M / kap) ** 2))
        p = next_prime(pmin)
        chi = chi_table(p)
        worst_a, worst_b = 0.0, 0.0
        na = nb = 0
        for m in range(k, M + 1):
            for kind in ["interval", "squares", "random", "random", "geom"]:
                if kind == "interval":
                    A = list(range(m))
                elif kind == "squares":
                    A = sorted(set((i * i) % p for i in range(1, m + 1)))
                elif kind == "geom":
                    A = sorted(set(pow(3, i, p) for i in range(m)))
                else:
                    A = rng.sample(range(p), m)
                mm = len(A)
                Fv = np.zeros(p, dtype=np.int64)
                for a in A:
                    Fv += np.roll(chi, -a)
                srt = np.sort(Fv)
                cs_top = np.concatenate([[0], np.cumsum(srt[::-1])])
                cs_bot = np.concatenate([[0], np.cumsum(srt)])
                n0 = math.ceil(tau * p / mm)
                for n in sorted(set([max(n0, mm), n0 + 1, 2 * n0, p // 2, p])):
                    if n < mm or n > p or mm * n < tau * p:
                        continue
                    S = max(int(cs_top[n]), -int(cs_bot[n]))
                    # |S| <= (1 - kappa/(tau M)) m n
                    ok = Fr(S) <= (1 - kap / (tau * M)) * mm * n
                    chk("F3_thm43a", ok, {"p": p, "k": k, "kappa": str(kap), "A": kind, "m": mm, "n": n, "S": S})
                    na += 1
                    worst_a = max(worst_a, S / (mm * n) / float(1 - kap / (tau * M)))
        # part (b): m > M, n >= 12 sqrt p
        for m in [M + 1, 2 * M, 5 * M, int(math.sqrt(p))]:
            for kind in ["interval", "squares", "random"]:
                if kind == "interval":
                    A = list(range(m))
                elif kind == "squares":
                    A = sorted(set((i * i) % p for i in range(1, m + 1)))
                else:
                    A = rng.sample(range(p), m)
                mm = len(A)
                if mm <= M:
                    continue
                Fv = np.zeros(p, dtype=np.int64)
                for a in A:
                    Fv += np.roll(chi, -a)
                srt = np.sort(Fv)
                cs_top = np.concatenate([[0], np.cumsum(srt[::-1])])
                cs_bot = np.concatenate([[0], np.cumsum(srt)])
                n12 = math.ceil(12 * math.sqrt(p))
                for n in sorted(set([max(n12, math.ceil(tau * p / mm)), 2 * n12, p // 3])):
                    if n > p or mm * n < tau * p or n < mm:
                        continue
                    S = max(int(cs_top[n]), -int(cs_bot[n]))
                    ok = Fr(S) <= Fr(841, 1000) * mm * n
                    chk("F3_thm43b", ok, {"p": p, "A": kind, "m": mm, "n": n, "S": S})
                    nb += 1
                    worst_b = max(worst_b, S / (mm * n))
        res43.append({"k": k, "kappa": str(kap), "tau": float(tau), "M": M, "p": p,
                      "checks_a": na, "max_ratio_a": worst_a, "checks_b": nb, "max_bias_b": worst_b})
    out["F3_thm43"] = res43
    # F4 Prop 4.1 witnesses: |A| = k, B = {e_b = 0}: bias ~ 1 at |A||B| ~ k p / 2^k
    wit = []
    p = 100003
    chi = chi_table(p)
    for k in [2, 3, 4, 5, 6]:
        A = rng.sample(range(p), k)
        e = np.zeros(p, dtype=np.int64)
        for a in A:
            e += (np.roll(chi, -a) == -1)
        B = np.nonzero(e == 0)[0]
        Fv = np.zeros(p, dtype=np.int64)
        for a in A:
            Fv += np.roll(chi, -a)
        S = int(Fv[B].sum())
        mn = k * len(B)
        wit.append({"k": k, "|B|": int(len(B)), "mn/p": mn / p, "k/2^k": k / 2 ** k, "S/(mn)": S / mn})
        chk("F4_prop41_witness", S >= mn - k and abs(mn / p - k / 2 ** k) < 0.01 * k, wit[-1])
    out["F4_prop41"] = wit
    # F5 the LP-consistency of complete bicliques: n0 = d/m satisfies (star)_e iff m <= 3d/2
    for d in [50, 500, 5000]:
        for m in range(1, min(2 * d, 400)):
            n0 = Fr(d, m)
            ok = all(n0 * (e + 1) * (m - Fr(3 * e, 2)) <= (e + 1) * (d - e) for e in range(0, (m - 1) // 2 + 1))
            chk("F5_biclique_lp_consistent", ok or m > Fr(3 * d, 2), {"d": d, "m": m})
    note("F_weil", out)
    RES["timing"]["F"] = round(time.time() - t, 1)


# ============================================================================ Section G
COS2 = {2: -2, 3: -1, 4: 0, 6: 1}   # 2*cos(2*pi/k)


def abs2_times2(N, k):
    """2|sum_j N_j zeta^j|^2 as an exact integer (k in 2,3,4,6)."""
    tot = 0
    c2 = {0: 2}
    for j in range(1, k):
        ang = 2 * math.pi * j / k
        c2[j] = int(round(2 * math.cos(ang)))
    for i in range(k):
        for j in range(k):
            tot += N[i] * N[j] * c2[(i - j) % k]
    return tot


def modulus_bound_ok(S2x2, Mtot, m, n, r, dk, k):
    """exact: |S|^2 <= M^2 (1 - 2 theta(1-theta)(1-cos 2pi/k)), theta = (1-u(w))/2,
    w = max((dk+r)/(mn), 1/4).  S2x2 = 2|S|^2 (integer)."""
    P, Q = dk + r, m * n
    if 4 * P < Q:
        P, Q = 1, 4
    cn = Fr(COS2[k], 2)                       # cos(2pi/k)
    wq = Fr(P, Q)
    Rp = 12 * wq - 3 * wq * wq
    one_minus_c = 1 - cn
    # RHS = M^2 [1 - (1-c)(2 - 6w + w^2 + w sqrt(Rp))/4]
    G = Fr(Mtot * Mtot) * (1 - one_minus_c * (2 - 6 * wq + wq * wq) / 4) - Fr(S2x2, 2)
    alpha = Fr(Mtot * Mtot) * one_minus_c * wq / 4
    return G >= 0 and alpha * alpha * Rp <= G * G


def bad_count_ok(bad, m, n, r, dk):
    """#omega-bad >= theta m n  <=>  1 - 2 bad/(mn) <= u(max(w_k,1/4))."""
    return thm14_ok(m * n - 2 * bad, m, n, r, dk)


def compositions(M, k):
    if k == 1:
        yield (M,)
        return
    for i in range(M + 1):
        for rest in compositions(M - i, k - 1):
            yield (i,) + rest


def cor52_ok(S2x2, m, n, p, k, lam):
    """exact: |S|^2 <= (mn)^2 (1 - 2 th1(1-th1)(1-cos 2pi/k)), th1 = (1-u(min(1,w1)))/2."""
    dk = (p - 1) // k
    X = (1 + lam) * dk
    m1 = math.isqrt(int(X))
    while Fr(m1 * m1) < X:
        m1 += 1
    w1 = Fr(dk + m1) / X
    if w1 >= 1:
        return Fr(S2x2, 2) <= (m * n) ** 2
    Q = m * n
    cn = Fr(COS2[k], 2)
    wq = w1
    Rp = 12 * wq - 3 * wq * wq
    omc = 1 - cn
    G = Fr(Q * Q) * (1 - omc * (2 - 6 * wq + wq * wq) / 4) - Fr(S2x2, 2)
    alpha = Fr(Q * Q) * omc * wq / 4
    return G >= 0 and alpha * alpha * Rp <= G * G


def section_G():
    t = time.time()
    out = {}
    # G0 vertex lemma: max |sum N_w w| over {N>=0, sum N = M, N_w <= (1-th)M} = M|1-th+th zeta|
    nv = 0
    for k, Mmax in [(3, 14), (4, 12), (6, 9)]:
        for Mt in range(1, Mmax + 1):
            for thn in range(0, Mt // 2 + 1):   # theta*M = thn
                best = 0
                for N in compositions(Mt, k):
                    if max(N) > Mt - thn:
                        continue
                    best = max(best, abs2_times2(N, k))
                target = abs2_times2([Mt - thn, thn] + [0] * (k - 2), k)
                chk("G0_vertex_lemma", best == target, {"k": k, "M": Mt, "thM": thn})
                nv += 1
    out["G0_vertex_cases"] = nv
    # G1 exhaustive over all (A,B) at p = 7, 13 (A containing 0; translation invariance);
    #    p = 19: sampled A, exhaustive worst B for the bad-pair count
    for p, ks in [(7, [3, 6]), (13, [3, 4, 6]), (19, [3, 6, 9])]:
        L, g = dlog_table(p)
        for k in ks:
            dk = (p - 1) // k
            cls = np.full((p, p), -1, dtype=np.int64)
            for a in range(p):
                for b in range(p):
                    y = (a + b) % p
                    if y:
                        cls[a, b] = L[y] % k
            nsub = 1 << (p - 1)
            if p == 19:
                # all A containing 0 with m <= d_k + 1 (the range of Theorem 5.1)
                idxs = [sum(1 << (i - 1) for i in T) for sz in range(0, dk + 1)
                        for T in itertools.combinations(range(1, p), sz)]
            else:
                idxs = range(nsub)
            minbad = {}
            maxS2 = {}
            nAB = 0
            allB = None
            exact_mod = (p <= 13 and k in COS2)
            if exact_mod:
                allB = np.array([[(bm >> b) & 1 for b in range(p)] for bm in range(1 << p)], dtype=np.int64)
                nBv = allB.sum(axis=1)
                c2 = [2] + [int(round(2 * math.cos(2 * math.pi * j / k))) for j in range(1, k)]
            for mask in idxs:
                A = [0] + [i + 1 for i in range(p - 1) if (mask >> i) & 1]
                m = len(A)
                cnt = np.zeros((p, k), dtype=np.int64)
                for a in A:
                    row = cls[a]
                    for b in range(p):
                        if row[b] >= 0:
                            cnt[b, row[b]] += 1
                dl = np.zeros(p, dtype=bool)
                for a in A:
                    dl[(-a) % p] = True
                if m <= dk + 1:
                    for wj in range(k):
                        badb = (m - dl.astype(np.int64)) - cnt[:, wj]
                        ins = np.sort(badb[dl])
                        outs = np.sort(badb[~dl])
                        ci = np.concatenate([[0], np.cumsum(ins)])
                        co = np.concatenate([[0], np.cumsum(outs)])
                        for r in range(0, len(ins) + 1):
                            for no in range(0, len(outs) + 1):
                                n = r + no
                                if n == 0:
                                    continue
                                key = (m, n, r)
                                v = int(ci[r] + co[no])
                                if key not in minbad or v < minbad[key]:
                                    minbad[key] = v
                if exact_mod:
                    Nsum = allB @ cnt
                    rB = allB[:, dl].sum(axis=1)
                    S2 = np.zeros(len(allB), dtype=np.int64)
                    for i in range(k):
                        for j in range(k):
                            S2 += Nsum[:, i] * Nsum[:, j] * c2[(i - j) % k]
                    nAB += len(allB)
                    key = nBv * 100 + rB
                    order = np.argsort(key, kind="stable")
                    ks_ = key[order]
                    uniq, start = np.unique(ks_, return_index=True)
                    mx = np.maximum.reduceat(S2[order], start)
                    for kk, v in zip(uniq.tolist(), mx.tolist()):
                        tup = (m, kk // 100, kk % 100)
                        if tup[1] == 0:
                            continue
                        if tup not in maxS2 or v > maxS2[tup]:
                            maxS2[tup] = v
            nb = 0
            for (m, n, r), bad in minbad.items():
                if dk + r < m * n:
                    chk("G1_theorem51_bad_pairs", bad_count_ok(bad, m, n, r, dk),
                        {"p": p, "k": k, "m": m, "n": n, "r": r, "bad": bad})
                    nb += 1
            nm = 0
            worstm = 0.0
            for (m, n, r), S2 in maxS2.items():
                Mt = m * n - r
                if m <= dk + 1 and dk + r < m * n:
                    chk("G1_theorem51_modulus", modulus_bound_ok(S2, Mt, m, n, r, dk, k),
                        {"p": p, "k": k, "m": m, "n": n, "r": r, "2|S|^2": S2})
                    nm += 1
                    w = max((dk + r) / (m * n), 0.25)
                    th = (1 - u_float(w)) / 2
                    bound = Mt * math.sqrt(1 - 2 * th * (1 - th) * (1 - math.cos(2 * math.pi / k)))
                    if bound > 0:
                        worstm = max(worstm, math.sqrt(S2 / 2) / bound)
            # Corollary 5.2 over all (A,B) (any m), lambda at the critical values
            n52 = 0
            worst52 = 0.0
            if exact_mod:
                mabs = {}
                for (m, n, r), S2 in maxS2.items():
                    mabs[(m, n)] = max(mabs.get((m, n), 0), S2)
                lam_cap = min(Fr(3), Fr(p - 1, k) - 1)
                for (m, n), S2 in mabs.items():
                    lmax = min(lam_cap, Fr(m * n, dk) - 1)
                    if lmax <= 0:
                        continue
                    cands = {lmax}
                    j = 1
                    while Fr(j * j, dk) - 1 <= lmax:
                        lj = Fr(j * j, dk) - 1
                        if lj > 0:
                            cands.add(lj)
                        j += 1
                    for lam in cands:
                        chk("G2_cor52_exhaustive", cor52_ok(S2, m, n, p, k, lam),
                            {"p": p, "k": k, "m": m, "n": n, "lam": str(lam), "2|S|^2": S2})
                        n52 += 1
            out[f"p{p}_k{k}"] = {"A_sets": len(idxs), "(A,B)_pairs_exhaustive": nAB,
                                 "bad_checks_(m,n,r)": nb, "modulus_checks_(m,n,r)": nm,
                                 "max_|S|/bound_thm51": worstm, "cor52_checks": n52}
    # G2b random A at larger p, B from direction-greedy selections (exact |S|^2 for k in 3,4,6)
    rng = random.Random(17)
    nr = 0
    for p, ks in [(31, [3, 6]), (37, [3, 4, 6]), (61, [3, 4, 6]), (73, [3, 4, 6])]:
        L, g = dlog_table(p)
        for k in ks:
            dk = (p - 1) // k
            for trial in range(40):
                m = rng.randint(1, dk + 1)
                A = rng.sample(range(p), m)
                cnt = np.zeros((p, k), dtype=np.int64)
                for a in A:
                    for b in range(p):
                        y = (a + b) % p
                        if y:
                            cnt[b, L[y] % k] += 1
                dl = np.zeros(p, dtype=bool)
                for a in A:
                    dl[(-a) % p] = True
                for wj in range(k):
                    badb = (m - dl.astype(np.int64)) - cnt[:, wj]
                    ins_ = np.sort(badb[dl])
                    outs_ = np.sort(badb[~dl])
                    ci = np.concatenate([[0], np.cumsum(ins_)])
                    co = np.concatenate([[0], np.cumsum(outs_)])
                    for r in range(0, len(ins_) + 1):
                        for no in range(0, len(outs_) + 1):
                            n = r + no
                            if n and dk + r < m * n:
                                chk("G2b_thm51_bad_pairs_random_A_all_B", bad_count_ok(int(ci[r] + co[no]), m, n, r, dk),
                                    {"p": p, "k": k, "A": A, "n": n, "r": r})
                z = np.exp(2j * np.pi * np.arange(k) / k)
                Fb = cnt @ z
                for ang in np.linspace(0, 2 * np.pi, 24, endpoint=False):
                    proj = (Fb * np.exp(-1j * ang)).real
                    order = np.argsort(-proj)
                    for n in range(1, p + 1, max(1, p // 12)):
                        B = order[:n]
                        N = cnt[B].sum(axis=0).tolist()
                        r = int(dl[B].sum())
                        S2 = abs2_times2(N, k)
                        if dk + r < m * n:
                            chk("G2b_thm51_modulus_random", modulus_bound_ok(S2, m * n - r, m, n, r, dk, k),
                                {"p": p, "k": k, "A": A, "n": n})
                            nr += 1
                        lmax = min(Fr(3), Fr(p - 1, k) - 1, Fr(m * n, dk) - 1)
                        if lmax > 0:
                            chk("G2b_cor52_random", cor52_ok(S2, m, n, p, k, lmax), {"p": p, "k": k, "A": A, "n": n})
                            nr += 1
    out["G2b_random_checks"] = nr
    # G3 sharpness example A = {0}, B = H u (part of gH)
    sh = {}
    for p in [1009, 10009]:
        L, g = dlog_table(p)
        for k in [3, 4, 6]:
            if (p - 1) % k:
                continue
            dk = (p - 1) // k
            H = [x for x in range(1, p) if L[x] % k == 0]
            gH = [x for x in range(1, p) if L[x] % k == 1]
            for lam in [Fr(1, 10), Fr(1, 2), Fr(1)]:
                extra = int(lam * dk)
                B = H + gH[:extra]
                N = [0] * k
                for b in B:
                    N[L[b] % k] += 1
                S2 = abs2_times2(N, k)
                exact = math.sqrt(S2 / 2) / len(B)
                lamr = extra / dk
                pred = math.sqrt(1 + 2 * lamr * math.cos(2 * math.pi / k) + lamr ** 2) / (1 + lamr)
                chk("G3_sharpness_example", abs(exact - pred) < 1e-12, {"p": p, "k": k, "lam": str(lam)})
                sh[f"p{p}_k{k}_lam{lam}"] = {"bias": exact, "saving": 1 - exact,
                                             "lam(1-cos)": float(lamr) * (1 - math.cos(2 * math.pi / k))}
        # threshold example A={0}, B=H: S = |A||B| = d_k exactly
        for k in [3, 4, 6]:
            if (p - 1) % k:
                continue
            H = [x for x in range(1, p) if L[x] % k == 0]
            chk("G3_threshold_example", len(H) == (p - 1) // k)
    out["G3_sharpness"] = sh
    note("G_order_k", out)
    RES["timing"]["G"] = round(time.time() - t, 1)


# ============================================================================ Section H
def section_H():
    out = {}
    tab = {}
    for kap in [0.01, 0.05, 0.1, 0.25, 0.5, 1, 1.5]:
        w = 1 / (1 + 2 * kap)
        tab[kap] = {"eta_star": 1 - u_float(w), "(1-sqrt w)^2": (1 - math.sqrt(w)) ** 2,
                    "eta/k^2": (1 - u_float(w)) / kap ** 2, "chung_saving_1-sqrt(2w)": 1 - math.sqrt(2 * w)}
    out["table"] = tab
    quoted = {0.01: 1.290e-4, 0.05: 2.844e-3, 0.1: 9.838e-3, 0.25: 0.04234, 0.5: 0.10436, 1: 0.20924, 1.5: 0.28647}
    for kap, v in quoted.items():
        chk("H1_eta_star_table", abs(tab[kap]["eta_star"] - v) < 6e-4 * max(v, 1e-3) + 5e-7, {"kappa": kap, "calc": tab[kap]["eta_star"], "quoted": v})
    quoted2 = {0.01: 0.971e-4, 0.05: 2.166e-3, 0.1: 7.591e-3, 0.25: 0.03367, 0.5: 0.08579, 1: 0.17863, 1.5: 0.25}
    for kap, v in quoted2.items():
        chk("H1_old_saving_table", abs(tab[kap]["(1-sqrt w)^2"] - v) < 6e-4 * max(v, 1e-3) + 5e-7, {"kappa": kap})
    # eta_star > (1-sqrt w)^2 on a grid, ratio -> 4/3
    ok = all(1 - u_float(w) > (1 - math.sqrt(w)) ** 2 for w in [i / 5001 for i in range(1, 5001)])
    chk("H2_eta_star_beats_old", ok)
    # where the elementary Chung bound sqrt(2w) beats u(w): w < 2 - sqrt 3, kappa > (1+sqrt3)/2
    w0 = 2 - math.sqrt(3)
    out["chung_crossing_w"] = w0
    out["chung_crossing_kappa"] = (1 / w0 - 1) / 2
    chk("H3_chung_crossing", abs(u_float(w0) - math.sqrt(2 * w0)) < 1e-12)
    # exact: at kappa = 3/2 (w = 1/4): sqrt(2w) = sqrt(1/2) < u(1/4)   <=>  (2sqrt(1/2)+1/4)^2 < 45/16
    # i.e. 2 + sqrt(2)/... do it with enclosures
    chk("H3_chung_beats_at_3/2", isqrt_frac_hi(Fr(1, 2)) < u_lo(Fr(1, 4)))
    # the small-p non-vacuity thresholds quoted for Theorem D are not in scope; the 2.14 constant:
    up = (12 - 6 * 0.25) / (2 * math.sqrt(12 * 0.25 - 3 * 0.0625))
    out["u'(1/4)"] = (up - 1) / 2
    chk("H4_uprime_quarter", abs((up - 1) / 2 - 1.0652) < 1e-4 and 2 * (up - 1) / 2 <= 2.14)
    note("H_numbers", out)


# ============================================================================ main
def main():
    secs = sys.argv[1:] or list("ABCDEFGH")
    for s in secs:
        t = time.time()
        globals()["section_" + s]()
        print(f"section {s} done in {time.time() - t:.1f}s; failures so far {len(RES['failures'])}", flush=True)
        dump(secs)


def jsonable(o):
    if isinstance(o, dict):
        return {str(k): jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [jsonable(v) for v in o]
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, Fr):
        return str(o)
    return o


def dump(secs):
    RES["sections_run"] = secs
    RES["total_checks"] = sum(RES["checks"].values())
    RES["n_failures"] = len(RES["failures"])
    RES["wall_seconds"] = round(time.time() - T0, 1)
    scratch = "/private/tmp/claude-501/-Users-shawwalters-Desktop-paleygraph/1b733267-6e28-4860-bf0c-a4c3dcfdbfb5/scratchpad/referee3"
    out = OUT if secs == list("ABCDEFGH") else os.path.join(scratch, "partial_" + "".join(secs) + ".json")
    with open(out, "w") as f:
        json.dump(jsonable(RES), f, indent=1)


if __name__ == "__main__":
    main()
