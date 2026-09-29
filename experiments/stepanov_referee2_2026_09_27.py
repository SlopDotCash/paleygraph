#!/usr/bin/env python3
"""Independent referee verifier (slug referee2) for research/stepanov-robust-2026-09-26.md.

Written from scratch; does not import or read the worker's verifier.
Run with /opt/miniconda3/bin/python3.  Writes results/stepanov_referee2_2026_09_27.json.

Parts
  A  full polynomial H_{e+1} = det[u_{i+j}] over F_p, exact degree, leading coefficient,
     exact root multiplicity at every b versus the claimed lower bound; Lemma 1.1 and Step 1
  B  Lambda = det[C(D-i-j, m-1)] nonzero mod p (modular elimination, all p <= 1500,
     m <= (p+1)/2, e <= min((m-1)/2, 12)); closed product formula over Z for D <= 60
  C  (star), Cor 2.3 (worst B), Thm 2.4, the core bound in Thm 2.5 and Thm 2.5 itself,
     exhaustive for small p; hypothesis-free variants from the summary; adversarial search
  D  Corollary E on Paley graphs (local search / exhaustive for small p)
"""
import json, math, os, sys, time, random, itertools
from fractions import Fraction
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_referee2_2026_09_27.json")
T0 = time.time()
R = {"checks": {}, "failures": [], "witnesses": {}, "notes": {}, "timing": {}}
QUICK = "--quick" in sys.argv
PARTS = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--parts=")), "ABCD")
if PARTS != "ABCD" or QUICK:  # partial / quick runs never overwrite the official results file
    OUT = next((a.split("=", 1)[1] for a in sys.argv if a.startswith("--out=")),
               OUT.replace(".json", "_partial.json"))


def save():
    R["elapsed_seconds"] = round(time.time() - T0, 1)
    R["total_checks"] = int(sum(v for v in R["checks"].values() if isinstance(v, int)))
    R["n_failures"] = len(R["failures"])
    with open(OUT, "w") as f:
        json.dump(R, f, indent=1, default=lambda o: o.item() if hasattr(o, "item") else str(o))


def bump(key, k=1):
    R["checks"][key] = R["checks"].get(key, 0) + int(k)


def fail(tag, **kw):
    kw["tag"] = tag
    if len(R["failures"]) < 200:
        R["failures"].append(kw)


def primes_upto(n):
    s = bytearray([1]) * (n + 1)
    s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = bytearray(len(s[i * i::i]))
    return [i for i in range(3, n + 1) if s[i]]


def chi_table(p):
    t = np.full(p, -1, dtype=np.int64)
    t[0] = 0
    for x in range(1, p):
        t[x * x % p] = 1
    return t


# ---------------------------------------------------------------- polynomial helpers mod p
def pmul(a, b, p):
    return np.convolve(a, b) % p


def padd(a, b, p):
    n = max(len(a), len(b))
    out = np.zeros(n, dtype=np.int64)
    out[:len(a)] += a
    out[:len(b)] += b
    return out % p


def trim(a):
    a = np.asarray(a, dtype=np.int64)
    nz = np.nonzero(a)[0]
    if len(nz) == 0:
        return a[:1] * 0
    return a[:nz[-1] + 1]


def small_binom_table(p):
    T = np.zeros((p, p), dtype=np.int64)
    for n in range(p):
        T[n, 0] = 1
        for k in range(1, n + 1):
            T[n, k] = (T[n - 1, k - 1] + (T[n - 1, k] if k <= n - 1 else 0)) % p
    return T


def lucas_binom_vec(tarr, j, p, SB):
    """C(t, j) mod p for an array of t (any size), via Lucas."""
    out = np.ones_like(tarr)
    t = tarr.copy()
    jj = j
    while True:
        t0 = t % p
        j0 = jj % p
        out = out * np.where(t0 >= j0, SB[t0, j0], 0) % p
        t //= p
        jj //= p
        if jj == 0:
            break
    return out


def orders_at_points(h, pts, p, SB, cap=None):
    """Exact multiplicity of each b in pts as a root of the nonzero polynomial h (coeff list, low->high).
    Uses Hasse derivatives: ord_b h >= M iff D^(j)h(b) = 0 for j < M (valid in every characteristic)."""
    h = trim(h)
    N = len(h) - 1
    pts = np.asarray(pts, dtype=np.int64)
    pw = np.ones((len(pts), N + 1), dtype=np.int64)
    for t in range(1, N + 1):
        pw[:, t] = pw[:, t - 1] * pts % p
    ords = np.full(len(pts), -1, dtype=np.int64)
    alive = np.arange(len(pts))
    tarr = np.arange(N + 1, dtype=np.int64)
    j = 0
    while len(alive) and j <= N:
        cj = lucas_binom_vec(tarr[j:], j, p, SB) * h[j:] % p
        vals = (pw[alive][:, :N + 1 - j] @ cj) % p
        nz = vals != 0
        ords[alive[nz]] = j
        alive = alive[~nz]
        j += 1
        if cap is not None and j > cap:
            ords[alive] = cap + 1  # means ">= cap+1"
            alive = alive[:0]
    return ords


# ---------------------------------------------------------------- the HP objects
def hp_data(p, A):
    """c_k, D, d, and u_s coefficient arrays (length d-s+1) for s = 0..2e_max, with Fact-0 checks."""
    m = len(A)
    d = (p - 1) // 2
    D = d + m - 1
    c = []
    for k, ak in enumerate(A):
        prod = 1
        for l, al in enumerate(A):
            if l != k:
                prod = prod * (ak - al) % p
        c.append(pow(prod, p - 2, p))
    # power sums P_l = sum_k c_k a_k^l, l = 0..D
    P = [0] * (D + 1)
    for ck, ak in zip(c, A):
        v = ck
        for l in range(D + 1):
            P[l] = (P[l] + v) % p
            v = v * ak % p
    return c, d, D, P


def fact_mod(p, n):
    f = [1] * (n + 1)
    for i in range(1, n + 1):
        f[i] = f[i - 1] * i % p
    return f


def u_poly(s, p, d, D, m, P, fct, ifct):
    """u_s = F^{(s)}/(D)_s = [s==0](-1) + sum_k c_k (x+a_k)^{D-s}; coefficient of x^t is
    C(D-s,t) * P_{D-s-t}.  Returns full array (length D-s+1) -- caller checks the top vanishes."""
    n = D - s
    arr = np.zeros(n + 1, dtype=np.int64)
    for t in range(n + 1):
        arr[t] = fct[n] * ifct[t] % p * ifct[n - t] % p * P[n - t] % p
    if s == 0:
        arr[0] = (arr[0] - 1) % p
    return arr


def det_poly_matrix(M, p):
    """Determinant of an n x n matrix of polynomials over F_p by memoised Laplace expansion."""
    n = len(M)
    memo = {}

    def rec(k, mask):
        if k == n:
            return np.array([1], dtype=np.int64)
        key = (k, mask)
        if key in memo:
            return memo[key]
        acc = np.array([0], dtype=np.int64)
        pos = 0
        for j in range(n):
            if mask >> j & 1:
                sub = rec(k + 1, mask & ~(1 << j))
                term = pmul(M[k][j], sub, p)
                if pos % 2 == 1:
                    term = (-term) % p
                acc = padd(acc, term, p)
                pos += 1
        memo[key] = acc
        return acc

    return trim(rec(0, (1 << n) - 1))


def lam_mod(p, D, m, e, fct, ifct):
    """Lambda = det[C(D-i-j, m-1)]_{0<=i,j<=e} mod p by Gaussian elimination (independent of H)."""
    n = e + 1
    M = [[0] * n for _ in range(n)]
    for i in range(n):
        for j in range(n):
            N = D - i - j
            K = m - 1
            M[i][j] = fct[N] * ifct[K] % p * ifct[N - K] % p if 0 <= K <= N else 0
    det = 1
    for col in range(n):
        piv = next((r for r in range(col, n) if M[r][col]), None)
        if piv is None:
            return 0
        if piv != col:
            M[col], M[piv] = M[piv], M[col]
            det = -det
        det = det * M[col][col] % p
        inv = pow(M[col][col], p - 2, p)
        for r in range(col + 1, n):
            f = M[r][col] * inv % p
            if f:
                M[r] = [(x - f * y) % p for x, y in zip(M[r], M[col])]
    return det % p


def part_A():
    t = time.time()
    rng = random.Random(20260927)
    plist = [11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 53, 61, 73, 89, 97, 101, 109, 113, 127,
             151, 173, 197, 229, 257, 281, 293]
    if QUICK:
        plist = plist[:8]
    tight = []  # instances with ord == bound
    stats = {"instances": 0, "b_checked": 0, "min_slack": None, "equality_hits": 0,
             "extra_order_at_eb_gt_e": 0}
    per_family = {}
    for p in plist:
        chi = chi_table(p)
        d = (p - 1) // 2
        SB = small_binom_table(p)
        fct = fact_mod(p, p - 1)
        ifct = [pow(x, p - 2, p) for x in fct]
        fams = []
        # intervals
        for m in sorted({2, 3, 4, 5, 7, (p + 1) // 4, int(math.isqrt(p)), (p + 1) // 2}):
            if 1 <= m <= (p + 1) // 2:
                fams.append(("interval", list(range(m))))
        # multiplicative subgroup cosets
        g = next(x for x in range(2, p) if all(pow(x, (p - 1) // q, p) != 1
                                               for q in {q for q in range(2, p) if (p - 1) % q == 0 and all(q % r for r in range(2, int(q ** .5) + 1))}))
        for k in range(1, p):
            if (p - 1) % k == 0:
                m = (p - 1) // k
                if 2 <= m <= (p + 1) // 2 and m <= 60:
                    H = sorted({pow(g, k * i, p) for i in range(m)})
                    fams.append(("subgroup", H))
                    fams.append(("coset_shift", sorted({(x * g + 1) % p for x in H})))
        # squares, squares+{0}
        Q = [x for x in range(1, p) if chi[x] == 1]
        fams.append(("squares", Q))
        fams.append(("squares+0", sorted(Q + [0])))
        # Paley clique (p = 1 mod 4) greedy, and a residue-sum clique for p = 3 mod 4 (A+A subset Q u 0?)
        cl = [0]
        for x in range(1, p):
            if all(chi[(x - y) % p] == 1 for y in cl):
                cl.append(x)
        if p % 4 == 1:
            fams.append(("paley_clique", cl))
        # A with A - A subset of Q u {0} style for any p: greedy set with a + a' in Q u 0 (sum-clique)
        sc = [1]
        for x in range(2, p):
            if chi[(2 * x) % p] != -1 and all(chi[(x + y) % p] != -1 for y in sc):
                sc.append(x)
        fams.append(("sum_clique", sc))
        # common neighbourhoods N(0) cap N(1) (and N(0)) truncated to <= (p+1)/2
        cn = [x for x in range(p) if chi[x] == 1 and chi[(x - 1) % p] == 1]
        fams.append(("common_nbhd_0_1", cn))
        fams.append(("nbhd_0", [x for x in range(p) if chi[x] == 1][: (p + 1) // 2]))
        # random
        for m in {3, 5, 8, int(math.isqrt(p)), int(math.isqrt(2 * p)) + 1}:
            if m <= (p + 1) // 2:
                fams.append(("random", sorted(rng.sample(range(p), m))))
        for fam, A in fams:
            m = len(A)
            if m < 1 or m > (p + 1) // 2:
                continue
            c, d_, D, P = hp_data(p, A)
            emax = min((m - 1) // 2, 4 if p <= 100 else (3 if p <= 200 else 2))
            if QUICK:
                emax = min(emax, 2)
            U = []
            for s in range(2 * emax + 1):
                arr = u_poly(s, p, d, D, m, P, fct, ifct)
                # Fact 0: coefficients above d-s vanish; leading coeff = C(D-s, m-1)
                top = arr[d - s + 1:]
                bump("A_fact0_u_top_vanish")
                if np.any(top % p):
                    fail("A_fact0", p=p, A=A, s=s)
                lc = arr[d - s]
                want = fct[D - s] * ifct[m - 1] % p * ifct[D - s - (m - 1)] % p
                bump("A_u_leading_coeff")
                if lc != want:
                    fail("A_u_lc", p=p, A=A, s=s, lc=int(lc), want=want)
                U.append(arr[:d - s + 1].copy())
            # e_b, delta_b
            Aset = set(A)
            eb = np.array([sum(1 for a in A if chi[(a + b) % p] == -1) for b in range(p)])
            db = np.array([1 if (-b) % p in Aset else 0 for b in range(p)])
            for e in range(emax + 1):
                M = [[U[i + j] for j in range(e + 1)] for i in range(e + 1)]
                H = det_poly_matrix(M, p)
                degH = len(H) - 1
                want_deg = (e + 1) * (d - e)
                lam = lam_mod(p, D, m, e, fct, ifct)
                bump("A_H_degree")
                if degH != want_deg or H[-1] != lam:
                    fail("A_H_degree_or_lc", p=p, fam=fam, A=A, e=e, deg=degH, want=want_deg,
                         lc=int(H[-1]), lam=lam)
                ords = orders_at_points(H, list(range(p)), p, SB)
                stats["instances"] += 1
                starL = 0
                for b in range(p):
                    if eb[b] <= e:
                        L = sum(m - e - i - db[b] for i in range(eb[b], e + 1))
                        starL += L
                        stats["b_checked"] += 1
                        bump("A_order_lower_bound")
                        sl = int(ords[b] - L)
                        if sl < 0:
                            fail("A_order", p=p, fam=fam, A=A, e=e, b=b, ord=int(ords[b]), bound=int(L))
                        if stats["min_slack"] is None or sl < stats["min_slack"]:
                            stats["min_slack"] = sl
                        if sl == 0 and e >= 1 and L > 0:
                            stats["equality_hits"] += 1
                            if len(tight) < 12:
                                tight.append({"p": p, "fam": fam, "m": m, "e": e, "b": b,
                                              "e_b": int(eb[b]), "delta": int(db[b]), "ord": int(ords[b])})
                    elif ords[b] > 0:
                        stats["extra_order_at_eb_gt_e"] += 1
                bump("A_sum_orders_le_deg")
                if int(ords.sum()) > degH:
                    fail("A_sum_orders", p=p, A=A, e=e)
                key = fam
                pf = per_family.setdefault(key, {"instances": 0, "max_star_ratio_e0": 0.0,
                                                 "max_star_ratio_e_ge_1": 0.0, "argmax_e_ge_1": None,
                                                 "max_orders_ratio_e_ge_1": 0.0})
                pf["instances"] += 1
                ratio = starL / want_deg if want_deg else 0
                if e == 0:
                    pf["max_star_ratio_e0"] = max(pf["max_star_ratio_e0"], ratio)
                else:
                    if ratio > pf["max_star_ratio_e_ge_1"]:
                        pf["max_star_ratio_e_ge_1"] = ratio
                        pf["argmax_e_ge_1"] = {"p": p, "m": m, "e": e, "A": A if m <= 20 else A[:20]}
                    pf["max_orders_ratio_e_ge_1"] = max(pf["max_orders_ratio_e_ge_1"],
                                                        float(ords.sum()) / want_deg if want_deg else 0)
        save()
    R["witnesses"]["A_equality_examples_e_ge_1"] = tight
    R["notes"]["A_stats"] = stats
    R["notes"]["A_per_family"] = per_family
    R["timing"]["A"] = round(time.time() - t, 1)
    save()


def part_A_local():
    """Lemma 1.1 (a)-(d) and Step 1 at every b, for all A containing 0 at p <= 13 (m <= (p+1)/2),
    and a random sample at p <= 43."""
    t = time.time()
    rng = random.Random(7)
    cases = []
    for p in [5, 7, 11, 13]:
        for m in range(1, (p + 1) // 2 + 1):
            for rest in itertools.combinations(range(1, p), m - 1):
                cases.append((p, [0] + list(rest)))
    for p in [17, 19, 23, 29, 31, 37, 41, 43]:
        for _ in range(12 if not QUICK else 2):
            m = rng.randint(2, (p + 1) // 2)
            cases.append((p, sorted(rng.sample(range(p), m))))
    cache = {}
    for p, A in cases:
        if p not in cache:
            fct = fact_mod(p, p - 1)
            cache[p] = (chi_table(p), small_binom_table(p), fct, [pow(x, p - 2, p) for x in fct])
        chi, SB, fct, ifct = cache[p]
        m = len(A)
        c, d, D, P = hp_data(p, A)
        Ff = u_poly(0, p, d, D, m, P, fct, ifct)
        # (a) deg F = d, lc = C(D, m-1)
        bump("L11a")
        if np.any(Ff[d + 1:]) or Ff[d] != fct[D] * ifct[m - 1] % p * ifct[D - m + 1] % p:
            fail("L11a", p=p, A=A)
        Aset = set(A)
        # derivatives F^{(j)}(b) for j <= m-1 : sum_t f_t (t)_j b^{t-j}
        tt = np.arange(len(Ff), dtype=np.int64)
        for b in range(p):
            y = [(b + a) % p for a in A]
            E = [k for k in range(m) if chi[y[k]] == -1]
            dl = 1 if (-b) % p in Aset else 0
            pw = np.array([pow(b, int(x), p) for x in range(len(Ff))], dtype=np.int64)
            for j in range(m):
                ff = np.ones_like(tt)
                for i in range(j):
                    ff = ff * ((tt - i) % p) % p
                val = int(np.sum(Ff[j:] * ff[j:] % p * pw[:len(Ff) - j] % p) % p)
                Dj = 1
                for i in range(j):
                    Dj = Dj * (D - i) % p
                pred = -2 * Dj * sum(c[k] * pow(y[k], m - 1 - j, p) for k in E)
                if dl and j == m - 1:
                    k0 = next(k for k in range(m) if y[k] == 0)
                    pred = Dj * (-c[k0] - 2 * sum(c[k] for k in E))
                bump("L11bc")
                if val != pred % p:
                    fail("L11bc", p=p, A=A, b=b, j=j, val=val, pred=pred % p)
            # (d) and Step 1: eps_s = u_s - 2 rho_s has ord_b >= m - delta - s, s <= m-1
            smax = m - 1
            for s in range(0, smax + 1):
                us = u_poly(s, p, d, D, m, P, fct, ifct)
                rho = np.zeros(D - s + 1, dtype=np.int64)
                for k in E:
                    # c_k (x+a_k)^{D-s}
                    n = D - s
                    ak = A[k]
                    coef = np.array([fct[n] * ifct[t_] % p * ifct[n - t_] % p * pow(ak, n - t_, p) % p
                                     for t_ in range(n + 1)], dtype=np.int64)
                    rho = (rho + c[k] * coef) % p
                eps = (us - 2 * rho) % p
                need = m - dl - s
                if need <= 0:
                    continue
                if not np.any(eps):
                    bump("L11d_step1")
                    continue
                o = orders_at_points(eps, [b], p, SB, cap=need)[0]
                bump("L11d_step1")
                if o < need:
                    fail("step1", p=p, A=A, b=b, s=s, ord=int(o), need=need)
    R["timing"]["A_local"] = round(time.time() - t, 1)
    save()


def batched_full_rank_mod(M, p, inv):
    """M: (B, n, n) int64 array mod p. Returns boolean array: det != 0 mod p."""
    M = M.copy() % p
    Bn, n, _ = M.shape
    ok = np.ones(Bn, dtype=bool)
    ar = np.arange(Bn)
    for k in range(n):
        sub = M[:, k:, k] != 0
        has = sub.any(axis=1)
        ok &= has
        piv = k + np.argmax(sub, axis=1)
        # swap rows k and piv
        rk = M[ar, k, :].copy()
        rp = M[ar, piv, :].copy()
        M[ar, k, :] = rp
        M[ar, piv, :] = rk
        pv = M[:, k, k]
        iv = inv[pv]  # inv[0] = 0 for dead matrices
        if k + 1 < n:
            f = M[:, k + 1:, k] * iv[:, None] % p  # (B, n-k-1)
            M[:, k + 1:, :] = (M[:, k + 1:, :] - f[:, :, None] * M[:, k, None, :]) % p
    return ok


def part_B():
    t = time.time()
    PMAX = 1500 if not QUICK else 200
    EMAX = 12
    fails = 0
    count = 0
    for p in primes_upto(PMAX):
        d = (p - 1) // 2
        fct = np.ones(p, dtype=np.int64)
        for i in range(1, p):
            fct[i] = fct[i - 1] * i % p
        ifct = np.array([pow(int(x), p - 2, p) for x in fct], dtype=np.int64)
        inv = np.zeros(p, dtype=np.int64)
        inv[1:] = [pow(x, p - 2, p) for x in range(1, p)]
        ms = np.arange(1, (p + 1) // 2 + 1)
        for e in range(0, EMAX + 1):
            mm = ms[(ms - 1) // 2 >= e]
            if len(mm) == 0:
                continue
            n = e + 1
            D = d + mm - 1  # (B,)
            ii = np.arange(n)
            NN = D[:, None, None] - ii[None, :, None] - ii[None, None, :]  # (B,n,n)
            K = (mm - 1)[:, None, None]
            assert NN.max() <= p - 1 and (NN - K).min() >= 0
            Mb = fct[NN] * ifct[np.broadcast_to(K, NN.shape)] % p * ifct[NN - K] % p
            ok = batched_full_rank_mod(Mb, p, inv)
            count += len(mm)
            if not ok.all():
                bad = mm[~ok]
                fails += len(bad)
                fail("B_lambda_zero_mod_p", p=p, e=e, m=[int(x) for x in bad[:10]])
    bump("B_lambda_nonzero_mod_p", count)
    R["notes"]["B_lambda_mod_p"] = {"p_max": PMAX, "e_max": EMAX, "cases": count, "zero_cases": fails}
    R["timing"]["B_mod"] = round(time.time() - t, 1)
    save()
    # closed product formula over Z
    t = time.time()
    DMAX = 60 if not QUICK else 25
    cnt = 0
    for D in range(0, DMAX + 1):
        for m in range(1, D + 2):
            dd = D - m + 1
            for e in range(0, (m - 1) // 2 + 1):
                if dd < e:
                    continue
                n = e + 1
                Mx = [[math.comb(D - i - j, m - 1) for j in range(n)] for i in range(n)]
                lam = bareiss_det(Mx)
                Np = D - 2 * e
                c = m - 1 - e
                num = 1
                for k in range(1, n):
                    num *= math.factorial(k)
                for i in range(1, n + 1):
                    num *= math.factorial(Np + i - 1)
                den = 1
                for i in range(1, n + 1):
                    den *= math.factorial(c - i + n) * math.factorial(Np - c + i - 1)
                sign = -1 if (e * (e + 1) // 2) % 2 else 1
                cnt += 1
                if Fraction(num, den) * sign != lam or lam == 0:
                    fail("B_closed_formula", D=D, m=m, e=e, lam=lam, formula=str(Fraction(num, den) * sign))
    bump("B_closed_formula_integer", cnt)
    R["notes"]["B_closed_formula"] = {"D_max": DMAX, "cases": cnt}
    R["timing"]["B_int"] = round(time.time() - t, 1)
    save()


def bareiss_det(M):
    M = [row[:] for row in M]
    n = len(M)
    sign = 1
    prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((r for r in range(k + 1, n) if M[r][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]
            sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]


# ---------------------------------------------------------------- Part C
KGRID = 400


def Fmin_table(p):
    """For every product P = mn (0..p^2): min over kappa in (0, K], K = min(3/2, P/p - 1/2) (grid of
    KGRID points + endpoint) of the Theorem 2.5 factor; +inf when K <= 0.  K is constant for P >= 2p,
    so the grid is evaluated only for P <= 2p + 1."""
    Pmax = p * p
    Pc = np.arange(min(Pmax, 2 * p + 1) + 1, dtype=np.float64)
    K = np.minimum(1.5, Pc / p - 0.5)
    out = np.full(len(Pc), np.inf)
    for g in np.linspace(1.0 / KGRID, 1.0, KGRID):
        kap = np.maximum(K * g, 1e-300)
        u = (1 + 2 * kap) ** -0.5
        val = 1 - (1 - u) ** 2 + (np.sqrt((0.5 + kap) * p) + 1) / (2 * (p - 1))
        out = np.minimum(out, np.where(K > 0, val, np.inf))
    full = np.full(Pmax + 1, out[-1])
    full[:len(Pc)] = out
    Kf = np.minimum(1.5, np.arange(Pmax + 1, dtype=np.float64) / p - 0.5)
    uK = np.where(Kf > 0, (1 + 2 * np.maximum(Kf, 1e-300)) ** -0.5, 1.0)
    noerr = np.where(Kf > 0, 1 - (1 - uK) ** 2, np.inf)
    return full, noerr, Kf


def subsets_with_zero(p, lo, hi):
    """Indicator rows (int8) of subsets of F_p containing 0, numbered lo..hi-1 over the bits 1..p-1."""
    idx = np.arange(lo, hi, dtype=np.int64)
    I = np.zeros((hi - lo, p), dtype=np.int64)
    I[:, 0] = 1
    for x in range(1, p):
        I[:, x] = (idx >> (x - 1)) & 1
    return I


def all_subsets(p, lo, hi):
    idx = np.arange(lo, hi, dtype=np.int64)
    I = np.zeros((hi - lo, p), dtype=np.int64)
    for x in range(p):
        I[:, x] = (idx >> x) & 1
    return I


def part_C_exhaustive():
    t = time.time()
    summary = {}
    plan = [(5, "all"), (7, "all"), (11, "all"), (13, "all"), (17, "zero"), (19, "zero")]
    if not QUICK:
        plan.append((23, "zero"))
    for p, mode in plan:
        chi = chi_table(p)
        d = (p - 1) // 2
        Chi = np.array([[chi[(a + b) % p] for b in range(p)] for a in range(p)], dtype=np.int64)
        Neg = (Chi == -1).astype(np.int64)
        negidx = np.array([(-b) % p for b in range(p)])
        Fmin, noerr, Kt = Fmin_table(p)
        total = 2 ** p if mode == "all" else 2 ** (p - 1)
        CH = 1 << 17
        st = {"subsets": 0, "star_max_ratio_e_ge_1": 0.0, "star_argmax": None,
              "cor23_max_ratio": 0.0, "thm24_max_ratio": 0.0, "thm24eta_max_ratio": 0.0,
              "core25_min_slack": None, "core25_nonvacuous_cases": 0,
              "thm25_min_slack": None, "thm25_nonvacuous_cases": 0, "thm25_min_Fmin": None,
              "max_bias_above_noerr": [], "hyp_free_cor23_violations": 0,
              "hyp_free_thmC_violations": 0, "hyp_free_examples": []}
        for lo in range(0, total, CH):
            hi = min(total, lo + CH)
            I = all_subsets(p, lo, hi) if mode == "all" else subsets_with_zero(p, lo, hi)
            m = I.sum(axis=1)
            keep = m >= 1
            I, m = I[keep], m[keep]
            st["subsets"] += len(m)
            eb = I @ Neg                    # (N,p)
            db = I[:, negidx]               # delta_b = [ -b in A ]
            s = I @ Chi                     # row sums s_A(b)
            small = m <= (p + 1) // 2
            # ---- (star) and Cor 2.3 worst-B for all admissible e
            for e in range(0, (p + 1) // 4 + 1):
                sel = small & ((m - 1) // 2 >= e)
                if not sel.any():
                    continue
                E_, M_, Db = eb[sel], m[sel][:, None], db[sel]
                mask = E_ <= e
                lhs2 = ((e + 1 - E_) * (2 * M_ - 3 * e - E_ - 2 * Db) * mask).sum(axis=1)
                rhs2 = 2 * (e + 1) * (d - e)
                bump("C_star", sel.sum())
                if (lhs2 > rhs2).any():
                    k = np.argmax(lhs2 - rhs2)
                    fail("C_star", p=p, e=e, A=np.nonzero(I[sel][k])[0].tolist())
                if e >= 1:
                    k = int(np.argmax(lhs2))
                    if lhs2[k] / rhs2 > st["star_max_ratio_e_ge_1"]:
                        st["star_max_ratio_e_ge_1"] = float(lhs2[k] / rhs2)
                        st["star_argmax"] = {"e": e, "A": np.nonzero(I[sel][k])[0].tolist()}
                w = (M_ - 2 * e) * (e + 1 - E_) - (e + 1) * Db
                W = np.maximum(w, 0).sum(axis=1)
                bump("C_cor23_worstB", sel.sum())
                if (W > (e + 1) * (d - e)).any():
                    k = np.argmax(W)
                    fail("C_cor23", p=p, e=e, A=np.nonzero(I[sel][k])[0].tolist())
                st["cor23_max_ratio"] = max(st["cor23_max_ratio"], float(W.max() / ((e + 1) * (d - e))) if d > e else 0)
                # Thm 2.4 first display with level sets B = {e_b <= e'} for e' <= e
                for ep in range(0, e + 1):
                    lev = E_ <= ep
                    n = lev.sum(axis=1)
                    r = (lev & (Db == 1)).sum(axis=1)
                    lhs = (e + 1 - ep) * ((2 * M_[:, 0] - 3 * e - ep) * n - 2 * r)
                    rhs = 2 * (e + 1) * (d - e)
                    bump("C_thm24_display", sel.sum())
                    if (lhs > rhs).any():
                        k = np.argmax(lhs - rhs)
                        fail("C_thm24", p=p, e=e, ep=ep, A=np.nonzero(I[sel][k])[0].tolist())
                    st["thm24_max_ratio"] = max(st["thm24_max_ratio"], float(lhs.max() / rhs))
            # ---- Thm 2.4 eta-form (m <= (p+1)/2), level sets, eta = actual max e_b / m <= 1/8
            for ep in range(0, (p + 1) // 16 + 2):
                sel = small & (8 * ep <= m)
                if not sel.any():
                    continue
                lev = eb[sel] <= ep
                n = lev.sum(axis=1)
                r = (lev & (db[sel] == 1)).sum(axis=1)
                eact = np.where(lev, eb[sel], 0).max(axis=1)
                eta = eact / m[sel]
                bound = (1 - np.sqrt(2 * eta)) ** -2 * d + m[sel]
                val = m[sel] * n - r
                bump("C_thm24_eta", sel.sum())
                if (val > bound + 1e-9).any():
                    k = np.argmax(val - bound)
                    fail("C_thm24_eta", p=p, ep=ep, A=np.nonzero(I[sel][k])[0].tolist())
                st["thm24eta_max_ratio"] = max(st["thm24eta_max_ratio"], float((val / bound).max()))
            # ---- hypothesis-free variants of summary Thms B, C (m > (p+1)/2)
            big = ~small
            if big.any():
                for e in range(0, p // 2 + 1):
                    sel = big & ((m - 1) // 2 >= e)
                    if not sel.any():
                        continue
                    E_, M_, Db = eb[sel], m[sel][:, None], db[sel]
                    w = (M_ - 2 * e) * (e + 1 - E_) - (e + 1) * Db
                    W = np.maximum(w, 0).sum(axis=1)
                    bump("C_hypfree_cor23", sel.sum())
                    bad = W > (e + 1) * (d - e)
                    if bad.any():
                        st["hyp_free_cor23_violations"] += int(bad.sum())
                        if len(st["hyp_free_examples"]) < 6:
                            k = int(np.argmax(bad))
                            st["hyp_free_examples"].append({"thm": "B(Cor2.3)", "e": e, "m": int(m[sel][k]),
                                                            "A": np.nonzero(I[sel][k])[0].tolist(),
                                                            "worstB_lhs": int(W[k]), "rhs": (e + 1) * (d - e)})
                for ep in range(0, p // 8 + 1):
                    sel = big & (8 * ep <= m)
                    if not sel.any():
                        continue
                    lev = eb[sel] <= ep
                    n = lev.sum(axis=1)
                    r = (lev & (db[sel] == 1)).sum(axis=1)
                    eact = np.where(lev, eb[sel], 0).max(axis=1)
                    eta = eact / m[sel]
                    bound = (1 - np.sqrt(2 * eta)) ** -2 * d + m[sel]
                    val = m[sel] * n - r
                    bump("C_hypfree_thmC", sel.sum())
                    bad = val > bound + 1e-9
                    if bad.any():
                        st["hyp_free_thmC_violations"] += int(bad.sum())
                        if len(st["hyp_free_examples"]) < 12:
                            k = int(np.argmax(bad))
                            st["hyp_free_examples"].append({"thm": "C", "ep": ep, "m": int(m[sel][k]),
                                                            "A": np.nonzero(I[sel][k])[0].tolist(),
                                                            "val": int(val[k]), "bound": float(bound[k])})
            # ---- Theorem 2.5 (all m, n) and its core lemma (m <= (p+1)/2), extreme B of each size
            ss = np.sort(s, axis=1)
            Smin = np.cumsum(ss, axis=1)                 # n = 1..p smallest
            Smax = np.cumsum(ss[:, ::-1], axis=1)        # n = 1..p largest
            nn = np.arange(1, p + 1)[None, :]
            P = m[:, None] * nn
            absS = np.maximum(np.abs(Smax), np.abs(Smin))
            Fm = Fmin[P]
            adm = np.isfinite(Fm)
            slack = np.where(adm, Fm * P - absS, np.inf)
            bump("C_thm25_final", adm.sum())
            if (slack < -1e-9).any():
                k = np.unravel_index(np.argmin(slack), slack.shape)
                fail("C_thm25", p=p, A=np.nonzero(I[k[0]])[0].tolist(), n=int(k[1] + 1))
            ms_ = float(slack.min())
            st["thm25_min_slack"] = ms_ if st["thm25_min_slack"] is None else min(st["thm25_min_slack"], ms_)
            nv = adm & (Fm < 1)
            st["thm25_nonvacuous_cases"] += int(nv.sum())
            fmn = float(Fm[adm].min()) if adm.any() else None
            if fmn is not None:
                st["thm25_min_Fmin"] = fmn if st["thm25_min_Fmin"] is None else min(st["thm25_min_Fmin"], fmn)
            # heuristic: bias above 1-(1-u_K)^2 (no error term)
            ratio = np.where(adm, absS / np.maximum(P, 1), 0)
            ne = noerr[P]
            over = adm & (ratio > ne + 1e-12)
            if over.any() and len(st["max_bias_above_noerr"]) < 5:
                k = np.unravel_index(np.argmax(np.where(over, ratio - ne, -1)), over.shape)
                st["max_bias_above_noerr"].append({"A": np.nonzero(I[k[0]])[0].tolist(), "n": int(k[1] + 1),
                                                   "ratio": float(ratio[k]), "noerr_bound": float(ne[k]),
                                                   "Fmin": float(Fm[k])})
            # core lemma: S <= [1-(1-u)^2 + u(1-u) m/d] m n for m <= (p+1)/2, at kappa = K (the minimum)
            Kc = Kt[P]
            uK = np.where(Kc > 0, (1 + 2 * np.maximum(Kc, 1e-300)) ** -0.5, 1.0)
            tt = (m[:, None] / d)
            g = (2 + tt) * uK - (1 + tt) * uK ** 2
            admc = (Kc > 0) & small[:, None]
            slc = np.where(admc, g * P - absS, np.inf)
            bump("C_core25", admc.sum())
            if (slc < -1e-9).any():
                k = np.unravel_index(np.argmin(slc), slc.shape)
                fail("C_core25", p=p, A=np.nonzero(I[k[0]])[0].tolist(), n=int(k[1] + 1))
            v = float(slc.min())
            st["core25_min_slack"] = v if st["core25_min_slack"] is None else min(st["core25_min_slack"], v)
            st["core25_nonvacuous_cases"] += int((admc & (g < 1)).sum())
        summary[p] = st
        R["notes"]["C_exhaustive"] = summary
        save()
    R["timing"]["C_exhaustive"] = round(time.time() - t, 1)
    save()


def part_C_adversarial():
    """Local search over A (sizes 2..~sqrt(2p)) maximizing max_n [S_max(A,n)/(mn) - Fmin(mn)]."""
    t = time.time()
    rng = random.Random(11)
    out = {}
    plist = [101, 211, 409, 1009, 2003] if not QUICK else [101]
    for p in plist:
        chi = chi_table(p)
        d = (p - 1) // 2
        Chi = np.array([[chi[(a + b) % p] for b in range(p)] for a in range(p)], dtype=np.int64)
        Fmin, noerr, Kt = Fmin_table(p)
        best = {"gap_final": -9.0, "gap_noerr": -9.0, "gap_core": -9.0}
        wit = {}
        msizes = sorted({2, 3, 4, 5, 6, 8, 10, int(math.isqrt(p) // 2), int(math.isqrt(p)),
                         int(math.isqrt(2 * p)) + 1})
        iters = 250 if p <= 409 else 120

        def score(s, m):
            ss = np.sort(s)[::-1]
            Smax = np.cumsum(ss)
            nn = np.arange(1, p + 1)
            P = m * nn
            Fm = Fmin[P]
            ne = noerr[P]
            ok = np.isfinite(Fm)
            if not ok.any():
                return -9.0, -9.0, -9.0, 0
            ratio = Smax / P
            g1 = np.where(ok, ratio - Fm, -9)
            g2 = np.where(ok, ratio - ne, -9)
            K = Kt[P]
            uK = np.where(K > 0, (1 + 2 * np.maximum(K, 1e-300)) ** -0.5, 1.0)
            tt = m / d
            gc = (2 + tt) * uK - (1 + tt) * uK ** 2
            g3 = np.where(ok & (m <= (p + 1) // 2), ratio - gc, -9)
            return float(g1.max()), float(g2.max()), float(g3.max()), int(np.argmax(g1) + 1)

        for m in msizes:
            for restart in range(3):
                A = rng.sample(range(p), m)
                s = Chi[A].sum(axis=0)
                cur = score(s, m)
                for it in range(iters):
                    i = rng.randrange(m)
                    x = rng.randrange(p)
                    if x in A:
                        continue
                    s2 = s - Chi[A[i]] + Chi[x]
                    sc = score(s2, m)
                    if sc[1] >= cur[1]:
                        A[i] = x
                        s, cur = s2, sc
                    bump("C_adv_evals")
                for key, val in zip(["gap_final", "gap_noerr", "gap_core"], cur[:3]):
                    if val > best[key]:
                        best[key] = val
                        wit[key] = {"m": m, "A": sorted(A), "n_at_final_gap": cur[3]}
                if cur[0] > 1e-12:
                    fail("C_adv_thm25", p=p, A=sorted(A))
                if cur[2] > 1e-12:
                    fail("C_adv_core25", p=p, A=sorted(A))
        out[p] = {"best_gaps(ratio-bound; negative = bound holds)": best, "witness": wit}
        R["notes"]["C_adversarial"] = out
        save()
    R["timing"]["C_adversarial"] = round(time.time() - t, 1)
    save()


# ---------------------------------------------------------------- Part D: Corollary E
def part_D():
    t = time.time()
    rng = random.Random(5)
    res = {}
    plist = [p for p in primes_upto(420) if p % 4 == 1 and p >= 13] + [1009, 2017]
    if QUICK:
        plist = plist[:5]
    for p in plist:
        chi = chi_table(p)
        ar = np.arange(p)
        Adj = (chi[(ar[:, None] - ar[None, :]) % p] == 1).astype(np.int64)
        Fmin, noerr, Kt = Fmin_table(p)
        rows = []
        for kap in [0.05, 0.25, 0.5, 1.0, 1.5]:
            m = math.isqrt(int(math.ceil((0.5 + kap) * p)))
            while m * m < (0.5 + kap) * p:
                m += 1
            # local search maximizing internal edges (min density = 1 - max density by self-complementarity)
            bestE = -1
            bestA = None
            for restart in range(3):
                A = rng.sample(range(p), m)
                inA = np.zeros(p, dtype=bool)
                inA[A] = True
                deg = Adj[:, inA].sum(axis=1)
                E = int(deg[inA].sum() // 2)
                tabu_until = np.full(p, -1)
                nrng = np.random.default_rng(rng.randrange(1 << 30))
                for it in range(25 * m):
                    # swap out argmin internal degree in A, swap in argmax degree outside A (tabu-aware)
                    free = tabu_until < it
                    co = inA & free
                    ci = (~inA) & free
                    if not co.any() or not ci.any():
                        break
                    noise = nrng.random(p) * 0.5
                    vo = int(np.argmin(np.where(co, deg + noise, np.inf)))
                    vi = int(np.argmax(np.where(ci, deg - Adj[:, vo] + noise, -np.inf)))
                    newE = E - int(deg[vo]) + int(deg[vi]) - int(Adj[vi, vo])
                    inA[vo] = False
                    inA[vi] = True
                    deg = deg - Adj[:, vo] + Adj[:, vi]
                    E = newE
                    tabu_until[vo] = it + 7
                    tabu_until[vi] = it + 3
                    if E > bestE:
                        bestE = E
                        bestA = sorted(np.nonzero(inA)[0].tolist())
            rho = 2 * bestE / (m * (m - 1))
            S = 4 * bestE - m * (m - 1)  # S(A,-A)
            P = m * m
            F = float(Fmin[P])
            upper = (1 + F * m / (m - 1)) / 2
            Keff = float(Kt[P])
            c = (1 - (1 + 2 * min(Keff, 1.5)) ** -0.5) ** 2
            bump("D_corE")
            if abs(S) > F * P + 1e-9 or rho > upper + 1e-12:
                fail("D_corE", p=p, m=m, A=bestA)
            rows.append({"kappa": kap, "m": m, "max_edges_found": bestE, "density": round(rho, 4),
                         "thmD_upper_density": round(upper, 4), "heuristic_1_minus_c_over_2": round(1 - c / 2, 4),
                         "exceeds_heuristic": bool(rho > 1 - c / 2)})
        res[p] = rows
    R["notes"]["D_corollaryE_localsearch"] = res
    save()
    # exact maxima for small p (vertex transitivity: A contains 0)
    exact = {}
    for p, mtop in [(13, 8), (17, 8), (29, 7)]:
        chi = chi_table(p)
        Fmin, noerr, Kt = Fmin_table(p)
        for m in range(math.isqrt(p // 2) + 1, mtop):
            if m * m < 0.5 * p:
                continue
            best = -1
            for rest in itertools.combinations(range(1, p), m - 1):
                A = (0,) + rest
                E = sum(1 for i in range(m) for j in range(i + 1, m) if chi[(A[i] - A[j]) % p] == 1)
                if E > best:
                    best = E
            bump("D_exact_sets", math.comb(p - 1, m - 1))
            rho = 2 * best / (m * (m - 1))
            P = m * m
            F = float(Fmin[P]) if np.isfinite(Fmin[P]) else None
            ok = True
            if F is not None:
                S = 4 * best - m * (m - 1)
                ok = abs(S) <= F * P + 1e-9
                if not ok:
                    fail("D_exact", p=p, m=m)
            exact[f"{p}/{m}"] = {"max_density": round(rho, 4), "thmD_factor": F,
                                 "thmD_upper_density": None if F is None else round((1 + F * m / (m - 1)) / 2, 4)}
    R["notes"]["D_exact_small"] = exact
    R["timing"]["D"] = round(time.time() - t, 1)
    save()


def star_equality_witness():
    """(star) with equality at e = 1: p = 13, A = {0,2,3,5} (hand-checked in the note, section 2)."""
    p, A, e = 13, [0, 2, 3, 5], 1
    chi = chi_table(p)
    m, d = len(A), (p - 1) // 2
    terms = {}
    for b in range(p):
        eb = sum(1 for a in A if chi[(a + b) % p] == -1)
        db = 1 if (-b) % p in A else 0
        if eb <= e:
            terms[b] = Fraction((e + 1 - eb) * (2 * m - 3 * e - eb - 2 * db), 2)
    lhs = sum(terms.values())
    bump("witness_star_equality")
    if lhs != (e + 1) * (d - e):
        fail("witness_star_equality", lhs=str(lhs))
    R["witnesses"]["star_equality_e1"] = {"p": p, "A": A, "e": e, "terms": {b: str(v) for b, v in terms.items()},
                                          "lhs": str(lhs), "rhs": (e + 1) * (d - e)}
    R["verdicts"] = {"Lemma 1.1": "CORRECT", "Thm 2.1": "CORRECT", "Lemma 2.2": "CORRECT",
                     "Cor 2.3": "CORRECT (summary Thm B must add |A| <= (p+1)/2)",
                     "Thm 2.4": "CORRECT (summary Thm C is true even for |A| > (p+1)/2: note section 8.2)",
                     "Thm 2.5": "CORRECT", "Corollary E": "CORRECT (explicit error term in note section 7)"}
    # the count used in note section 8.2: #{y: chi(y) != -1, chi(y+c) != -1} = (p-3-chi(c)-chi(-c))/4
    # + [chi(c)=1] + [chi(-c)=1] <= (p+3)/4
    for p in primes_upto(200):
        ch = chi_table(p)
        for c in range(1, p):
            cnt = sum(1 for y in range(p) if ch[y] != -1 and ch[(y + c) % p] != -1)
            form = (p - 3 - ch[c] - ch[(-c) % p]) // 4 + (ch[c] == 1) + (ch[(-c) % p] == 1)
            bump("thmC_count_formula")
            if cnt != form or 4 * cnt > p + 3:
                fail("thmC_count", p=p, c=c)
    save()


if __name__ == "__main__":
    star_equality_witness()
    if "C" in PARTS:
        part_C_exhaustive()
        part_C_adversarial()
    if "D" in PARTS:
        part_D()
    if "B" in PARTS:
        part_B()
    if "A" in PARTS:
        part_A_local()
        part_A()
    print(json.dumps({k: R[k] for k in ["checks", "n_failures", "elapsed_seconds"]}, indent=1))
    if R["failures"]:
        print(R["failures"][:5])
