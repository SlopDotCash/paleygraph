#!/usr/bin/env python3
"""Exact verifier for research/stepanov-robust-2026-09-26.md (worker `robust`).

Checks (all exact integer / rational arithmetic, mod-p polynomial algebra):
  A. HP baseline: deg F = d, and the derivative formula
       F^{(j)}(b) = -2 (D)_j sum_{k in E(b)} c_k (b+a_k)^{m-1-j}   (0<=j<=m-1, b not in -A;
       0<=j<=m-2 and the j=m-1 correction for b in -A).
  B. Cheap robust bound (HP on random subsets), the identity-level inequality
       sum_b C(m-1-e_b,t-1)(m-e_b-[b in -A]) <= C(m,t) d   for 1<=t<=m,
     and the lower bound min_t C(m,t)/C(m-1-e,t-1) >= e*exp(1-e/m).
  C. Theorem 1 (Hankel-minor Stepanov inequality):
       H_{e+1}(x) = det[u_{i+j}(x)]_{0<=i,j<=e},  u_s = F^{(s)}/(D)_s,
     (C1) leading coefficient det[C(D-i-j,m-1)] equals Krattenthaler (3.12) product and is
          nonzero mod p; (C2) full-polynomial degree check on a sample; (C3) vanishing orders
          ord_b H_{e+1} >= sum_{i=e_b}^{e} (m-e-i-[b in -A]) for every b with e_b<=e;
     (C4) the weighted inequality (star).  Exhaustive over all A containing 0 for small p,
     random and structured A for larger p.
  D. Consequences: bias bound S <= mn - r - 2(e+1)[n-(d-e+r)/(m-2e)], RHP with
     f(eta)=(1-sqrt(2 eta))^-2, g(m)=m, and the kappa-theorem.
  E. No-go construction for subset-HP (abstract set system), exact integers.
  F. Prime sensitivity: the normalisation (D)_j vanishes mod p over F_{p^2} for A=F_p.
Writes results/stepanov_robust_2026_09_26.json.
"""
import itertools, json, math, os, random, sys, time
from fractions import Fraction
import numpy as np

T0 = time.time(); C0 = time.process_time()
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "stepanov_robust_2026_09_26.json")
R = {"checks": {}, "witnesses": {}, "notes": {}}
FAIL = []

def bump(key, n=1):
    R["checks"][key] = R["checks"].get(key, 0) + n

def fail(key, info):
    FAIL.append((key, info))
    R["witnesses"].setdefault("FAILURES", []).append({"key": key, "info": str(info)[:500]})

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

def ff(x, j):  # falling factorial (x)_j as integer
    r = 1
    for i in range(j):
        r *= (x - i)
    return r

class Inst:
    """Everything attached to (p, A)."""
    def __init__(self, p, A, idx=2):
        self.p = p; self.A = [a % p for a in A]; self.m = m = len(A)
        assert (p - 1) % idx == 0
        self.d = d = (p - 1) // idx; self.D = D = d + m - 1
        assert len(set(self.A)) == m and 1 <= m and D <= p - 1
        self.c = []
        for k in range(m):
            pr = 1
            for l in range(m):
                if l != k:
                    pr = pr * (self.A[k] - self.A[l]) % p
            self.c.append(pow(pr, p - 2, p))
        self.chi = [[leg(a + b, p) for a in self.A] for b in range(p)]
        if idx == 2:
            self.eb = [sum(1 for v in row if v == -1) for row in self.chi]
        else:  # bad partner: a+b nonzero and (a+b)^d != 1
            self.eb = [sum(1 for a in self.A if (a + b) % p and pow((a + b) % p, d, p) != 1)
                       for b in range(p)]
        self.negA = set((-a) % p for a in self.A)
        # binomials C(D,j) mod p (D<p so Lucas-free)
        self.CD = [math.comb(D, j) % p for j in range(D + 1)]

    def F_coeffs(self):
        p, D, m = self.p, self.D, self.m
        F = [0] * (D + 1)
        for k in range(m):
            ak = self.A[k]; pw = 1
            pws = [1] * (D + 1)
            for i in range(1, D + 1):
                pws[i] = pws[i - 1] * ak % p
            for t in range(D + 1):
                F[t] = (F[t] + self.c[k] * self.CD[t] * pws[D - t]) % p
        F[0] = (F[0] - 1) % p
        while F and F[-1] == 0:
            F.pop()
        return F

    def hasse(self, b, J):
        """tau_j(b) = F^{(j)}(b)/j! for j=0..J (exact mod p)."""
        p, D = self.p, self.D
        ys = [(b + a) % p for a in self.A]
        out = []
        for j in range(J + 1):
            if j > D:
                out.append(0); continue
            s = 0
            for k, y in enumerate(ys):
                if y == 0:
                    if D - j == 0:
                        s += self.c[k]
                    continue
                s += self.c[k] * pow(y, D - j, p)
            s = s * self.CD[j] % p
            if j == 0:
                s = (s - 1) % p
            out.append(s)
        return out

def inv(x, p):
    return pow(x % p, p - 2, p)

# ---------------------------------------------------------------- A. HP baseline
def check_HP_derivatives(I):
    p, m, D, d = I.p, I.m, I.D, I.d
    F = I.F_coeffs()
    if len(F) - 1 != d:
        fail("A.degF", (p, I.A, len(F) - 1)); return
    bump("A.degF_exact")
    for b in range(p):
        tau = I.hasse(b, m - 1)
        E = [k for k in range(m) if I.chi[b][k] == -1]
        k0 = [k for k in range(m) if (b + I.A[k]) % p == 0]
        for j in range(m):
            Fj = tau[j] * math.factorial(j) % p           # F^{(j)}(b)
            pred = 0
            for k in E:
                pred += I.c[k] * pow((b + I.A[k]) % p, m - 1 - j, p)
            pred = (-2 * ff(D, j) * pred) % p
            if k0 and j == m - 1:
                sE = sum(I.c[k] for k in E)
                pred = (ff(D, m - 1) * (-I.c[k0[0]] - 2 * sE)) % p
            if Fj != pred:
                fail("A.deriv", (p, I.A, b, j, Fj, pred))
            bump("A.deriv_formula")
    # HP itself on the complete biclique B = {e_b = 0}
    B0 = [b for b in range(p) if I.eb[b] == 0]
    r0 = sum(1 for b in B0 if b in I.negA)
    if m * len(B0) - r0 > d:
        fail("A.HP", (p, I.A))
    bump("A.HP_inequality")

# ---------------------------------------------------------------- B. cheap bound
def check_cheap(I):
    p, m, d = I.p, I.m, I.d
    for t in range(1, m + 1):
        lhs = 0
        for b in range(p):
            eb = I.eb[b]
            if m - 1 - eb >= t - 1:
                lhs += math.comb(m - 1 - eb, t - 1) * (m - eb - (1 if b in I.negA else 0))
        if lhs > math.comb(m, t) * d:
            fail("B.cheap_identity", (p, I.A, t, lhs, math.comb(m, t) * d))
        bump("B.cheap_identity")

def cheap_constant(m, e):
    best = None
    for t in range(1, m - e + 1):
        v = Fraction(math.comb(m, t), math.comb(m - 1 - e, t - 1))
        if best is None or v < best:
            best = v
    return best

# ---------------------------------------------------------------- C. Hankel minors
def series_u(I, b, e, T):
    """Truncated Taylor series (length T) at b of u_s, s=0..2e, as numpy int64 arrays."""
    p, D = I.p, I.D
    tau = I.hasse(b, 2 * e + T)
    us = []
    for s in range(2 * e + 1):
        fs = ff(D, s) % p
        base = math.factorial(s) * inv(fs, p) % p
        arr = np.zeros(T, dtype=np.int64)
        for i in range(T):
            if s + i < len(tau):
                arr[i] = tau[s + i] * math.comb(s + i, i) % p * base % p
        us.append(arr)
    return us

def series_det(us, e, p, T):
    n = e + 1
    tot = np.zeros(T, dtype=np.int64)
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sgn = -sgn
        t = np.zeros(T, dtype=np.int64); t[0] = 1
        for i in range(n):
            t = np.convolve(t, us[i + perm[i]])[:T] % p
        tot = (tot + sgn * t) % p
    return tot

def valuation(arr):
    nz = np.nonzero(arr)[0]
    return int(nz[0]) if len(nz) else None

def ord_lower(m, e, eb, neg):
    return sum(max(0, m - e - i - (1 if neg else 0)) for i in range(eb, e + 1))

def check_theorem1(I, emax=None, full_poly=False):
    """(C3)+(C4) for all admissible e; returns list of per-e slack records."""
    p, m, d = I.p, I.m, I.d
    recs = []
    top = (m - 1) // 2 if emax is None else min(emax, (m - 1) // 2)
    for e in range(0, top + 1):
        T = (e + 1) * m + 2
        lhs = Fraction(0)
        sumord = 0
        for b in range(p):
            eb = I.eb[b]
            if eb > e:
                continue
            neg = b in I.negA
            us = series_u(I, b, e, T)
            H = series_det(us, e, p, T)
            v = valuation(H)
            need = ord_lower(m, e, eb, neg)
            if v is not None and v < need:
                fail("C3.order", (p, I.A, e, b, eb, v, need))
            bump("C3.order_lower_bound")
            sumord += (T if v is None else v)
            lhs += Fraction((e + 1 - eb) * (2 * m - 3 * e - eb - (2 if neg else 0)), 2)
        rhs = (e + 1) * (d - e)
        if lhs > rhs:
            fail("C4.star", (p, I.A, e, lhs, rhs))
        bump("C4.star_inequality")
        recs.append((float(lhs / rhs) if rhs else 0.0, e))
    return recs

def poly_mul(a, b, p):
    return (np.convolve(np.array(a, dtype=object), np.array(b, dtype=object)) % p).tolist()

def full_H(I, e):
    """Exact polynomial H_{e+1}(x) mod p (object arithmetic)."""
    p, D = I.p, I.D
    F = I.F_coeffs()
    us = []
    for s in range(2 * e + 1):
        iv = inv(ff(D, s) % p, p)
        us.append([F[t] * ff(t, s) % p * iv % p for t in range(s, len(F))])
    n = e + 1
    L = (e + 1) * (I.d - e) + 1
    tot = [0] * L
    for perm in itertools.permutations(range(n)):
        sgn = 1
        for i in range(n):
            for j in range(i + 1, n):
                if perm[i] > perm[j]:
                    sgn = -sgn
        t = [1]
        for i in range(n):
            t = poly_mul(t, us[i + perm[i]], p)
        for k, v in enumerate(t):
            tot[k] = (tot[k] + sgn * v) % p
    while tot and tot[-1] == 0:
        tot.pop()
    return tot

def poly_ord(poly, b, p):
    a = list(poly); k = 0
    while a:
        n = len(a) - 1; q = [0] * n; r = a[n]
        for i in range(n - 1, -1, -1):
            q[i] = r; r = (a[i] + r * b) % p
        if r != 0:
            return k
        k += 1; a = q
        while a and a[-1] == 0:
            a.pop()
    return None

def kratt_product(Nprime, c, n):
    """Krattenthaler (3.12), q=1, L_i=c-i (i=1..n), A=Nprime: det[C(A, c-i+j)]_{1..n}."""
    num = 1; den = 1
    for i in range(1, n + 1):
        for j in range(i + 1, n + 1):
            num *= (j - i)          # L_i - L_j
    for i in range(1, n + 1):
        num *= math.factorial(Nprime + i - 1)
        den *= math.factorial(c - i + n) * math.factorial(Nprime - c + i - 1)
    assert num % den == 0
    return num // den

def int_det(M):
    """Fraction-free Bareiss determinant of an integer matrix."""
    M = [row[:] for row in M]; n = len(M); sign = 1; prev = 1
    for k in range(n - 1):
        if M[k][k] == 0:
            sw = next((i for i in range(k + 1, n) if M[i][k] != 0), None)
            if sw is None:
                return 0
            M[k], M[sw] = M[sw], M[k]; sign = -sign
        for i in range(k + 1, n):
            for j in range(k + 1, n):
                M[i][j] = (M[i][j] * M[k][k] - M[i][k] * M[k][j]) // prev
        prev = M[k][k]
    return sign * M[n - 1][n - 1]

def mod_det(M, p):
    M = [row[:] for row in M]; n = len(M); det = 1
    for k in range(n):
        piv = next((i for i in range(k, n) if M[i][k] % p), None)
        if piv is None:
            return 0
        if piv != k:
            M[k], M[piv] = M[piv], M[k]; det = -det
        det = det * M[k][k] % p
        iv = pow(M[k][k], p - 2, p)
        for i in range(k + 1, n):
            f = M[i][k] * iv % p
            for j in range(k, n):
                M[i][j] = (M[i][j] - f * M[k][j]) % p
    return det % p

def check_leading_coefficient():
    # (C1a) exact integer identity det[C(D-i-j,m-1)] = (-1)^{e(e+1)/2} * Kratt(D-2e, m-1-e, e+1)
    cnt = 0
    for D in range(2, 41):
        for m in range(1, min(D + 2, 26)):
            K = m - 1
            for e in range(0, (m - 1) // 2 + 1):
                if D - 2 * e < K or 2 * e > D:   # need D-2e >= m-1 (i.e. d >= e) as in our setting
                    continue
                M = [[math.comb(D - i - j, K) for j in range(e + 1)] for i in range(e + 1)]
                lhs = int_det(M)
                rhs = (-1) ** (e * (e + 1) // 2) * kratt_product(D - 2 * e, K - e, e + 1)
                if lhs != rhs:
                    fail("C1.kratt_identity", (D, m, e, lhs, rhs))
                cnt += 1
    bump("C1.leading_coeff_identity", cnt)
    # (C1b) independent check: det[C(D-i-j,m-1)] mod p by modular elimination, nonzero
    cnt = 0
    for p in primes_upto(300)[1:]:
        d = (p - 1) // 2
        for m in range(1, min(40, (p + 1) // 2) + 1):
            D = d + m - 1
            for e in range(0, (m - 1) // 2 + 1):
                M = [[math.comb(D - i - j, m - 1) % p for j in range(e + 1)] for i in range(e + 1)]
                if mod_det(M, p) == 0:
                    fail("C1.nonzero_mod_p", (p, m, e))
                cnt += 1
    bump("C1.leading_coeff_nonzero_mod_p", cnt)

# ---------------------------------------------------------------- D. consequences
def S_value(I, B):
    return sum(I.chi[b][k] for b in B for k in range(I.m))

def check_bias(I, B, tag):
    p, m, d = I.p, I.m, I.d
    n = len(B)
    if n == 0:
        return
    r = sum(1 for b in B if b in I.negA)
    Nneg = sum(I.eb[b] for b in B)
    S = S_value(I, B)
    assert S == m * n - r - 2 * Nneg
    for e in range(0, (m - 1) // 2 + 1):
        # (m-2e)((e+1)n - Nneg) <= (e+1)(d-e+r)
        if (m - 2 * e) * ((e + 1) * n - Nneg) > (e + 1) * (d - e + r):
            fail("D.weighted_N", (tag, p, I.A, e))
        bump("D.weighted_Nneg_inequality")
        bound = Fraction(m * n - r) - 2 * (e + 1) * (Fraction(n) - Fraction(d - e + r, m - 2 * e))
        if S > bound:
            fail("D.bias_bound", (tag, p, I.A, e, S, bound))
        bump("D.bias_bound")
    # theta-form, both signs (non-residue dilation gives -S)
    for num in range(0, 26):
        th = Fraction(num, 100)
        e = int(th * m)  # floor
        if e > (m - 1) // 2:
            continue
        if (1 - 2 * th) * m * n >= d + m:
            bnd = m * n - 2 * th * (m * n - Fraction(d + m) / (1 - 2 * th))
            if abs(S) > bnd:
                fail("D.theta_form", (tag, p, I.A, float(th), S, float(bnd)))
            bump("D.theta_form")

def check_RHP(I):
    p, m, d = I.p, I.m, I.d
    for ep in range(0, (m - 1) // 2 + 1):
        B = [b for b in range(p) if I.eb[b] <= ep]
        n = len(B); r = sum(1 for b in B if b in I.negA)
        lhs = m * n - r
        for e in range(ep, (m - 1) // 2 + 1):
            rhs = Fraction((e + 1) * m * d, (e + 1 - ep) * (m - 2 * e)) + Fraction(2 * e * r, m - 2 * e)
            if lhs > rhs:
                fail("D.RHP_integer_form", (p, I.A, ep, e, lhs, rhs))
            bump("D.RHP_integer_form")
        eta = ep / m
        if eta <= 1 / 8:
            f = (1 - math.sqrt(2 * eta)) ** -2
            if lhs > f * d + m + 1e-9:
                fail("D.RHP_sqrt_form", (p, I.A, ep, lhs, f * d + m))
            bump("D.RHP_sqrt_form")

def thm25_factor(p, kappa):
    u = (1 + 2 * kappa) ** -0.5
    return 1 - (1 - u) ** 2 + (math.sqrt((0.5 + kappa) * p) + 1) / (2 * (p - 1))

def check_theorem25(p, trials, rng):
    """For A random/structured and every n, B = the n points maximising (resp. minimising)
    F_A(b) is the worst case of S for that A and n; check Theorem 2.5 on all of them."""
    worst = (0, None)
    cands = [("pair01", [0, 1])] + [("random", rng.sample(range(p), rng.randint(2, min(12, (p + 1) // 2)))) for _ in range(trials)]
    cands += [(t, A) for (t, A) in structured_sets(p) if len(A) <= 12]
    for tag, A in cands:
        I = Inst(p, A)
        Fv = [sum(I.chi[b]) for b in range(p)]
        order_hi = sorted(range(p), key=lambda b: -Fv[b])
        order_lo = sorted(range(p), key=lambda b: Fv[b])
        m = I.m
        for order in (order_hi, order_lo):
            S = 0
            for n in range(1, p + 1):
                S += Fv[order[n - 1]]
                mn = m * n
                if mn < 0.5 * p + 1e-9:
                    continue
                kappa = min(mn / p - 0.5, 1.5)
                if kappa <= 0:
                    continue
                fac = thm25_factor(p, kappa)
                ratio = abs(S) / mn
                if ratio > fac + 1e-12:
                    fail("D.theorem25", (p, tag, A, n, S, mn, fac))
                bump("D.theorem25_extremal_B")
                gap = fac - ratio
                if fac < 1 and (worst[1] is None or gap < worst[0]):
                    worst = (gap, {"p": p, "A": list(A), "n": n, "S": S, "mn": mn,
                                   "kappa": round(kappa, 4), "bound_factor": round(fac, 5), "ratio": round(ratio, 5)})
    return worst

def sharpness_pair(p):
    """Remark 2.7: A={0,1}, B={e_b=0}: |B|=(p+3)/4 (p=1 mod 4), S=|A||B|-2 = d; padding B with
    points of F_A=0 gives S/mn = d/mn."""
    I = Inst(p, [0, 1])
    B0 = [b for b in range(p) if I.eb[b] == 0]
    S0 = S_value(I, B0)
    zeros = [b for b in range(p) if sum(I.chi[b]) == 0]
    out = {"p": p, "|B0|": len(B0), "(p+3)/4": (p + 3) / 4, "S(B0)": S0, "mn": 2 * len(B0),
           "ratio": S0 / (2 * len(B0))}
    pads = []
    for kappa in (0.05, 0.1, 0.25):
        n = math.ceil((0.5 + kappa) * p / 2)
        if n - len(B0) <= len(zeros):
            B = B0 + zeros[:n - len(B0)]
            S = S_value(I, B)
            pads.append({"kappa": kappa, "n": n, "ratio": S / (2 * n), "1/(1+2kappa)": 1 / (1 + 2 * kappa),
                         "thm25_factor": thm25_factor(p, kappa)})
            if S / (2 * n) > thm25_factor(p, kappa):
                fail("D.sharpness_vs_thm25", (p, kappa))
            bump("D.sharpness_padding")
    out["padded"] = pads
    return out

# ---------------------------------------------------------------- E. no-go construction
def nogo_system(m, k=1):
    """Odd m. Spike: w0 points missing exactly one set (each k); base: every j-subset,
    j=(m-1)/2, with multiplicity mu; z empty points."""
    assert m % 2 == 1
    j = (m - 1) // 2
    w0 = k * math.comb(m, j)
    mu = m * (m - 3) * k
    d = w0 * (m - 1) ** 2 // 2
    assert w0 * (m - 1) ** 2 % 2 == 0
    npts = m * w0 + mu * math.comb(m, j)
    z = 2 * d + 1 - npts
    return dict(m=m, j=j, w0=w0, mu=mu, d=d, z=z, npts=npts)

def nogo_check(sysd):
    m, j, w0, mu, d, z = (sysd[x] for x in ("m", "j", "w0", "mu", "d", "z"))
    ok = z >= 0
    viol = []
    I = []
    for t in range(1, m + 1):
        It = w0 * (m - t) + (mu * math.comb(m - t, j - t) if t <= j else 0)
        I.append(It)
        if t == 1 and It != d:
            ok = False; viol.append(("size", It, d))
        if t * It > d:
            ok = False; viol.append((t, It, Fraction(d, t)))
    covered_m_minus_1 = m * w0
    return ok, viol, covered_m_minus_1, I

def nogo_bruteforce_small():
    """Materialise a smaller-m instance (m=7, not satisfying all t) to cross-check the
    symmetric intersection formula against explicit enumeration."""
    m = 7; sysd = nogo_system(m, 1)
    pts = []
    for k in range(m):
        pts += [frozenset(set(range(m)) - {k})] * sysd["w0"]
    for T in itertools.combinations(range(m), sysd["j"]):
        pts += [frozenset(T)] * sysd["mu"]
    pts += [frozenset()] * sysd["z"]
    _, _, _, Ifor = nogo_check(sysd)
    for t in range(1, m + 1):
        for T in itertools.combinations(range(m), t):
            cnt = sum(1 for P in pts if set(T) <= P)
            if cnt != Ifor[t - 1]:
                fail("E.formula", (T, cnt, Ifor[t - 1]))
            bump("E.symmetric_formula_vs_enumeration")
    if len(pts) != 2 * sysd["d"] + 1:
        fail("E.universe", len(pts))

# ---------------------------------------------------------------- driver
def structured_sets(p):
    out = []
    g = next(x for x in range(2, p) if all(pow(x, (p - 1) // q, p) != 1
             for q in set(f for f in range(2, p) if (p - 1) % f == 0 and all(f % r for r in range(2, int(f ** .5) + 1)))))
    for m in range(2, min(14, (p + 1) // 2) + 1):
        out.append(("interval", list(range(m))))
        if (p - 1) % m == 0:
            h = pow(g, (p - 1) // m, p)
            out.append(("subgroup", sorted({pow(h, i, p) for i in range(m)})))
        # greedy near-clique inside Q: points x with chi(x - y)=+1 for many chosen y
        Q = [x for x in range(1, p) if leg(x, p) == 1]
        A = [0]
        for x in Q:
            if len(A) >= m:
                break
            if all(leg(x - y, p) == 1 for y in A):
                A.append(x)
        if len(A) == m:
            out.append(("clique", A))
        out.append(("squares", sorted({(i * i) % p for i in range(1, 2 * m)})[:m]))
    return out

def main():
    random.seed(20260926)
    # A,B,C (per-instance) -------------------------------------------------
    slack = []
    # exhaustive: every A containing 0, p in {5,7,11,13}
    for p in [5, 7, 11, 13]:
        for m in range(1, (p + 1) // 2 + 1):
            for rest in itertools.combinations(range(1, p), m - 1):
                I = Inst(p, (0,) + rest)
                check_HP_derivatives(I)
                if m >= 2:
                    check_cheap(I)
                recs = check_theorem1(I)
                for (s, e) in recs:
                    slack.append((s, p, (0,) + rest, e))
                check_RHP(I)
                bump("exhaustive_sets")
    R["notes"]["exhaustive_primes"] = [5, 7, 11, 13]
    # p=17: exhaustive for m<=6
    p = 17
    for m in range(1, 7):
        for rest in itertools.combinations(range(1, p), m - 1):
            I = Inst(p, (0,) + rest)
            recs = check_theorem1(I)
            for (s, e) in recs:
                slack.append((s, p, (0,) + rest, e))
            check_RHP(I)
            bump("exhaustive_sets_p17_m_le6")
    print("exhaustive done", round(time.time() - T0, 1), "s", flush=True)
    # random + structured for larger p
    for p in [x for x in primes_upto(211) if x > 17]:
        insts = [("random", random.sample(range(p), random.randint(2, min(16, (p + 1) // 2)))) for _ in range(3)]
        insts += structured_sets(p)
        for tag, A in insts:
            I = Inst(p, A)
            if p <= 80:
                check_HP_derivatives(I)
            check_cheap(I)
            recs = check_theorem1(I, emax=3)
            for (s, e) in recs:
                slack.append((s, p, tuple(A), e))
            check_RHP(I)
            # bias checks on several B: extremal level sets and random subsets
            for ep in range(0, 3):
                B = [b for b in range(p) if I.eb[b] <= ep]
                check_bias(I, B, tag + f"_levelset{ep}")
            for _ in range(2):
                B = random.sample(range(p), random.randint(1, p - 1))
                check_bias(I, B, tag + "_randomB")
    print("random/structured done", round(time.time() - T0, 1), "s", flush=True)
    # full polynomial cross-check on a sample
    for p in [29, 37, 41, 53, 61]:
        for tag, A in structured_sets(p)[:6]:
            I = Inst(p, A)
            for e in range(0, min(2, (I.m - 1) // 2) + 1):
                H = full_H(I, e)
                if len(H) - 1 != (e + 1) * (I.d - e):
                    fail("C2.degree", (p, A, e, len(H) - 1))
                bump("C2.full_poly_degree")
                for b in range(p):
                    if I.eb[b] <= e:
                        o = poly_ord(H, b, p)
                        if o is None or o < ord_lower(I.m, e, I.eb[b], b in I.negA):
                            fail("C2.full_poly_order", (p, A, e, b, o))
                        bump("C2.full_poly_order")
    # index-k generalisation (Remark 5c): cubic and quartic residues
    for idx, plist in ((3, [31, 37, 43, 61, 67, 73]), (4, [29, 37, 41, 53, 61])):
        for p in plist:
            for tag, A in [("random", random.sample(range(p), random.randint(2, 9))) for _ in range(3)] + \
                          [("interval", list(range(7)))]:
                I = Inst(p, A, idx=idx)
                check_theorem1(I, emax=2)
                bump(f"C5.index{idx}_instances")
    check_leading_coefficient()
    print("leading coeff done", round(time.time() - T0, 1), "s", flush=True)
    # Theorem 2.5 on extremal B, and sharpness of the 1/2 threshold (Remark 2.7)
    w25 = []
    for p in [11, 13, 29, 61, 101, 197, 401, 1009, 2003]:
        w25.append(check_theorem25(p, 4 if p < 1000 else 2, random))
    R["witnesses"]["theorem25_smallest_gaps"] = [w[1] for w in w25]
    R["witnesses"]["sharpness_pairs"] = [sharpness_pair(p) for p in [101, 401, 1009]]
    for sp in R["witnesses"]["sharpness_pairs"]:
        if sp["|B0|"] != (sp["p"] + 3) // 4 or sp["S(B0)"] != (sp["p"] - 1) // 2:
            fail("D.sharpness_pair", sp)
        bump("D.sharpness_pair")
    print("thm 2.5 done", round(time.time() - T0, 1), "s", flush=True)
    slack.sort(key=lambda x: -x[0])
    R["witnesses"]["theorem1_tightest"] = [
        {"ratio_lhs_over_rhs": round(s, 6), "p": p, "A": list(A), "e": e} for (s, p, A, e) in slack[:15]]
    # cheap-bound constant lower bound
    for m in range(2, 120):
        for e in range(1, (m - 1) // 2 + 1):
            cc = cheap_constant(m, e)
            if float(cc) < e * math.exp(1 - e / m) * (1 - 1e-12):
                fail("B.cheap_lower", (m, e, float(cc)))
            bump("B.cheap_constant_lower_bound")
    R["witnesses"]["cheap_constants"] = {f"m={m},e={e}": float(cheap_constant(m, e))
                                         for (m, e) in [(40, 1), (40, 2), (40, 4), (100, 1), (100, 5), (100, 10), (1000, 1), (1000, 10)]}
    # E. no-go construction
    nogo = []
    for m in range(5, 61, 2):
        sysd = nogo_system(m, 1)
        ok, viol, cov, _ = nogo_check(sysd)
        nogo.append({"m": m, "all_constraints_hold": ok, "first_violation": str(viol[:1]),
                     "m*count/d": float(Fraction(m * cov, sysd["d"])), "d_digits": len(str(sysd["d"]))})
        bump("E.nogo_systems_checked")
        if not ok:
            fail("E.nogo", (m, viol[:2]))
    R["witnesses"]["nogo_systems"] = nogo
    nogo_bruteforce_small()
    # explicit prime instance of the no-go at m=43
    m = 43; k = 1
    while True:
        sysd = nogo_system(m, k)
        import sympy
        if sympy.isprime(2 * sysd["d"] + 1):
            break
        k += 1
    ok, viol, cov, _ = nogo_check(sysd)
    R["witnesses"]["nogo_prime_instance"] = {"m": m, "k": k, "d": str(sysd["d"]), "p=2d+1_prime": True,
                                             "constraints_hold": ok, "points_in_exactly_m-1_sets": str(cov),
                                             "m*count/d": float(Fraction(m * cov, sysd["d"]))}
    if not ok:
        fail("E.nogo_prime", viol[:2])
    bump("E.nogo_prime_instance")
    # F. prime sensitivity over F_{p^2}: A=F_p, m=p, D=(p^2-1)/2+p-1; (D)_j == 0 mod p for some j<=m-1
    Fq = []
    for p in [3, 5, 7, 11, 13]:
        q = p * p; D = (q - 1) // 2 + p - 1
        jz = next(j for j in range(1, p) if ff(D, j) % p == 0)
        Fq.append({"p": p, "D": D, "first_j_with_(D)_j=0_mod_p": jz, "m-1": p - 1})
        bump("F.subfield_normalisation_breaks")
    R["witnesses"]["F_p2"] = Fq
    R["elapsed_seconds"] = round(time.time() - T0, 1)
    R["cpu_seconds"] = round(time.process_time() - C0, 1)
    R["total_checks"] = sum(R["checks"].values())
    R["failures"] = len(FAIL)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(R, fh, indent=1)
    print(json.dumps({"total_checks": R["total_checks"], "failures": len(FAIL),
                      "elapsed": R["elapsed_seconds"]}, indent=1))
    if FAIL:
        print("FAILURES:", FAIL[:10]); sys.exit(1)

if __name__ == "__main__":
    main()
