#!/usr/bin/env python3
"""
Verifier for the sigma 'biclique' direction (2026-09-05).

Complete bicliques of the Paley relation: pairs A, B ⊆ F_p with A + B ⊆ Q ∪ {0} (Q = nonzero squares).
For A ⊆ F_p put B(A) = {b : a + b ∈ Q ∪ {0} for all a ∈ A}.  Quantities:
  omega(p)      exact clique number of the Paley graph G_p (p ≡ 1 mod 4);
  M_k(p)        max_{|A| = k} |B(A)|  (profile of complete bicliques);
  b(p)          max min(|A|, |B(A)|)   (balanced biclique number);
  P(p)          max |A||B| over complete bicliques with |A|, |B| ≥ 2.

What this script does (standard library + numpy/sympy only; target < 10 minutes):
  1. recomputes omega(p) exactly in Python (bitset branch and bound with greedy colouring, cliques normalised
     to contain 0 and 1) for every prime p ≡ 1 (mod 4) up to OMEGA_PY_LIMIT, and checks the clique witnesses of
     the C++ helper (experiments/sigma_biclique_2026_09_05.cpp, table stored in
     results/sigma_biclique_2026_09_05_search.json) for every prime in that table; also checks the 21 prime
     entries of Brouwer's table (Shearer/Exoo) for p ≤ 197, and computes the sum-clique number s(p) exactly;
  2. recomputes M_k(p) exactly in Python for small k / p and checks every biclique witness of the table;
  3. verifies |B({0,1})| = (p+3)/4 (so P(p) ≥ (p+3)/2) and re-runs, in Python, the Hanson–Petridis-reduced
     search showing P(p) = (p+3)/2 for p ≤ P3_PY_LIMIT; brute-force (no reduction) confirmation for p ≤ 100;
  4. checks the elementary inequalities used in the note on every extremal biclique of the table:
     the Weil-count bound for |B(A)|, the Hanson–Petridis bound with the clique refinement, the Weil-amplified
     moment inequality, and the multiplicative-energy chain of §5 (T ≥ n^4/|SS|, T ≥ T_diag, and the exact
     comparison of the two lower bounds);
  5. recomputes the geometric-progression profile gp_k(p) and the Theorem-6.2-type lower bound.
Writes results/sigma_biclique_2026_09_05.json with all counts and witnesses.
"""
import json, os, sys, time, math
from fractions import Fraction
from itertools import combinations
from sympy import primerange

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
OUT = os.path.join(ROOT, 'results', 'sigma_biclique_2026_09_05.json')
SEARCH = os.path.join(ROOT, 'results', 'sigma_biclique_2026_09_05_search.json')

OMEGA_PY_LIMIT = 1300     # exact Python recomputation of omega up to here
PROFILE_LIMITS = {2: 3000, 3: 1000, 4: 1000, 5: 400, 6: 200}   # exact Python recomputation of M_k up to p ≤ limit
P3_PY_LIMIT = 500         # Python re-run of the HP-reduced product check
BRUTE_PRODUCT_LIMIT = 100 # brute force (no reduction) product maximum
BALANCED_PY_LIMIT = 120   # exact Python balanced biclique number
GP_LIMIT, GP_KMAX = 1000, 10

T0 = time.time()
CHECKS = {'count': 0, 'failures': []}
def check(name, cond, witness=None):
    CHECKS['count'] += 1
    if not cond:
        CHECKS['failures'].append({'name': name, 'witness': witness})
        print('FAIL', name, witness, flush=True)

# ----------------------------------------------------------------------------------------------------------
# basic tables
# ----------------------------------------------------------------------------------------------------------
def chi_table(p):
    chi = [-1] * p; chi[0] = 0
    for x in range(1, p): chi[x * x % p] = 1
    return chi

def tables(p):
    chi = chi_table(p)
    adj = [0] * p; Np = [0] * p
    for x in range(p):
        a = 0; n = 0
        for y in range(p):
            if chi[(x - y) % p] == 1: a |= 1 << y
            s = (x + y) % p
            if s == 0 or chi[s] == 1: n |= 1 << y
        adj[x] = a; Np[x] = n
    return chi, adj, Np

def bits(m):
    out = []
    while m:
        v = (m & -m).bit_length() - 1; out.append(v); m &= m - 1
    return out

def B_of(A, Np, p):
    m = (1 << p) - 1
    for a in A: m &= Np[a]
    return m

def least_nonresidue(chi):
    n0 = 2
    while chi[n0] != -1: n0 += 1
    return n0

# ----------------------------------------------------------------------------------------------------------
# 1. exact clique number (BBMC with greedy colouring), cliques normalised to contain 0 and 1
# ----------------------------------------------------------------------------------------------------------
def omega_exact(p, chi, adj):
    best = [2, [0, 1]]; nodes = [0]
    P = 0
    for x in range(2, p):
        if chi[x] == 1 and chi[x - 1] == 1: P |= 1 << x
    sys.setrecursionlimit(10000)
    def expand(P, cur):
        nodes[0] += 1
        order = []; color = []; U = P; col = 0
        while U:
            col += 1; Qc = U
            while Qc:
                v = (Qc & -Qc).bit_length() - 1
                Qc &= ~(1 << v); U &= ~(1 << v); Qc &= ~adj[v]
                order.append(v); color.append(col)
        for i in range(len(order) - 1, -1, -1):
            if len(cur) + color[i] <= best[0]: return
            v = order[i]; NP = P & adj[v]
            cur.append(v)
            if NP == 0:
                if len(cur) > best[0]: best[0] = len(cur); best[1] = list(cur)
            else: expand(NP, cur)
            cur.pop(); P &= ~(1 << v)
    expand(P, [0, 1])
    return best[0], best[1], nodes[0]

def is_clique(C, chi, p):
    return all(chi[(a - b) % p] == 1 for a, b in combinations(C, 2))

# sum-clique number s(p) = max |A| with A + A ⊆ Q ∪ {0}  (then (A, A) is a complete biclique, so s(p) ≤ b(p))
def is_sumclique(A, chi, p):
    return all(((a + b) % p == 0 or chi[(a + b) % p] == 1) for a in A for b in A)

def sumclique_exact(p, chi):
    n0 = least_nonresidue(chi)
    V = [x for x in range(p) if x == 0 or chi[2 * x % p] == 1]
    sadj = [0] * p
    for x in V:
        m = 0
        for y in V:
            if y != x and ((x + y) % p == 0 or chi[(x + y) % p] == 1): m |= 1 << y
        sadj[x] = m
    v0 = 1 if chi[2] == 1 else n0      # scaling by squares: a sum-clique with a nonzero element contains v0 wlog
    best = [1, [v0]]; nodes = [0]
    def expand(P, cur):
        nodes[0] += 1
        order = []; color = []; U = P; col = 0
        while U:
            col += 1; Qc = U
            while Qc:
                v = (Qc & -Qc).bit_length() - 1
                Qc &= ~(1 << v); U &= ~(1 << v); Qc &= ~sadj[v]
                order.append(v); color.append(col)
        for i in range(len(order) - 1, -1, -1):
            if len(cur) + color[i] <= best[0]: return
            v = order[i]; NP = P & sadj[v]
            cur.append(v)
            if NP == 0:
                if len(cur) > best[0]: best[0] = len(cur); best[1] = list(cur)
            else: expand(NP, cur)
            cur.pop(); P &= ~(1 << v)
    expand(sadj[v0], [v0])
    return best[0], best[1], nodes[0]

def sumclique_bruteforce(p, chi):
    V = [x for x in range(p) if x == 0 or chi[2 * x % p] == 1]
    best = 0
    for k in range(1, 12):
        if any(is_sumclique(A, chi, p) for A in combinations(V, k)): best = k
        else: break
    return best

# ----------------------------------------------------------------------------------------------------------
# 2. exact profile M_k(p) = max_{|A|=k} |B(A)|  (A normalised: 0 ∈ A and 1 ∈ A or n0 ∈ A)
# ----------------------------------------------------------------------------------------------------------
def profile_exact(p, chi, Np, k, lb=0):
    n0 = least_nonresidue(chi)
    best = [lb, None]; nodes = [0]
    def dfs(j, last, B, e, A):
        nodes[0] += 1
        cnt = B.bit_count()
        if j == k:
            if cnt > best[0]: best[0] = cnt; best[1] = list(A)
            return
        if cnt <= best[0]: return
        r = k - j
        cand = []
        for a in range(last + 1, p):
            if a == e: continue
            v = (B & Np[a]).bit_count()
            if v > best[0]: cand.append((a, v))
        if len(cand) < r: return
        vs = sorted((v for _, v in cand), reverse=True)
        if vs[r - 1] <= best[0]: return
        for a, v in cand:
            if v <= best[0]: continue
            A.append(a); dfs(j + 1, a, B & Np[a], e, A); A.pop()
    for e in (1, n0):
        B = Np[0] & Np[e]
        if k == 2:
            if B.bit_count() > best[0]: best[0] = B.bit_count(); best[1] = [0, e]
            continue
        dfs(2, 0, B, e, [0, e])
    return best[0], best[1], nodes[0]

# balanced biclique number b(p) = max min(|A|, |B(A)|), exact
def balanced_exact(p, chi, Np, lb=0):
    n0 = least_nonresidue(chi)
    best = [lb, None]
    def dfs(j, last, B, e, A):
        cnt = B.bit_count()
        if min(j, cnt) > best[0]: best[0] = min(j, cnt); best[1] = list(A)
        if cnt <= best[0]: return
        cand = []
        for a in range(last + 1, p):
            if a == e: continue
            v = (B & Np[a]).bit_count()
            if v > best[0]: cand.append((a, v))
        if not cand: return
        vs = sorted((v for _, v in cand), reverse=True)
        bound = max(min(j + r, vs[r - 1]) for r in range(1, len(vs) + 1))
        if bound <= best[0]: return
        for a, v in cand:
            if v <= best[0]: continue
            A.append(a); dfs(j + 1, a, B & Np[a], e, A); A.pop()
    for e in (1, n0):
        dfs(2, 0, Np[0] & Np[e], e, [0, e])
    return best[0], best[1]

# ----------------------------------------------------------------------------------------------------------
# 3. product maximum
# ----------------------------------------------------------------------------------------------------------
def product_bruteforce(p, chi, Np):
    """max |A||B(A)| over |A| ≥ 2, |B(A)| ≥ 2, A normalised; exact branch and bound without the HP reduction."""
    n0 = least_nonresidue(chi)
    best = [0, None]
    def dfs(j, last, B, e, A):
        cnt = B.bit_count()
        if cnt >= 2 and j * cnt > best[0]: best[0] = j * cnt; best[1] = list(A)
        if cnt < 2: return
        cand = []
        for a in range(last + 1, p):
            if a == e: continue
            v = (B & Np[a]).bit_count()
            if v >= 2: cand.append((a, v))
        if not cand: return
        vs = sorted((v for _, v in cand), reverse=True)
        bound = max((j + r) * vs[r - 1] for r in range(1, len(vs) + 1))
        if bound <= best[0]: return
        for a, v in cand:
            cb = max((j + r) * min(v, vs[r - 1]) for r in range(1, len(vs) + 1))
            if cb <= best[0]: continue
            A.append(a); dfs(j + 1, a, B & Np[a], e, A); A.pop()
    for e in (1, n0):
        dfs(2, 0, Np[0] & Np[e], e, [0, e])
    return best[0], best[1]

def product3_check(p, chi, Np):
    """Hanson–Petridis-reduced search.  A complete biclique with |A| ≤ |B| and |A||B| ≥ (p+5)/2 has
    |B ∩ (-A)| ≥ 3 (HP Theorem 1.2 with d = (p-1)/2), C = A ∩ (-B) is a clique, and after normalising three of
    its elements to 0, 1, c ∈ A (so 0, -1, -c ∈ B) one has A ⊆ CA = Np[0]∩Np[-1]∩Np[-c], B ⊆ Np[0]∩Np[1]∩Np[c],
    |B| ≥ ceil(sqrt((p+5)/2)), |A| ≤ |B|.  Returns (violator_found, best_product_found, triangles, nodes)."""
    thr = (p + 3) // 2
    minB = 1
    while minB * minB < (p + 5) // 2: minB += 1
    best = [thr, None]; nodes = [0]; ntri = 0
    def dfs(idx, B, j, CAlist, A):
        nodes[0] += 1
        cnt = B.bit_count()
        if cnt < minB or j > cnt: return
        if cnt >= 3 and j * cnt > best[0]: best[0] = j * cnt; best[1] = list(A)
        cand = []
        for i in range(idx, len(CAlist)):
            v = (B & Np[CAlist[i]]).bit_count()
            if v >= minB: cand.append((i, v))
        if not cand: return
        vs = sorted((v for _, v in cand), reverse=True)
        bound = max((j + r) * vs[r - 1] for r in range(1, len(vs) + 1))
        if bound <= best[0]: return
        for i, v in cand:
            cb = max((j + r) * min(v, vs[r - 1]) for r in range(1, len(vs) + 1))
            if cb <= best[0]: continue
            A.append(CAlist[i]); dfs(i + 1, B & Np[CAlist[i]], j + 1, CAlist, A); A.pop()
    for c in range(2, p):
        if not (chi[c] == 1 and chi[c - 1] == 1): continue
        ntri += 1
        CA = Np[0] & Np[p - 1] & Np[p - c]
        CA &= ~(1 << 0); CA &= ~(1 << 1); CA &= ~(1 << c)
        CB = Np[0] & Np[1] & Np[c]
        dfs(0, CB, 3, bits(CA), [0, 1, c])
    return best[0] > thr, best[0], best[1], ntri, nodes[0]

# ----------------------------------------------------------------------------------------------------------
# 4. inequalities on extremal bicliques
# ----------------------------------------------------------------------------------------------------------
def weil_count_bound(p, k):
    """Upper bound (real number) on #{b : chi(b + a) = 1 for all a ∈ A}, |A| = k, from Weil:
       2^{-k} [ p - C(k,2) + sum_{j>=3} C(k,j) (j-1) sqrt(p) ]  (the b ∈ -A terms of the product sum are ≥ 0)."""
    s = p - math.comb(k, 2) + sum(math.comb(k, j) * (j - 1) for j in range(3, k + 1)) * math.sqrt(p)
    return s / 2 ** k

def le_with_sqrt(lhs_int, c0, c1, p):
    """exact test lhs ≤ c0 + c1*sqrt(p) for integers lhs, c0 and a nonnegative rational c1"""
    d = Fraction(lhs_int) - Fraction(c0)
    if d <= 0: return True
    return d * d <= Fraction(c1) * Fraction(c1) * p

def moment_data(p, chi, A, B, kmax=3):
    """Exact Weil-amplified moment check for A + B ⊆ Q (B1 = B minus (-A) so that chi(a+b) = 1 for all pairs):
       |A||B'|^{2k} ≤ S_k := sum_x F(x)^{2k}, F(x) = sum_{b∈B'} chi(x+b), and
       S_k ≤ p N_k + (|B'|^{2k} - N_k) (2k-1) sqrt(p), N_k = (2k)! [x^{2k}] cosh(x)^{|B'|} ≤ (2k-1)!! |B'|^k."""
    Bp = [b for b in B if all(chi[(a + b) % p] == 1 for a in A)]
    n = len(Bp)
    F = [sum(chi[(x + b) % p] for b in Bp) for x in range(p)]
    out = []
    # cosh(x)^n truncated: coefficients c_m of x^m / m!  (only even m): number of n-tuples... use exact ints:
    # N_k = number of 2k-tuples over B' with all multiplicities even = (2k)! [x^{2k}] cosh(x)^n
    # compute via exponential generating function with Fractions
    from math import factorial
    egf = [Fraction(1)] + [Fraction(0)] * (2 * kmax)
    cosh = [Fraction(1 if m % 2 == 0 else 0, factorial(m)) for m in range(2 * kmax + 1)]
    for _ in range(n):
        new = [Fraction(0)] * (2 * kmax + 1)
        for i, ci in enumerate(egf):
            if ci == 0: continue
            for j in range(0, 2 * kmax + 1 - i):
                new[i + j] += ci * cosh[j]
        egf = new
    for k in range(1, kmax + 1):
        Sk = sum(f ** (2 * k) for f in F)
        Nk = egf[2 * k] * factorial(2 * k)
        assert Nk.denominator == 1
        Nk = int(Nk)
        dfact = 1
        for t in range(1, 2 * k, 2): dfact *= t
        lhs_ok = len(A) * n ** (2 * k) <= Sk
        weil_ok = le_with_sqrt(Sk - p * Nk, 0, (n ** (2 * k) - Nk) * (2 * k - 1), p)
        matching_ok = Nk <= dfact * n ** k
        out.append({'k': k, 'S_k': Sk, 'N_k': Nk, 'A_n2k': len(A) * n ** (2 * k), 'lhs_ok': lhs_ok, 'weil_ok': weil_ok,
                    'matching_ok': matching_ok, 'Bprime_size': n})
    return out

def energy_data(p, A, B):
    """Restrict to B1 = B minus (-A) so that S = A + B1 is contained in Q.  n = |A||B1|, r(s) = representation
    counts, T = number of solutions of (a+b)(a2+b2) = (a3+b3)(a4+b4) = sum_lambda rho(lambda)^2 with
    rho(lambda) = sum_t r(lambda t) r(t), sum_lambda rho(lambda) = n^2, rho(1) = E^+(A,B1); the ratio set S/S lies
    in Q so T >= n^4/|S/S| >= 2n^4/(p-1); T_diag = 2E^2 - sum r^4 (solutions with {s,t} = {s2,t2} as multisets,
    weighted) <= T, and E^2 <= T_diag."""
    Bp = [b for b in B if (-b) % p not in set(A)]
    r = [0] * p
    for a in A:
        for b in Bp: r[(a + b) % p] += 1
    assert r[0] == 0
    S = [s for s in range(1, p) if r[s]]
    n = len(A) * len(Bp)
    rho = [0] * p
    for s in S:
        for t in S: rho[s * pow(t, p - 2, p) % p] += r[s] * r[t]
    T = sum(x * x for x in rho)
    ratio_size = sum(1 for x in rho if x)
    Eplus = rho[1]
    r4 = sum(x ** 4 for x in r)
    Tdiag = 2 * Eplus * Eplus - r4
    maxrho = max(rho); argmax = [l for l in range(p) if rho[l] == maxrho]
    Qb = Fraction(2 * n ** 4, p - 1)
    chi = chi_table(p)
    return {'n': n, 'Bprime_size': len(Bp), 'T': T, 'ratio_set_size': ratio_size, 'E_plus': Eplus, 'T_diag': Tdiag,
            'T_over_Qbound': float(Fraction(T * (p - 1), 2 * n ** 4)) if n else None,
            'T_minus_Tdiag_over_Qbound': float(Fraction((T - Tdiag) * (p - 1), 2 * n ** 4)) if n else None,
            'S_size': len(S), 'sum_rho_is_n2': sum(rho) == n * n, 'Eplus_eq_sum_r2': Eplus == sum(x * x for x in r),
            'Q_bound': Qb, 'Q_bound_le_T': Qb <= T, 'ratio_bound_le_T': (Fraction(n ** 4, ratio_size) <= T) if ratio_size else True,
            'Tdiag_le_T': Tdiag <= T, 'E2_le_Tdiag': Eplus * Eplus <= Tdiag,
            'Q_bound_le_Tdiag': Qb <= Tdiag, 'max_rho': maxrho, 'argmax_rho': argmax[:5],
            'ratio_set_in_Q': all(chi[l] == 1 for l in range(1, p) if rho[l])}

def ap_pair_energy(p, N):
    """A = B = {1..N} in F_p (2N < p): exact E^+ and T; the obstruction U(N,N) >= T >= E^2 >= (2N^3/3 - N^2)^2."""
    A = list(range(1, N + 1))
    r = [0] * p
    for a in A:
        for b in A: r[(a + b) % p] += 1
    S = [s for s in range(1, p) if r[s]]
    rho = [0] * p
    for s in S:
        for t in S: rho[s * pow(t, p - 2, p) % p] += r[s] * r[t]
    T = sum(x * x for x in rho); E = rho[1]
    return {'N': N, 'n': N * N, 'E_plus': E, 'T': T, 'E_formula': (N - 1) * N * (2 * N - 1) // 3 + N * N,
            'E_ge_2N3over3_minus_N2': 3 * E >= 2 * N ** 3 - 3 * N * N, 'T_ge_E2': T >= E * E}

def hp_bound_check(p, A, B, omega):
    """Hanson–Petridis Thm 1.2 with d=(p-1)/2: |A||B| ≤ (p-1)/2 + |B ∩ (-A)|, and the refinement |B ∩ (-A)| ≤ omega."""
    r = sum(1 for b in B if (-b) % p in set(A))
    return {'r': r, 'hp_ok': len(A) * len(B) <= (p - 1) // 2 + r, 'r_le_omega': r <= omega}

# ----------------------------------------------------------------------------------------------------------
# 5. geometric progressions
# ----------------------------------------------------------------------------------------------------------
def gp_profile(p, chi, Np, kmax):
    n0 = least_nonresidue(chi)
    full = (1 << p) - 1
    res = {}
    for k in range(2, kmax + 1):
        best = (-1, None, None)
        for r in range(2, p):
            U = []; x = 1; ok = True
            for i in range(k):
                if x in U: ok = False; break
                U.append(x); x = x * r % p
            if not ok: continue
            for c in (1, n0):
                B = full
                for u in U: B &= Np[c * u % p]
                v = B.bit_count()
                if v > best[0]: best = (v, r, c)
        res[k] = best
    return res

def primitive_root(p):
    from sympy import primitive_root as pr
    return pr(p)

# ----------------------------------------------------------------------------------------------------------
def main():
    search = json.load(open(SEARCH)) if os.path.exists(SEARCH) else None
    results = {'checks': CHECKS, 'omega_python': {}, 'profile_python': {}, 'product': {}, 'balanced_python': {},
               'inequalities': {}, 'gp': {}, 'table': {}, 'notes': []}
    brouwer = {5: 2, 13: 3, 17: 3, 29: 4, 37: 4, 41: 5, 53: 5, 61: 5, 73: 5, 89: 5, 97: 6, 101: 5, 109: 6, 113: 7,
               137: 7, 149: 7, 157: 7, 173: 8, 181: 7, 193: 7, 197: 8}   # Brouwer's table (Shearer; 173, 197 Exoo)
    primes = [p for p in primerange(5, 3001) if p % 4 == 1]
    ctab = {int(k): v for k, v in search['table'].items()} if search else {}

    # ---- 1. omega ----
    t = time.time()
    for p in primes:
        chi, adj, Np = tables(p) if p <= OMEGA_PY_LIMIT else (chi_table(p), None, None)
        if p <= OMEGA_PY_LIMIT:
            om, C, nd = omega_exact(p, chi, adj)
            check('python clique witness is a clique', is_clique(C, chi, p) and len(C) == om, {'p': p, 'C': C})
            results['omega_python'][p] = {'omega': om, 'clique': C, 'nodes': nd}
            if p in brouwer: check('omega matches Brouwer/Shearer table', om == brouwer[p], {'p': p, 'omega': om})
            if p in ctab: check('omega Python == omega C++', om == ctab[p]['omega'], {'p': p, 'py': om, 'c': ctab[p]['omega']})
        if p <= OMEGA_PY_LIMIT:
            s, SA, snd = sumclique_exact(p, chi)
            results['omega_python'][p]['sumclique'] = s; results['omega_python'][p]['sumclique_A'] = SA
            check('python sum-clique witness is a sum-clique', is_sumclique(SA, chi, p) and len(SA) == s, {'p': p, 'A': SA})
            check('s(p) ≤ (1+sqrt(2p-1))/2', s <= (1 + math.sqrt(2 * p - 1)) / 2, {'p': p, 's': s})
            check('HP-type: s^2 ≤ (p-1)/2 + |A ∩ -A|', s * s <= (p - 1) // 2 + sum(1 for a in SA if (-a) % p in set(SA)), {'p': p, 'A': SA})
            if p <= 61: check('s(p) brute force agrees', sumclique_bruteforce(p, chi) == s, {'p': p, 's': s})
            if p in ctab and ctab[p].get('sumclique') is not None:
                check('sumclique Python == C++', s == ctab[p]['sumclique'], {'p': p, 'py': s, 'c': ctab[p]['sumclique']})
        if p in ctab:
            C = ctab[p]['clique']
            check('C++ clique witness is a clique of size omega', is_clique(C, chi, p) and len(C) == ctab[p]['omega'], {'p': p})
            if ctab[p].get('sumclique_A'):
                check('C++ sum-clique witness is a sum-clique of size s', is_sumclique(ctab[p]['sumclique_A'], chi, p) and len(ctab[p]['sumclique_A']) == ctab[p]['sumclique'], {'p': p})
                if ctab[p]['balanced']['exact'] and ctab[p]['balanced']['b'] is not None:
                    check('s(p) ≤ b(p) (exact b)', ctab[p]['sumclique'] <= ctab[p]['balanced']['b'], {'p': p})
            check('omega ≤ Hanson–Petridis (sqrt(2p-1)+1)/2', ctab[p]['omega'] <= (math.sqrt(2 * p - 1) + 1) / 2, {'p': p})
    results['timing_omega'] = time.time() - t
    print('omega done', round(time.time() - t, 1), 's', flush=True)

    # ---- 2. profile and biclique witnesses ----
    t = time.time()
    for p in primes:
        if p > 1000 and p not in ctab: continue
        chi, adj, Np = tables(p)
        # every witness in the C table
        if p in ctab:
            for key in ('profile', 'balanced', 'greedy', 'gp'):
                for w in ctab[p].get(key + '_witnesses', []):
                    A = w['A']; B = B_of(A, Np, p); m = B.bit_count()
                    check('witness |B(A)| as recorded (%s)' % key, m == w['B_size'] and len(set(A)) == len(A), {'p': p, 'A': A, 'rec': w['B_size'], 'got': m})
                    check('witness is a complete biclique', all((chi[(a + b) % p] == 1 or (a + b) % p == 0) for a in A for b in bits(B)), {'p': p, 'A': A})
        # construction |B({0,1})| = (p+3)/4
        m01 = (Np[0] & Np[1]).bit_count()
        check('|B({0,1})| = (p+3)/4', 4 * m01 == p + 3, {'p': p, 'm': m01})
        # exact Python profile for small k / p
        prof = {}
        for k, lim in PROFILE_LIMITS.items():
            if p <= lim:
                M, A, nd = profile_exact(p, chi, Np, k)
                prof[k] = {'M': M, 'A': A, 'nodes': nd}
                if p in ctab and str(k) in ctab[p]['M'] and ctab[p]['M'][str(k)]['exact']:
                    check('M_k Python == M_k C++', M == ctab[p]['M'][str(k)]['M'], {'p': p, 'k': k, 'py': M, 'c': ctab[p]['M'][str(k)]['M']})
                # Weil-count bound |B(A)| ≤ k + 2^{-k}[p - C(k,2) + sum_{j≥3} C(k,j)(j-1) sqrt p]
                c0 = Fraction(p - math.comb(k, 2), 2 ** k) + k
                c1 = Fraction(sum(math.comb(k, j) * (j - 1) for j in range(3, k + 1)), 2 ** k)
                check('M_k ≤ Weil count bound', le_with_sqrt(M, c0, c1, p), {'p': p, 'k': k, 'M': M})
        results['profile_python'][p] = prof
        if p <= 1000 and 2 in prof: check('M_2 = (p+3)/4', 4 * prof[2]['M'] == p + 3, {'p': p})
    results['timing_profile'] = time.time() - t
    print('profile done', round(time.time() - t, 1), 's', flush=True)

    # ---- 3. product maximum ----
    t = time.time()
    for p in primes:
        if p > P3_PY_LIMIT: break
        chi, adj, Np = tables(p)
        found, bestv, A, ntri, nd = product3_check(p, chi, Np)
        # P(p) = max((p+3)/2, best product of a biclique containing a normalised triangle pair) -- exact by the HP reduction
        Pp = max((p + 3) // 2, bestv)
        results['product'][p] = {'hp_reduced_violator': found, 'best': bestv, 'P': Pp, 'A': A, 'triangles': ntri, 'nodes': nd}
        if found:
            B = bits(B_of(A, Np, p)); r = sum(1 for b in B if (-b) % p in set(A))
            check('violator witness: product as recorded, |B ∩ (-A)| ≥ 3, and HP equality-type', len(A) * len(B) == bestv and r >= 3 and bestv <= (p - 1) // 2 + r, {'p': p, 'A': A, 'B': B})
        if p in ctab:
            check('product3 C++ agrees (flag and value)', ctab[p]['product3']['violator_found'] == found and ctab[p]['product3']['best'] == max(bestv, (p + 3) // 2), {'p': p, 'py': bestv, 'c': ctab[p]['product3']['best']})
        if p <= BRUTE_PRODUCT_LIMIT:
            pr, A = product_bruteforce(p, chi, Np)
            results['product'][p]['bruteforce_max'] = pr; results['product'][p]['bruteforce_A'] = A
            check('brute-force product max = max((p+3)/2, HP-reduced best)', pr == Pp, {'p': p, 'pr': pr, 'A': A, 'P': Pp})
        if p <= BALANCED_PY_LIMIT:
            b, A = balanced_exact(p, chi, Np)
            results['balanced_python'][p] = {'b': b, 'A': A, 'B': bits(B_of(A, Np, p))}
            if p in ctab and ctab[p]['balanced']['exact']:
                check('balanced Python == balanced C++', b == ctab[p]['balanced']['b'], {'p': p, 'py': b, 'c': ctab[p]['balanced']['b']})
            if p in ctab: check('b(p) ≥ omega(p)', b >= ctab[p]['omega'], {'p': p})
    results['timing_product'] = time.time() - t
    print('product done', round(time.time() - t, 1), 's', flush=True)

    # ---- 4. inequalities on extremal bicliques of the C table ----
    t = time.time()
    ineq = {}
    for p in primes:
        if p not in ctab or p > 1000: continue
        chi, adj, Np = tables(p)
        om = ctab[p]['omega']
        rows = []
        cand_sets = [(w['A'], 'M_%d' % w['k']) for w in ctab[p]['profile_witnesses']]
        cand_sets.append((ctab[p]['clique'], 'clique'))
        if ctab[p]['balanced'].get('A'): cand_sets.append((ctab[p]['balanced']['A'], 'balanced'))
        for A, label in cand_sets:
            B = bits(B_of(A, Np, p))
            hp = hp_bound_check(p, A, B, om)
            check('HP bound on extremal biclique', hp['hp_ok'], {'p': p, 'A': A})
            check('|B ∩ (-A)| ≤ omega', hp['r_le_omega'], {'p': p, 'A': A})
            en = energy_data(p, A, B)
            check('sum_lambda rho = n^2 and rho(1) = E^+', en['sum_rho_is_n2'] and en['Eplus_eq_sum_r2'], {'p': p, 'A': A})
            check('ratio set (A+B)/(A+B) ⊆ Q', en['ratio_set_in_Q'], {'p': p, 'A': A})
            check('T ≥ n^4/|S/S|', en['ratio_bound_le_T'], {'p': p, 'A': A})
            check('T ≥ 2n^4/(p-1)', en['Q_bound_le_T'], {'p': p, 'A': A})
            check('T ≥ T_diag ≥ E^2', en['Tdiag_le_T'] and en['E2_le_Tdiag'], {'p': p, 'A': A})
            n = en['n']
            # the ratio-set lower bound 2n^4/(p-1) is implied by the diagonal count whenever 2n^2 ≤ p-1
            if 2 * n * n <= p - 1:
                check('2n^4/(p-1) ≤ T_diag when 2n^2 ≤ p-1', en['Q_bound_le_Tdiag'], {'p': p, 'A': A, 'n': n})
            check('rho(1) = E^+ ≥ n', en['E_plus'] >= n, {'p': p, 'A': A})
            mo = moment_data(p, chi, A, B, kmax=3) if len(A) * len(B) <= 400 else []
            for m in mo:
                check('Weil moment: |A||B|^{2k} ≤ S_k', m['lhs_ok'], {'p': p, 'A': A, 'k': m['k']})
                check('Weil moment: S_k ≤ pN_k + (n^{2k}-N_k)(2k-1)sqrt p', m['weil_ok'], {'p': p, 'A': A, 'k': m['k']})
                check('N_k ≤ (2k-1)!! n^k', m['matching_ok'], {'p': p, 'k': m['k']})
            rows.append({'label': label, 'A': A, 'B_size': len(B), 'hp': hp,
                         'energy': {k: (str(v) if isinstance(v, Fraction) else v) for k, v in en.items()},
                         'moments': mo})
        ineq[p] = rows
    ap = []
    for (pp, N) in [(1009, 10), (1009, 20), (2003, 30), (3001, 38)]:
        d = ap_pair_energy(pp, N); ap.append(d)
        check('AP pair: E^+ formula and E ≥ 2N^3/3 - N^2', d['E_plus'] == d['E_formula'] and d['E_ge_2N3over3_minus_N2'], d)
        check('AP pair: T ≥ E^2', d['T_ge_E2'], d)
    results['ap_pair_obstruction'] = ap
    results['inequalities'] = ineq
    results['timing_inequalities'] = time.time() - t
    print('inequalities done', round(time.time() - t, 1), 's', flush=True)

    # ---- 5. geometric progressions ----
    t = time.time()
    for p in primes:
        if p > GP_LIMIT: break
        chi, adj, Np = tables(p)
        g = gp_profile(p, chi, Np, GP_KMAX)
        results['gp'][p] = {k: {'gp': v[0], 'r': v[1], 'c': v[2]} for k, v in g.items()}
        if p in ctab:
            for k in range(2, GP_KMAX + 1):
                if str(k) in ctab[p]['gp']:
                    check('gp_k Python == gp_k C++', g[k][0] == ctab[p]['gp'][str(k)]['gp'], {'p': p, 'k': k})
                if str(k) in ctab[p]['M'] and ctab[p]['M'][str(k)]['exact']:
                    check('gp_k ≤ M_k', g[k][0] <= ctab[p]['M'][str(k)]['M'], {'p': p, 'k': k})
        # Theorem 6.2-type lower bound for U = {r^i : 0 ≤ i < k}, r primitive root
        r = int(primitive_root(p)); full = (1 << p) - 1
        for k in range(4, 8):
            U = [pow(r, i, p) for i in range(k)]; B = full
            for u in U: B &= Np[u]
            # N(U) counts b with chi(u+b) = 1 for all u; B(U) additionally allows b = -u; so |B(U)| ≥ |N(U)|
            NU = sum(1 for b in range(p) if all(chi[(u + b) % p] == 1 for u in U))
            lb = (p - math.sqrt(p) * ((k - 2) * 2 ** (k - 1) + 1)) / 2 ** k - k / 2
            check('GP lower bound |N(U_k)| ≥ 2^{-k}[p - sqrt p((k-2)2^{k-1}+1)] - k/2', NU >= lb - 1e-9, {'p': p, 'k': k, 'NU': NU, 'lb': lb})
            check('|B(U_k)| ≥ |N(U_k)|', B.bit_count() >= NU, {'p': p, 'k': k})
    results['timing_gp'] = time.time() - t
    print('gp done', round(time.time() - t, 1), 's', flush=True)

    # ---- Weil exclusion: for which min-sides m does Prop 2.5 alone exclude |A||B| > (p+3)/2 ? ----
    weil_excl = {}
    for p in primes:
        ok = []
        for m in range(3, 41):
            c0 = Fraction(p - math.comb(m, 2), 2 ** m) + m
            c1 = Fraction(sum(math.comb(m, j) * (j - 1) for j in range(3, m + 1)), 2 ** m)
            # m * (c0 + c1 sqrt p) < (p+3)/2  <=>  m*c0 + m*c1*sqrt(p) < (p+3)/2
            lhs = Fraction(p + 3, 2) - m * c0
            excl = lhs > 0 and (m * c1) ** 2 * p < lhs * lhs
            ok.append(excl)
        first_fail = next((m for m, e in zip(range(3, 41), ok) if not e), None)
        weil_excl[p] = {'m_star': (first_fail - 1) if first_fail else 40, 'excluded_m': [m for m, e in zip(range(3, 41), ok) if e]}
        if p >= 116: check('Weil exclusion covers m=3 for p ≥ 116', 3 in weil_excl[p]['excluded_m'], {'p': p})
        if p >= 213: check('Weil exclusion covers m=4 for p ≥ 213', 4 in weil_excl[p]['excluded_m'], {'p': p})
        if p == 113: check('Weil exclusion fails at m=3, p=113', 3 not in weil_excl[p]['excluded_m'], {'p': p})
    results['weil_exclusion'] = weil_excl

    # ---- table with comparisons ----
    for p in primes:
        row = {'p': p, 'sqrt_p_over_2': round(math.sqrt(p / 2), 3), 'p_over_2': p / 2,
               'graham_ringrose_log_p_logloglog_p': round(math.log(p) * math.log(math.log(math.log(p))), 3),
               'log2_p': round(math.log2(p), 3)}
        if p in ctab:
            # root fix 2026-09-05: table rows may lack optional keys when a C++ run did not finish
            c = ctab[p]
            row.update({'omega': c.get('omega'), 'B_of_clique': c.get('B_of_clique'), 'sumclique': c.get('sumclique'),
                        'M': c.get('M'), 'balanced': c.get('balanced'),
                        'greedy_balanced': (c.get('greedy') or {}).get('b'),
                        'product3_violator': (c.get('product3') or {}).get('violator_found'), 'gp': c.get('gp')})
        if p in results['omega_python']: row['omega_python'] = results['omega_python'][p]['omega']
        results['table'][p] = row

    results['summary'] = {'checks': CHECKS['count'], 'failures': len(CHECKS['failures']), 'time_s': time.time() - T0,
                          'omega_python_limit': OMEGA_PY_LIMIT, 'profile_limits': PROFILE_LIMITS,
                          'p3_python_limit': P3_PY_LIMIT, 'brute_product_limit': BRUTE_PRODUCT_LIMIT,
                          'balanced_python_limit': BALANCED_PY_LIMIT, 'search_table_present': search is not None}
    json.dump(results, open(OUT, 'w'), indent=1, default=str)
    print('checks', CHECKS['count'], 'failures', len(CHECKS['failures']), 'time', round(time.time() - T0, 1), 's')

if __name__ == '__main__':
    main()
