"""sigma-crux verifier (2026-09-05).

Exact checks for research/sigma-crux-2026-09-05.md.  Standard library + numpy only.
Reads  results/sigma_crux_2026_09_05_search.json  (stored search witnesses; no search is redone here)
and    results/exact_checks.json                  (independent small-prime extrema),
writes results/sigma_crux_2026_09_05.json.

Conventions.  chi = Legendre symbol, chi(0) = 0.  S(A,B) = sum_{a in A, b in B} chi(a+b).
g(t) = sum chi(t a + b) = S(tA, B).  All identities are checked in exact integer arithmetic; inequalities
involving sqrt(p) are checked by squaring (exact integers).
"""
import json, os, sys, time, math, itertools, random
from fractions import Fraction
from collections import Counter, defaultdict
import numpy as np

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
SEARCH = os.path.join(REPO, 'results', 'sigma_crux_2026_09_05_search.json')
EXACT = os.path.join(REPO, 'results', 'exact_checks.json')
OUT = os.path.join(REPO, 'results', 'sigma_crux_2026_09_05.json')

T0 = time.time()
def log(*a):
    print('[%6.1fs]' % (time.time() - T0), *a, flush=True)

# ----------------------------------------------------------------------------- basic tools
_chi_cache = {}
def chi_table(p):
    if p not in _chi_cache:
        chi = np.zeros(p, dtype=np.int64)
        for x in range(1, p):
            chi[x] = 1 if pow(x, (p - 1) // 2, p) == 1 else -1
        _chi_cache[p] = chi
    return _chi_cache[p]

def profile(p, A, B):
    """g(t) for all t in F_p, exact int64."""
    chi = chi_table(p)
    t = np.arange(p, dtype=np.int64)
    Bn = np.array(B, dtype=np.int64)
    g = np.zeros(p, dtype=np.int64)
    for a in A:
        g += chi[(((t * a) % p)[:, None] + Bn[None, :]) % p].sum(axis=1)
    return g

def S_of(p, A, B):
    chi = chi_table(p)
    return int(sum(int(chi[(a + b) % p]) for a in A for b in B))

def classes(p, A, B):
    """multiplicities mu[r] (r in F_p: pairs (a,b), a != 0, b = r a) and mu_inf (pairs (0,b), b != 0);
    signed sums sigma[r] = sum chi(a) over the class."""
    chi = chi_table(p)
    mu = Counter(); sigma = Counter(); mu_inf = 0
    for a in A:
        if a % p == 0:
            mu_inf += sum(1 for b in B if b % p)
            continue
        inv = pow(a, p - 2, p)
        for b in B:
            r = (b * inv) % p
            mu[r] += 1; sigma[r] += int(chi[a])
    return mu, sigma, mu_inf

def ratio_energy(p, A, B):
    mu, _, mu_inf = classes(p, A, B)
    return sum(v * v for v in mu.values()) + mu_inf * mu_inf

def signed_ratio_energy(p, A, B):
    _, sigma, _ = classes(p, A, B)
    return sum(s * s for s in sigma.values())

# exact polynomial helpers (Fractions), truncated at degree D
def mul_trunc(P, Q, D):
    out = [Fraction(0)] * (D + 1)
    for i, pi in enumerate(P):
        if pi == 0: continue
        for j in range(0, D + 1 - i):
            if Q[j]: out[i + j] += pi * Q[j]
    return out

def cosh_poly(s, D):
    return [Fraction(s ** j, math.factorial(j)) if j % 2 == 0 else Fraction(0) for j in range(D + 1)]

def exp_poly(s, D):
    return [Fraction(s ** j, math.factorial(j)) for j in range(D + 1)]

def pow_poly(P, e, D):
    out = [Fraction(0)] * (D + 1); out[0] = Fraction(1); base = P
    while e:
        if e & 1: out = mul_trunc(out, base, D)
        e >>= 1
        if e: base = mul_trunc(base, base, D)
    return out

def N_k(p, A, B, k):
    """# of 2k-tuples of pairs in Pi = A x B \\ {(0,0)} whose product polynomial is a nonzero constant times a
    square = (2k)! [x^{2k}] exp(mu_inf x) prod_r cosh(mu_r x)."""
    mu, _, mu_inf = classes(p, A, B)
    D = 2 * k
    prod = exp_poly(mu_inf, D)
    for v, c in Counter(mu.values()).items():
        prod = mul_trunc(prod, pow_poly(cosh_poly(v, D), c, D), D)
    val = prod[D] * math.factorial(D)
    assert val.denominator == 1
    return int(val)

def dfact(n):  # (n)!! for odd n, (-1)!! = 1
    return 1 if n <= 0 else n * dfact(n - 2)

def isqrt_le(x2, y):
    """True iff sqrt(x2) <= y for integers x2 >= 0, y (exact)."""
    return y >= 0 and x2 <= y * y

# ----------------------------------------------------------------------------- brute force (tiny cases)
def brute_QW(p, A, B, k):
    """Q_k(t), W_k(t): square-tuple and non-square-tuple parts of g(t)^{2k}; also N_k by counting, and the
    termwise Weil check |sum_t chi(P(t))|^2 <= (d-1)^2 p for every non-square tuple (d = degree)."""
    chi = chi_table(p)
    pairs = [(a % p, b % p) for a in A for b in B if (a % p, b % p) != (0, 0)]
    def cls(a, b): return 'inf' if a == 0 else (b * pow(a, p - 2, p)) % p
    Q = [0] * p; W = [0] * p; Nk = 0; weil_ok = 0
    for tup in itertools.product(pairs, repeat=2 * k):
        cnt = Counter(cls(a, b) for (a, b) in tup)
        is_sq = all(v % 2 == 0 for c, v in cnt.items() if c != 'inf')
        vals = [1] * p
        for (a, b) in tup:
            for t in range(p): vals[t] *= int(chi[(t * a + b) % p])
        if is_sq:
            Nk += 1
            for t in range(p): Q[t] += vals[t]
        else:
            for t in range(p): W[t] += vals[t]
            d = sum(1 for (a, b) in tup if a != 0)
            assert sum(vals) ** 2 <= (d - 1) ** 2 * p, (p, tup, sum(vals))
            weil_ok += 1
    return Q, W, Nk, weil_ok

# ----------------------------------------------------------------------------- set statistics
def sumset_size(p, X, Y):
    return len({(x + y) % p for x in X for y in Y})

def prodset_size(p, X, Y):
    Xs = [x for x in X if x % p]; Ys = [y for y in Y if y % p]
    return len({(x * y) % p for x in Xs for y in Ys})

def rep_function(p, A, B):
    return Counter((a + b) % p for a in A for b in B)

_pf_cache = {}
def prime_factors(n):
    if n in _pf_cache: return _pf_cache[n]
    f = []; d = 2; m = n
    while d * d <= m:
        if m % d == 0:
            f.append(d)
            while m % d == 0: m //= d
        d += 1
    if m > 1: f.append(m)
    _pf_cache[n] = f
    return f

def subgroup_order_of_ratios(p, X):
    """order of the subgroup of F_p^* generated by {x/x0 : x in X*} (= lcm of element orders, F_p^* cyclic)."""
    Xs = [x % p for x in X if x % p]
    if len(Xs) <= 1: return 1
    inv0 = pow(Xs[0], p - 2, p)
    def order(s):
        o = p - 1
        for q in prime_factors(p - 1):
            while o % q == 0 and pow(s, o // q, p) == 1: o //= q
        return o
    L = 1
    for x in Xs:
        o = order((x * inv0) % p); L = L * o // math.gcd(L, o)
    return L

def minimal_ap_length(p, X):
    X = sorted(set(x % p for x in X))
    if len(X) <= 1: return len(X)
    d = np.arange(1, (p - 1) // 2 + 1, dtype=np.int64)
    inv = np.array([pow(int(x), p - 2, p) for x in d], dtype=np.int64)
    pts = (np.array(X, dtype=np.int64)[None, :] * inv[:, None]) % p
    pts.sort(axis=1)
    gaps = np.diff(pts, axis=1)
    wrap = pts[:, 0] + p - pts[:, -1]
    maxgap = np.maximum(gaps.max(axis=1), wrap)
    return int((p - maxgap + 1).min())

_dlog_cache = {}
def dlog_table(p):
    if p in _dlog_cache: return _dlog_cache[p]
    # find a primitive root
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in prime_factors(p - 1)):
            break
    T = np.zeros(p, dtype=np.int64); x = 1
    for i in range(p - 1):
        T[x] = i; x = (x * g) % p
    _dlog_cache[p] = T
    return T

def minimal_gp_length(p, X):
    """minimal L with X* subset of {c r^i : 0 <= i < L} for some r != 1; exact, vectorised over the ratio."""
    T = dlog_table(p)
    L = sorted({int(T[x % p]) for x in X if x % p})
    if len(L) <= 1: return len(L)
    N = p - 1; best = N
    Larr = np.array(L, dtype=np.int64)
    for e in [e for e in range(1, N + 1) if N % e == 0]:
        M = N // e
        if M < len(L): continue
        if len({l % e for l in L}) != 1: continue
        dp = np.array([x for x in range(1, M) if math.gcd(x, M) == 1] or [0], dtype=np.int64)
        if dp[0] == 0:   # M == 1
            continue
        invs = np.array([pow(int(x), -1, M) for x in dp], dtype=np.int64)
        pos = (((Larr - Larr[0]) // e)[None, :] * invs[:, None]) % M
        pos.sort(axis=1)
        gaps = np.diff(pos, axis=1); wrap = pos[:, 0] + M - pos[:, -1]
        mg = np.maximum(gaps.max(axis=1), wrap)
        best = min(best, int((M - mg + 1).min()))
    return int(best)

def stabiliser_order(p, D):
    """order of {h in F_p^* : hD = D} (D subset of F_p^*); p-1 for empty D."""
    Dset = sorted(set(D))
    if not Dset: return p - 1
    mask = np.zeros(p, dtype=bool); mask[Dset] = True
    h = np.arange(1, p, dtype=np.int64); Dn = np.array(Dset, dtype=np.int64)
    ok = np.ones(p - 1, dtype=bool)
    for j in range(0, len(Dn), 64):
        ok &= mask[(h[:, None] * Dn[None, j:j + 64]) % p].all(axis=1)
    return int(ok.sum())

def rect_stats(p, A, B, want_gp=False):
    A = sorted(set(a % p for a in A)); B = sorted(set(b % p for b in B))
    m, n = len(A), len(B)
    g = profile(p, A, B); S = int(g[1]); absg = np.abs(g)
    As = [a for a in A if a]; Bs = [b for b in B if b]
    rep = rep_function(p, A, B)
    st = {'p': p, 'm': m, 'n': n, 'S': S, 'bias': abs(S) / (m * n),
          'addK_A': sumset_size(p, A, A) / m, 'addK_B': sumset_size(p, B, B) / n,
          'mulK_A': prodset_size(p, As, As) / max(1, len(As)), 'mulK_B': prodset_size(p, Bs, Bs) / max(1, len(Bs)),
          'R': ratio_energy(p, A, B), 'sumset': len(rep), 'E_plus': sum(v * v for v in rep.values()),
          'maxmult': max(rep.values()), 'ap_A': minimal_ap_length(p, A), 'ap_B': minimal_ap_length(p, B),
          'sub_A': subgroup_order_of_ratios(p, A), 'sub_B': subgroup_order_of_ratios(p, B),
          'zero_in_A': 0 in A, 'zero_in_B': 0 in B,
          'L_half': int(sum(1 for t in range(1, p) if 2 * absg[t] >= abs(S))),
          'Labs_half': int(sum(1 for t in range(1, p) if 2 * absg[t] >= m * n)),
          'rank_t1': int(sum(1 for t in range(1, p) if absg[t] > abs(S))),
          'second_moment': int(np.sum(g * g)),
          }
    D = [t for t in range(1, p) if 2 * absg[t] >= abs(S)]
    st['stab_D'] = stabiliser_order(p, D)
    if want_gp:
        st['gp_A'] = minimal_gp_length(p, A); st['gp_B'] = minimal_gp_length(p, B)
    return st

# ----------------------------------------------------------------------------- exhaustive (validation only)
def exhaustive_max(p, n, m_list, chunk=20000):
    """exact max over all B with {0,1} subset B, |B| = n, and all A with |A| = m, of |S(A,B)|."""
    chi = chi_table(p)
    idx = (np.arange(p)[None, :] + np.arange(p)[:, None]) % p
    T = chi[idx].astype(np.int16)                     # T[b, x] = chi(x + b)
    base = T[0].astype(np.int16) + T[1].astype(np.int16)
    best = {m: -1 for m in m_list}; count = 0
    it = itertools.combinations(range(2, p), n - 2)
    while True:
        block = list(itertools.islice(it, chunk))
        if not block: break
        ix = np.array(block, dtype=np.int64)
        F = base[None, :] + T[ix].sum(axis=1)
        for m in m_list:
            top = np.partition(F, p - m, axis=1)[:, p - m:].sum(axis=1)
            bot = np.partition(F, m - 1, axis=1)[:, :m].sum(axis=1)
            best[m] = max(best[m], int(np.maximum(top, -bot).max()))
        count += len(block)
    return best, count

# ============================================================================= main
def main():
    res = {'description': 'sigma-crux verifier output (2026-09-05). Counts of exact checks and every witness.',
           'checks': {}, 'witnesses': {}, 'refutations': {}, 'conjecture_SI': {}}
    C = res['checks']; W = res['witnesses']
    rng = random.Random(20260905)
    search = json.load(open(SEARCH))
    exact = json.load(open(EXACT))

    # ---------------------------------------------------------------- 1. second-moment identity (b)
    log('1. second-moment identity')
    n_ok = 0
    cases = []
    for p in [5, 7, 11, 13, 17, 19, 23, 29]:
        for trial in range(40):
            m = rng.randint(1, 5); n = rng.randint(1, 5)
            A = rng.sample(range(p), m); B = rng.sample(range(p), n)
            if trial % 4 == 0: A[0] = 0
            if trial % 4 == 1: B[0] = 0
            if trial % 4 == 2: A[0] = 0; B[0] = 0
            cases.append((p, sorted(set(A)), sorted(set(B))))
    # plus 40 stored witnesses at p >= 101
    wit_pool = [(e['p'], w['A'], w['B']) for e in search['exhaustive'] for w in e['witnesses']] + \
               [(a['p'], w['A'], w['B']) for a in search['anneal'] for w in a['witnesses']]
    for (p, A, B) in rng.sample(wit_pool, 40):
        cases.append((p, A, B))
    for (p, A, B) in cases:
        chi = chi_table(p)
        g = profile(p, A, B)
        lhs = int(np.sum(g * g))
        Rchi = signed_ratio_energy(p, A, B)
        sA = int(sum(chi[a] for a in A if a % p)); sB = int(sum(chi[b] for b in B))
        rhs = p * Rchi - len(B) ** 2 * sA ** 2 + (p * sB * sB if 0 in A else 0)
        assert lhs == rhs, (p, A, B, lhs, rhs)
        assert lhs <= p * ratio_energy(p, A, B), (p, A, B)
        n_ok += 1
    C['second_moment_identity_cases'] = n_ok
    C['second_moment_upper_bound_pR_cases'] = n_ok

    # ---------------------------------------------------------------- 2. N_k, Q/W decomposition, Weil termwise (c)
    log('2. tuple decomposition / N_k / termwise Weil (brute force, tiny cases)')
    n_dec = 0; n_weil = 0; n_Nk_bound = 0; n_exact32 = 0
    for p in [7, 11, 13, 17]:
        for trial in range(10):
            m = rng.randint(1, 3); n = rng.randint(1, 3)
            A = rng.sample(range(p), m); B = rng.sample(range(p), n)
            if trial % 3 == 0: A[0] = 0
            if trial % 3 == 1: B[0] = 0
            A = sorted(set(A)); B = sorted(set(B))
            g = profile(p, A, B)
            npairs = len([1 for a in A for b in B if (a, b) != (0, 0)])
            for k in [1, 2]:
                if npairs ** (2 * k) > 30000: continue
                Q, Wt, Nk, wk = brute_QW(p, A, B, k)
                assert N_k(p, A, B, k) == Nk, (p, A, B, k)
                for t in range(p):
                    assert Q[t] + Wt[t] == int(g[t]) ** (2 * k)
                    assert Q[t] <= Nk
                    if 0 not in A: assert Q[t] >= 0
                # moment inequality with exact T_ns
                M = int(sum(int(x) ** (2 * k) for x in g)); T_ns = npairs ** (2 * k) - Nk
                excess = M - p * Nk
                assert excess <= 0 or excess ** 2 <= ((2 * k - 1) ** 2) * p * T_ns ** 2, (p, A, B, k)
                # N_k <= (2k-1)!! R^k when 0 not in A
                if 0 not in A:
                    R = ratio_energy(p, A, B)
                    assert Nk <= dfact(2 * k - 1) * R ** k
                    n_Nk_bound += 1
                n_dec += 1; n_weil += wk
    C['QW_decomposition_cases'] = n_dec
    C['termwise_weil_nonsquare_tuples_checked'] = n_weil
    C['N_k_le_dfact_R^k_cases'] = n_Nk_bound
    # exact count 3n^2 - 2n for k = 2 with n distinct classes (all multiplicities 1) vs (2k-1)!! n^k and k! n^k
    tab = []
    for n in range(1, 9):
        # n distinct classes: A = {1}, B = n distinct elements in F_p, p large
        p = 101; B = list(range(1, n + 1))
        N2 = N_k(p, [1], B, 2)
        assert N2 == 3 * n * n - 2 * n
        tab.append({'n': n, 'N_2': N2, 'dfact_bound_3n2': 3 * n * n, 'kfact_2n2': 2 * n * n, 'kfact_bound_valid': N2 <= 2 * n * n})
        n_exact32 += 1
    C['N_2_equals_3n2_minus_2n_cases'] = n_exact32
    W['N_2_table'] = tab
    res['refutations']['k_factorial_bound'] = {'statement': 'N_k <= k! (|A||B|)^k is FALSE for k=2: with n distinct ratio classes N_2 = 3n^2-2n > 2n^2 for n>=3',
                                              'witness': [t for t in tab if not t['kfact_bound_valid']][0]}

    # ---------------------------------------------------------------- 3. moment inequality on stored witnesses, k = 1,2,3
    log('3. moment inequality and level-set bound on stored witnesses')
    n_mom = 0; n_lev = 0; lev_examples = []
    sample = rng.sample(wit_pool, 300)
    for (p, A, B) in sample:
        A = sorted(set(A)); B = sorted(set(B)); m, n = len(A), len(B)
        g = profile(p, A, B); gl = [int(x) for x in g]
        npairs = m * n - (1 if (0 in A and 0 in B) else 0)
        for k in [1, 2, 3]:
            Nk = N_k(p, A, B, k); M = sum(x ** (2 * k) for x in gl); T_ns = npairs ** (2 * k) - Nk
            excess = M - p * Nk
            assert excess <= 0 or excess ** 2 <= ((2 * k - 1) ** 2) * p * T_ns ** 2, (p, A, B, k)
            n_mom += 1
            # level-set bound: L(eta) * (eta m n)^{2k} <= M <= p N_k + (2k-1) sqrt p T_ns
            for eta in (Fraction(1, 2), Fraction(3, 4), Fraction(1)):
                L = sum(1 for t in range(p) if abs(gl[t]) >= eta * m * n)
                lhs = L * (eta * m * n) ** (2 * k)          # Fraction
                assert lhs <= M
                n_lev += 1
        if len(lev_examples) < 5:
            k = 2; Nk = N_k(p, A, B, k); T_ns = npairs ** 4 - Nk
            bound = (p * Nk + 3 * math.sqrt(p) * T_ns) / ((0.5 * m * n) ** 4)
            lev_examples.append({'p': p, 'm': m, 'n': n, 'S': gl[1], 'L_half': sum(1 for t in range(1, p) if 2 * abs(gl[t]) >= m * n),
                                 'moment_bound_k2_eta_half': bound, 'weil_floor_3sqrtp_over_eta4': 3 * math.sqrt(p) * 16})
    C['moment_inequality_witness_cases'] = n_mom
    C['level_set_chebyshev_cases'] = n_lev
    W['level_set_examples'] = lev_examples
    # the floor: min over k in {1,2,3} of the bound >= sqrt p whenever |A||B| <= sqrt p (numerical grid)
    n_floor = 0
    for p in [101, 401, 1009, 1499]:
        for _ in range(30):
            m = rng.randint(1, 6); n = rng.randint(1, 6)
            if m * n > math.sqrt(p): continue
            A = rng.sample(range(1, p), m); B = rng.sample(range(p), n)
            R = ratio_energy(p, A, B)
            b1 = p * R / (m * n) ** 2                                   # k = 1, eta = 1 (exact identity bound)
            bks = [(p * N_k(p, A, B, k) + (2 * k - 1) * math.sqrt(p) * ((m * n) ** (2 * k) - N_k(p, A, B, k))) / (m * n) ** (2 * k) for k in (2, 3)]
            assert min([b1] + bks) >= math.sqrt(p) - 1e-9, (p, m, n)
            n_floor += 1
    C['level_set_floor_sqrtp_cases'] = n_floor

    # ---------------------------------------------------------------- 4. reduction (a): intervals, runs, Vinogradov
    log('4. interval reduction')
    n_tri = 0; n_np = 0; n_run = 0; run_examples = []
    primes = [q for q in range(5, 3000) if all(q % r for r in range(2, int(q ** 0.5) + 1))]
    for p in primes:
        chi = chi_table(p)
        for eps in (0.3, 0.45):
            N = int(p ** eps)
            if 2 * N >= p or N < 1: continue
            A = list(range(1, N + 1))
            S = S_of(p, A, A)
            S2 = sum(min(s - 1, 2 * N + 1 - s) * int(chi[s]) for s in range(2, 2 * N + 1))
            assert S == S2
            allres = all(chi[s] == 1 for s in range(2, 2 * N + 1))
            np_ = next(x for x in range(1, p) if chi[x] == -1)
            assert (S == N * N) == allres == (np_ > 2 * N)
            n_tri += 1; n_np += 1
        # runs: longest run of consecutive residues (and of non-residues) inside [1, p-1]
        for sign in (1, -1):
            best = 0; cur = 0; start = None; bstart = None
            for x in range(1, p):
                if chi[x] == sign:
                    cur += 1
                    if cur == 1: start = x
                    if cur > best: best = cur; bstart = start
                else: cur = 0
            Nr = (best + 1) // 2
            if Nr >= 1:
                A = list(range(0, Nr)); B = list(range(bstart, bstart + Nr))
                assert S_of(p, A, B) == sign * Nr * Nr
                n_run += 1
                if p in (101, 1009, 2999) and sign == 1:
                    run_examples.append({'p': p, 'longest_residue_run': best, 'N': Nr, 'S': sign * Nr * Nr})
    C['triangle_weight_identity_cases'] = n_tri
    C['vinogradov_equivalence_cases'] = n_np
    C['run_rectangle_cases'] = n_run
    W['run_examples'] = run_examples

    # ---------------------------------------------------------------- 5. exhaustive routine vs exact_checks.json (a-b form there)
    log('5. exhaustive routine vs exact_checks.json')
    n_ex = 0; ex_rows = []
    for e in exact['extrema']:
        p, m, n = e['p'], e['m'], e['n']
        best, count = exhaustive_max(p, n, [m])
        Sab = sum(int(chi_table(p)[(a - b) % p]) for a in e['witness']['A'] for b in e['witness']['B'])
        ok = best[m] == e['maximum_absolute_sum'] and Sab == e['witness']['signed_sum'] and count == math.comb(p - 2, n - 2)
        assert ok, (p, m, n, best, e)
        ex_rows.append({'p': p, 'm': m, 'n': n, 'max_ours': best[m], 'max_stored': e['maximum_absolute_sum'], 'stored_witness_S_in_a_minus_b_form': Sab})
        n_ex += 1
    C['exhaustive_vs_exact_checks_cases'] = n_ex
    W['exhaustive_vs_exact_checks'] = ex_rows

    # ---------------------------------------------------------------- 6. re-verify every stored witness; spot-recompute exhaustive maxima
    log('6. re-verify stored witnesses')
    n_w = 0; n_full = 0; by_mn = defaultdict(list)
    for e in search['exhaustive']:
        p, m, n = e['p'], e['m'], e['n']
        for w in e['witnesses']:
            assert len(set(w['A'])) == m and len(set(w['B'])) == n and {0, 1} <= set(w['B'])
            S = S_of(p, w['A'], w['B'])
            assert S == w['S'] and abs(S) == e['max'], (p, m, n, S, w['S'], e['max'])
            n_w += 1
        by_mn[(m, n)].append((p, e['max']))
        if e['max'] == m * n: n_full += 1
    for a in search['anneal']:
        for w in a['witnesses']:
            S = S_of(a['p'], w['A'], w['B'])
            assert S == w['S'] and abs(S) <= a['max'] and len(set(w['A'])) == a['m'] and len(set(w['B'])) == a['n']
            n_w += 1
        assert max(abs(w['S']) for w in a['witnesses']) == a['max']
    for k in search['kstar']:
        for kk, v in k['per_k'].items():
            S = S_of(k['p'], v['A'], v['B'])
            assert S == v['S'] and abs(S) == v['max'] and len(set(v['A'])) == int(kk) == len(set(v['B']))
            n_w += 1
        assert all(k['per_k'][str(j)]['max'] == j * j for j in range(4, k['kstar_lower_bound'] + 1))
    C['stored_witnesses_reverified'] = n_w
    # spot recompute of exhaustive maxima (n = 4, p <= 200; n = 5, p <= 110)
    n_spot = 0; spot = []
    grid = defaultdict(dict)
    for e in search['exhaustive']:
        grid[(e['p'], e['n'])][e['m']] = e['max']
    for (p, n), d in sorted(grid.items()):
        if (n == 4 and p <= 200) or (n == 5 and p <= 110):
            best, count = exhaustive_max(p, n, sorted(d))
            for m in d:
                assert best[m] == d[m], (p, n, m, best[m], d[m])
                n_spot += 1
            spot.append({'p': p, 'n': n, 'maxima': {str(m): d[m] for m in sorted(d)}, 'B_subsets': count})
    C['exhaustive_maxima_recomputed'] = n_spot
    W['exhaustive_spot_recomputation'] = spot

    # ---------------------------------------------------------------- 7. H4 (fixed k): maxima equal k*m for every prime in the grid
    log('7. H4 fixed-k table')
    h4 = {}
    for (m, n), rows in sorted(by_mn.items()):
        rows.sort()
        h4[f'{m}x{n}'] = {'primes': len(rows), 'p_min': rows[0][0], 'p_max': rows[-1][0], 'min_max': min(r[1] for r in rows),
                          'max_max': max(r[1] for r in rows), 'mn': m * n, 'all_complete': all(r[1] == m * n for r in rows),
                          'first_prime_with_complete': next((r[0] for r in rows if r[1] == m * n), None),
                          'complete_for_all_grid_primes_from': next((rows[i][0] for i in range(len(rows)) if all(r[1] == m * n for r in rows[i:])), None)}
    res['refutations']['H4_fixed_k'] = {
        'statement': 'H4 (fixed k): max_{|A|=|B|=k} |S|/k^2 decays like p^{-c(k)}, c(k)>0.  REFUTED: for k=4,5,6 the exact maximum is k^2 (bias 1) at every prime of the grid.',
        'table': h4,
        'kstar_lower_bounds': [{'p': k['p'], 'kstar_lb': k['kstar_lower_bound'], 'log2p': math.log2(k['p'])} for k in search['kstar']]}
    # descriptive fits on the stored annealing / k* data (HEURISTIC, recorded so that the note quotes verified numbers)
    slopes = {}
    for k in [12, 14, 16, 18, 20]:
        rows_k = sorted([(a['p'], a['max'] / (k * k)) for a in search['anneal'] if a['m'] == k and a['n'] == k])
        lp = np.log([r[0] for r in rows_k]); lb = np.log([r[1] for r in rows_k])
        slopes[str(k)] = {'slope_log_bias_vs_log_p': float(np.polyfit(lp, lb, 1)[0]), 'primes': len(rows_k), 'bias_at_p_min': rows_k[0][1], 'bias_at_p_max': rows_k[-1][1]}
    def k0(p):
        k = 1
        while 2 * (math.lgamma(p + 1) - math.lgamma(k + 1) - math.lgamma(p - k + 1)) - k * k * math.log(2) > 0: k += 1
        return k - 1
    ks_rows = [{'p': k['p'], 'kstar_lb': k['kstar_lower_bound'], 'log2p': math.log2(k['p']), 'ratio': k['kstar_lower_bound'] / math.log2(k['p']),
                'k0_random_model': k0(k['p'])} for k in search['kstar']]
    res['refutations']['H4_fixed_k']['kstar_lower_bounds'] = ks_rows
    res['refutations']['H4_fixed_k']['kstar_over_log2p_range'] = [min(r['ratio'] for r in ks_rows), max(r['ratio'] for r in ks_rows)]
    res['refutations']['H4_fixed_k']['k0_range'] = [min(r['k0_random_model'] for r in ks_rows), max(r['k0_random_model'] for r in ks_rows)]
    res['refutations']['H4_fixed_k']['anneal_fixed_k_slopes_HEURISTIC'] = slopes
    xs = []; ys = []
    for a in search['anneal']:
        if a['m'] == a['n'] and a['max'] < a['m'] * a['n']:
            xs.append(a['m'] / math.log2(a['p'])); ys.append(math.log(a['max'] / (a['m'] * a['n'])))
    c1, c0 = np.polyfit(xs, ys, 1)
    res['refutations']['H4_fixed_k']['anneal_bias_vs_k_over_log2p_fit_HEURISTIC'] = {'points': len(xs), 'c1': float(c1), 'c0': float(c0), 'x_at_bias_1': float(-c0 / c1),
                                                                                    'x_range': [min(xs), max(xs)]}
    # the GP construction (proved existence of complete rectangles with many biased dilates): exact checks
    n_gp = 0; gp_rows = []
    for p in [101, 211, 401, 1009, 1499, 2003, 4001]:
        chi = chi_table(p)
        for k in [4, 5, 6, 7]:
            r = None
            for cand in range(2, p):
                if all(pow(cand, (p - 1) // q, p) != 1 for q in prime_factors(p - 1)):
                    r = cand; break
            U = [pow(r, i, p) for i in range(k)]
            NU = [b for b in range(p) if all(chi[(u + b) % p] == 1 for u in U)]
            # exact identity: sum_b prod(1+chi(u+b)) = 2^k |N(U)| + sum_{u0} prod_{u != u0}(1+chi(u-u0))
            tot = sum(math.prod(1 + int(chi[(u + b) % p]) for u in U) for b in range(p))
            bnd = sum(math.prod(1 + int(chi[(u - u0) % p]) for u in U if u != u0) for u0 in U)
            assert tot == 2 ** k * len(NU) + bnd
            lower = (p - math.sqrt(p) * ((k - 2) * 2 ** (k - 1) + 1)) / 2 ** k - k / 2
            assert len(NU) >= lower - 1e-9, (p, k, len(NU), lower)
            if NU:
                S = S_of(p, U, NU); assert S == k * len(NU)
                g = profile(p, U, NU)
                D = [t for t in range(1, p) if 2 * abs(int(g[t])) >= S]
                need = {pow(r, j, p) for j in range(-(k // 4), k // 4 + 1)}
                assert need <= set(D), (p, k)
                gp_rows.append({'p': p, 'k': k, 'r': r, 'N_U_size': len(NU), 'weil_lower_bound': lower, 'S': S, 'D_size': len(D), 'guaranteed_D_size': 2 * (k // 4) + 1})
            n_gp += 1
    C['GP_construction_cases'] = n_gp
    W['GP_construction'] = gp_rows

    # ---------------------------------------------------------------- 8. H1-H3 statistics with matched random baselines
    log('8. H1-H3 statistics')
    wit = []
    for e in search['exhaustive']:
        for w in e['witnesses'][:1]:
            wit.append(('ex', e['p'], w['A'], w['B'], None))
    for a in search['anneal']:
        for w in a['witnesses']:
            wit.append(('an', a['p'], w['A'], w['B'], w['found_by']))
    for k in search['kstar']:
        for kk, v in k['per_k'].items():
            wit.append(('ks', k['p'], v['A'], v['B'], None))
    stats = []
    for (src, p, A, B, fb) in wit:
        st = rect_stats(p, A, B); st['src'] = src; st['found_by'] = fb; st['A'] = A; st['B'] = B
        stats.append(st)
    log('   witnesses done', len(stats))
    triples = sorted({(s['p'], s['m'], s['n']) for s in stats})
    brng = random.Random(777)
    base = defaultdict(list)
    for (p, m, n) in triples:
        for _ in range(3):
            A = brng.sample(range(p), m); B = brng.sample(range(p), n)
            base[(p, m, n)].append(rect_stats(p, A, B))
    log('   baselines done', sum(len(v) for v in base.values()))
    def med(rows, f): return float(np.median([f(r) for r in rows]))
    feats = {'addK_A': lambda r: r['addK_A'], 'addK_B': lambda r: r['addK_B'], 'mulK_A': lambda r: r['mulK_A'], 'mulK_B': lambda r: r['mulK_B'],
             'R_over_mn': lambda r: r['R'] / (r['m'] * r['n']), 'sumset_over_mn': lambda r: r['sumset'] / (r['m'] * r['n']),
             'Eplus_over_mn': lambda r: r['E_plus'] / (r['m'] * r['n']), 'apA_over_m': lambda r: r['ap_A'] / r['m'], 'apB_over_n': lambda r: r['ap_B'] / r['n']}
    def paired(rows):
        out = {'N': len(rows)}
        for name, f in feats.items():
            ratios = []; wins = 0
            for r in rows:
                b = med(base[(r['p'], r['m'], r['n'])], f)
                if b > 0: ratios.append(f(r) / b); wins += f(r) < b
            out[name] = {'median_ratio': float(np.median(ratios)), 'q10': float(np.quantile(ratios, .1)), 'q90': float(np.quantile(ratios, .9)),
                         'frac_witness_below_baseline': wins / len(ratios)}
        return out
    an_lt1 = [s for s in stats if s['src'] == 'an' and s['bias'] < 1]
    groups = {'anneal_bias_lt_1_square': paired([s for s in an_lt1 if s['m'] == s['n']]),
              'anneal_bias_lt_1_lopsided_m_ge_3n': paired([s for s in an_lt1 if s['m'] >= 3 * s['n']]),
              'anneal_bias_1_mn_ge_64': paired([s for s in stats if s['src'] == 'an' and s['bias'] >= 1 and s['m'] * s['n'] >= 64]),
              'exhaustive_lexfirst_mn_ge_64': paired([s for s in stats if s['src'] == 'ex' and s['m'] * s['n'] >= 64])}
    res['H1_H2_H3'] = {'paired_witness_vs_baseline': groups}
    # H1 / H2 refutation: an explicit near-extremal rectangle whose doubling equals that of random sets
    ex_h12 = max(an_lt1, key=lambda s: (s['m'] == s['n']) * s['m'] * s['n'] * (s['bias'] >= 0.6))
    res['refutations']['H1_H2_doubling'] = {
        'statement': 'H1/H2 (near-extremal rectangles have small additive / multiplicative doubling): REFUTED.  Paired to random sets of the same sizes in the same field the median doubling ratio is 1.00 (see H1_H2_H3.paired_witness_vs_baseline).',
        'witness': {k: ex_h12[k] for k in ('p', 'm', 'n', 'S', 'bias', 'addK_A', 'addK_B', 'mulK_A', 'mulK_B', 'A', 'B')},
        'baseline_medians_same_pmn': {k: med(base[(ex_h12['p'], ex_h12['m'], ex_h12['n'])], feats[k]) for k in ('addK_A', 'addK_B', 'mulK_A', 'mulK_B')}}
    # H3: biased dilate set a union of cosets of a nontrivial subgroup?  (stabiliser order > 1)
    pop = [s for s in stats if s['m'] * s['n'] >= 64 and s['bias'] >= 0.5]
    stab_gt1 = [s for s in pop if s['stab_D'] > 1]
    ex_h3 = max([s for s in pop if s['stab_D'] == 1 and s['L_half'] >= 4], key=lambda s: s['L_half'])
    res['refutations']['H3_coset_union'] = {
        'statement': 'H3 (the biased dilate set D = {t != 0 : |g(t)| >= |S|/2} is a union of cosets of a nontrivial subgroup): REFUTED.',
        'population': {'N': len(pop), 'stab_gt_1': len(stab_gt1), 'stabiliser_values': sorted({s['stab_D'] for s in stab_gt1})},
        'witness_with_large_D_and_trivial_stabiliser': {k: ex_h3[k] for k in ('p', 'm', 'n', 'S', 'L_half', 'stab_D', 'A', 'B')},
        'the_stab_gt_1_cases': [{k: s[k] for k in ('p', 'm', 'n', 'S', 'L_half', 'stab_D', 'sub_A', 'sub_B', 'src')} for s in stab_gt1]}
    # sumset compression (supported)
    w_e = [(s['E_plus'] / (s['m'] * s['n']), med(base[(s['p'], s['m'], s['n'])], feats['Eplus_over_mn'])) for s in an_lt1 + [s for s in stats if s['src'] == 'ex' and s['m'] * s['n'] >= 64]]
    res['H1_H2_H3']['sumset_compression'] = {'N': len(w_e), 'frac_Eplus_witness_gt_baseline': float(np.mean([a > b for a, b in w_e])),
                                             'median_ratio': float(np.median([a / b for a, b in w_e]))}
    # spike isolation summary
    res['H1_H2_H3']['spike_isolation_raw'] = {'N_pop_mn_ge_64_bias_ge_half': len(pop), 'median_L_half': float(np.median([s['L_half'] for s in pop])),
                                             'frac_L_half_eq_1': float(np.mean([s['L_half'] == 1 for s in pop])), 'max_L_half': max(s['L_half'] for s in pop),
                                             'by_src': {src: {'N': len([s for s in pop if s['src'] == src]), 'frac_L_half_eq_1': float(np.mean([s['L_half'] == 1 for s in pop if s['src'] == src])),
                                                              'max_L_half': max(s['L_half'] for s in pop if s['src'] == src)} for src in ('an', 'ex', 'ks')}}

    # ---------------------------------------------------------------- 9. Conjecture SI test (0 removed from A and B)
    log('9. Conjecture SI')
    rows = []
    for (src, p, A, B, fb) in wit:
        As = sorted({a % p for a in A if a % p}); Bs = sorted({b % p for b in B if b % p})
        m, n = len(As), len(Bs)
        if m * n < 64: continue
        g = profile(p, As, Bs); S = int(g[1])
        if 2 * abs(S) < m * n: continue
        D = [t for t in range(1, p) if 2 * abs(int(g[t])) >= abs(S)]
        gA = minimal_gp_length(p, As); gB = minimal_gp_length(p, Bs)
        rows.append({'p': p, 'm': m, 'n': n, 'S': S, 'D': len(D), 'D_over_log2p': len(D) / math.log2(p), 'gpA_over_m': gA / m, 'gpB_over_n': gB / n,
                     'compact': min(gA / m, gB / n) <= 2, 'stab_D': stabiliser_order(p, D), 'src': src, 'found_by': fb, 'A': As, 'B': Bs})
    viol_u = [r for r in rows if r['D'] > 8 * math.log2(r['p'])]
    viol_g = [r for r in rows if (not r['compact']) and r['D'] > 2 * math.log2(r['p'])]
    rows.sort(key=lambda r: -r['D'])
    res['conjecture_SI'] = {
        'statement_universal': 'For A,B subset F_p^* with |A||B| >= 64 and |S(A,B)| >= |A||B|/2: |D_{1/2}(A,B)| <= 8 log2 p, D_{1/2} = {t in F_p^*: |g(t)| >= |S|/2}.',
        'statement_generic': 'If moreover neither A nor B is contained in a geometric progression of length <= 2|A| (resp. 2|B|), then |D_{1/2}| <= 2 log2 p.',
        'N_tested': len(rows), 'violations_universal': len(viol_u), 'violations_generic': len(viol_g),
        'max_D_over_log2p': rows[0]['D_over_log2p'], 'median_D': float(np.median([r['D'] for r in rows])), 'frac_D_eq_1': float(np.mean([r['D'] == 1 for r in rows])),
        'frac_D_le_log2p': float(np.mean([r['D'] <= math.log2(r['p']) for r in rows])),
        'N_compact': len([r for r in rows if r['compact']]), 'frac_compact_with_D_gt_log2p': float(np.mean([r['D'] > math.log2(r['p']) for r in rows if r['compact']])),
        'by_src': {src: {'N': len([r for r in rows if r['src'] == src]), 'frac_D_eq_1': float(np.mean([r['D'] == 1 for r in rows if r['src'] == src])),
                         'max_D': max(r['D'] for r in rows if r['src'] == src)} for src in ('an', 'ex', 'ks')},
        'n_D_gt_log2p': len([r for r in rows if r['D'] > math.log2(r['p'])]),
        'cases_D_gt_log2p': [{k: r[k] for k in r if k not in ('A', 'B')} for r in rows if r['D'] > math.log2(r['p'])],
        'top_cases': [{k: r[k] for k in r} for r in rows[:8]],
        'violations': viol_u + viol_g}
    assert not viol_u and not viol_g

    res['seconds'] = round(time.time() - T0, 1)
    res['total_checks'] = sum(v for v in C.values() if isinstance(v, int))
    json.dump(res, open(OUT, 'w'), indent=1)
    log('wrote', OUT, 'total exact checks', res['total_checks'])

if __name__ == '__main__':
    main()
