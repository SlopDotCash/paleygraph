#!/usr/bin/env python3
"""
sigma pass, direction lit2026 (2026-09-05): exact small-prime tests of the key
lemmas behind the 2024-2026 literature surveyed in
research/sigma-lit2026-2026-09-05.md.

No source found in the survey claims a proof of the two-set Paley conjecture,
a clique bound p^{1/2-c}, or a sub-sqrt(p) biclique bound, so there is no
"square-root-breaking" lemma to refute.  Instead this script tests, by exact
integer arithmetic, every lemma that the survey identifies as load-bearing in
the genuinely new 2025-2026 results, plus the elementary link between the
Paley RIP and the clique number, plus the best cited clique bound:

 T1  Hanson-Petridis Theorem 1.2 (Proc. LMS 2021; local copy
     sources/sigma-hanson-petridis-1905.09134.txt):
         A+B <= mu_d u {0}  ==>  |A||B| <= d + |B n (-A)|.
     Exhaustive over all A (WLOG 0 in A by translation) with B maximal.
 T2  Kalmynin (arXiv:2504.10202v2) Lemma 4: for every d-critical pair
     (equality in T1) the Hanson-Petridis polynomial factors as
         HP(x;A,d) = C prod_{b in B} (x-b)^{alpha-eps(b)}.
 T3  Kalmynin Theorem 3 (Sarkozy's conjecture, all p): Q_p != A+B with
     |A|,|B|>1; Theorem 2: mu_d = A+B, |A|,|B|>1 ==> |A|=|B|=sqrt(d);
     Theorem 1 (Lev-Sonn): A-A = mu_d u {0} ==> d in {2,6}.
 T4  Yip-Yoo (arXiv:2607.25711): proper H <= F_p^*, |H|>=7 ==> H != A (+)^ A
     (restricted sumset); exceptional sizes only 1,3,6.
 T5  Satake (arXiv:2405.08608v2) Lemma 28: (K,delta)-RIP of the Paley ETF
     ==> |sum_{u,v in U} chi(u-v)| <= delta sqrt(p) |U| for |U|<=K; and the
     elementary corollary omega(G_p) <= delta_K sqrt(p) + 1 whenever
     omega(G_p) <= K.  Gram matrix built exactly from Satake 2011.02907 Def. 3.
 T6  Hanson-Petridis Corollary 1.5: omega(G_p) <= (sqrt(2p-1)+1)/2, checked
     against exact clique numbers for p = 1 mod 4, p <= P_CLIQUE.
 T7  Zhang (arXiv:2506.02434v6) statement: sum_{a<=(p-1)/2} (a/p) > 0 for
     p = 3 mod 4 (Dirichlet's theorem; the statement, not the proof, is tested).
 T8  Wang-Shen-Kobzar Theorem 3.3 at t=1: L_1(G_p) = theta(G_p) = sqrt(p),
     via the exact Lovasz/Hoffman eigenvalue formula for Paley graphs.

Standard library + numpy only.  Writes results/sigma_lit2026_2026_09_05.json.
"""
import json, math, os, sys, time, itertools
import numpy as np

T0 = time.time()
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(os.path.dirname(HERE), 'results', 'sigma_lit2026_2026_09_05.json')

def primes_upto(n):
    s = bytearray([1]) * (n + 1); s[0] = s[1] = 0
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(n + 1) if s[i]]

def legendre_table(p):
    chi = [0] * p
    for x in range(1, p):
        chi[x] = -1
    for x in range(1, p):
        chi[(x * x) % p] = 1
    return chi

def primitive_root(p):
    if p == 2: return 1
    fac = []
    m = p - 1; q = 2
    while q * q <= m:
        if m % q == 0:
            fac.append(q)
            while m % q == 0: m //= q
        q += 1
    if m > 1: fac.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fac): return g
    raise RuntimeError

def subgroup(p, d):
    g = primitive_root(p)
    h = pow(g, (p - 1) // d, p)
    H = set(); x = 1
    for _ in range(d):
        H.add(x); x = (x * h) % p
    assert len(H) == d
    return H

def mask(S):
    m = 0
    for s in S: m |= 1 << s
    return m

def bits(m):
    out = []
    while m:
        b = m & -m; out.append(b.bit_length() - 1); m ^= b
    return out

def popcount(m): return bin(m).count('1')

def shift_mask(m, c, p):
    """mask of S + c (mod p) for S given by mask m."""
    r = 0
    for s in bits(m): r |= 1 << ((s + c) % p)
    return r

def neg_mask(m, p):
    r = 0
    for s in bits(m): r |= 1 << ((-s) % p)
    return r

results = {'date': '2026-09-05', 'direction': 'lit2026', 'checks': {}, 'witnesses': {}, 'timing_s': {}}
checks = results['checks']; wit = results['witnesses']

# ---------------------------------------------------------------- T1 + T2 + T3
def dfs_pairs(p, S, want_sum_equal=None):
    """Enumerate all A containing 0 with B_max(A) = intersect_{a in A} (S - a)
    of size >= 2 (plus the |B_max| = 1 boundary), depth-first, elements
    increasing.  Yields (Amask, Bmask).  S is a mask."""
    Sm = S
    stack = [(1, Sm, 1)]          # A = {0}, B_max = S - 0 = S, next element 1
    # B_max for A: intersection over a in A of (S - a)
    shifted = [shift_mask(Sm, -a, p) for a in range(p)]
    while stack:
        A, B, nxt = stack.pop()
        yield A, B
        if popcount(B) <= 1:
            continue
        for a in range(nxt, p):
            B2 = B & shifted[a]
            if popcount(B2) >= 1:
                stack.append((A | (1 << a), B2, a + 1))

def sumset_mask(A, B, p):
    r = 0
    for a in bits(A): r |= shift_mask(B, a, p)
    return r

def hp_polynomial(p, A, d):
    """Kalmynin's HP(x;A,d) = sum_a c_a (x+a)^{d+alpha-1} - 1, c_a = 1/prod_{a'!=a}(a-a'),
    as list of coefficients (low to high) mod p."""
    alpha = len(A)
    n = d + alpha - 1
    coef = [0] * (n + 1)
    # binomials mod p (n < p is NOT guaranteed; use exact integers)
    binom = [math.comb(n, k) for k in range(n + 1)]
    for a in A:
        den = 1
        for a2 in A:
            if a2 != a: den = den * ((a - a2) % p) % p
        c = pow(den, p - 2, p)
        # (x+a)^n = sum_k C(n,k) a^{n-k} x^k
        apow = [1] * (n + 1)
        for k in range(1, n + 1): apow[k] = apow[k - 1] * a % p
        for k in range(n + 1):
            coef[k] = (coef[k] + c * (binom[k] % p) * apow[n - k]) % p
    coef[0] = (coef[0] - 1) % p
    return coef

def poly_mul(f, g, p):
    r = [0] * (len(f) + len(g) - 1)
    for i, fi in enumerate(f):
        if fi == 0: continue
        for j, gj in enumerate(g):
            r[i + j] = (r[i + j] + fi * gj) % p
    return r

def strip(f):
    while len(f) > 1 and f[-1] == 0: f = f[:-1]
    return f

def run_T1_T2_T3(P_ALL_D=31, P_QCASE=53, time_cap=150.0):
    t0 = time.time()
    n_hp = 0; hp_max_excess = -10**9; critical = []; hp_viol = []
    n_lemma4 = 0; lemma4_fail = []
    n_sark = 0; sark_wit = []
    n_thm2 = 0; thm2_wit = []; thm2_bad = []
    n_thm1 = 0; thm1_wit = []; thm1_bad = []
    n_binom_zero = 0
    primes_done = []
    for p in primes_upto(P_QCASE):
        if p < 5: continue
        if time.time() - t0 > time_cap: break
        ds = [d for d in range(2, p - 1) if (p - 1) % d == 0] if p <= P_ALL_D else [(p - 1) // 2]
        for d in ds:
            H = subgroup(p, d)
            Hm = mask(H); S = Hm | 1          # mu_d u {0}
            # ---- T1/T2/T3(Thm 2) over pairs with A+B <= S (sumset problems use S=mu_d)
            for A, B in dfs_pairs(p, S):
                Ab = bits(A); Bb = bits(B)
                if len(Bb) == 0: continue
                # T1: HP Theorem 1.2 at maximal B (value is strictly increasing in B)
                negA = neg_mask(A, p)
                val = len(Ab) * len(Bb) - popcount(B & negA) - d
                n_hp += 1
                if val > hp_max_excess: hp_max_excess = val
                if val > 0:
                    hp_viol.append({'p': p, 'd': d, 'A': Ab, 'B': Bb, 'excess': val})
                if val == 0 and len(Ab) > 1 and len(Bb) > 1:
                    # d-critical pair (A, B): test Kalmynin Lemma 4 exactly
                    alpha = len(Ab)
                    hp = strip(hp_polynomial(p, Ab, d))
                    prod = [1]
                    for b in Bb:
                        eps = 1 if ((-b) % p) in set(Ab) else 0
                        for _ in range(alpha - eps):
                            prod = poly_mul(prod, [(-b) % p, 1], p)
                    prod = strip(prod)
                    C = math.comb(d + alpha - 1, alpha - 1) % p
                    if C == 0: n_binom_zero += 1
                    target = strip([(C * c) % p for c in prod])
                    ok = (hp == target)
                    n_lemma4 += 1
                    critical.append({'p': p, 'd': d, 'A': Ab, 'B': Bb, 'lemma4_ok': ok, 'binom_mod_p': C})
                    if not ok:
                        lemma4_fail.append({'p': p, 'd': d, 'A': Ab, 'B': Bb, 'HP': hp, 'target': target})
            # ---- T3 Theorem 2 (mu_d = A+B) and Sarkozy (d=(p-1)/2): pairs with A+B <= mu_d (no 0)
            for A, B in dfs_pairs(p, Hm):
                Ab = bits(A); Bb = bits(B)
                if len(Ab) < 2 or len(Bb) < 2: continue
                if sumset_mask(A, B, p) == Hm:
                    # every sub-B with A+B = mu_d
                    subs = []
                    if len(Bb) <= 14:
                        for k in range(2, len(Bb) + 1):
                            for Bs in itertools.combinations(Bb, k):
                                if sumset_mask(A, mask(Bs), p) == Hm: subs.append(list(Bs))
                    else:
                        subs = [Bb]
                    for Bs in subs:
                        n_thm2 += 1
                        rec = {'p': p, 'd': d, 'A': Ab, 'B': Bs, '|A|': len(Ab), '|B|': len(Bs)}
                        thm2_wit.append(rec)
                        if not (len(Ab) == len(Bs) and len(Ab) ** 2 == d): thm2_bad.append(rec)
                        if d == (p - 1) // 2:
                            sark_wit.append(rec)
                if d == (p - 1) // 2: n_sark += 1
            # ---- T3 Theorem 1 (Lev-Sonn): A - A = mu_d u {0}; A a clique of Cay(F_p, mu_d) with 0
            if (-1) % p in H:
                # DFS over cliques containing 0
                stack = [(1, 1)]
                while stack:
                    A, nxt = stack.pop()
                    n_thm1 += 1
                    Ab = bits(A)
                    if len(Ab) >= 2:
                        diff = 0
                        for a in Ab:
                            diff |= shift_mask(neg_mask(A, p), a, p)
                        if diff == S:
                            thm1_wit.append({'p': p, 'd': d, 'A': Ab})
                            if d not in (2, 6): thm1_bad.append({'p': p, 'd': d, 'A': Ab})
                    for a in range(nxt, p):
                        if all(((a - b) % p) in H for b in Ab):
                            stack.append((A | (1 << a), a + 1))
        primes_done.append(p)
    checks['T1_hanson_petridis_pairs_checked'] = n_hp
    checks['T1_max_excess_over_bound'] = hp_max_excess
    checks['T2_lemma4_dcritical_pairs_tested'] = n_lemma4
    checks['T2_lemma4_binomial_vanishing_cases'] = n_binom_zero
    checks['T3_sarkozy_A_nodes_checked'] = n_sark
    checks['T3_thm2_decompositions_found'] = n_thm2
    checks['T3_thm1_clique_nodes_checked'] = n_thm1
    checks['T1_T3_primes_all_d'] = [p for p in primes_done if p <= P_ALL_D]
    checks['T1_T3_primes_Q_case'] = primes_done
    wit['T1_violations'] = hp_viol
    wit['T2_lemma4_failures'] = lemma4_fail
    wit['T2_dcritical_pairs_sample'] = critical[:40]
    wit['T3_sarkozy_counterexamples'] = sark_wit
    wit['T3_thm2_decompositions'] = thm2_wit
    wit['T3_thm2_violations'] = thm2_bad
    wit['T3_thm1_decompositions'] = thm1_wit
    wit['T3_thm1_violations'] = thm1_bad
    results['timing_s']['T1_T2_T3'] = round(time.time() - t0, 2)

# ---------------------------------------------------------------- T4
def run_T4(P_MAX=43, time_cap=90.0):
    t0 = time.time()
    n_nodes = 0; found = []; bad = []; primes_done = []
    for p in primes_upto(P_MAX):
        if p < 5: continue
        if time.time() - t0 > time_cap: break
        for d in range(1, p - 1):
            if (p - 1) % d: continue
            H = subgroup(p, d); Hm = mask(H)
            # DFS over A (elements increasing) with A (+)^ A <= H
            stack = [(0, 0, 0)]       # (Amask, restricted-sumset mask, next)
            while stack:
                A, R, nxt = stack.pop()
                n_nodes += 1
                if R == Hm and popcount(A) >= 2:
                    rec = {'p': p, '|H|': d, 'A': bits(A)}
                    found.append(rec)
                    if d >= 7: bad.append(rec)
                for a in range(nxt, p):
                    add = shift_mask(A, a, p)        # {x + a : x in A}
                    if add & ~Hm: continue
                    stack.append((A | (1 << a), R | add, a + 1))
        primes_done.append(p)
    checks['T4_restricted_sumset_nodes_checked'] = n_nodes
    checks['T4_primes'] = primes_done
    checks['T4_sizes_realised'] = sorted(set(r['|H|'] for r in found))
    wit['T4_subgroup_equals_restricted_sumset'] = found[:60]
    wit['T4_violations_size_ge_7'] = bad
    results['timing_s']['T4'] = round(time.time() - t0, 2)

# ---------------------------------------------------------------- T5
def run_T5(PS=(13, 17, 29, 37, 41), KS=(2, 3, 4), extra=((13, 5), (17, 5), (29, 5))):
    t0 = time.time()
    rows = []; n = 0; viol = []
    combos = [(p, K) for p in PS for K in KS] + list(extra)
    for p, K in combos:
        chi = legendre_table(p)
        # Satake Def. 3, p = 1 mod 4 (r = 0): <phi_i,phi_j> = chi(a_i - a_j)/sqrt(p) for i != j <= p.
        M = np.array([[chi[(i - j) % p] for j in range(p)] for i in range(p)], dtype=float)
        delta_sqrtp = 0.0; lhs_max = 0; best = None
        clique_max = 0
        for U in itertools.combinations(range(p), K):
            sub = M[np.ix_(U, U)]
            nrm = float(np.max(np.abs(np.linalg.eigvalsh(sub))))
            s = int(abs(sub.sum()))
            if nrm > delta_sqrtp: delta_sqrtp = nrm
            if s > lhs_max: lhs_max = s; best = U
            if all(chi[(u - v) % p] == 1 for u, v in itertools.combinations(U, 2)): clique_max = K
            n += 1
        # Lemma 28: |sum chi(u-v)| <= delta sqrt(p) |U|,  delta sqrt(p) = max spectral norm
        ratio = lhs_max / (delta_sqrtp * K)
        rows.append({'p': p, 'K': K, 'delta_K_sqrt_p': round(delta_sqrtp, 6), 'delta_K': round(delta_sqrtp / math.sqrt(p), 6),
                     'max_|sum chi(u-v)|': lhs_max, 'ratio': round(ratio, 6), 'has_clique_of_size_K': bool(clique_max == K),
                     'clique_bound_delta_sqrtp_plus_1': round(delta_sqrtp + 1, 6)})
        if ratio > 1 + 1e-9: viol.append(rows[-1])
    checks['T5_subsets_checked'] = n
    checks['T5_rows'] = rows
    wit['T5_lemma28_violations'] = viol
    results['timing_s']['T5'] = round(time.time() - t0, 2)

# ---------------------------------------------------------------- T6
def max_clique(adj):
    """Exact maximum clique (vertex list) of a graph given as neighbour bitmasks:
    branch and bound with a greedy colouring bound (Tomita-style)."""
    n = len(adj)
    best = [[]]
    def colour_order(cand):
        order = []; col = 0; rem = cand
        while rem:
            col += 1; q = rem
            while q:
                v = (q & -q).bit_length() - 1
                q &= ~(1 << v); q &= ~adj[v]
                rem &= ~(1 << v); order.append((v, col))
        return order
    def expand(cur, cand):
        order = colour_order(cand)
        for v, c in reversed(order):
            if len(cur) + c <= len(best[0]): return
            nc = cand & adj[v]
            if nc: expand(cur + [v], nc)
            elif len(cur) + 1 > len(best[0]): best[0] = cur + [v]
            cand &= ~(1 << v)
    expand([], (1 << n) - 1)
    return best[0]

def max_clique_size_plain(adj):
    """Independent second algorithm: plain Carraghan-Pardalos recursion with the
    trivial bound |cur| + |cand| (no colouring)."""
    n = len(adj); best = [0]
    def rec(size, cand):
        if size + popcount(cand) <= best[0]: return
        if cand == 0:
            best[0] = max(best[0], size); return
        while cand:
            if size + popcount(cand) <= best[0]: return
            v = (cand & -cand).bit_length() - 1
            cand &= ~(1 << v)
            rec(size + 1, cand & adj[v])
    rec(0, (1 << n) - 1)
    return best[0]

def run_T6(P_CLIQUE=421, time_cap=150.0):
    t0 = time.time()
    rows = []; viol = []; n = 0; bad_witness = []; disagree = []
    for p in primes_upto(P_CLIQUE):
        if p % 4 != 1: continue
        if time.time() - t0 > time_cap: break
        chi = legendre_table(p)
        Q = [x for x in range(1, p) if chi[x] == 1]
        # every maximum clique can be moved (translation, multiplication by Q) to contain 0 and 1,
        # so omega = 2 + omega(G_p[Q n (Q+1)]).
        N = [x for x in Q if chi[(x - 1) % p] == 1]
        adj = [0] * len(N)
        for i, u in enumerate(N):
            for j, v in enumerate(N):
                if i != j and chi[(u - v) % p] == 1: adj[i] |= 1 << j
        sub = max_clique(adj) if N else []
        clique = [0, 1] + [N[i] for i in sub]
        omega = len(clique)
        # explicit verification of the witness
        ok = all(chi[(u - v) % p] == 1 for u, v in itertools.combinations(clique, 2))
        if not ok: bad_witness.append({'p': p, 'clique': clique})
        # independent second exact algorithm
        omega2 = 2 + (max_clique_size_plain(adj) if N else 0)
        if omega2 != omega: disagree.append({'p': p, 'omega_colouring_bb': omega, 'omega_plain_bb': omega2})
        hp = (math.sqrt(2 * p - 1) + 1) / 2
        rows.append({'p': p, 'omega': omega, 'omega_second_algorithm': omega2, 'max_clique_witness': clique, 'witness_verified': ok,
                     'HP_bound': round(hp, 4), 'floor_HP': math.floor(hp),
                     'sqrt_p_over_2_plus_1': round(math.sqrt(p / 2) + 1, 4), 'ratio_omega_over_sqrtp': round(omega / math.sqrt(p), 4)})
        n += 1
        if omega > hp + 1e-12: viol.append(rows[-1])
    checks['T6_primes_checked'] = n
    checks['T6_rows'] = rows
    checks['T6_max_ratio_omega_over_sqrtp'] = max(r['ratio_omega_over_sqrtp'] for r in rows)
    checks['T6_omega_table'] = {r['p']: r['omega'] for r in rows}
    wit['T6_HP_violations'] = viol
    wit['T6_unverified_witnesses'] = bad_witness
    wit['T6_algorithm_disagreements'] = disagree
    results['timing_s']['T6'] = round(time.time() - t0, 2)

# ---------------------------------------------------------------- T7
def run_T7(P_MAX=20000):
    t0 = time.time(); n = 0; bad = []; minA = None
    for p in primes_upto(P_MAX):
        if p % 4 != 3 or p < 3: continue
        chi = legendre_table(p)
        A = sum(chi[a] for a in range(1, (p - 1) // 2 + 1))
        n += 1
        if A <= 0: bad.append({'p': p, 'A': A})
        if minA is None or A < minA[1]: minA = (p, A)
    checks['T7_primes_3mod4_checked'] = n
    checks['T7_min_A_p'] = {'p': minA[0], 'A': minA[1]}
    wit['T7_nonpositive_A_p'] = bad
    results['timing_s']['T7'] = round(time.time() - t0, 2)

# ---------------------------------------------------------------- T8
def run_T8(PS=(5, 13, 17, 29, 37, 41, 53, 61, 73, 89, 97, 101)):
    t0 = time.time(); rows = []; bad = []
    for p in PS:
        chi = legendre_table(p)
        Aadj = np.array([[1.0 if chi[(i - j) % p] == 1 else 0.0 for j in range(p)] for i in range(p)])
        ev = np.linalg.eigvalsh(Aadj)
        lmax, lmin = float(ev[-1]), float(ev[0])
        theta = -p * lmin / (lmax - lmin)         # Lovasz theta for vertex-transitive graphs
        rows.append({'p': p, 'lambda_max': round(lmax, 6), 'lambda_min': round(lmin, 6), 'theta': round(theta, 6), 'sqrt_p': round(math.sqrt(p), 6)})
        if abs(theta - math.sqrt(p)) > 1e-6: bad.append(rows[-1])
    checks['T8_rows'] = rows
    wit['T8_theta_not_sqrt_p'] = bad
    results['timing_s']['T8'] = round(time.time() - t0, 2)

if __name__ == '__main__':
    quick = '--quick' in sys.argv
    run_T1_T2_T3(P_ALL_D=23 if quick else 31, P_QCASE=31 if quick else 53)
    print('T1-T3 done', results['timing_s'], flush=True)
    run_T4(P_MAX=29 if quick else 43)
    print('T4 done', results['timing_s'], flush=True)
    run_T5()
    print('T5 done', results['timing_s'], flush=True)
    run_T6(P_CLIQUE=200 if quick else 421)
    print('T6 done', results['timing_s'], flush=True)
    run_T7(P_MAX=3000 if quick else 20000)
    run_T8()
    results['total_time_s'] = round(time.time() - T0, 2)
    results['summary'] = {
        'T1_hanson_petridis_holds': len(wit['T1_violations']) == 0,
        'T2_kalmynin_lemma4_holds': len(wit['T2_lemma4_failures']) == 0,
        'T3_sarkozy_holds': len(wit['T3_sarkozy_counterexamples']) == 0,
        'T3_thm2_holds': len(wit['T3_thm2_violations']) == 0,
        'T3_thm1_holds': len(wit['T3_thm1_violations']) == 0,
        'T4_yip_yoo_holds': len(wit['T4_violations_size_ge_7']) == 0,
        'T5_satake_lemma28_holds': len(wit['T5_lemma28_violations']) == 0,
        'T6_hp_clique_bound_holds': len(wit['T6_HP_violations']) == 0,
        'T6_clique_witnesses_verified_and_algorithms_agree': len(wit['T6_unverified_witnesses']) == 0 and len(wit['T6_algorithm_disagreements']) == 0,
        'T7_zhang_statement_holds': len(wit['T7_nonpositive_A_p']) == 0,
        'T8_theta_equals_sqrt_p': len(wit['T8_theta_not_sqrt_p']) == 0,
    }
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as f: json.dump(results, f, indent=1)
    print(json.dumps(results['summary'], indent=1))
    print('checks:', {k: v for k, v in checks.items() if not isinstance(v, list)})
    print('total time', results['total_time_s'], 's ->', OUT)
