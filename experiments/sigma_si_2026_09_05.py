#!/usr/bin/env python3
"""Verifier for research/sigma-si-2026-09-05.md (sigma pass, direction `si`: Conjecture SI).

Reads the precomputed witnesses in results/sigma_si_2026_09_05_search.json, re-verifies every stored witness by exact
integer computation (numpy int arithmetic; the smallest refutation witness also in pure Python), re-runs the cheap
from-scratch scans (all subgroups of all primes p <= 4000; the 1243 stored crux rectangles; the exact (1,2)/(1,3)
formulas), and writes results/sigma_si_2026_09_05.json with the counts of checks and every witness.
Standard library + numpy + sympy only.  Runtime about 1-2 minutes.

Conventions: p odd prime, chi Legendre symbol with chi(0)=0, A,B subsets of F_p^*, S(A,B) = sum chi(a+b),
g(t) = S(tA,B), D = D_{1/2}(A,B) = {t in F_p^* : 2|g(t)| >= |g(1)|}, R(A,B) = #{(a,b,a',b'): ab' = a'b} (ratio energy),
z2 = S^2/R (signal-to-noise), T_B(x) = sum_b chi(x+b), M_H = max_{c != 0}|T_H(c)|, H^+ = H cap squares.
"""
import json, os, sys, time
from math import log, log2, sqrt, erfc, gcd
from itertools import combinations
import numpy as np
from sympy import primitive_root, divisors

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, '..'))
SEARCH = os.path.join(ROOT, 'results', 'sigma_si_2026_09_05_search.json')
OUT = os.path.join(ROOT, 'results', 'sigma_si_2026_09_05.json')

# ----------------------------------------------------------------------------------------------------------- helpers
def chi_table(p):
    chi = -np.ones(p, dtype=np.int8); chi[0] = 0
    x = np.arange(1, (p - 1) // 2 + 1, dtype=np.int64); chi[(x * x) % p] = 1
    return chi

def T_profile(chi, B, p):
    T = np.zeros(p, dtype=np.int32)
    for b in B: T += np.roll(chi, -int(b) % p).astype(np.int32)
    return T

def g_profile(chi, A, B, p, T=None):
    if T is None: T = T_profile(chi, B, p)
    t = np.arange(p, dtype=np.int64); g = np.zeros(p, dtype=np.int64)
    for a in A: g += T[(t * (int(a) % p)) % p]
    return g

def D_half(g):
    thr = abs(int(g[1])); return (np.nonzero(2 * np.abs(g[1:]) >= thr)[0] + 1).tolist()

def ratio_energy(A, B, p, chi=None):
    from collections import Counter
    mu, sig = Counter(), Counter()
    for a in A:
        ia = pow(int(a), p - 2, p); ca = int(chi[int(a) % p]) if chi is not None else 1
        for b in B:
            r = (int(b) * ia) % p; mu[r] += 1; sig[r] += ca
    return int(sum(v * v for v in mu.values())), int(sum(v * v for v in sig.values()))

def legendre_py(a, p):
    a %= p
    if a == 0: return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1

def mult_order(r, p):
    o, x = 1, r % p
    while x != 1: x = x * r % p; o += 1
    return o

def gp_length(A, p):
    """Shortest geometric progression {c r^i : 0 <= i < L} in F_p^* containing A (exhaustive over ratios)."""
    A = sorted(int(a) % p for a in A); k = len(A)
    if k == 1: return 1
    gr = int(primitive_root(p)); dlog = {}; x = 1
    for i in range(p - 1): dlog[x] = i; x = x * gr % p
    logs = sorted(dlog[a] for a in A); n = p - 1; best = p
    for e in range(1, p - 1):
        d = gcd(e, n)
        if any((l - logs[0]) % d for l in logs): continue
        m = n // d; einv = pow(e // d, -1, m) if m > 1 else 0
        pos = sorted(((l - logs[0]) // d * einv) % m for l in logs)
        gaps = [(pos[(i + 1) % k] - pos[i]) % m for i in range(k)]
        L = m if k == m else m - max(gaps) + 1
        best = min(best, L)
        if best == k: break
    return best

def primes_upto(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]: s[i * i::i] = False
    return np.nonzero(s)[0].tolist()

phibar2 = lambda z: erfc(z / sqrt(2))

class Checker:
    def __init__(self): self.n = 0; self.fail = []
    def check(self, ok, tag, info=None):
        self.n += 1
        if not ok: self.fail.append({'tag': tag, 'info': info})
        return ok

C = Checker()
res = {'checks': {}, 'witnesses': {}}
T0 = time.time()
S = json.load(open(SEARCH))

# ------------------------------------------------------------------- 1. refutation witnesses (random-like families)
sec = {'witnesses': []}
for w in S['refutation_random']:
    p, A, B = w['p'], w['A'], w['B']; m, n = len(A), len(B)
    chi = chi_table(p)
    g = g_profile(chi, A, B, p); Sv = int(g[1]); D = D_half(g)
    R, Rchi = ratio_energy(A, B, p, chi)
    C.check(all(1 <= a < p for a in A) and all(1 <= b < p for b in B) and len(set(A)) == m and len(set(B)) == n, 'ref_random_subsets', w['family'])
    C.check(Sv == w['S'] and len(D) == w['D_size'] and R == w['R'], 'ref_random_recompute', (w['family'], p, Sv, len(D), R))
    C.check(m * n >= 64 and 2 * abs(Sv) >= m * n, 'ref_random_hypotheses', (w['family'], p))
    viol1 = len(D) > 8 * log2(p)
    rec = {'family': w['family'], 'p': p, 'm': m, 'n': n, 'S': Sv, 'bias': abs(Sv) / (m * n), 'R': R, 'z2': Sv * Sv / R, 'z2_over_lnp': Sv * Sv / R / log(p),
           'D_size': len(D), 'eight_log2p': 8 * log2(p), 'two_log2p': 2 * log2(p), 'clause1_violated': viol1,
           'random_model_1+(p-2)*2Phibar(z/2)': 1 + (p - 2) * phibar2(abs(Sv) / sqrt(R) / 2), 'A': A, 'B': B}
    if m == 8 and n == 8 and p <= 10000:
        gA, gB = gp_length(A, p), gp_length(B, p)
        rec['gp_len_A'], rec['gp_len_B'] = gA, gB
        rec['clause2_applicable'] = gA > 2 * m and gB > 2 * n
        rec['clause2_violated'] = rec['clause2_applicable'] and len(D) > 2 * log2(p)
    sec['witnesses'].append(rec)
    if p <= 5000:  # independent pure-Python recomputation of S and of |D| for the p = 1009 8x8 witness
        Spy = sum(legendre_py(a + b, p) for a in A for b in B)
        C.check(Spy == Sv, 'ref_random_pure_python_S', (w['family'], p))
        if p == 1009 and m == 8:
            Dpy = sum(1 for t in range(1, p) if 2 * abs(sum(legendre_py(t * a + b, p) for a in A for b in B)) >= abs(Spy))
            C.check(Dpy == len(D), 'ref_random_pure_python_D', (p, Dpy, len(D)))
sec['n_clause1_violations'] = sum(r['clause1_violated'] for r in sec['witnesses'])
sec['n_clause2_violations'] = sum(r.get('clause2_violated', False) for r in sec['witnesses'])
res['witnesses']['refutation_random'] = sec
# smallest-p hill-climbed 1x64 family: re-verify every stored B
sm = S['smallest_p_1x64_half']; nv = 0; first = None
for e in sm:
    p, B = e['p'], e['B']; chi = chi_table(p); T = T_profile(chi, B, p)
    ok = int(T[1]) == 32 and len(B) == 64 and len(set(B)) == 64 and 0 not in B
    d = int(np.count_nonzero(2 * np.abs(T[1:]) >= 32))
    C.check(ok and d == e['D_size'], 'smallest_p_family_recompute', (p, d, e['D_size']))
    if d > 8 * log2(p):
        nv += 1
        if first is None: first = {'p': p, 'A': [1], 'B': B, 'S': 32, 'D_size': d, 'eight_log2p': 8 * log2(p)}
res['witnesses']['smallest_p_1x64_half'] = {'primes_scanned': [e['p'] for e in sm], 'n_violations': nv, 'first_violation': first,
                                             'all_p_ge_first_violate': all(e['D_size'] > 8 * log2(e['p']) for e in sm if first and e['p'] >= first['p'])}
res['checks']['1_refutation_random'] = C.n; print('1. refutation witnesses', C.n, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ------------------------------------------------------------- 2. subgroups: from-scratch scan p <= 4000 and rectangles
n0 = C.n
PMAX = S['subgroup_scan']['PMAX']
viol = []; maxratio = (0, None); rects = []; sub_viol = []; smallest_sub = None; nsub = 0
for p in primes_upto(PMAX):
    if p < 11: continue
    chi = chi_table(p); gr = int(primitive_root(p)); L8 = 8 * log2(p)
    for N in divisors(p - 1):
        if N < 2 or N > 2 * (p + 1) ** 0.5 + 1: continue
        r = pow(gr, (p - 1) // N, p); H = [1]
        for _ in range(N - 1): H.append(H[-1] * r % p)
        T = T_profile(chi, H, p); absT = np.abs(T[1:]); M = int(absT.max()); c = int(np.argmax(absT)) + 1; nsub += 1
        # coset covariance T_H(ch) = chi(h) T_H(c) (Prop. 2.1 of the structured note), spot check on the maximiser
        C.check(all(int(T[(c * h) % p]) == int(chi[h]) * int(T[c]) for h in H[:5]), 'coset_covariance', (p, N))
        # elementary bound M_H^2 <= p - N (second moment on one coset)
        C.check(M * M <= p - N, 'M_H_sq_le_p_minus_N', (p, N, M))
        if N > L8:
            if M / N > maxratio[0]: maxratio = (M / N, {'p': p, 'N': N, 'M_H': M, 'c': c})
            if 2 * M >= N: viol.append({'p': p, 'N': N, 'M_H': M, 'c': c})
        if 2 * M >= N:
            Hp = [h for h in H if chi[h] == 1]
            if len(Hp) * N >= 64:
                cnt = int(np.count_nonzero(2 * absT >= M)); C.check(cnt % N == 0, 'D_union_of_cosets', (p, N))
                A = sorted((c * h) % p for h in Hp)
                if p <= 400:  # direct recomputation of D from the definition for small p
                    g = g_profile(chi, A, H, p, T); C.check(abs(int(g[1])) == M * len(Hp) and len(D_half(g)) == cnt, 'subgroup_rect_D_direct', (p, N))
                rec = {'p': p, 'N': N, 'M_H': M, 'c': c, 'm': len(Hp), 'S': M * len(Hp), 'R': N * len(Hp) ** 2, 'z2_over_lnp': M * M / N / log(p),
                       'D_size': cnt, 'D_cosets': cnt // N, 'D_over_log2p': cnt / log2(p)}
                rects.append(rec)
                if cnt > L8:
                    sub_viol.append(rec)
                    if smallest_sub is None: smallest_sub = dict(rec, A=A, H=sorted(H))
stored = {(r['p'], r['N']): r for r in S['subgroup_rectangles']['rectangles']}
C.check(len(rects) == len(stored) and all((r['p'], r['N']) in stored and stored[(r['p'], r['N'])]['D_size'] == r['D_size'] for r in rects), 'subgroup_rectangles_match_stored', (len(rects), len(stored)))
C.check(len(viol) == 0, 'no_subgroup_with_M_ge_half_above_8log2p', viol[:5])
# pure-Python re-verification of the smallest structured witness (p = 97)
if smallest_sub:
    p, A, H = smallest_sub['p'], smallest_sub['A'], smallest_sub['H']
    Spy = sum(legendre_py(a + b, p) for a in A for b in H)
    Dpy = [t for t in range(1, p) if 2 * abs(sum(legendre_py(t * a + b, p) for a in A for b in H)) >= abs(Spy)]
    C.check(Spy == smallest_sub['S'] and len(Dpy) == smallest_sub['D_size'], 'smallest_subgroup_witness_pure_python', (p, Spy, len(Dpy)))
    C.check(all(legendre_py(h, p) == 1 for h in H) and sorted(set((h1 * h2) % p for h1 in H for h2 in H)) == H, 'p97_H_is_subgroup_of_squares', p)
    smallest_sub['D'] = Dpy
res['witnesses']['subgroups'] = {'PMAX': PMAX, 'n_subgroups_scanned': nsub, 'violations_M_ge_half_above_8log2p': viol,
                                 'max_M_over_N_above_8log2p': maxratio, 'max_N_over_log2p_with_M_ge_half': S['subgroup_scan']['max_N_over_log2p_with_M_ge_half'],
                                 'n_subgroup_rectangles_bias_ge_half_mn_ge_64': len(rects), 'n_literal_SI_violations': len(sub_viol),
                                 'smallest_violation': smallest_sub, 'largest_ratio_violation': max(sub_viol, key=lambda r: r['D_over_log2p']) if sub_viol else None,
                                 'max_z2_over_lnp_among_violations': max(r['z2_over_lnp'] for r in sub_viol) if sub_viol else None,
                                 'max_z2_over_lnp_all_subgroup_rectangles': max(r['z2_over_lnp'] for r in rects)}
res['checks']['2_subgroups'] = C.n - n0; print('2. subgroups', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ------------------------------------------------------------------------ 3. stored crux rectangles: SNR versus |D|
n0 = C.n
chis = {}
rows = []
for r in S['stored_crux']['rectangles']:
    p, A, B = r['p'], r['A'], r['B']; m, n = len(A), len(B)
    if p not in chis: chis[p] = chi_table(p)
    chi = chis[p]
    g = g_profile(chi, A, B, p); Sv = int(g[1]); D = D_half(g); R, Rchi = ratio_energy(A, B, p, chi)
    C.check(Sv == r['S'] and len(D) == r['D_size'] and R == r['R'], 'stored_crux_recompute', (p, m, n))
    C.check(m * n >= 64 and 2 * abs(Sv) >= m * n and 0 not in A and 0 not in B, 'stored_crux_hypotheses', (p, m, n))
    # exact second-moment identity (Theorem 2.1 of the crux note, 0 not in A): sum_t g(t)^2 = p R_chi - n^2 s_A^2
    sA = int(sum(int(chi[a]) for a in A))
    C.check(int(np.sum(g.astype(object) ** 2)) == p * Rchi - n * n * sA * sA, 'second_moment_identity', (p, m, n))
    # Theorem 6.1 bound |D| <= 4 p R / S^2 and Chung |S| <= sqrt(mn(p-n))
    C.check(len(D) * Sv * Sv <= 4 * p * R, 'D_le_4pR_over_S2', (p, m, n))
    C.check(Sv * Sv <= m * n * (p - n) and Sv * Sv <= m * n * (p - m), 'chung', (p, m, n))
    rows.append((p, m, n, Sv, R, len(D), Sv * Sv / R / log(p)))
bins = [(0, 2), (2, 4), (4, 8), (8, 12), (12, 16), (16, 1e9)]; tab = []
for lo, hi in bins:
    sub = [x for x in rows if lo <= x[6] < hi]
    if sub:
        w = max(sub, key=lambda x: x[5] / log2(x[0]))
        tab.append({'z2_over_lnp_bin': [lo, hi], 'count': len(sub), 'max_D_over_log2p': w[5] / log2(w[0]), 'witness_p_m_n_S_R_D': w[:6],
                    'frac_D_eq_1': sum(x[5] == 1 for x in sub) / len(sub), 'n_gt_log2p': sum(x[5] > log2(x[0]) for x in sub)})
res['witnesses']['stored_crux'] = {'n_rectangles': len(rows), 'snr_table': tab,
                                   'max_D_over_log2p_z2_ge_8lnp': max([x[5] / log2(x[0]) for x in rows if x[6] >= 8], default=None),
                                   'max_D_over_log2p_z2_ge_4lnp': max([x[5] / log2(x[0]) for x in rows if x[6] >= 4], default=None),
                                   'p1009_example': S['model'].get('p1009_example')}
res['checks']['3_stored_crux'] = C.n - n0; print('3. stored crux rectangles', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ------------------------------------------------------------------------------------------- 4. clique rectangles
n0 = C.n
cl = S['cliques']; clout = {}
for p_s, v in cl['per_p'].items():
    p = int(p_s); chi = chi_table(p); k = v['omega']
    for key in ('best', 'best_z8'):
        w = v.get(key)
        if not w: continue
        A, B = w['A'], w['B']
        C.check(all(chi[(a - a2) % p] == 1 for a in A for a2 in A if a != a2) and len(A) == k and 0 not in A, 'clique_property', (p, key))
        C.check(sorted((-a) % p for a in A) == sorted(B), 'clique_B_is_minus_A', (p, key))
        g = g_profile(chi, A, B, p); Sv = int(g[1]); D = D_half(g)
        C.check(Sv == k * (k - 1) == w['S'] and len(D) == w['D_size'], 'clique_recompute', (p, key, Sv, len(D)))
        j0 = (k * (k - 1)) // (2 * (2 * k - 1)); Aset = set(A); Dset = set(D)
        sym = [s for s in range(1, p) if sum(((s * a) % p) in Aset for a in A) >= k - j0]
        C.check(all(s in Dset for s in sym) and len(sym) == w['sym_j0_size'], 'clique_Sym_j0_subset_D', (p, key))
        R, _ = ratio_energy(A, B, p); C.check(R == w['R'], 'clique_R', (p, key))
    # omega(p) cross-check for small p by brute force over cliques containing 0,1 (p <= 61)
    if p <= 61:
        cand = [x for x in range(2, p) if chi[x] == 1 and chi[(x - 1) % p] == 1]
        best = 2
        for s in range(1, len(cand) + 1):
            found = any(all(chi[(a - b) % p] == 1 for a, b in combinations(cc, 2)) for cc in combinations(cand, s))
            if found: best = 2 + s
            else: break
        C.check(best == k, 'omega_bruteforce', (p, best, k))
    if k >= 8 and p <= 400:  # SI clause 2 status of the extremal clique rectangle
        gA = gp_length(v['best']['A'], p); v['best']['gp_len_A'] = gA
        v['best']['clause2_violated'] = gA > 2 * k and v['best']['D_size'] > 2 * log2(p)
    clout[p] = {'omega': k, 'n_max_cliques_containing_01': v['n_max_cliques_containing_01'], 'n_rectangles': v['n_rectangles'], 'max_D': v['max_D'],
                'max_D_over_log2p': v['max_D_over_log2p'], 'median_D': v['median_D'], 'frac_D_eq_1': v['frac_D_eq_1'], 'z2_over_lnp_of_max': v['best']['z2_over_lnp'],
                'sym_j0_of_max': v['best']['sym_j0_size'], 'n_rectangles_z8': v['n_rectangles_z8'], 'max_D_z8': v['max_D_z8'],
                'best': {kk: v['best'][kk] for kk in ('A', 'B', 'S', 'R', 'D_size', 'stab', 'sym_j0_size', 'c', 'gp_len_A', 'clause2_violated') if kk in v['best']}}
res['witnesses']['cliques'] = {'per_p': clout, 'sym_subset_checks_precomputed': cl['checks_sym_subset_D'], 'sym_subset_fails_precomputed': cl['fails_sym_subset_D'],
                               'max_D_over_log2p_all': max(v['max_D_over_log2p'] for v in clout.values()),
                               'max_D_z8_over_log2p': max([v['max_D_z8'] / log2(p) for p, v in clout.items() if v['max_D_z8']], default=None),
                               'n_rectangles_total': sum(v['n_rectangles'] for v in clout.values()),
                               'clause2_violations_k_ge_8': [(p, v['best']['D_size'], v['best'].get('gp_len_A')) for p, v in clout.items() if v['best'].get('clause2_violated')]}
res['checks']['4_cliques'] = C.n - n0; print('4. cliques', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ---------------------------------------------------------------------------------------- 5. GP-family rectangles
n0 = C.n
def NU(chi, U, p):
    mask = np.ones(p, dtype=bool); mask[0] = False
    for u in U: mask &= (np.roll(chi, -u) == 1)
    return np.nonzero(mask)[0].tolist()
gpout = {}
for grp in ('p_le_1500_all_ratios', 'large_primes_sampled_ratios'):
    gpout[grp] = {}
    for p_s, v in S['gp_family'][grp]['per_p'].items():
        p = int(p_s); chi = chi_table(p); rec = {}
        for key in ('all', 'z4', 'z8', 'z16'):
            w = v.get(key)
            if not w: continue
            r, L, k = w['r'], w['L'], w['k']; U = [pow(r, i, p) for i in range(L)]
            C.check(len(set(U)) == L and w['A'] == U[:k], 'gp_A_is_prefix_of_U', (p, key))
            B = NU(chi, U, p); C.check(sorted(B) == sorted(w['B']) and len(B) * k >= 64, 'gp_B_is_N(U)', (p, key))
            g = g_profile(chi, w['A'], B, p); Sv = int(g[1]); D = D_half(g)
            C.check(Sv == k * len(B) == w['S'] and len(D) == w['D_size'], 'gp_recompute', (p, key, Sv, len(D)))
            inner = [pow(r, j, p) % p for j in range(-(k // 4), L - k + 1 + k // 4)]
            Dset = set(D); C.check(all(x in Dset for x in inner) and len(set(inner)) == w['LB'], 'gp_inner_set_subset_D', (p, key))
            R, _ = ratio_energy(w['A'], B, p); C.check(R == w['R'], 'gp_R', (p, key))
            rec[key] = {kk: w[kk] for kk in ('r', 'L', 'k', 'nB', 'S', 'R', 'z2_over_lnp', 'D_size', 'D_over_log2p', 'LB')}
        gpout[grp][p] = rec
res['witnesses']['gp_family'] = {'n_rectangles_p_le_1500': S['gp_family']['p_le_1500_all_ratios']['n_rectangles'], 'per_p': gpout}
for key in ('all', 'z4', 'z8', 'z16'):
    vals = [(v[key]['D_over_log2p'], p) for grp in gpout for p, v in gpout[grp].items() if key in v]
    res['witnesses']['gp_family']['max_D_over_log2p_' + key] = max(vals) if vals else None
    d_lb = [v[key]['D_size'] - v[key]['LB'] for grp in gpout for p, v in gpout[grp].items() if key in v]
    res['witnesses']['gp_family']['D_minus_LB_' + key] = {'min': min(d_lb), 'max': max(d_lb)} if d_lb else None
res['checks']['5_gp_family'] = C.n - n0; print('5. GP family', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# -------------------------------------------------- 6. tiny rectangles: exact formulas (1,2), (1,3) and exhaustive maxima
n0 = C.n
small = S['small_exact']; C.check(small['formula_failures'] == 0, 'precomputed_formula_failures', small['formula_failures'])
form = {'(1,2)': 0, '(1,3)': 0}
for p in (101, 211, 401):
    chi = chi_table(p); x = np.arange(p)
    mx2 = 0
    for b in range(2, p):
        T = chi[(x + 1) % p].astype(np.int32) + chi[(x + b) % p].astype(np.int32)
        d = int(np.count_nonzero(np.abs(T[1:]) >= 1)); mx2 = max(mx2, d)
        C.check(d == (p + 1) // 2 - (1 if chi[b] == 1 else 0), 'formula_1x2', (p, b)); form['(1,2)'] += 1
    st = [e for e in small['exhaustive'] if e['m'] == 1 and e['n'] == 2 and e['p'] == p][0]
    C.check(st['max_D'] == mx2, 'exhaustive_1x2_max', (p, mx2, st['max_D']))
    if p <= 211:
        mx3 = 0; T1 = chi[(x + 1) % p].astype(np.int32)
        for b in range(2, p):
            Tb = T1 + chi[(x + b) % p].astype(np.int32)
            for c in range(b + 1, p):
                T = Tb + chi[(x + c) % p].astype(np.int32)
                if not np.any(np.abs(T[1:]) == 3): continue
                d = int(np.count_nonzero(2 * np.abs(T[1:]) >= 3)); mx3 = max(mx3, d)
                Bt = (1, b, c); ssum = 0
                for i in range(3):
                    j, kk = [t for t in range(3) if t != i]; ssum += int(chi[((Bt[j] - Bt[i]) * (Bt[kk] - Bt[i])) % p])
                C.check((p + ssum) % 4 == 0 and d == (p + ssum) // 4 - (1 if abs(int(T[0])) == 3 else 0), 'formula_1x3', (p, b, c)); form['(1,3)'] += 1
        st = [e for e in small['exhaustive'] if e['m'] == 1 and e['n'] == 3 and e['p'] == p][0]
        C.check(st['max_D'] == mx3, 'exhaustive_1x3_max', (p, mx3, st['max_D']))
# exhaustive (2,2),(3,3) maxima re-check at p = 31 and the stored witnesses of every exhaustive entry
for e in small['exhaustive']:
    p, m, n, w = e['p'], e['m'], e['n'], e['witness']
    if w is None: continue
    chi = chi_table(p); B = w['B']
    if m == 1:  # the stored witness records B and |T_B(a)|; pick an a attaining it
        T = T_profile(chi, B, p); ta = w.get('abs_T_a', 3)
        A = [int(np.nonzero(np.abs(T[1:]) == ta)[0][0]) + 1]
    else:
        A = w['A']
    g = g_profile(chi, A, B, p); Sv = int(g[1]); D = D_half(g)
    C.check(len(D) == e['max_D'] and 2 * abs(Sv) >= m * n, 'small_exhaustive_witness', (p, m, n, len(D), e['max_D']))
for p, mn_list in ((23, ((2, 2), (3, 3))), (31, ((2, 2), (2, 3))), (29, ((3, 3),))):
    chi = chi_table(p)
    for (m, n) in mn_list:
        mx = 0
        for Arest in combinations(range(2, p), m - 1):
            A = (1,) + Arest
            for B in combinations(range(1, p), n):
                g = g_profile(chi, A, B, p); Sv = int(g[1])
                if 2 * abs(Sv) < m * n: continue
                mx = max(mx, len(D_half(g)))
        st = [e for e in small['exhaustive'] if e['m'] == m and e['n'] == n and e['p'] == p][0]
        C.check(st['max_D'] == mx, 'exhaustive_mxn_max', (p, m, n, mx, st['max_D']))
summary = {}
for (m, n) in ((1, 2), (1, 3), (1, 4), (1, 5), (2, 2), (2, 3), (3, 3)):
    es = [e for e in small['exhaustive'] if e['m'] == m and e['n'] == n]
    if es:
        summary[f'({m},{n})'] = {'n_primes': len(es), 'p_range': [min(e['p'] for e in es), max(e['p'] for e in es)],
                                 'min_maxD_over_p': min(e['max_D_over_p'] for e in es), 'max_maxD_over_p': max(e['max_D_over_p'] for e in es),
                                 'example_largest_p': max(es, key=lambda e: e['p'])}
res['witnesses']['small_exact'] = {'formula_checks_here': form, 'formula_checks_precomputed': small['formula_checks'], 'summary': summary}
res['checks']['6_small_exact'] = C.n - n0; print('6. tiny rectangles', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ---------------------------------------------------------- 7. identities on random pairs + adversarial generic-A rectangles
n0 = C.n
rng = np.random.default_rng(5)
for p in (101, 401, 1009):
    chi = chi_table(p)
    for _ in range(30):
        m, n = int(rng.integers(2, 12)), int(rng.integers(2, 12))
        A = sorted(int(a) for a in rng.choice(np.arange(1, p), m, replace=False)); B = sorted(int(b) for b in rng.choice(np.arange(1, p), n, replace=False))
        g = g_profile(chi, A, B, p); R, Rchi = ratio_energy(A, B, p, chi); sA = int(sum(int(chi[a]) for a in A))
        C.check(int(np.sum(g.astype(object) ** 2)) == p * Rchi - n * n * sA * sA, 'second_moment_identity_random', (p, m, n))
        C.check(Rchi <= R <= m * n * min(m, n), 'R_bounds', (p, m, n))
        # coset structure: D_{A,cH} = c D_{A,H} for a subgroup H (size of D invariant)
    for N in [d for d in divisors(p - 1) if 4 <= d <= 40][:3]:
        gr = int(primitive_root(p)); r = pow(gr, (p - 1) // N, p); H = [pow(r, i, p) for i in range(N)]
        A = sorted(int(a) for a in rng.choice(np.arange(1, p), 6, replace=False)); c = int(rng.integers(2, p))
        gH = g_profile(chi, A, H, p); gcH = g_profile(chi, A, [(c * h) % p for h in H], p)
        C.check(all(int(gcH[t]) == int(chi[c]) * int(gH[(t * pow(c, p - 2, p)) % p]) for t in range(p)), 'coset_dilation_covariance', (p, N))
        C.check(all(int(gH[(t * h) % p]) == int(chi[h]) * int(gH[t]) for t in range(1, p, 37) for h in H), 'g_H_covariance', (p, N))
adv = []
for w in S['adversarial_generic']:
    p = w['p']; chi = chi_table(p); A, B = w['A'], w['B']
    g = g_profile(chi, A, B, p); Sv = int(g[1]); D = D_half(g); R, _ = ratio_energy(A, B, p)
    C.check(Sv == len(A) * len(B) == w['S'] and len(D) == w['D_size'] and R == w['R'], 'adversarial_recompute', p)
    adv.append({k: w[k] for k in ('p', 'L', 'k', 'extra', 'S', 'R', 'z2_over_lnp', 'D_size', 'D_over_log2p', 'gpA')})
res['witnesses']['adversarial_generic_high_snr'] = adv
res['checks']['7_identities_adversarial'] = C.n - n0; print('7. identities/adversarial', C.n - n0, 'checks', len(C.fail), 'failures', '%.0fs' % (time.time() - T0))

# ------------------------------------------------------------------------------ 8. corrected conjecture: combined table
res['witnesses']['snr_table_all_sources'] = S['snr_table_all_sources']
res['witnesses']['model_stored'] = S['model'].get('model_stored'); res['witnesses']['Vprime_over_R'] = S['model'].get('Vprime_over_R')
res['checks']['total'] = C.n; res['failures'] = C.fail; res['n_failures'] = len(C.fail); res['seconds'] = time.time() - T0
json.dump(res, open(OUT, 'w'), indent=1)
print('TOTAL checks', C.n, 'failures', len(C.fail), '%.0fs' % res['seconds'], '->', OUT)
if C.fail: print(json.dumps(C.fail[:10], indent=1, default=str))
