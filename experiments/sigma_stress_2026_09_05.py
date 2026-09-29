#!/usr/bin/env python3
"""Verifier for research/sigma-stress-2026-09-05.md (sigma pass, direction `stress`).

Re-derives, by exact integer arithmetic, every witness of the note:
  * the moment route  T_ns = sum_{t in F_p^*} g(t)^{2k} - T_sq  (tuple note Prop. 1.2) against brute-force tuple
    enumeration (own code) and against the independent D-decomposition of experiments/sigma_tuple_2026_09_05.py;
  * the subgroup reduction  T_ns(cH+, H; k) = |H+|^{2k} Delta_k(H)  and the closed form of Delta_2 for H in Q;
  * the refutation of Conjecture T(k): A = B = Q (ratio ~ (2k-1)!! sqrt p), A = B = H (ratio ~ sqrt h), gp rectangles;
  * the exponent analysis (T'(k), exponent k - 1/2) on all stored families;
  * the refutation of Conjecture SI*: designed 1 x n rectangles (stored B re-verified; the smallest one in pure Python),
    the strong-form variant, the GP-type, subgroup and annealing scans.
Reads results/sigma_stress_2026_09_05_search.json (precomputed heavy runs) and writes results/sigma_stress_2026_09_05.json.
Standard library + numpy only (sympy not needed).  Runtime: a few minutes.
"""
import itertools, json, math, os, random, sys, time
from collections import Counter
from fractions import Fraction
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
SEARCH = os.path.join(ROOT, 'results', 'sigma_stress_2026_09_05_search.json')
OUT = os.path.join(ROOT, 'results', 'sigma_stress_2026_09_05.json')
sys.path.insert(0, HERE)

T0 = time.time()
CHECKS = {}
FAILURES = []
WIT = {}


def check(name, ok, witness=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not ok:
        FAILURES.append({'check': name, 'witness': witness})
    return ok


def el(t):
    return f'{time.time() - T0:.0f}s'


# ------------------------------------------------------------------------------------------------ exact core
def legendre_table(p):
    chi = -np.ones(p, dtype=np.int8); chi[0] = 0
    x = np.arange(1, (p - 1) // 2 + 1, dtype=np.int64); chi[(x * x) % p] = 1
    return chi


def legendre_py(a, p):
    a %= p
    if a == 0:
        return 0
    return 1 if pow(a, (p - 1) // 2, p) == 1 else -1


def class_data(p, A, B, chi):
    nu = Counter(); sg = Counter()
    for a in A:
        a %= p; assert a != 0
        ia = pow(a, p - 2, p); ca = int(chi[a])
        for b in B:
            r = (-(b % p) * ia) % p; nu[r] += 1; sg[r] += ca
    ratios = sorted(nu)
    return ratios, [nu[r] for r in ratios], [sg[r] for r in ratios]


def profile_from_classes(p, ratios, sigma, chi):
    t = np.arange(1, p, dtype=np.int64); g = np.zeros(p - 1, dtype=np.int64)
    for r, s in zip(ratios, sigma):
        if s:
            g += s * chi[(t - r) % p].astype(np.int64)
    return g


def profile_direct(p, A, B, chi):
    t = np.arange(1, p, dtype=np.int64); g = np.zeros(p - 1, dtype=np.int64)
    for a in A:
        for b in B:
            g += chi[(t * a + b) % p].astype(np.int64)
    return g


def moment_exact(g, k):
    vals, cnts = np.unique(g, return_counts=True)
    return sum(int(c) * int(v) ** (2 * k) for v, c in zip(vals, cnts))


def even_poly(sigmas, k):
    c = [Fraction(0)] * (k + 1); c[0] = Fraction(1)
    for s in sigmas:
        if s == 0:
            continue
        f = [Fraction(s ** (2 * j), math.factorial(2 * j)) for j in range(k + 1)]
        new = [Fraction(0)] * (k + 1)
        for i in range(k + 1):
            if c[i] == 0:
                continue
            for j in range(k + 1 - i):
                new[i + j] += c[i] * f[j]
        c = new
    return c


def Nk_words(L, k):
    v = even_poly([1] * L, k)[k] * math.factorial(2 * k)
    assert v.denominator == 1
    return int(v)


def T_sq_exact(p, ratios, sigma, k):
    f2k = math.factorial(2 * k)
    full = even_poly(sigma, k)[k] * f2k; assert full.denominator == 1
    Lnz = sum(1 for r in ratios if r != 0)
    total = (p - 1 - Lnz) * int(full)
    for i, r in enumerate(ratios):
        if r == 0:
            continue
        q = even_poly(sigma[:i] + sigma[i + 1:], k)[k] * f2k; assert q.denominator == 1
        total += int(q)
    return total


def T_ns_exact(p, A, B, k, chi):
    ratios, nu, sigma = class_data(p, A, B, chi)
    g = profile_from_classes(p, ratios, sigma, chi)
    M = moment_exact(g, k); Tsq = T_sq_exact(p, ratios, sigma, k)
    return {'T_ns': M - Tsq, 'T_sq': Tsq, 'moment': M, 'g': g, 'R_nu': sum(x * x for x in nu),
            'R_chi': sum(x * x for x in sigma), 'L': len(ratios), 'S': int(g[0])}


def rho_T(Tns, p, m, n, k):
    return abs(Tns) * min(m, n) ** k / (math.sqrt(p) * (m * n) ** (2 * k))


def rho_prime(Tns, p, m, n, k):
    return abs(Tns) * min(m, n) ** (k - 0.5) / (math.sqrt(p) * (m * n) ** (2 * k))


def violates_Tk_exact(Tns, p, m, n, k, C):
    """exact integer test of  |T_ns| min^k > C sqrt(p) (mn)^{2k}  for rational C."""
    lhs = abs(Tns) * min(m, n) ** k; rhs = Fraction(C) * (m * n) ** (2 * k)
    return Fraction(lhs * lhs) > rhs * rhs * p


def brute_force_T_ns(p, A, B, k, chi):
    pairs = [(a % p, b % p) for a in A for b in B]; T_sq = T_ns = 0
    for tau in itertools.product(pairs, repeat=2 * k):
        w = 0
        for t in range(1, p):
            pr = 1
            for (a, b) in tau:
                pr = pr * (t * a + b) % p
            w += int(chi[pr])
        cnt = Counter((-b * pow(a, p - 2, p)) % p for (a, b) in tau)
        if all(c % 2 == 0 for c in cnt.values()):
            T_sq += w
        else:
            T_ns += w
    return T_sq, T_ns


def primitive_root(p):
    q = p - 1; fac = []; d = 2
    while d * d <= q:
        if q % d == 0:
            fac.append(d)
            while q % d == 0:
                q //= d
        d += 1
    if q > 1:
        fac.append(q)
    for g in range(2, p):
        if all(pow(g, (p - 1) // f, p) != 1 for f in fac):
            return g
    raise RuntimeError


def subgroup(p, h, gen=None):
    assert (p - 1) % h == 0
    gen = gen or primitive_root(p)
    x = pow(gen, (p - 1) // h, p); H = []; y = 1
    for _ in range(h):
        H.append(y); y = y * x % p
    return sorted(H)


def T_H_profile(p, H, chi):
    s = np.arange(p, dtype=np.int64); T = np.zeros(p, dtype=np.int64)
    for h in H:
        T += chi[(s + h) % p].astype(np.int64)
    return T


def coset_values(p, H, chi, gen):
    """T_H(g^j) for j < d (one representative per coset), exact ints."""
    h = len(H); d = (p - 1) // h; Hn = np.array(H, dtype=np.int64)
    reps = np.empty(d, dtype=np.int64); x = 1
    for j in range(d):
        reps[j] = x; x = x * gen % p
    vals = np.empty(d, dtype=np.int64); CH = max(1, (2 * 10 ** 7) // h)
    for s in range(0, d, CH):
        R = reps[s:s + CH]
        vals[s:s + CH] = chi[(R[:, None] + Hn[None, :]) % p].astype(np.int32).sum(axis=1)
    return vals


def Delta_from_vals(p, h, vals, k):
    uv, cnt = np.unique(vals, return_counts=True)
    M = h * sum(int(c) * int(v) ** (2 * k) for v, c in zip(uv, cnt))
    return M - (p - 1 - h) * Nk_words(h, k) - h * Nk_words(h - 1, k)


def Delta_k(p, H, chi, k, T=None):
    h = len(H)
    if T is None:
        T = T_H_profile(p, H, chi)
    return moment_exact(T[1:], k) - (p - 1 - h) * Nk_words(h, k) - h * Nk_words(h - 1, k)


def primes_upto(N):
    s = np.ones(N + 1, dtype=bool); s[:2] = False
    for i in range(2, int(N ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = False
    return np.nonzero(s)[0].tolist()


def loglog_slope(xs, ys):
    X = np.log(np.array(xs, dtype=float)); Y = np.log(np.array(ys, dtype=float))
    return float(((X - X.mean()) * (Y - Y.mean())).sum() / ((X - X.mean()) ** 2).sum())


# ------------------------------------------------------------------------------------------------ SI* helpers
def T_of(chi, B, p):
    x = np.arange(p, dtype=np.int64); T = np.zeros(p, dtype=np.int64)
    for b in B:
        T += chi[(x + b) % p]
    return T


def ratio_energy(A, B, p):
    mu = Counter()
    for a in A:
        ia = pow(int(a), p - 2, p)
        for b in B:
            mu[(int(b) * ia) % p] += 1
    return sum(v * v for v in mu.values())


def rect_profile(p, A, B, chi):
    TB = T_of(chi, B, p); t = np.arange(p, dtype=np.int64); g = np.zeros(p, dtype=np.int64)
    for a in A:
        g += TB[(t * int(a)) % p]
    return g


def D_half_size(g):
    S = int(g[1]); return int(np.count_nonzero(2 * np.abs(g[1:]) >= abs(S))), S


# ================================================================================================ main
def main():
    S = json.load(open(SEARCH))
    rng = random.Random(2026)
    out = {'sections': {}}

    # ---------------------------------------------------------------- 1. validation of the moment route
    n0 = sum(CHECKS.values())
    for p in (7, 11, 13):
        chi = legendre_table(p)
        for (m, n, k) in ((2, 2, 2), (2, 3, 2), (3, 3, 2), (2, 2, 3), (3, 2, 3)):
            A = rng.sample(range(1, p), m); B = rng.sample(range(0, p), n)
            bsq, bns = brute_force_T_ns(p, A, B, k, chi); r = T_ns_exact(p, A, B, k, chi)
            check('bruteforce_vs_moment_route', (bsq, bns) == (r['T_sq'], r['T_ns']), {'p': p, 'A': A, 'B': B, 'k': k})
            check('profile_classes_vs_direct', bool((profile_direct(p, A, B, chi) == r['g']).all()), {'p': p, 'A': A, 'B': B})
    try:
        import sigma_tuple_2026_09_05 as TUP
        have_tup = True
    except Exception as e:  # pragma: no cover
        have_tup = False; WIT['tuple_verifier_import_error'] = repr(e)
    if have_tup:
        for p in (101, 211, 401):
            chi = legendre_table(p); specs = []
            for (m, n, k) in ((6, 6, 2), (6, 4, 2), (4, 4, 3), (5, 5, 3)):
                specs.append((rng.sample(range(1, p), m), rng.sample(range(0, p), n), k))
            for h in [d for d in range(2, 21) if (p - 1) % (2 * d) == 0]:
                H = subgroup(p, h); specs.append((H, H, 2))
                if h <= 8:
                    specs.append((H, H, 3))
                c = rng.randrange(2, p); specs.append((H, [c * x % p for x in H], 2))
            for A, B, k in specs:
                an = TUP.analyze(p, A, B, k, chi, tag='stress', D_stats=False); r = T_ns_exact(p, A, B, k, chi)
                check('tuple_verifier_Ddecomposition_vs_moment_route',
                      an['T_ns'] == r['T_ns'] and an['T_sq'] == r['T_sq'] and an['moment'] == r['moment'], {'p': p, 'k': k, 'm': len(A), 'n': len(B)})
        check('tuple_verifier_internal_checks_pass', len(TUP.FAILURES) == 0, TUP.FAILURES[:3])
    # subgroup reduction and closed form
    for p in (101, 151, 211, 307, 401, 1009):
        chi = legendre_table(p); gen = primitive_root(p)
        for h in [d for d in range(2, 30) if (p - 1) % d == 0]:
            H = subgroup(p, h, gen); T = T_H_profile(p, H, chi)
            for k in (2, 3):
                D = Delta_k(p, H, chi, k, T)
                check('coset_moment_equals_full_moment', Delta_from_vals(p, h, coset_values(p, H, chi, gen), k) == D, {'p': p, 'h': h, 'k': k})
                Hp = [x for x in H if chi[x] == 1]; c = rng.randrange(1, p); A = [c * x % p for x in Hp]
                r = T_ns_exact(p, A, H, k, chi)
                check('subgroup_reduction_Tns_eq_m2k_Delta', r['T_ns'] == len(A) ** (2 * k) * D, {'p': p, 'h': h, 'k': k, 'c': c})
            if all(chi[x] == 1 for x in H) and 4 <= h <= 16:
                s = np.arange(1, p, dtype=np.int64); tot = 0
                for Dset in itertools.combinations(H, 4):
                    pr = np.ones(p - 1, dtype=np.int64)
                    for r_ in Dset:
                        pr *= chi[(s + r_) % p]
                    tot += int(pr.sum())
                Tm1 = int(T[p - 1])
                check('closed_form_Delta2_subgroup', Delta_k(p, H, chi, 2, T) == 24 * tot - 2 * h * (h - 1) * (6 * h - 11) - 6 * h * Tm1 * Tm1, {'p': p, 'h': h})
    out['sections']['1_validation'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('1. validation', out['sections']['1_validation'], flush=True)

    # ---------------------------------------------------------------- 2. T(k): A = B = Q
    n0 = sum(CHECKS.values()); qrows = []
    for p in (61, 101, 151, 211, 401, 1009, 2003, 4001, 10007):
        chi = legendre_table(p); Q = [x for x in range(1, p) if chi[x] == 1]; h = len(Q); T = T_H_profile(p, Q, chi)
        # exact profile: T_Q(s) = -(1 + chi(s))/2 for s != 0
        check('T_Q_profile_formula', all(int(T[s]) == -(1 + int(chi[s])) // 2 for s in range(1, p)), {'p': p})
        for k in (2, 3):
            D = Delta_k(p, Q, chi, k, T); Tns = h ** (2 * k) * D
            if p <= 401:
                check('Q_Tns_moment_route_matches', T_ns_exact(p, Q, Q, k, chi)['T_ns'] == Tns, {'p': p, 'k': k})
            # exact value: sum_{t != 0} g^{2k} = h^{2k+1}; T_sq = h^{2k}[(p-1-h) N_k(h) + h N_k(h-1)]
            check('Q_exact_Tns_formula', Tns == h ** (2 * k + 1) - h ** (2 * k) * ((p - 1 - h) * Nk_words(h, k) + h * Nk_words(h - 1, k)), {'p': p, 'k': k})
            poly = -2 * h * (h - 1) * (3 * h - 2) if k == 2 else -h * (30 * h ** 3 - 105 * h ** 2 + 137 * h - 62)
            check('Q_Delta_k_closed_polynomial', D == poly, {'p': p, 'k': k, 'Delta': D, 'poly': poly})
            C = 7.8 if k == 2 else 36.0
            viol = violates_Tk_exact(Tns, p, h, h, k, C)
            check('Q_violates_Tk_with_observed_constant', viol or p < 100, {'p': p, 'k': k})
            qrows.append({'p': p, 'h': h, 'k': k, 'Delta_k': D, 'T_ns': Tns, 'rho_T': rho_T(Tns, p, h, h, k),
                          'rho_T_over_sqrtp': rho_T(Tns, p, h, h, k) / math.sqrt(p), 'rho_prime': rho_prime(Tns, p, h, h, k),
                          'violates_Tk_C_obs': bool(viol), 'C_obs': C})
    # ratio/sqrt p -> (2k-1)!!: exact statement  |T_ns| min^k / (sqrt p (mn)^{2k}) = (2k-1)!! sqrt p (1 - O(1/h))
    for r in qrows:
        dfac = 3 if r['k'] == 2 else 15
        check('Q_rho_over_sqrtp_within_5pct_of_(2k-1)!!_for_p>=1009', abs(r['rho_T_over_sqrtp'] / dfac - 1) < 0.05 or r['p'] < 1009, r)
    WIT['T_k_refutation_A_eq_B_eq_Q'] = qrows
    out['sections']['2_Tk_Q'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('2. A=B=Q', out['sections']['2_Tk_Q'], flush=True)

    # ---------------------------------------------------------------- 3. T(k): subgroups A = B = H
    n0 = sum(CHECKS.values())
    # 3a. full scan p <= 3000 from scratch (exact, full profile)
    scan = []
    for p in primes_upto(3000):
        if p < 7:
            continue
        chi = legendre_table(p); gen = primitive_root(p); q = (p - 1) // 2
        for h in range(4, 2 * int(math.isqrt(p)) + 1):
            if q % h:
                continue
            H = subgroup(p, h, gen); T = T_H_profile(p, H, chi); rec = {'p': p, 'h': h}
            for k in (2, 3):
                D = Delta_k(p, H, chi, k, T); rec[f'Delta{k}'] = D
                rec[f'rho{k}'] = abs(D) / (math.sqrt(p) * h ** k); rec[f'rhop{k}'] = abs(D) / (math.sqrt(p) * h ** (k + 0.5))
            scan.append(rec)
    stored_scan = [r for r in S['subgroup_scan_3000'] if r['h'] >= 4]   # the stored scan also holds h = 2, 3 (no |D| = 4 terms)
    check('subgroup_scan_3000_row_count', len(scan) == len(stored_scan), {'n': len(scan), 'stored': len(stored_scan)})
    sd = {(r['p'], r['h']): r for r in stored_scan}
    for r in scan:
        s_ = sd.get((r['p'], r['h']))
        check('subgroup_scan_3000_Delta_matches_stored', s_ is not None and s_['Delta2'] == r['Delta2'] and s_['Delta3'] == r['Delta3'], {'p': r['p'], 'h': r['h']})
    summ = {}
    for k in (2, 3):
        C = 7.8 if k == 2 else 36.0
        nviol = sum(1 for r in scan if violates_Tk_exact(r[f'Delta{k}'] * r['h'] ** (2 * k), r['p'], r['h'], r['h'], k, C))
        worst = max(scan, key=lambda r: r[f'rho{k}'])
        bands = {}
        for lo, hi in ((4, 7), (8, 15), (16, 31), (32, 63), (64, 127)):
            b = [r for r in scan if lo <= r['h'] <= hi]
            if b:
                rs = sorted(r[f'rho{k}'] for r in b); rp = sorted(r[f'rhop{k}'] for r in b)
                bands[f'{lo}-{hi}'] = {'n': len(b), 'median_rho': rs[len(rs) // 2], 'max_rho': rs[-1], 'median_rho_prime': rp[len(rp) // 2], 'max_rho_prime': rp[-1],
                                       'frac_Delta_negative': sum(1 for r in b if r[f'Delta{k}'] < 0) / len(b)}
        xs = []; ys = []; yp = []
        for h in sorted(set(r['h'] for r in scan if r['h'] >= 8)):
            b = [r for r in scan if r['h'] == h]
            if len(b) >= 5:
                xs.append(h); ys.append(sorted(r[f'rho{k}'] for r in b)[len(b) // 2]); yp.append(sorted(r[f'rhop{k}'] for r in b)[len(b) // 2])
        summ[f'k={k}'] = {'n_subgroups': len(scan), 'n_violating_Tk_with_C_obs': nviol, 'C_obs': C,
                          'worst': {x: worst[x] for x in ('p', 'h', f'Delta{k}', f'rho{k}', f'rhop{k}')}, 'bands': bands,
                          'slope_median_rho_vs_h': loglog_slope(xs, ys) if len(xs) >= 3 else None,
                          'slope_median_rho_prime_vs_h': loglog_slope(xs, yp) if len(xs) >= 3 else None}
        check('subgroup_scan_rho_grows_with_h_slope_gt_0.3', summ[f'k={k}']['slope_median_rho_vs_h'] > 0.3, summ[f'k={k}'])
        check('subgroup_scan_worst_rho_exceeds_C_obs_by_10x', worst[f'rho{k}'] > 10 * C, worst)
    WIT['T_k_refutation_subgroups_p_le_3000'] = summ
    # 3b. stored p <= 20000 rows: re-verify a random sample of 250 exactly (full profile), all 136 targeted rows (coset method)
    st20 = S['subgroup_scan_20000']
    sample = rng.sample(st20, min(250, len(st20)))
    for r in sample:
        p, h = r['p'], r['h']; chi = legendre_table(p); H = subgroup(p, h); T = T_H_profile(p, H, chi)
        check('subgroup_scan_20000_sample_exact', Delta_k(p, H, chi, 2, T) == r['Delta2'] and Delta_k(p, H, chi, 3, T) == r['Delta3'], {'p': p, 'h': h})
    for r in S['subgroup_targeted']:
        p, h = r['p'], r['h']; chi = legendre_table(p); gen = primitive_root(p); H = subgroup(p, h, gen); vals = coset_values(p, H, chi, gen)
        ok = Delta_from_vals(p, h, vals, 2) == r['Delta2'] and Delta_from_vals(p, h, vals, 3) == r['Delta3']
        check('subgroup_targeted_exact', ok, {'p': p, 'h': h})
        check('subgroup_targeted_maxT_matches', int(np.abs(vals).max()) == r['maxT'], {'p': p, 'h': h})
        if p <= 200000:
            check('subgroup_targeted_coset_equals_full_profile', Delta_k(p, H, chi, 2) == r['Delta2'], {'p': p, 'h': h})
    # 3c. the named witnesses
    named = []
    for (p, h) in ((2833, 59), (13183, 169), (38011, 181), (1073153, 1024), (14109313, 724), (5969921, 512)):
        chi = legendre_table(p); gen = primitive_root(p); H = subgroup(p, h, gen); vals = coset_values(p, H, chi, gen)
        rec = {'p': p, 'h': h, 'M_H': int(np.abs(vals).max()), 'd': (p - 1) // h}
        for k in (2, 3):
            D = Delta_from_vals(p, h, vals, k); C = 7.8 if k == 2 else 36.0
            rec[f'Delta{k}'] = D; rec[f'rho{k}'] = abs(D) / (math.sqrt(p) * h ** k); rec[f'rhop{k}'] = abs(D) / (math.sqrt(p) * h ** (k + 0.5))
            rec[f'violates_Tk_C_obs_k{k}'] = bool(violates_Tk_exact(D * h ** (2 * k), p, h, h, k, C))
            check('named_subgroup_witness_violates_Tk', rec[f'violates_Tk_C_obs_k{k}'], rec)
        rec['coset_values_sorted'] = sorted(int(v) for v in vals) if len(vals) <= 100 else None
        named.append(rec)
    WIT['T_k_refutation_subgroups_named'] = named
    big = [r for r in st20 + S['subgroup_targeted'] if r['h'] >= 8]
    tg = S['subgroup_targeted']
    WIT['T_k_subgroup_growth'] = {}
    for k in (2, 3):
        xs = []; ys = []
        for h in sorted(set(r['h'] for r in tg)):
            b = [r[f'rho{k}'] for r in tg if r['h'] == h]; xs.append(h); ys.append(float(np.mean(b)))
        WIT['T_k_subgroup_growth'][f'k={k}'] = {'targeted_mean_rho_by_h': dict(zip(map(str, xs), ys)), 'slope_mean_rho_vs_h_targeted': loglog_slope(xs, ys),
                                                 'max_rho_all': max(r[f'rho{k}'] for r in big), 'max_rho_prime_all': max(r[f'rhop{k}'] for r in big),
                                                 'rms_rho_prime_h_le_sqrtp': float(math.sqrt(np.mean([r[f'rhop{k}'] ** 2 for r in big if r['h'] ** 2 <= r['p']]))),
                                                 'n_rows': len(big)}
    out['sections']['3_Tk_subgroups'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('3. subgroups', out['sections']['3_Tk_subgroups'], flush=True)

    # ---------------------------------------------------------------- 4. other families (stored rows re-verified exactly)
    n0 = sum(CHECKS.values()); fam = S['families']
    fam_summary = {}
    for r in fam:
        if r['A'] is None or r['B'] is None:
            continue
        p = r['p']
        if p > 40009 and r['tag'] != 'Q_Q':
            continue
        chi = legendre_table(p); res = T_ns_exact(p, r['A'], r['B'], r['k'], chi)
        check('family_row_exact', res['T_ns'] == r['T_ns'] and res['T_sq'] == r['T_sq'], {'tag': r['tag'], 'p': p, 'k': r['k']})
    for tag in sorted(set(r['tag'] for r in fam)):
        for k in (2, 3):
            rs = [r for r in fam if r['tag'] == tag and r['k'] == k]
            if not rs:
                continue
            w = max(rs, key=lambda r: r['rho']); wp = max(rs, key=lambda r: r['rhop'])
            fam_summary[f'{tag}|k={k}'] = {'n': len(rs), 'max_rho': w['rho'], 'worst_rho': {x: w[x] for x in ('p', 'm', 'n', 'S', 'T_ns')},
                                          'max_rho_prime': wp['rhop'], 'worst_rho_prime': {x: wp[x] for x in ('p', 'm', 'n', 'S', 'T_ns')},
                                          'n_violating_Tk_C_obs': sum(1 for r in rs if r['rho'] > (7.8 if k == 2 else 36.0))}
    # gp: growth of rho with K (m = n = K)
    for k in (2, 3):
        gp = [r for r in fam if r['tag'] == 'gp' and r['k'] == k and r['m'] == r['n']]
        xs = []; ys = []
        for K in sorted(set(r['m'] for r in gp)):
            b = [r['rho'] for r in gp if r['m'] == K]
            if len(b) >= 3:
                xs.append(K); ys.append(float(np.sqrt(np.mean(np.array(b) ** 2))))
        fam_summary[f'gp_rms_rho_by_K|k={k}'] = {'K': xs, 'rms_rho': ys, 'slope': loglog_slope(xs, ys) if len(xs) >= 3 else None}
    WIT['families'] = fam_summary
    out['sections']['4_families'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('4. families', out['sections']['4_families'], flush=True)

    # ---------------------------------------------------------------- 5. SI*: designed 1 x n witnesses
    n0 = sum(CHECKS.values()); si = []
    for w in S['si_designed']:
        p, n, B = w['p'], w['n'], w['B']; chi = legendre_table(p); L2 = math.log2(p); lnp = math.log(p)
        check('si_B_size_and_range', len(B) == n and len(set(B)) == n and all(0 < b < p for b in B), {'p': p})
        T = T_of(chi, B, p); Sv = int(T[1]); D, _ = D_half_size(T)
        R = n  # singleton A: all ratios distinct
        check('si_R_equals_n', ratio_energy([1], B, p) == n, {'p': p})
        hyp = (2 * Sv >= n) and (Sv * Sv >= 8 * n * lnp)
        check('si_stored_S_and_D_match', Sv == w['S'] and D == w['D'], {'p': p, 'S': Sv, 'D': D, 'stored': (w['S'], w['D'])})
        check('si_hypothesis_holds', hyp, {'p': p, 'S': Sv, 'n': n})
        check('si_D_exceeds_2log2p', D > 2 * L2, {'p': p, 'D': D, '2log2p': 2 * L2})
        # exact hypothesis in integers: S^2 >= 8 n ln p  <=  S^2 >= 8 n * (ln p rounded up at 1e-9) -- use float lnp with margin check
        si.append({'p': p, 'n': n, 'S': Sv, 'bias': Sv / n, 'R': R, 'z2_over_lnp': Sv * Sv / n / lnp, 'D': D, 'D_over_log2p': D / L2, '2log2p': 2 * L2,
                   'kind': w.get('kind', 'relative'), 'B': B if n <= 400 else None, 'B_stored_in_search_file': True})
    # pure-Python recomputation of the smallest witness
    w0 = min([w for w in S['si_designed'] if w.get('kind', 'relative') == 'relative'], key=lambda w: w['p'])
    p, n, B = w0['p'], w0['n'], w0['B']
    Tpy = [sum(legendre_py(x + b, p) for b in B) for x in range(p)]
    Spy = Tpy[1]; Dpy = sum(1 for x in range(1, p) if 2 * abs(Tpy[x]) >= abs(Spy))
    check('si_smallest_witness_pure_python', Spy == w0['S'] and Dpy == w0['D'] and 2 * Spy >= n and Spy * Spy >= 8 * n * math.log(p) and Dpy > 2 * math.log2(p), {'p': p, 'S': Spy, 'D': Dpy})
    WIT['SI_star_refutation_designed'] = si
    # strong form
    st = []
    for w in S['si_strong']:
        p, n, B = w['p'], w['n'], w['B']; chi = legendre_table(p); L2 = math.log2(p); lnp = math.log(p)
        T = T_of(chi, B, p); Sv = int(T[1]); D, _ = D_half_size(T)
        check('si_strong_bias_one', Sv == n and all(chi[(1 + b) % p] == 1 for b in B), {'p': p})
        check('si_strong_stored_D_matches', D == w['D'], {'p': p, 'D': D, 'stored': w['D']})
        st.append({'p': p, 'n': n, 'z2_over_lnp': n / lnp, 'D_abs': D, 'D_over_log2p': D / L2, 'exceeds_2log2p': D > 2 * L2})
    WIT['SI_star_strong_form'] = st
    out['sections']['5_SI_designed'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('5. SI designed', out['sections']['5_SI_designed'], flush=True)

    # ---------------------------------------------------------------- 6. SI*: GP-type, subgroup scan, annealing
    n0 = sum(CHECKS.values())
    gpw = []
    for p_s, v in S['si_gp'].items():
        p = int(p_s); w = v['best_z8']
        if w is None:
            continue
        chi = legendre_table(p); r, L, k = w['r'], w['L'], w['k']
        U = [pow(r, i, p) for i in range(L)]
        mask = np.ones(p, dtype=bool); mask[0] = False
        for u in U:
            mask &= (np.roll(chi, -u) == 1)
        B = np.nonzero(mask)[0].tolist(); A = U[:k]
        g = rect_profile(p, A, B, chi); D, Sv = D_half_size(g); R = ratio_energy(A, B, p)
        check('si_gp_recompute', len(B) == w['nB'] and Sv == w['S'] and D == w['D'] and R == w['R'], {'p': p})
        check('si_gp_high_snr_D_le_2log2p', Sv * Sv >= 8 * R * math.log(p) and D <= 2 * math.log2(p), {'p': p, 'D': D})
        gpw.append({'p': p, 'r': r, 'L': L, 'k': k, 'nB': len(B), 'S': Sv, 'R': R, 'z2_over_lnp': Sv * Sv / R / math.log(p), 'D': D, 'D_over_log2p': D / math.log2(p), 'n_rectangles_searched': v['n_rectangles']})
    WIT['SI_star_gp_family_large_p'] = gpw
    # subgroup SI* scan from scratch, p <= 20000
    rows = []
    for p in primes_upto(20000):
        if p < 11:
            continue
        chi = legendre_table(p); gen = primitive_root(p); lnp = math.log(p); L2 = math.log2(p)
        for h in range(2, 2 * int(math.isqrt(p)) + 2):
            if (p - 1) % h:
                continue
            H = subgroup(p, h, gen); vals = coset_values(p, H, chi, gen); M = int(np.abs(vals).max())
            if 2 * M < h:
                continue
            inQ = bool(chi[H[1]] == 1); hp = h if inQ else h // 2
            Sv = M * hp; R = h * hp * hp; D = h * int(np.count_nonzero(2 * np.abs(vals) >= M))
            rows.append({'p': p, 'h': h, 'M_H': M, 'in_Q': inQ, 'S': Sv, 'R': R, 'z2_over_lnp': Sv * Sv / R / lnp, 'D': D, 'D_over_log2p': D / L2, 'h_over_lnp': h / lnp})
    stored_rows = S['si_subgroups_20000']
    check('si_subgroup_scan_count', len(rows) == len(stored_rows), {'n': len(rows), 'stored': len(stored_rows)})
    mx = max(rows, key=lambda r: r['z2_over_lnp'])
    check('si_subgroup_none_in_SI_star_range', mx['z2_over_lnp'] < 8, mx)
    # verify one stored row exactly (the largest-h one) via the rectangle itself
    hb = max(rows, key=lambda r: r['h_over_lnp']); p = hb['p']; chi = legendre_table(p); gen = primitive_root(p); H = subgroup(p, hb['h'], gen)
    vals = coset_values(p, H, chi, gen); j = int(np.argmax(np.abs(vals))); c = pow(gen, j, p)
    Hp = [x for x in H if chi[x] == 1]; A = [c * x % p for x in Hp]
    g = rect_profile(p, A, H, chi); D, Sv = D_half_size(g)
    check('si_subgroup_rectangle_direct', abs(Sv) == hb['S'] and D == hb['D'] and ratio_energy(A, H, p) == hb['R'], {'p': p, 'h': hb['h'], 'S': Sv, 'D': D})
    WIT['SI_star_subgroup_scan_20000'] = {'n_biased_subgroup_rectangles': len(rows), 'max_z2_over_lnp': mx,
                                          'max_h_over_lnp_with_M_ge_half': hb,
                                          'n_with_D_gt_2log2p': sum(1 for r in rows if r['D'] > 2 * math.log2(r['p'])),
                                          'max_z2_over_lnp_with_D_gt_2log2p': max([r for r in rows if r['D'] > 2 * math.log2(r['p'])], key=lambda r: r['z2_over_lnp'])}
    # annealing
    an = {}
    for p_s, v in S['si_anneal'].items():
        p = int(p_s); chi = legendre_table(p); an[p] = {}
        for key, w in v.items():
            if w is None:
                an[p][key] = None; continue
            g = rect_profile(p, w['A'], w['B'], chi); D, Sv = D_half_size(g); R = ratio_energy(w['A'], w['B'], p)
            m, n = len(w['A']), len(w['B'])
            check('si_anneal_recompute', D == w['D'] and Sv == w['S'] and R == w['R'] and 2 * Sv >= m * n and Sv * Sv >= 8 * R * math.log(p), {'p': p, 'key': key})
            an[p][key] = {'D': D, 'D_over_log2p': D / math.log2(p), 'S': Sv, 'R': R, 'z2_over_lnp': Sv * Sv / R / math.log(p)}
    WIT['SI_star_annealing_small_rectangles'] = an
    out['sections']['6_SI_other'] = {'checks': sum(CHECKS.values()) - n0, 'seconds': round(time.time() - T0, 1)}
    print('6. SI other', out['sections']['6_SI_other'], flush=True)

    out['checks'] = CHECKS; out['total_checks'] = sum(CHECKS.values()); out['n_failures'] = len(FAILURES); out['failures'] = FAILURES[:20]
    out['witnesses'] = WIT; out['seconds'] = round(time.time() - T0, 1)
    json.dump(out, open(OUT, 'w'), indent=1, default=float)
    print(f'total checks {out["total_checks"]}, failures {len(FAILURES)}, {out["seconds"]} s')
    for f in FAILURES[:10]:
        print('FAILURE', f)


if __name__ == '__main__':
    main()
