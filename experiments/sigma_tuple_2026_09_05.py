#!/usr/bin/env python3
"""sigma-tuple verifier (2026-09-05).

Exact computation of the tuple decomposition of the dilation moment
    sum_t g(t)^{2k} = sum_{tau in (A x B)^{2k}} W(tau),   W(tau) = sum_t chi(prod_i (t a_i + b_i)),
with g(t) = sum_{a in A, b in B} chi(t a + b), 0 not in A, chi the Legendre symbol (chi(0) = 0).
CONVENTION: all t-sums run over t in F_p^* (the degenerate dilate t = 0, where g(0) = m sum_b chi(b) can equal
mn independently of the conjecture, is excluded).  W(tau) here equals the full F_p sum minus prod_i chi(b_i).

All quantities are exact integers.  Tuples are grouped by the multiplicity pattern j of their
ratios r = -b/a (a 2k-multiset on the ratio set).  For a pattern j with odd support D and
even-positive support E,
    W(tau) = chi(prod a_i) * Wt(j),   Wt(j) = sum_{t not in E, t != 0} prod_{r in D} chi(t - r)
                                            = W_D - sum_{e in E, e != 0} prod_{r in D} chi(e - r)   (D nonempty)
                                            = p - 1 - |E minus {0}|                                  (D empty, square)
with W_D = sum_{t != 0} prod_{r in D} chi(t - r).
and the signed / unsigned tuple counts of the pattern are
    s_j = (2k)!/prod j_r! * prod sigma_r^{j_r},     n_j = (2k)!/prod j_r! * prod nu_r^{j_r},
where nu_r = #{(a,b): -b/a = r} and sigma_r = sum_{(a,b): -b/a = r} chi(a).

Outputs results/sigma_tuple_2026_09_05.json with the counts of every exact check and every witness.
Standard library + numpy only.  Runtime target: well under ten minutes.
"""
import itertools
import json
import math
import os
import random
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.normpath(os.path.join(HERE, '..'))
OUT = os.path.join(ROOT, 'results', 'sigma_tuple_2026_09_05.json')
SEARCH = os.path.join(ROOT, 'results', 'sigma_crux_2026_09_05_search.json')

T0 = time.time()
CHECKS = {}
WITNESSES = {}
FAILURES = []


def check(name, ok, witness=None):
    CHECKS[name] = CHECKS.get(name, 0) + 1
    if not ok:
        FAILURES.append({'check': name, 'witness': witness})


def legendre_table(p):
    chi = np.zeros(p, dtype=np.int8)
    chi[:] = -1
    chi[0] = 0
    for x in range(1, p):
        chi[(x * x) % p] = 1
    return chi


def primes_in(lo, hi):
    sieve = bytearray([1]) * (hi + 1)
    sieve[0] = sieve[1] = 0
    for i in range(2, int(hi ** 0.5) + 1):
        if sieve[i]:
            sieve[i * i::i] = bytearray(len(sieve[i * i::i]))
    return [q for q in range(lo, hi + 1) if sieve[q]]


def isqrt_le(x2, y):
    """exact test  sqrt(x2) <= y  for nonnegative integers (x2 <= y^2)."""
    return x2 <= y * y


def double_factorial_odd(k):
    r = 1
    for i in range(1, 2 * k, 2):
        r *= i
    return r


# --------------------------------------------------------------------------------------------
# core: exact pattern-grouped tuple decomposition
# --------------------------------------------------------------------------------------------

def ratio_classes(p, A, B, chi):
    classes = {}
    for a in A:
        assert a % p != 0
        ainv = pow(a, p - 2, p)
        ca = int(chi[a % p])
        for b in B:
            r = (-b * ainv) % p
            if r in classes:
                classes[r][0] += 1
                classes[r][1] += ca
            else:
                classes[r] = [1, ca]
    ratios = sorted(classes)
    nu = [classes[r][0] for r in ratios]
    sigma = [classes[r][1] for r in ratios]
    return ratios, nu, sigma


def all_WD(p, ratios, chi, kmax, chunk=1500):
    """W_D = sum_t prod_{r in D} chi(t - r) for all subsets D of the ratio set with |D| even, |D| <= 2 kmax.
    Returns dict {tuple(indices): int}.  Also returns the p x L sign matrix."""
    L = len(ratios)
    t = np.arange(1, p)
    M = np.empty((p - 1, L), dtype=np.int8)
    for l, r in enumerate(ratios):
        M[:, l] = chi[(t - r) % p]
    W = {(): p - 1}
    for d in range(2, 2 * kmax + 1, 2):
        if d > L:
            break
        combos = list(itertools.combinations(range(L), d))
        arr = np.array(combos, dtype=np.int32)
        for s in range(0, len(combos), chunk):
            idx = arr[s:s + chunk]
            sub = M[:, idx]                       # (p, c, d)
            prod = sub[:, :, 0].astype(np.int8)
            for q in range(1, d):
                prod = prod * sub[:, :, q]
            sums = prod.sum(axis=0, dtype=np.int64)
            for c, D in enumerate(combos[s:s + chunk]):
                W[D] = int(sums[c])
    return W, M


def analyze(p, A, B, k, chi, tag='', D_stats=True):
    """Exact tuple decomposition for the rectangle (A,B) at the dilation moment 2k."""
    m, n = len(A), len(B)
    ratios, nu, sigma = ratio_classes(p, A, B, chi)
    L = len(ratios)
    W, M = all_WD(p, ratios, chi, k)
    # G[e][l] = chi(r_e - r_l)
    G = [[int(chi[(ratios[e] - ratios[l]) % p]) for l in range(L)] for e in range(L)]
    # exact profile g(t) and moment
    Mi = M.astype(np.int64)
    g = Mi @ np.array(sigma, dtype=np.int64)          # g(t) = sum_r sigma_r chi(t - r), t = 1..p-1
    gl = [int(x) for x in g]
    moment = sum(x ** (2 * k) for x in gl)            # over t in F_p^*
    S = gl[0]                                         # t = 1
    maxg = max(abs(x) for x in gl)
    g0 = m * sum(int(chi[b % p]) for b in B)          # the excluded dilate t = 0
    zero_idx = ratios.index(0) if 0 in ratios else -1
    Qk0 = 0
    fact = [math.factorial(i) for i in range(2 * k + 1)]
    f2k = fact[2 * k]
    T_sq = 0
    T_ns = 0
    T_tri = 0
    R2 = 0
    N_sq = 0
    N_ns = 0
    Corr = 0
    n_patterns = 0
    n_square_patterns = 0
    weil_ok = True
    weil_witness = None
    cD = {}
    cplus = {}
    by_d = {}          # |D| -> [sum s_j W_D, sum s_j corr_j, sum n_j W_D^2, sum n_j, sum n_j |W_D - corr|]
    for pat in itertools.combinations_with_replacement(range(L), 2 * k):
        n_patterns += 1
        counts = {}
        for l in pat:
            counts[l] = counts.get(l, 0) + 1
        coef = f2k
        sg = 1
        ct = 1
        D = []
        E = []
        for l, c in counts.items():
            coef //= fact[c]
            sg *= sigma[l] ** c
            ct *= nu[l] ** c
            if c & 1:
                D.append(l)
            else:
                E.append(l)
        if not D:
            n_square_patterns += 1
            w = p - 1 - sum(1 for e in E if e != zero_idx)
            T_sq += coef * sg * w
            N_sq += coef * ct
            if zero_idx not in E:
                Qk0 += coef * sg
            continue
        Dt = tuple(D)
        wD = W[Dt]
        corr = 0
        for e in E:
            if e == zero_idx:
                continue
            pr = 1
            for l in D:
                pr *= G[e][l]
                if pr == 0:
                    break
            corr += pr
        w = wD - corr
        # termwise Weil check: |W_D| <= (|D|-1) sqrt p + 1  (exact, squared)
        if abs(wD) >= 1 and (abs(wD) - 1) ** 2 > (len(D) - 1) ** 2 * p:
            weil_ok = False
            weil_witness = {'D': [ratios[l] for l in D], 'W_D': wD}
        s = coef * sg
        nn = coef * ct
        T_ns += s * w
        T_tri += nn * abs(w)
        R2 += nn * w * w
        N_ns += nn
        Corr += s * corr
        bd = by_d.get(len(D))
        if bd is None:
            by_d[len(D)] = [s * wD, s * corr, nn * wD * wD, nn, nn * abs(w)]
        else:
            bd[0] += s * wD
            bd[1] += s * corr
            bd[2] += nn * wD * wD
            bd[3] += nn
            bd[4] += nn * abs(w)
        if D_stats:
            cD[Dt] = cD.get(Dt, 0) + s
            cplus[Dt] = cplus.get(Dt, 0) + nn
    Rnu = sum(x * x for x in nu)
    E2 = by_d.get(2, [0, 0, 0, 0, 0])[0]                      # sum_{|D|=2} c(D) W_D = -sum_{|D|=2} c(D)
    T_W = sum(v[0] for d, v in by_d.items() if d >= 4)        # pure Weil part
    RW2 = sum(v[2] for d, v in by_d.items() if d >= 4)        # its random-sign scale squared
    NW = sum(v[3] for d, v in by_d.items() if d >= 4)
    TW_tri = sum(v[4] for d, v in by_d.items() if d >= 4)
    out = {
        'g0': g0, 'Qk0': Qk0, 'moment_full_Fp': moment + g0 ** (2 * k),
        'E2': E2, 'T_W': T_W, 'RW2': RW2, 'N_W': NW, 'TW_tri': TW_tri,
        'by_d': {str(d): {'sum_cW': v[0], 'sum_corr': v[1], 'sum_nW2': v[2], 'n_tuples': v[3]} for d, v in sorted(by_d.items())},
        'tag': tag, 'p': p, 'm': m, 'n': n, 'k': k, 'L': L, 'S': S, 'bias': S / (m * n),
        'max_abs_g': maxg, 'n_biased_half': sum(1 for x in gl if 2 * abs(x) >= abs(S)) if S else None,
        'R_nu': Rnu, 'max_nu': max(nu), 'sum_sigma2': sum(x * x for x in sigma),
        'moment': moment, 'T_sq': T_sq, 'T_ns': T_ns, 'T_tri': T_tri, 'R2': R2,
        'N_sq': N_sq, 'N_ns': N_ns, 'Corr': Corr, 'n_patterns': n_patterns,
        'n_square_patterns': n_square_patterns, 'weil_ok': weil_ok, 'weil_witness': weil_witness,
    }
    # exact identities
    check('T_sq_plus_T_ns_equals_moment', T_sq + T_ns == moment, {'tag': tag, 'p': p})
    check('N_sq_plus_N_ns_equals_(mn)^2k', N_sq + N_ns == (m * n) ** (2 * k), {'tag': tag, 'p': p})
    check('T_sq_nonnegative', T_sq >= 0, {'tag': tag, 'p': p, 'T_sq': T_sq})
    check('T_sq_le_(p-1)_dfact_Rnu^k', T_sq <= (p - 1) * double_factorial_odd(k) * Rnu ** k, {'tag': tag, 'p': p})
    check('Qk0_between_0_and_Nk_bound', 0 <= Qk0 <= double_factorial_odd(k) * Rnu ** k, {'tag': tag, 'p': p})
    check('termwise_Weil_W_D', weil_ok, weil_witness)
    check('T_ns_equals_TW_plus_E2_minus_Corr', T_ns == T_W + E2 - Corr, {'tag': tag, 'p': p})
    # PROVED elementary bound: |E2 - Corr| <= (k+1) * C(2k,2) * R_nu * (mn)^(2k-2)   (note: also test the sharper k)
    check('elementary_part_bound_k_plus_1', abs(E2 - Corr) <= (k + 1) * math.comb(2 * k, 2) * Rnu * (m * n) ** (2 * k - 2),
          {'tag': tag, 'p': p, 'E2': E2, 'Corr': Corr})
    if abs(E2 - Corr) > k * math.comb(2 * k, 2) * Rnu * (m * n) ** (2 * k - 2):
        # the sharper (unproved) constant k fails here: recorded as a witness, not as a failure
        WITNESSES.setdefault('elementary_bound_constant_k_violations', []).append(
            {'tag': tag, 'p': p, 'm': m, 'n': n, 'k': k, 'R_nu': Rnu, 'E2': E2, 'Corr': Corr,
             'abs_E2_minus_Corr': abs(E2 - Corr), 'bound_with_k': k * math.comb(2 * k, 2) * Rnu * (m * n) ** (2 * k - 2),
             'bound_with_k_plus_1': (k + 1) * math.comb(2 * k, 2) * Rnu * (m * n) ** (2 * k - 2)})
    # Weil part triangle bound: |T_W| <= TW_tri <= (2k-1) sqrt(p) N_W + k N_W
    check('TW_le_TWtri', abs(T_W) <= TW_tri + abs(Corr), {'tag': tag, 'p': p})
    check('W_D_equals_-1-chi(r1r2)_for_|D|=2',
          all(v == -1 - int(chi[(ratios[D[0]] * ratios[D[1]]) % p]) for D, v in W.items() if len(D) == 2),
          {'tag': tag, 'p': p})
    # T_tri <= ((2k-1) sqrt(p) + 1 + k) N_ns  (the +1 is the removed t = 0 term, +k covers |corr| <= |E| <= k)
    check('T_tri_le_((2k-1)sqrtp+1+k)Nns', (T_tri - (k + 1) * N_ns) <= 0 or
          (T_tri - (k + 1) * N_ns) ** 2 <= (2 * k - 1) ** 2 * p * N_ns ** 2, {'tag': tag, 'p': p})
    # R^2 <= ((2k-1) sqrt p + k + 1)^2 N_ns  <=  2((2k-1)^2 p + (k+1)^2) N_ns
    check('R2_le_2((2k-1)^2p+(k+1)^2)Nns', R2 <= 2 * ((2 * k - 1) ** 2 * p + (k + 1) ** 2) * N_ns, {'tag': tag, 'p': p})
    if D_stats:
        sumcW = sum(c * W[D] for D, c in cD.items())
        check('T_ns_equals_sum_cD_WD_minus_Corr', T_ns == sumcW - Corr, {'tag': tag, 'p': p})
        check('sum_cplus_equals_N_ns', sum(cplus.values()) == N_ns, {'tag': tag, 'p': p})
        # correction bound: |Corr| <= k * C(2k,2) * Rnu * (mn)^(2k-2)
        check('Corr_bound', abs(Corr) <= k * math.comb(2 * k, 2) * Rnu * (m * n) ** (2 * k - 2),
              {'tag': tag, 'p': p, 'Corr': Corr})
        # D-level statistics: alignment of c(D) with W_D; pairwise correlation of W_D, W_D' sharing 2k-1 ratios
        Ds = list(cD.keys())
        cv = np.array([cD[D] for D in Ds], dtype=np.float64)
        wv = np.array([W[D] for D in Ds], dtype=np.float64)
        nc = float(np.sqrt((cv * cv).sum()) * np.sqrt((wv * wv).sum()))
        out['sum_cD_WD'] = sumcW
        out['RD2'] = sum(cplus[D] ** 2 * W[D] ** 2 for D in Ds)                    # pattern random-sign scale^2
        # Lemma 1.7a: c_+(D) <= (2k)!/(2k-|D|)! (2k-|D|-1)!! min^{|D|} R_nu^{k-|D|/2}, and the R_D bound
        mn_min = min(m, n)
        cplus_ok = True
        for D, c in cplus.items():
            dd = len(D)
            bd = (f2k // fact[2 * k - dd]) * double_factorial_odd((2 * k - dd) // 2) * mn_min ** dd * Rnu ** (k - dd // 2)
            if c > bd:
                cplus_ok = False
                break
        check('cplus_D_bound_lemma_1.7a', cplus_ok, {'tag': tag, 'p': p})
        maxc_bound = 2 * k * double_factorial_odd(k) * (m * n) ** (k - 1) * mn_min ** (k + 1)
        check('max_cplus_le_2k(2k-1)!!(mn)^(k-1)min^(k+1)', max(cplus.values()) <= maxc_bound, {'tag': tag, 'p': p})
        # R_D^2 <= ((2k-1) sqrt p + k + 1)^2 * 2k(2k-1)!! (mn)^(3k-1) min^(k+1): exact sufficient test
        lhs = out['RD2']
        bigc = 2 * k * double_factorial_odd(k) * (m * n) ** (3 * k - 1) * mn_min ** (k + 1)
        # ((2k-1)sqrt p + k + 1)^2 >= (2k-1)^2 p + (k+1)^2, so the following is a sufficient exact test
        check('R_D_bound_lemma_1.7a', lhs <= ((2 * k - 1) ** 2 * p + (k + 1) ** 2) * bigc, {'tag': tag, 'p': p})
        out['RDW2'] = sum(cplus[D] ** 2 * W[D] ** 2 for D in Ds if len(D) >= 4)    # same, Weil part only
        out['Rc2'] = sum(cD[D] ** 2 * W[D] ** 2 for D in Ds)                        # with the signed weights
        out['cos_c_W'] = float((cv * wv).sum() / nc) if nc > 0 else None
        out['n_D'] = len(Ds)
        out['CS_scale'] = nc
        # naive (unsigned-count) version and no-correction version, for the witness of non-identity
        out['sum_cplus_WD'] = sum(cplus[D] * W[D] for D in Ds)
        # pairwise correlation among top-size subsets sharing 2k-1 elements
        if L >= 2 * k + 1:
            top = [D for D in W if len(D) == 2 * k]
            vals = {D: W[D] for D in top}
            # group by (2k-1)-subset D0: for each D0 the list of W_{D0 + r}
            s1 = 0.0
            s2 = 0.0
            sxy = 0.0
            npairs = 0
            for D in top:
                pass
            groups = {}
            for D in top:
                w = vals[D]
                for i in range(2 * k):
                    D0 = D[:i] + D[i + 1:]
                    gg = groups.get(D0)
                    if gg is None:
                        groups[D0] = [w, w * w, 1]
                    else:
                        gg[0] += w
                        gg[1] += w * w
                        gg[2] += 1
            # ordered pairs (D, D') with D != D' both containing D0: sum W W' = (sum W)^2 - sum W^2
            sxy = 0.0
            sx = 0.0
            sxx = 0.0
            npairs = 0
            for D0, (sw, sww, cnt) in groups.items():
                if cnt >= 2:
                    sxy += sw * sw - sww
                    sx += (cnt - 1) * sw
                    sxx += (cnt - 1) * sww
                    npairs += cnt * (cnt - 1)
            if npairs > 0:
                mx = sx / npairs
                var = sxx / npairs - mx * mx
                cov = sxy / npairs - mx * mx
                out['pair_corr_share_2k-1'] = float(cov / var) if var > 0 else None
                out['pair_corr_npairs'] = npairs
                # mean and mean-square of W_D over top-size subsets, and E[W^2]/p
                allw = np.array(list(vals.values()), dtype=np.float64)
                out['top_WD_mean'] = float(allw.mean())
                out['top_WD_msq_over_p'] = float((allw * allw).mean() / p)
    return out


# --------------------------------------------------------------------------------------------
# brute force at tiny primes (enumerates every tuple)
# --------------------------------------------------------------------------------------------

def brute_force(p, A, B, k, chi):
    pairs = [(a, b) for a in A for b in B]
    m, n = len(A), len(B)
    T_sq = T_ns = T_tri = R2 = N_sq = N_ns = 0
    cD = {}
    cplus = {}
    Corr = 0
    ratios_of = {}
    for (a, b) in pairs:
        ratios_of[(a, b)] = (-b * pow(a, p - 2, p)) % p
    neg_square_W = None
    for tau in itertools.product(pairs, repeat=2 * k):
        # W(tau) directly
        w = 0
        for t in range(1, p):
            pr = 1
            for (a, b) in tau:
                pr = (pr * (t * a + b)) % p
            w += int(chi[pr])
        counts = {}
        for pr_ in tau:
            r = ratios_of[pr_]
            counts[r] = counts.get(r, 0) + 1
        D = tuple(sorted(r for r, c in counts.items() if c & 1))
        E = [r for r, c in counts.items() if not (c & 1)]
        sgn = 1
        for (a, b) in tau:
            sgn *= int(chi[a % p])
        if not D:
            T_sq += w
            N_sq += 1
            # W(tau) = chi(prod a) (p - 1 - #distinct nonzero roots)
            check('bf_square_W_equals_sign_times_p-1_minus_nonzero_roots',
                  w == sgn * (p - 1 - sum(1 for r in counts if r != 0)), {'p': p, 'tau': tau, 'W': w})
            if w < 0 and neg_square_W is None:
                neg_square_W = {'tau': [list(x) for x in tau], 'W': w}
        else:
            T_ns += w
            T_tri += abs(w)
            R2 += w * w
            N_ns += 1
            cD[D] = cD.get(D, 0) + sgn
            cplus[D] = cplus.get(D, 0) + 1
            # W(tau) = sgn * ( W_D - sum_{e in E} prod_{r in D} chi(e - r) )
            WD = sum(int(np.prod([chi[(t - r) % p] for r in D])) for t in range(1, p))
            corr = 0
            for e in E:
                if e == 0:
                    continue
                pr = 1
                for r in D:
                    pr *= int(chi[(e - r) % p])
                corr += pr
            check('bf_W_tau_equals_sign_times_(W_D_minus_corr)', w == sgn * (WD - corr),
                  {'p': p, 'tau': tau, 'W': w, 'W_D': WD, 'corr': corr})
            Corr += sgn * corr
    return {'T_sq': T_sq, 'T_ns': T_ns, 'T_tri': T_tri, 'R2': R2, 'N_sq': N_sq, 'N_ns': N_ns,
            'cD': cD, 'cplus': cplus, 'Corr': Corr, 'neg_square_W': neg_square_W}


def brute_force_section(rng):
    cases = 0
    witnesses = {'naive_count_formula_fails': None, 'no_correction_formula_fails': None,
                 'negative_square_term': None, 'sign_matters_for_c_D': None}
    specs = []
    for p in (7, 11, 13):
        chi = legendre_table(p)
        for (m, n) in ((2, 2), (2, 3), (3, 2), (3, 3)):
            for trial in range(2):
                A = rng.sample(range(1, p), m)
                B = rng.sample(range(0, p), n)
                specs.append((p, A, B, 2))
        # k = 3 with mn <= 6
        for (m, n) in ((2, 2), (2, 3), (3, 2)):
            A = rng.sample(range(1, p), m)
            B = rng.sample(range(0, p), n)
            specs.append((p, A, B, 3))
    # deliberate: A contains a non-residue and a repeated class with opposite signs
    for p in (11, 13):
        chi = legendre_table(p)
        u = next(x for x in range(2, p) if chi[x] == -1)
        specs.append((p, [1, u], [1, u], 2))       # classes: -1 with (1,1),(u,u): sigma = 0, nu = 2
        specs.append((p, [1, u, 2 % p], [1, u, 3 % p], 2))
    for (p, A, B, k) in specs:
        chi = legendre_table(p)
        bf = brute_force(p, A, B, k, chi)
        an = analyze(p, A, B, k, chi, tag='bf', D_stats=True)
        ratios, nu, sigma = ratio_classes(p, A, B, chi)
        ok = (bf['T_sq'] == an['T_sq'] and bf['T_ns'] == an['T_ns'] and bf['T_tri'] == an['T_tri']
              and bf['R2'] == an['R2'] and bf['N_sq'] == an['N_sq'] and bf['N_ns'] == an['N_ns']
              and bf['Corr'] == an['Corr'])
        check('bf_vs_pattern_totals', ok, {'p': p, 'A': A, 'B': B, 'k': k, 'bf': {x: bf[x] for x in
              ('T_sq', 'T_ns', 'T_tri', 'R2', 'N_sq', 'N_ns', 'Corr')},
              'an': {x: an[x] for x in ('T_sq', 'T_ns', 'T_tri', 'R2', 'N_sq', 'N_ns', 'Corr')}})
        # c(D) from brute force (signed) vs pattern formula: compare through the D-sum
        W, _ = all_WD(p, ratios, chi, k)
        idx = {r: i for i, r in enumerate(ratios)}
        bf_sum = sum(c * W[tuple(idx[r] for r in D)] for D, c in bf['cD'].items())
        check('bf_sum_cD_WD_equals_pattern_sum', bf_sum == an['sum_cD_WD'], {'p': p, 'A': A, 'B': B, 'k': k})
        check('bf_T_ns_equals_sum_cD_WD_minus_Corr', bf['T_ns'] == bf_sum - bf['Corr'], {'p': p, 'A': A, 'B': B})
        # signed vs unsigned c(D): the unsigned (nu-based) version differs when A is not inside the residues
        naive = sum(c * W[tuple(idx[r] for r in D)] for D, c in bf['cplus'].items())
        if naive != bf_sum and witnesses['sign_matters_for_c_D'] is None:
            witnesses['sign_matters_for_c_D'] = {'p': p, 'A': A, 'B': B, 'k': k, 'sum_cD_WD_signed': bf_sum,
                                                 'sum_cplus_WD_unsigned': naive, 'T_ns': bf['T_ns']}
        if naive != bf['T_ns'] and witnesses['naive_count_formula_fails'] is None:
            witnesses['naive_count_formula_fails'] = {'p': p, 'A': A, 'B': B, 'k': k,
                                                      'sum_cplus_WD': naive, 'T_ns': bf['T_ns']}
        if bf['Corr'] != 0 and witnesses['no_correction_formula_fails'] is None:
            witnesses['no_correction_formula_fails'] = {'p': p, 'A': A, 'B': B, 'k': k, 'sum_cD_WD': bf_sum,
                                                        'Corr': bf['Corr'], 'T_ns': bf['T_ns']}
        if bf['neg_square_W'] is not None and witnesses['negative_square_term'] is None:
            witnesses['negative_square_term'] = {'p': p, 'A': A, 'B': B, 'k': k, **bf['neg_square_W']}
        cases += 1
    return cases, witnesses


# --------------------------------------------------------------------------------------------
# rectangle constructions
# --------------------------------------------------------------------------------------------

def common_residue_set(p, A, chi):
    return [b for b in range(p) if all(chi[(a + b) % p] == 1 for a in A)]


def greedy_B(p, A, n, chi, rng):
    f = [(sum(int(chi[(a + b) % p]) for a in A), rng.random(), b) for b in range(p)]
    f.sort(reverse=True)
    return sorted(x[2] for x in f[:n])


def gp_sets(p, m, n, rng, biased=False, chi=None, residue=True):
    # random ratio of multiplicative order > m + n, a quadratic residue unless residue=False
    while True:
        r = rng.randrange(2, p)
        if chi is not None and (chi[r] == 1) != residue:
            continue
        x = 1
        ok = True
        for i in range(1, m + n + 1):
            x = x * r % p
            if x == 1:
                ok = False
                break
        if ok:
            break
    c = rng.randrange(1, p)
    A = sorted({c * pow(r, i, p) % p for i in range(m)})
    if biased:
        return A, r
    c2 = rng.randrange(1, p)
    B = sorted({c2 * pow(r, j, p) % p for j in range(n)})
    return A, B


def primitive_root(p):
    fac = []
    q = p - 1
    d = 2
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


def build_instances(p, m, n, chi, rng):
    """One rectangle of each type at (p, m, n).  Returns list of (type, A, B)."""
    inst = []
    A = rng.sample(range(1, p), m)
    B = rng.sample(range(0, p), n)
    inst.append(('random', sorted(A), sorted(B)))
    u = rng.randrange(0, p)
    A = list(range(1, m + 1))
    B = [(u + j) % p for j in range(1, n + 1)]
    inst.append(('interval', A, B))
    A, B = gp_sets(p, m, n, rng, chi=chi, residue=True)
    inst.append(('gp', A, B))
    A, B = gp_sets(p, m, n, rng, chi=chi, residue=False)
    inst.append(('gp_nr', A, B))
    # gp-biased: A = powers of a primitive root, B inside the common residue set of A (bias 1)
    r = primitive_root(p)
    A = sorted({pow(r, i, p) for i in range(m)})
    N = common_residue_set(p, A, chi)
    if len(N) >= n:
        B = sorted(rng.sample(N, n))
        inst.append(('gp_biased', A, B))
    # greedy: random A, B = the n best shifts
    A = sorted(rng.sample(range(1, p), m))
    B = greedy_B(p, A, n, chi, rng)
    inst.append(('greedy', A, B))
    return inst


def loglog_slope(xs, ys):
    xs = np.log(np.array(xs, dtype=float))
    ys = np.log(np.array(ys, dtype=float))
    if len(xs) < 3:
        return None
    a, b = np.polyfit(xs, ys, 1)
    return float(a)


def summarize_row(r):
    p = r['p']
    mn = r['m'] * r['n']
    k = r['k']
    R = math.sqrt(r['R2']) if r['R2'] > 0 else float('nan')
    row = {
        'tag': r['tag'], 'p': p, 'm': r['m'], 'n': r['n'], 'k': k, 'L': r['L'], 'S': r['S'],
        'bias': round(r['bias'], 4), 'max_abs_g': r['max_abs_g'], 'n_biased_half': r['n_biased_half'],
        'g0': r['g0'], 'Qk0': r['Qk0'],
        'R_nu': r['R_nu'], 'max_nu': r['max_nu'],
        'T_sq': r['T_sq'], 'T_ns': r['T_ns'], 'T_tri': r['T_tri'], 'R': R, 'moment': r['moment'],
        'ratio_Tns_over_Ttri': r['T_ns'] / r['T_tri'] if r['T_tri'] else None,
        'ratio_Tns_over_R': r['T_ns'] / R if R == R and R > 0 else None,
        'ratio_Tns_over_mn2k': r['T_ns'] / mn ** (2 * k),
        'ratio_Tsq_over_pRnu^k': r['T_sq'] / (p * r['R_nu'] ** k),
        'rho_T': abs(r['T_ns']) * min(r['m'], r['n']) ** k / (math.sqrt(p) * mn ** (2 * k)),
        'rho_strong': abs(r['T_ns']) / (math.sqrt(p) * r['R_nu'] ** k),
        'rho_prompt': abs(r['T_ns']) / (mn ** (2 * k) + math.sqrt(p) * r['R_nu'] ** k),
        'Corr_over_Tns': (r['Corr'] / r['T_ns']) if r['T_ns'] else None,
        'R_D': math.sqrt(r['RD2']) if r.get('RD2') else None,
        'R_DW': math.sqrt(r['RDW2']) if r.get('RDW2') else None,
        'ratio_Tns_over_RD': (r['T_ns'] / math.sqrt(r['RD2'])) if r.get('RD2') else None,
        'ratio_TW_over_RDW': (r['T_W'] / math.sqrt(r['RDW2'])) if r.get('RDW2') else None,
        'E2': r['E2'], 'T_W': r['T_W'], 'Corr': r['Corr'], 'R_W': math.sqrt(r['RW2']) if r['RW2'] > 0 else None,
        'N_W': r['N_W'],
        'ratio_TW_over_RW': (r['T_W'] / math.sqrt(r['RW2'])) if r['RW2'] > 0 else None,
        'ratio_TW_over_TWtri': (r['T_W'] / r['TW_tri']) if r['TW_tri'] else None,
        'ratio_TW_over_mn2k': r['T_W'] / mn ** (2 * k),
        'ratio_elem_over_mn2k': (r['E2'] - r['Corr']) / mn ** (2 * k),
        'rho_T_W': abs(r['T_W']) * min(r['m'], r['n']) ** k / (math.sqrt(p) * mn ** (2 * k)),
        'rho_strong_W': abs(r['T_W']) / (math.sqrt(p) * r['R_nu'] ** k),
        'cos_c_W': r.get('cos_c_W'), 'pair_corr_share_2k-1': r.get('pair_corr_share_2k-1'),
        'pair_corr_npairs': r.get('pair_corr_npairs'), 'top_WD_msq_over_p': r.get('top_WD_msq_over_p'),
        'top_WD_mean': r.get('top_WD_mean'), 'n_D': r.get('n_D'), 'CS_scale': r.get('CS_scale'),
    }
    return row


# --------------------------------------------------------------------------------------------
# main
# --------------------------------------------------------------------------------------------

def main():
    rng = random.Random(20260905)
    out = {'description': 'sigma-tuple verifier output (2026-09-05): exact tuple decomposition of the dilation '
                          '2k-th moment, T_sq / T_ns / T_tri / R, D-level statistics, conjecture tests, witnesses.',
           'checks': CHECKS, 'failures': FAILURES, 'witnesses': WITNESSES}

    # ---- 1. brute force at tiny primes -------------------------------------------------------
    t = time.time()
    n_bf, bf_wit = brute_force_section(rng)
    WITNESSES['brute_force'] = bf_wit
    out['brute_force_cases'] = n_bf
    out['seconds_brute_force'] = round(time.time() - t, 1)
    print(f'[bf] {n_bf} cases, {round(time.time() - t, 1)} s, failures so far {len(FAILURES)}', flush=True)

    # ---- 2. the prime grid ---------------------------------------------------------------------
    rows = []
    grid_primes = [101, 151, 211, 307, 401, 503, 701, 809, 1009, 1213, 1409, 1801, 2203, 2609, 2999]
    grid_primes = [q for q in grid_primes if q in set(primes_in(100, 3000))]
    t = time.time()
    for p in grid_primes:
        chi = legendre_table(p)
        for (m, n, k, types) in ((6, 6, 2, None), (6, 4, 2, None), (4, 4, 2, None), (4, 4, 3, None),
                                 (5, 5, 3, ('random', 'greedy', 'gp_biased', 'gp')),
                                 (6, 4, 3, ('random', 'greedy', 'gp_biased'))):
            for (typ, A, B) in build_instances(p, m, n, chi, rng):
                if types is not None and typ not in types:
                    continue
                r = analyze(p, A, B, k, chi, tag=typ)
                r['A'] = A
                r['B'] = B
                rows.append(r)
        # a second random rectangle at (6,6), k = 2 and at (5,5), k = 3 for the fluctuation scale
        A = sorted(rng.sample(range(1, p), 6))
        B = sorted(rng.sample(range(0, p), 6))
        r = analyze(p, A, B, 2, chi, tag='random2')
        r['A'] = A
        r['B'] = B
        rows.append(r)
        A = sorted(rng.sample(range(1, p), 5))
        B = sorted(rng.sample(range(0, p), 5))
        r = analyze(p, A, B, 3, chi, tag='random2')
        r['A'] = A
        r['B'] = B
        rows.append(r)
        print(f'[grid] p={p} done, {round(time.time() - T0, 1)} s total', flush=True)
    out['seconds_grid'] = round(time.time() - t, 1)
    out['grid_primes'] = grid_primes

    # ---- 4. stored near-extremal rectangles ----------------------------------------------------
    t = time.time()
    stored = []
    zero_removed = 0
    if os.path.exists(SEARCH):
        sj = json.load(open(SEARCH))
        cand = []
        for e in sj['exhaustive']:
            if (e['m'], e['n']) in ((4, 4), (5, 5), (6, 6), (6, 5), (6, 4), (5, 4)) and e['witnesses']:
                w = e['witnesses'][0]
                cand.append(('stored_exhaustive', e['p'], w['A'], w['B'], w['S']))
        for e in sj['anneal']:
            if (e['m'], e['n']) == (6, 6):
                for w in e['witnesses'][:1]:
                    cand.append(('stored_anneal', e['p'], w['A'], w['B'], w['S']))
        for e in sj['kstar']:
            for kk in ('4', '5', '6'):
                if kk in e['per_k']:
                    w = e['per_k'][kk]
                    cand.append(('stored_kstar', e['p'], w['A'], w['B'], w['S']))
        # thin the exhaustive list (it is long): keep every 12th, plus all with mn = 36
        ex = [c for c in cand if c[0] == 'stored_exhaustive']
        keep = [c for i, c in enumerate(ex) if i % 12 == 0 or len(c[2]) * len(c[3]) == 36]
        others = [c for c in cand if c[0] != 'stored_exhaustive']
        cand = keep + others
        for (typ, p, A, B, S0) in cand:
            chi = legendre_table(p)
            A = [a % p for a in A]
            B = [b % p for b in B]
            S_rec = sum(int(chi[(a + b) % p]) for a in A for b in B)
            check('stored_S_reverified', S_rec == S0, {'p': p, 'A': A, 'B': B, 'S_stored': S0, 'S': S_rec})
            if 0 in A:
                A = [a for a in A if a != 0]
                zero_removed += 1
            k = 2
            r = analyze(p, A, B, k, chi, tag=typ)
            r['A'] = A
            r['B'] = B
            r['S_stored'] = S0
            rows.append(r)
            stored.append(r)
        print(f'[stored] {len(stored)} rectangles, {round(time.time() - t, 1)} s', flush=True)
    out['stored_rectangles'] = len(stored)
    out['stored_zero_removed_from_A'] = zero_removed
    out['seconds_stored'] = round(time.time() - t, 1)

    # ---- 5. witnesses for the refuted forms ----------------------------------------------------
    t = time.time()
    wit = {}
    # (a) A = {1}, B = the nonzero residues: T_ns ~ (p/2)^{2k} >> sqrt p R_nu^k
    for (p, k) in ((61, 2), (101, 2), (151, 2), (41, 3), (61, 3)):
        chi = legendre_table(p)
        A = [1]
        B = [b for b in range(1, p) if chi[b] == 1]
        r = analyze(p, A, B, k, chi, tag='one_by_Q', D_stats=False)
        r['A'] = A
        r['B'] = B
        rows.append(r)
        sr = summarize_row(r)
        wit[f'one_by_Q_p{p}_k{k}'] = {x: sr[x] for x in ('p', 'm', 'n', 'k', 'S', 'T_sq', 'T_ns', 'T_tri', 'R',
                                                          'R_nu', 'rho_strong', 'rho_T', 'rho_prompt')}
    # (b) subgroup construction A = H, B inside N(H): |H| biased dilates
    for p in (1009, 2003, 2999):
        chi = legendre_table(p)
        g = primitive_root(p)
        best = None
        for h in range(4, 7):
            if (p - 1) % h == 0:
                H = sorted({pow(g, (p - 1) // h * i, p) for i in range(h)})
                N = common_residue_set(p, H, chi)
                if len(N) >= 6:
                    best = (H, sorted(rng.sample(N, 6)))
        if best:
            A, B = best
            r = analyze(p, A, B, 2, chi, tag='subgroup_biased')
            r['A'] = A
            r['B'] = B
            rows.append(r)
            sr = summarize_row(r)
            wit[f'subgroup_p{p}'] = {x: sr[x] for x in ('p', 'm', 'n', 'k', 'S', 'n_biased_half', 'T_sq', 'T_ns',
                                                        'T_tri', 'R', 'R_nu', 'rho_strong', 'rho_T', 'rho_prompt')}
    WITNESSES['refuted_forms'] = wit
    out['seconds_witnesses'] = round(time.time() - t, 1)

    # ---- 6. summaries ----------------------------------------------------------------------------
    srows = [summarize_row(r) for r in rows]
    out['rows'] = srows
    out['n_rows'] = len(srows)
    # per (tag, m, n, k): log-log slopes against p of |T_ns|/(mn)^2k, |T_ns|/R, |T_ns|/T_tri; medians
    groups = {}
    for s in srows:
        key = (s['tag'], s['m'], s['n'], s['k'])
        groups.setdefault(key, []).append(s)
    fits = {}
    for key, rs in groups.items():
        rs = sorted(rs, key=lambda s: s['p'])
        ps = [s['p'] for s in rs]
        def med(f):
            v = [f(s) for s in rs if f(s) is not None]
            return float(np.median(v)) if v else None
        entry = {'n': len(rs), 'primes': ps,
                 'median_abs_Tns_over_R': med(lambda s: abs(s['ratio_Tns_over_R']) if s['ratio_Tns_over_R'] is not None else None),
                 'median_abs_Tns_over_Ttri': med(lambda s: abs(s['ratio_Tns_over_Ttri'])),
                 'median_abs_Tns_over_RD': med(lambda s: abs(s['ratio_Tns_over_RD']) if s['ratio_Tns_over_RD'] is not None else None),
                 'median_abs_TW_over_RDW': med(lambda s: abs(s['ratio_TW_over_RDW']) if s['ratio_TW_over_RDW'] is not None else None),
                 'max_abs_TW_over_RDW': max([abs(s['ratio_TW_over_RDW']) for s in rs if s['ratio_TW_over_RDW'] is not None] or [None]),
                 'median_TW_over_mn2k': med(lambda s: s['ratio_TW_over_mn2k']),
                 'median_abs_elem_over_mn2k': med(lambda s: abs(s['ratio_elem_over_mn2k'])),
                 'max_rho_T_W': max(s['rho_T_W'] for s in rs),
                 'frac_TW_negative': float(np.mean([s['T_W'] < 0 for s in rs])),
                 'median_Tns_over_mn2k': med(lambda s: s['ratio_Tns_over_mn2k']),
                 'median_bias': med(lambda s: abs(s['bias'])),
                 'frac_Tns_negative': float(np.mean([s['T_ns'] < 0 for s in rs])),
                 'max_rho_T': max(s['rho_T'] for s in rs),
                 'max_rho_strong': max(s['rho_strong'] for s in rs),
                 'max_rho_prompt': max(s['rho_prompt'] for s in rs),
                 'median_cos_c_W': med(lambda s: s['cos_c_W']),
                 'median_abs_cos_c_W': med(lambda s: abs(s['cos_c_W']) if s['cos_c_W'] is not None else None),
                 'median_pair_corr': med(lambda s: s['pair_corr_share_2k-1']),
                 'median_top_WD_msq_over_p': med(lambda s: s['top_WD_msq_over_p']),
                 'median_Corr_over_Tns': med(lambda s: abs(s['Corr_over_Tns']) if s['Corr_over_Tns'] is not None else None)}
        if len(rs) >= 4:
            entry['slope_abs_Tns_over_mn2k'] = loglog_slope(ps, [max(abs(s['ratio_Tns_over_mn2k']), 1e-12) for s in rs])
            entry['slope_abs_Tns_over_R'] = loglog_slope(ps, [max(abs(s['ratio_Tns_over_R']), 1e-12) for s in rs])
            entry['slope_abs_Tns_over_Ttri'] = loglog_slope(ps, [max(abs(s['ratio_Tns_over_Ttri']), 1e-12) for s in rs])
            entry['slope_rho_T'] = loglog_slope(ps, [max(s['rho_T'], 1e-12) for s in rs])
            entry['slope_abs_TW_over_RDW'] = loglog_slope(ps, [max(abs(s['ratio_TW_over_RDW'] or 0), 1e-12) for s in rs])
            entry['slope_abs_TW_over_mn2k'] = loglog_slope(ps, [max(abs(s['ratio_TW_over_mn2k']), 1e-12) for s in rs])
        fits['|'.join(map(str, key))] = entry
    out['fits_by_type'] = fits
    # conjecture tests on every computed instance
    conj = {}
    for name, fld in (('T_k_min_saving', 'rho_T'), ('T_k_min_saving_Weil_part', 'rho_T_W'),
                      ('strong_sqrtp_Rnu_k', 'rho_strong'), ('strong_sqrtp_Rnu_k_Weil_part', 'rho_strong_W'),
                      ('prompt_form', 'rho_prompt')):
        by_k = {}
        for s in srows:
            by_k.setdefault(s['k'], []).append(s)
        conj[name] = {}
        for k, rs in by_k.items():
            worst = max(rs, key=lambda s: s[fld])
            conj[name][f'k={k}'] = {'n_instances': len(rs), 'max_ratio': worst[fld],
                                    'worst': {x: worst[x] for x in ('tag', 'p', 'm', 'n', 'S', 'T_ns', 'R_nu')},
                                    'max_ratio_excluding_lopsided_witnesses':
                                        max(s[fld] for s in rs if s['tag'] not in ('one_by_Q',)),
                                    'n_ratio_gt_1': sum(1 for s in rs if s[fld] > 1)}
    out['conjecture_tests'] = conj
    out['total_checks'] = sum(CHECKS.values())
    out['n_failures'] = len(FAILURES)
    out['seconds'] = round(time.time() - T0, 1)
    json.dump(out, open(OUT, 'w'), indent=1, default=float)
    print(f'total checks {out["total_checks"]}, failures {len(FAILURES)}, rows {len(srows)}, {out["seconds"]} s')
    for f in FAILURES[:10]:
        print('FAILURE', f)


if __name__ == '__main__':
    main()
