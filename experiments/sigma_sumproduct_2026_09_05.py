#!/usr/bin/env python3
"""
sigma_sumproduct_2026_09_05.py -- verifier for research/sigma-sumproduct-2026-09-05.md

Exact-integer checks of the explicit Burgess--Chang chain (Theorem E and Lemmas A-D, G, G'),
the exponent bookkeeping (Corollary 3.1), the no-go Proposition F with witnesses, and
measurements of the expansion factor lambda and the ratio-energy deviation Phi for random,
interval, subgroup-coset and mixed sets.  Standard library + numpy only.

Writes results/sigma_sumproduct_2026_09_05.json.  Runs in well under ten minutes.
"""
import json, os, sys, time, math, random
from fractions import Fraction
from collections import Counter
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, '..', 'results', 'sigma_sumproduct_2026_09_05.json')
T0 = time.time()
RES = {'checks': 0, 'failures': [], 'witnesses': [], 'sections': {}}


def check(name, cond, info=None):
    RES['checks'] += 1
    if not cond:
        RES['failures'].append({'name': name, 'info': info})
        print('FAIL', name, info)


def legendre_table(p):
    """chi(a) for a in 0..p-1 as an int8 numpy array, chi(0)=0."""
    t = np.zeros(p, dtype=np.int64)
    for a in range(1, p):
        v = pow(a, (p - 1) // 2, p)
        t[a] = 1 if v == 1 else -1
    return t


def dfact(n):  # (2r-1)!! with n = 2r-1
    return math.prod(range(1, n + 1, 2)) if n >= 1 else 1


# ---------------------------------------------------------------- set constructors
def primitive_root(p):
    fac = []
    m = p - 1
    d = 2
    while d * d <= m:
        if m % d == 0:
            fac.append(d)
            while m % d == 0:
                m //= d
        d += 1
    if m > 1:
        fac.append(m)
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fac):
            return g
    raise ValueError


def subgroup_coset(p, d, shift_mult=1):
    """coset shift_mult * H where H is the subgroup of order d (d | p-1)."""
    g = primitive_root(p)
    h = pow(g, (p - 1) // d, p)
    H = []
    x = 1
    for _ in range(d):
        H.append(x * shift_mult % p)
        x = x * h % p
    return sorted(set(H))


def rand_set(p, n, rng):
    return sorted(rng.sample(range(p), n))


# ---------------------------------------------------------------- core quantities
def S_sum(chi, p, A, B):
    A = np.array(A); B = np.array(B)
    return int(chi[(A[:, None] + B[None, :]) % p].sum())


def fA_vals(chi, p, A, ys):
    A = np.array(A); ys = np.array(ys)
    return chi[(ys[:, None] + A[None, :]) % p].sum(axis=1)


def inv_table(p):
    inv = np.zeros(p, dtype=np.int64)
    for z in range(1, p):
        inv[z] = pow(z, p - 2, p)
    return inv


def nu_counts(p, A, Y, Z, inv, chunk_elems=4_000_000):
    """nu(u1,u2) for all (u1,u2): returns (keys, counts) of the nonzero entries, computed in chunks
    over (y,z) pairs so that memory stays bounded.  Exact integer counts."""
    A = np.array(A, dtype=np.int64); Y = np.array(Y, dtype=np.int64); Z = np.array(Z, dtype=np.int64)
    zinv = inv[Z]
    ys = np.repeat(Y, len(Z)); zi = np.tile(zinv, len(Y))  # all (y,z) pairs
    acc = np.zeros(p * p, dtype=np.int64)
    rows_per_chunk = max(1, chunk_elems // (len(A) * len(A)))
    for s in range(0, len(ys), rows_per_chunk):
        yb = ys[s:s + rows_per_chunk]; zb = zi[s:s + rows_per_chunk]
        W = ((yb[:, None] + A[None, :]) % p) * zb[:, None] % p  # rows x |A|
        keys = (W[:, :, None] * p + W[:, None, :]).ravel()
        acc += np.bincount(keys, minlength=p * p)
    nz = np.nonzero(acc)[0]
    return nz, acc[nz]


def g_table(chi, p, H):
    """g(u1,u2) = sum_{t=1..H} chi(u1+t)chi(u2+t) as a p x p int array."""
    u = np.arange(p)
    M = np.stack([chi[(u + t) % p] for t in range(1, H + 1)], axis=1)  # p x H
    return M @ M.T  # p x p


def add_energy(p, X):
    X = np.array(X)
    d = (X[:, None] - X[None, :]) % p
    c = np.bincount(d.ravel(), minlength=p)
    return int((c.astype(np.int64) ** 2).sum())


def mult_energy(p, X):
    X = [x for x in X if x % p != 0]
    X = np.array(X)
    pr = (X[:, None] * X[None, :]) % p
    c = np.bincount(pr.ravel(), minlength=p)
    return int((c.astype(np.int64) ** 2).sum())


def ratio_count_T(p, A, Y, Z, inv):
    """T(A,Y,Z) = #{(x,x',y,y',z,z') : (x+y)/z = (x'+y')/z'} = sum_w n(w)^2."""
    A = np.array(A); Y = np.array(Y); Z = np.array(Z)
    W = (((A[:, None] + Y[None, :]) % p)[:, :, None] * inv[Z][None, None, :]) % p
    c = np.bincount(W.ravel(), minlength=p)
    return int((c.astype(np.int64) ** 2).sum())


def max_pair_count(p, A, Y, Z, inv):
    """max_{x,x'} #{(y,y',z,z') : (x+y)/z = (x'+y')/z'}  = max_{x,x'} sum_w n_x(w) n_{x'}(w)."""
    Y = np.array(Y); Z = np.array(Z)
    rows = []
    for x in A:
        W = (((x + Y) % p)[:, None] * inv[Z][None, :]) % p
        rows.append(np.bincount(W.ravel(), minlength=p).astype(np.int64))
    R = np.stack(rows)
    G = R @ R.T
    return int(G.max())


# ================================================================ 5.1 Lemma D
def section_lemmaD():
    out = []
    for p in [101, 199, 293]:
        chi = legendre_table(p)
        for H in [3, 5, 8]:
            G = g_table(chi, p, H)
            for r in [1, 2, 3]:
                Wr = int((G.astype(object) ** (2 * r)).sum()) if p <= 101 else int((G.astype(np.int64) ** (2 * r)).sum())
                bound = dfact(2 * r - 1) * H ** r * p ** 2 + (2 * r - 1) ** 2 * H ** (2 * r) * p
                check(f'LemmaD p={p} H={H} r={r}', Wr <= bound, {'Wr': Wr, 'bound': bound})
                # cross-check via t-vector expansion for r=1 (exact formula) and r=2
                if r == 1:
                    exact = H * (p - 1) ** 2 + H * (H - 1)
                    check(f'LemmaD-exact-r1 p={p} H={H}', Wr == exact, {'Wr': Wr, 'exact': exact})
                if r == 2 and H <= 5:
                    tot = 0
                    u = np.arange(p)
                    for t1 in range(1, H + 1):
                        for t2 in range(1, H + 1):
                            for t3 in range(1, H + 1):
                                for t4 in range(1, H + 1):
                                    F = ((u + t1) * (u + t2) % p) * ((u + t3) * (u + t4) % p) % p
                                    w = int(chi[F].sum())
                                    tot += w * w
                    check(f'LemmaD-tvec-r2 p={p} H={H}', tot == Wr, {'tvec': tot, 'Wr': Wr})
                out.append({'p': p, 'H': H, 'r': r, 'W_r': Wr, 'bound': bound, 'ratio': Wr / bound})
    RES['sections']['5.1_lemmaD'] = out
    print('5.1 done', time.time() - T0)


# ================================================================ 5.2 Theorem E, Lemmas A,B,C,G,G'
def theoremE_instance(p, chi, inv, A, B, Z, H, r, label):
    I = list(range(1, H + 1))
    Bset = set(B)
    Y = sorted({(b - z * t) % p for b in B for z in Z for t in I})
    lam = len(Y) / len(B)
    LHS = abs(S_sum(chi, p, A, B))
    # Lemma A
    fB = np.abs(fA_vals(chi, p, A, B)).sum()
    box = 0
    fY = {}
    tot_sq = 0
    Yarr = np.array(Y)
    for z in Z:
        for t in I:
            vals = np.abs(fA_vals(chi, p, A, (Yarr + z * t) % p))
            box += int(vals.sum())
            tot_sq += int((vals ** 2).sum())
    check(f'LemmaA {label}', fB * len(Z) * H <= box, {'lhs': int(fB), 'rhs/ZH': box / (len(Z) * H)})
    # Lemma B identity: sum |f|^2 over box == sum nu * g
    uk, cnt = nu_counts(p, A, Y, Z, inv)
    G = g_table(chi, p, H)
    u1 = uk // p; u2 = uk % p
    nug = int((cnt * G[u1, u2]).sum())
    check(f'LemmaB-identity {label}', nug == tot_sq, {'nu.g': nug, 'sum|f|^2': tot_sq})
    Snu = int(cnt.sum()); Snu2 = int((cnt.astype(np.int64) ** 2).sum())
    check(f'Snu {label}', Snu == len(A) ** 2 * len(Y) * len(Z))
    # Lemma C on the actual data
    gabs = np.abs(G[u1, u2]).astype(np.float64)
    lhsC = float((cnt * gabs).sum())
    W_r_actual = float((np.abs(G).astype(np.float64) ** (2 * r)).sum())
    rhsC = Snu ** (1 - 1 / r) * Snu2 ** (1 / (2 * r)) * W_r_actual ** (1 / (2 * r))
    check(f'LemmaC {label}', lhsC <= rhsC * (1 + 1e-9), {'lhs': lhsC, 'rhs': rhsC})
    # Theorem E (first display, with Lemma D's bound)
    Wbound = dfact(2 * r - 1) * H ** r * p ** 2 + (2 * r - 1) ** 2 * H ** (2 * r) * p
    RHS2 = (len(A) ** 2 * len(Y) ** 2 / H) * (Snu2 / Snu ** 2) ** (1 / (2 * r)) * Wbound ** (1 / (2 * r))
    check(f'TheoremE {label}', LHS ** 2 <= RHS2 * (1 + 1e-9), {'LHS^2': LHS ** 2, 'RHS': RHS2})
    Phi = p ** 2 * Snu2 / Snu ** 2
    bracket = dfact(2 * r - 1) ** (1 / (2 * r)) * H ** -0.5 + (2 * r - 1) ** (1 / r) * p ** (-1 / (2 * r))
    RHS_E = lam * len(A) * len(B) * Phi ** (1 / (4 * r)) * bracket ** 0.5
    check(f'TheoremE-form2 {label}', LHS <= RHS_E * (1 + 1e-9), {'LHS': LHS, 'RHS': RHS_E})
    check(f'Phi>=1 {label}', Phi >= 1 - 1e-12)
    # Lemma G and G'
    Znz = [z for z in Z if z % p != 0]
    M = max_pair_count(p, A, Y, Znz, inv)
    EZ = mult_energy(p, Znz)
    maxE = max(mult_energy(p, [(x + y) % p for y in Y]) for x in A)
    check(f'LemmaG-1 {label}', Snu2 <= len(A) ** 3 * M, {'Snu2': Snu2, '|A|^3 M': len(A) ** 3 * M})
    check(f'LemmaG-2 {label}', M <= math.sqrt(EZ * maxE) * (1 + 1e-12), {'M': M, 'sqrt(EZ maxE)': math.sqrt(EZ * maxE)})
    T = ratio_count_T(p, A, Y, Znz, inv)
    check(f"LemmaG' {label}", Snu2 <= len(A) * T, {'Snu2': Snu2, '|A| T': len(A) * T})
    # Proposition F lower bound on |Y|
    Hp = sorted({z * t % p for z in Z for t in I})
    EB = add_energy(p, B); EH = add_energy(p, Hp)
    lowerY = len(B) ** 2 * len(Hp) ** 2 / math.sqrt(EB * EH)
    check(f'PropF-|Y| {label}', len(Y) >= lowerY - 1e-9, {'|Y|': len(Y), 'lower': lowerY})
    check(f'PropF-RHS {label}', RHS_E >= len(A) * len(B) * len(B) / math.sqrt(EB) - 1e-9,
          {'RHS_E': RHS_E, 'floor': len(A) * len(B) * len(B) / math.sqrt(EB)})
    return {'label': label, 'p': p, '|A|': len(A), '|B|': len(B), '|Z|': len(Z), 'H': H, 'r': r,
            '|Y|': len(Y), 'lambda': lam, '|S(A,B)|': LHS, 'trivial': len(A) * len(B),
            'Snu': Snu, 'Snu2': Snu2, 'Phi': Phi, 'RHS_E': RHS_E, 'RHS_E/trivial': RHS_E / (len(A) * len(B)),
            'E+(B)': EB, 'C=E+(B)/|B|^2': EB / len(B) ** 2, 'LemmaG_M': M, 'T': T,
            'Snu2/(|A|^2|B|^3)': Snu2 / (len(A) ** 2 * len(B) ** 3)}


def section_theoremE():
    rng = random.Random(20260905)
    out = []
    for p in [101, 197, 293]:
        chi = legendre_table(p); inv = inv_table(p)
        n = int(round(p ** 0.45))
        # structured: intervals, Z a short interval
        A = list(range(1, n + 1)); B = list(range(5, 5 + n))
        Z = list(range(1, 4)); H = 3
        for r in [1, 2]:
            out.append(theoremE_instance(p, chi, inv, A, B, Z, H, r, f'int-int p={p} r={r}'))
        # arbitrary A, interval B
        A = rand_set(p, n, rng)
        out.append(theoremE_instance(p, chi, inv, A, B, Z, H, 2, f'rand-int p={p} r=2'))
        # random both, small shift set
        B2 = rand_set(p, n, rng); Z2 = rand_set(p - 1, 3, rng); Z2 = [z + 1 for z in Z2]
        out.append(theoremE_instance(p, chi, inv, A, B2, Z2, 3, 2, f'rand-rand p={p} r=2'))
        # coset B, Z inside the coset (multiplicatively structured shift)
        d = max(dd for dd in range(2, p) if (p - 1) % dd == 0 and dd <= n + 5)
        Bc = subgroup_coset(p, d, 3)
        Zc = Bc[:3]
        out.append(theoremE_instance(p, chi, inv, A, Bc, Zc, 3, 2, f'rand-coset p={p} r=2'))
    RES['sections']['5.2_theoremE'] = out
    print('5.2 done', time.time() - T0)


# ================================================================ 5.3 exponent bookkeeping
def section_exponents():
    rows = {
        'Chang Prop1 (3,11/4)': (Fraction(3), Fraction(11, 4)),
        'AMRS via Lemma G (3,5/2)': (Fraction(3), Fraction(5, 2)),
        "AMRS via Lemma G' (5/2,3)": (Fraction(5, 2), Fraction(3)),
        'Alsetri-Shao rank<=2 (3,2)': (Fraction(3), Fraction(2)),
        'diagonal floor (2,3)': (Fraction(2), Fraction(3)),
        'absolute floor (2,2) [condition |A||B|>p^{1/2}]': (Fraction(2), Fraction(2)),
    }
    out = {}
    for k, (a, b) in rows.items():
        theta = Fraction(1) / (8 - a - b)
        out[k] = str(theta)
    check('theta(3,11/4)=4/9', out['Chang Prop1 (3,11/4)'] == '4/9')
    check('theta(3,5/2)=2/5', out['AMRS via Lemma G (3,5/2)'] == '2/5')
    check("theta(5/2,3)=2/5", out["AMRS via Lemma G' (5/2,3)"] == '2/5')
    check('theta(3,2)=1/3', out['Alsetri-Shao rank<=2 (3,2)'] == '1/3')
    check('theta(2,3)=1/3', out['diagonal floor (2,3)'] == '1/3')
    # Chang's own arithmetic: (4/9)*(9/8) = 1/2 exactly
    check('Chang 4/9*9/8=1/2', Fraction(4, 9) * Fraction(9, 8) == Fraction(1, 2))
    # (2r-1)!! <= r^r for r<=12 (used only to compare with Chang/SV constants)
    check('dfact<=r^r', all(dfact(2 * r - 1) <= r ** r for r in range(1, 13)))
    RES['sections']['5.3_exponents'] = out
    print('5.3 done', time.time() - T0)


# ================================================================ 5.4 Proposition F witnesses
def section_propF():
    rng = random.Random(1)
    out = []
    for p in [1009, 2003]:
        chi = legendre_table(p); inv = inv_table(p)
        for alpha in [0.40, 0.45]:
            n = int(round(p ** alpha))
            B = rand_set(p, n, rng); A = rand_set(p, n, rng)
            EB = add_energy(p, B); C = EB / n ** 2
            floor = n / math.sqrt(EB)  # RHS(E)/(|A||B|) >= this, for every Z,I,r
            rec = {'p': p, 'alpha': alpha, '|B|': n, 'B': B, 'E+(B)': EB, 'C': C,
                   'chain_floor_RHS_over_trivial': floor, 'shift_trials': []}
            # a few explicit (Z,I,r) to show lambda, Phi and the actual RHS(E)
            for (m, H, r) in [(1, 2, 1), (3, 3, 2), (6, 8, 3), (12, 16, 4)]:
                Z = [z + 1 for z in rand_set(p - 1, m, rng)]
                I = range(1, H + 1)
                Y = sorted({(b - z * t) % p for b in B for z in Z for t in I})
                lam = len(Y) / n
                Hp = sorted({z * t % p for z in Z for t in I})
                lowerY = n ** 2 * len(Hp) ** 2 / math.sqrt(EB * add_energy(p, Hp))
                check(f'PropF-|Y| p={p} a={alpha} m={m} H={H}', len(Y) >= lowerY - 1e-9)
                trial = {'|Z|': m, 'H': H, 'r': r, '|Y|': len(Y), 'lambda': lam, 'lambda_lower': lowerY / n}
                if len(Y) * len(Z) * n * n <= 6_000_000:
                    uk, cnt = nu_counts(p, A, Y, Z, inv)
                    Snu = int(cnt.sum()); Snu2 = int((cnt.astype(np.int64) ** 2).sum())
                    Phi = p ** 2 * Snu2 / Snu ** 2
                    bracket = dfact(2 * r - 1) ** (1 / (2 * r)) * H ** -0.5 + (2 * r - 1) ** (1 / r) * p ** (-1 / (2 * r))
                    RHS_E = lam * n * n * Phi ** (1 / (4 * r)) * bracket ** 0.5
                    S = abs(S_sum(chi, p, A, B))
                    check(f'TheoremE-large p={p} a={alpha} m={m}', S <= RHS_E * (1 + 1e-9))
                    check(f'PropF-RHS-large p={p} a={alpha} m={m}', RHS_E >= n * n * floor - 1e-9)
                    trial.update({'Phi': Phi, 'RHS_E/trivial': RHS_E / (n * n), '|S|/trivial': S / (n * n)})
                rec['shift_trials'].append(trial)
            out.append(rec)
            RES['witnesses'].append({'hypothesis': 'H_chain(alpha): some (Z,I,r) makes bound (E) <= p^{-delta}|A||B| for all A,B of size p^alpha',
                                     'status': 'REFUTED by Proposition F', 'p': p, 'alpha': alpha, '|B|': n,
                                     'B': B, 'E+(B)': EB, 'E+(B)/|B|^2': C,
                                     'chain_bound_over_trivial_is_at_least': floor})
    RES['sections']['5.4_propF'] = out
    print('5.4 done', time.time() - T0)


# ================================================================ 5.5 lambda / Phi measurements (E* test)
def measure(p, chi, inv, lab, A, B, Z, H, r=2, do_true=True):
    I = range(1, H + 1)
    Y = sorted({(b - z * t) % p for b in B for z in Z for t in I})
    lam = len(Y) / len(B)
    uk, cnt = nu_counts(p, A, B, Z, inv)  # idealised: Y := B
    Snu = int(cnt.sum()); Snu2 = int((cnt.astype(np.int64) ** 2).sum())
    Phi_ideal = p ** 2 * Snu2 / Snu ** 2
    rec = {'p': p, 'config': lab, '|A|': len(A), '|B|': len(B), '|Z|': len(Z), 'H': H,
           '|Y|': len(Y), 'lambda': lam, 'Snu_ideal': Snu, 'Snu2_ideal': Snu2, 'Phi_ideal(Y=B)': Phi_ideal,
           'Snu2/(|A|^2|B|^3)': Snu2 / (len(A) ** 2 * len(B) ** 3),
           'Snu2/max(Snu,Snu^2/p^2)': Snu2 / max(Snu, Snu ** 2 / p ** 2),
           'log_p(Phi_ideal)': math.log(Phi_ideal) / math.log(p),
           'Phi_ideal_over_AMRS_G': Phi_ideal / (p ** 2 * len(A) ** 3 * len(B) ** 2.5 / (len(A) ** 4 * len(B) ** 4)),
           'Phi_ideal_over_Chang': Phi_ideal / (p ** 2 * len(A) ** 3 * len(B) ** 2.75 / (len(A) ** 4 * len(B) ** 4))}
    if do_true and len(Y) * len(Z) * len(A) ** 2 <= 400_000_000:
        uk2, cnt2 = nu_counts(p, A, Y, Z, inv)
        Snu_t = int(cnt2.sum()); Snu2_t = int((cnt2.astype(np.int64) ** 2).sum())
        Phi_t = p ** 2 * Snu2_t / Snu_t ** 2
        bracket = dfact(2 * r - 1) ** (1 / (2 * r)) * H ** -0.5 + (2 * r - 1) ** (1 / r) * p ** (-1 / (2 * r))
        RHS_E = lam * len(A) * len(B) * Phi_t ** (1 / (4 * r)) * bracket ** 0.5
        S = abs(S_sum(chi, p, A, B))
        check(f'TheoremE-meas p={p} {lab}', S <= RHS_E * (1 + 1e-9), {'S': S, 'RHS': RHS_E})
        rec.update({'Phi_true': Phi_t, 'Snu2_true': Snu2_t, 'r': r, 'RHS_E/trivial': RHS_E / (len(A) * len(B)),
                    '|S|/trivial': S / (len(A) * len(B)), 'lambda^(4r)*Phi_true/p': lam ** (4 * r) * Phi_t / p})
    return rec


def section_measurements():
    rng = random.Random(7)
    out = []
    for p in [1009, 2003]:
        chi = legendre_table(p); inv = inv_table(p)
        for alpha in [0.60, 0.68]:
            n = int(round(p ** alpha))
            M, H = max(2, n // 10), max(2, n // 10)
            configs = []
            Bi = list(range(1, n + 1)); Zi = list(range(1, M + 1))
            configs.append(('interval-B / random-A', rand_set(p, n, rng), Bi, Zi))
            configs.append(('interval-B / interval-A', list(range(1, n + 1)), Bi, Zi))
            configs.append(('interval-B / AP-A(step 7)', [7 * k % p for k in range(1, n + 1)], Bi, Zi))
            Br = rand_set(p, n, rng); Zr = [z + 1 for z in rand_set(p - 1, M, rng)]
            configs.append(('random-B / random-A', rand_set(p, n, rng), Br, Zr))
            d = max(dd for dd in range(2, p) if (p - 1) % dd == 0 and dd <= n + 8)
            Bc = subgroup_coset(p, d, 5); Zc = Bc[:M]
            configs.append(('coset-B / random-A', rand_set(p, n, rng), Bc, Zc))
            configs.append(('coset-B / interval-A', list(range(1, n + 1)), Bc, Zc))
            Bm = sorted(set(list(range(1, n // 2 + 1)) + rand_set(p, n - n // 2, rng)))
            configs.append(('mixed-B / random-A', rand_set(p, n, rng), Bm, Zi))
            configs.append(('interval-B / coset-A', Bc, Bi, Zi))
            for (lab, A, B, Z) in configs:
                rec = measure(p, chi, inv, lab + f' a={alpha}', A, B, Z, H)
                out.append(rec)
                print(f"  {p} {lab:28s} n={n:4d} lam={rec['lambda']:7.2f} Phi_ideal={rec['Phi_ideal(Y=B)']:10.2f} "
                      f"Snu2/(A^2B^3)={rec['Snu2/(|A|^2|B|^3)']:.3f} Phi_true={rec.get('Phi_true', float('nan')):10.2f} "
                      f"RHS/triv={rec.get('RHS_E/trivial', float('nan')):8.2f}", flush=True)
    RES['sections']['5.5_measurements'] = out
    print('5.5 done', time.time() - T0, flush=True)


def section_scaling():
    """Size scaling of the idealised ratio energy Snu2 (Y := B, Z = first |B|/10 of B - shift) at fixed p."""
    rng = random.Random(11)
    p = 2003
    chi = legendre_table(p); inv = inv_table(p)
    out = []
    for kind in ['interval-B/interval-A', 'interval-B/random-A', 'random-B/random-A']:
        pts = []
        for n in [25, 50, 100, 200]:
            if kind == 'interval-B/interval-A':
                A = list(range(1, n + 1)); B = list(range(1, n + 1)); Z = list(range(1, max(2, n // 10) + 1))
            elif kind == 'interval-B/random-A':
                A = rand_set(p, n, rng); B = list(range(1, n + 1)); Z = list(range(1, max(2, n // 10) + 1))
            else:
                A = rand_set(p, n, rng); B = rand_set(p, n, rng); Z = [z + 1 for z in rand_set(p - 1, max(2, n // 10), rng)]
            uk, cnt = nu_counts(p, A, B, Z, inv)
            Snu = int(cnt.sum()); Snu2 = int((cnt.astype(np.int64) ** 2).sum())
            pts.append({'n': n, '|Z|': len(Z), 'Snu': Snu, 'Snu2': Snu2, 'Snu2/(n^2 n^3)': Snu2 / n ** 5,
                        'Snu2/(Snu^2/p^2)': Snu2 / (Snu ** 2 / p ** 2)})
        # empirical exponent of Snu2 in n between consecutive sizes (|Z| ~ n/10 scales too)
        ex = [math.log(pts[i + 1]['Snu2'] / pts[i]['Snu2']) / math.log(pts[i + 1]['n'] / pts[i]['n']) for i in range(len(pts) - 1)]
        out.append({'kind': kind, 'points': pts, 'local_exponents_of_Snu2_in_n': ex})
        print('  scaling', kind, ['%.2f' % e for e in ex], flush=True)
    RES['sections']['5.7_scaling'] = out
    print('5.7 done', time.time() - T0, flush=True)


def section_random_sweep():
    """Theorem E and Lemmas A-C,G,G' on many random configurations at p <= 300."""
    rng = random.Random(99)
    out = []
    for p in [101, 197, 293]:
        chi = legendre_table(p); inv = inv_table(p)
        for k in range(12):
            nA = rng.randint(3, 25); nB = rng.randint(3, 25); m = rng.randint(1, 4); H = rng.randint(1, 5); r = rng.randint(1, 3)
            kind = rng.choice(['rand', 'int', 'ap'])
            A = rand_set(p, nA, rng)
            if kind == 'rand':
                B = rand_set(p, nB, rng)
            elif kind == 'int':
                s = rng.randrange(p); B = [(s + i) % p for i in range(nB)]
            else:
                s = rng.randrange(p); q = rng.randrange(1, p); B = [(s + q * i) % p for i in range(nB)]
            Z = [z + 1 for z in rand_set(p - 1, m, rng)]
            out.append(theoremE_instance(p, chi, inv, A, B, Z, H, r, f'sweep p={p} k={k} {kind}'))
    RES['sections']['5.8_random_sweep'] = {'n_instances': len(out), 'min_RHS_over_trivial': min(o['RHS_E/trivial'] for o in out),
                                         'max_|S|_over_RHS': max(o['|S(A,B)|'] / o['RHS_E'] for o in out)}
    print('5.8 done', time.time() - T0, flush=True)


def section_exhaustive_lambda():
    """Exact minimum of lambda over ALL singleton shift sets Z={z} with I={1,2} for a random B: witness for Prop F."""
    rng = random.Random(5)
    out = []
    for p, alpha in [(1009, 0.40), (2003, 0.45)]:
        n = int(round(p ** alpha)); B = rand_set(p, n, rng)
        EB = add_energy(p, B)
        best = None
        for z in range(1, p):
            Y = {(b - z * t) % p for b in B for t in (1, 2)}
            lam = len(Y) / n
            if best is None or lam < best[0]:
                best = (lam, z, len(Y))
        # Prop F lower bound with H'={z,2z}: |H'|=2, E+(H')=... computed exactly
        Hp = [best[1] % p, 2 * best[1] % p]
        lower = n ** 2 * 4 / math.sqrt(EB * add_energy(p, Hp)) / n
        check(f'exhaustive-lambda p={p}', best[0] >= lower - 1e-12, {'min_lambda': best[0], 'lower': lower})
        rec = {'p': p, 'alpha': alpha, '|B|': n, 'B': B, 'E+(B)': EB, 'min_lambda_over_all_z_with_I={1,2}': best[0],
               'argmin_z': best[1], 'PropF_lower_bound_on_lambda': lower,
               'RHS(E)/trivial_is_at_least': best[0] * 2 ** -0.25}
        out.append(rec)
        RES['witnesses'].append({'hypothesis': 'some singleton shift set Z={z}, I={1,2} makes B nearly shift-invariant (lambda <= 1+o(1))',
                                 'status': 'REFUTED (exhaustive over all z)', **{k: v for k, v in rec.items() if k != 'B'}, 'B': B})
        print('  exhaustive', p, n, 'min lambda = %.4f' % best[0], flush=True)
    RES['sections']['5.9_exhaustive_lambda'] = out
    print('5.9 done', time.time() - T0, flush=True)


# ================================================================ 5.6 multiplicative energy of translated intervals
def section_translated_interval_energy():
    rng = random.Random(3)
    out = []
    for p in [1009, 2003]:
        N = int(round(math.sqrt(p) / 2))  # below Chang's (sqrt p -1)/2 threshold
        worst = 0; worst_x = None
        xs = rng.sample(range(p), 60) + [0, p - N // 2]
        for x in xs:
            E = mult_energy(p, [(x + t) % p for t in range(1, N + 1)])
            if E > worst:
                worst, worst_x = E, x
        rec = {'p': p, 'N': N, 'max_E_x_over_sampled_x': worst, 'argmax_x': worst_x,
               'E/(N^2 log p)': worst / (N ** 2 * math.log(p)), 'E/N^{11/4}': worst / N ** 2.75,
               'E/(N^2+N^4/p)': worst / (N ** 2 + N ** 4 / p)}
        out.append(rec)
        # Chang Prop 1 (n=1) with C = 2^{9/4}: E < C log p N^{11/4}
        check(f'ChangProp1 p={p}', worst < 2 ** 2.25 * math.log(p) * N ** 2.75, rec)
        # Alsetri-Shao form with a generous constant 8: E <= 8 (N^2 + N^4/p) log p (HEURISTIC constant)
        rec['AS_form_holds_with_C=8'] = bool(worst <= 8 * (N ** 2 + N ** 4 / p) * math.log(p))
    RES['sections']['5.6_translated_interval_energy'] = out
    print('5.6 done', time.time() - T0)


def main():
    section_lemmaD()
    section_theoremE()
    section_exponents()
    section_propF()
    section_measurements()
    section_scaling()
    section_random_sweep()
    section_exhaustive_lambda()
    section_translated_interval_energy()
    RES['status'] = 'ALL CHECKS PASSED' if not RES['failures'] else f"{len(RES['failures'])} FAILURES"
    RES['elapsed_seconds'] = time.time() - T0
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, 'w') as fh:
        json.dump(RES, fh, indent=1, default=lambda o: int(o) if isinstance(o, (np.integer,)) else float(o))
    print(RES['status'], 'checks =', RES['checks'], 'elapsed = %.1fs' % RES['elapsed_seconds'])


if __name__ == '__main__':
    main()
