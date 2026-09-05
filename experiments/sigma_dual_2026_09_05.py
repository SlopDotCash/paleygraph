#!/usr/bin/env python3
"""
sigma / dual (2026-09-05): verifier for research/sigma-dual-2026-09-05.md.

Direction: the multiplicative-Fourier dual of the shifted-subgroup sum
    T_H(c) = sum_{h in H} chi(h + c),   H <= F_p^*,  c != 0,
its expansion over cosets of H,
    G * T_H(c) = sum_{cosets xi H} chi(xi) Hhat(xi) Hchi(c xi),
    Hhat(xi) = sum_h e(xi h/p),   Hchi(eta) = sum_h chi(h) e(eta h/p),   G = sum_x chi(x) e(x/p),
the l2 obstruction, the coset-cancellation ratio
    rho(H,c) = sqrt(p) |T_H(c)| / D(H,c),   D(H,c)^2 = sum_{cosets} |Hhat(xi)|^2 |Hchi(c xi)|^2,
its exact expression through the signed additive energy R_H(c) (Theorem 3.1 of the note),
and the Gauss-period / cyclotomic-number / Jacobi-sum reformulations (Section 4 of the note).

Section 1 is exact integer arithmetic in Z[zeta_p] (normal form: coefficient of zeta^0 removed
using 1 + zeta + ... + zeta^{p-1} = 0), cross-validated against sympy's polynomial remainder
modulo the cyclotomic polynomial.  Section 3 uses double-precision FFTs for the periods and
verifies the resulting D(H,c)^2 against the exact rational formula at the witness shifts.

Writes results/sigma_dual_2026_09_05.json.  Standard library + numpy + sympy only.
"""
import json, os, sys, time, math
from math import isqrt
from fractions import Fraction
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "sigma_dual_2026_09_05.json")
T0 = time.time()

P_EXACT = 200      # exact Z[zeta_p] checks: all primes <= P_EXACT, all subgroups, all c
P_PRODUCT = 100    # exact period product relations
P_SCAN = 3000      # rho scan: all primes <= P_SCAN, all subgroups, all cosets of c

RES = {"checks": 0, "failures": [], "witnesses": {}, "sections": {}, "params":
       {"P_EXACT": P_EXACT, "P_PRODUCT": P_PRODUCT, "P_SCAN": P_SCAN}}

def check(cond, name, info=None):
    RES["checks"] += 1
    if not cond:
        RES["failures"].append({"name": name, "info": info})
        print("FAIL", name, info)

# ------------------------------------------------------------------ helpers
def primes_in(lo, hi):
    s = bytearray([1]) * (hi + 1); s[0] = s[1] = 0
    for i in range(2, isqrt(hi) + 1):
        if s[i]: s[i*i::i] = bytearray(len(s[i*i::i]))
    return [i for i in range(lo, hi + 1) if s[i]]

def factor(n):
    f = {}; d = 2
    while d * d <= n:
        while n % d == 0: f[d] = f.get(d, 0) + 1; n //= d
        d += 1
    if n > 1: f[n] = f.get(n, 0) + 1
    return f

def divisors(n):
    ds = [1]
    for q, e in factor(n).items():
        ds = [d * q**i for d in ds for i in range(e + 1)]
    return sorted(ds)

def primitive_root(p):
    fs = list(factor(p - 1))
    for g in range(2, p):
        if all(pow(g, (p - 1) // q, p) != 1 for q in fs): return g
    raise ValueError

def legendre(p):
    chi = np.zeros(p, dtype=np.int64)
    chi[1:] = -1
    for x in range(1, p): chi[x * x % p] = 1
    return chi

def dlog_table(p, g):
    logs = np.zeros(p, dtype=np.int64); x = 1
    for k in range(p - 1):
        logs[x] = k; x = x * g % p
    return logs

def subgroup(p, g, d):
    """H = <g^d>, order f = (p-1)/d, as an int64 array."""
    f = (p - 1) // d
    return np.array([pow(g, d * k, p) for k in range(f)], dtype=np.int64)

# ------------------------------------------------------------------ exact Z[zeta_p]
def normal(v):
    """Normal form in Z[zeta_p] w.r.t. the basis zeta^1..zeta^{p-1}: subtract v[0]*(1,...,1)."""
    v = np.asarray(v, dtype=np.int64)
    return v - v[0]

def sympy_crosscheck():
    """Validate the normal form + cyclic-convolution product against sympy at p in {7,11,13}."""
    import sympy
    x = sympy.symbols("x")
    rng = np.random.default_rng(20260905)
    sec = {"cases": 0}
    for p in (7, 11, 13):
        Phi = sympy.Poly(sympy.cyclotomic_poly(p, x), x)
        for _ in range(4):
            a = rng.integers(-5, 6, size=p); b = rng.integers(-5, 6, size=p)
            # cyclic convolution = product in Z[zeta_p]
            prod = np.zeros(p, dtype=np.int64)
            for i in range(p):
                prod[(i + np.arange(p)) % p] += a[i] * b
            mine = normal(prod)
            # sympy: reduce the product modulo Phi_p; basis 1..x^{p-2}
            Pa = sympy.Poly(list(reversed([int(t) for t in a])), x)
            Pb = sympy.Poly(list(reversed([int(t) for t in b])), x)
            rem = (Pa * Pb).rem(Phi)
            coeffs = [0] * (p - 1)
            for (k,), cf in rem.terms():
                coeffs[k] = int(cf)
            # convert mine (basis zeta^1..zeta^{p-1}) to basis 1..zeta^{p-2}: zeta^{p-1} = -sum_{k<p-1} zeta^k
            conv = [int(mine[k]) - int(mine[p - 1]) for k in range(p - 1)]
            check(conv == coeffs, "sympy_crosscheck", {"p": p})
            sec["cases"] += 1
    RES["sections"]["B_sympy_crosscheck"] = sec

# ------------------------------------------------------------------ Section 1: exact identities
def section1():
    sec = {"primes": 0, "subgroups": 0, "identity_checks": 0, "c_in_minus_H_checks": 0, "c_zero_checks": 0,
           "l2_checks": 0, "product_relation_checks": 0, "cyclotomic_alt_checks": 0,
           "inversion_symmetry_checks": 0}
    for p in primes_in(3, P_EXACT):
        chi = legendre(p); g = primitive_root(p); logs = dlog_table(p, g)
        Gvec = chi.copy()                       # G = sum_x chi(x) zeta^x
        sec["primes"] += 1
        for d in divisors(p - 1):
            f = (p - 1) // d; H = subgroup(p, g, d); chiH = chi[H]
            iota = 1 if d % 2 == 0 else 0       # H subset of QR  <=>  d even
            reps = np.array([pow(g, j, p) for j in range(d)], dtype=np.int64)
            chireps = chi[reps]
            sec["subgroups"] += 1
            # T_H(c) for all c, directly
            T = np.array([chi[(H + c) % p].sum() for c in range(p)], dtype=np.int64)
            xi_h = (reps[:, None] * H[None, :]) % p                 # d x f  : exponents of Hhat(xi)
            w_base = (chireps[:, None, None] * chiH[None, None, :]) * np.ones((1, f, 1), dtype=np.int64)
            for c in range(p):
                idx = (xi_h[:, :, None] + ((c * xi_h) % p)[:, None, :]) % p
                rhs = np.bincount(idx.ravel(), weights=w_base.ravel().astype(float), minlength=p)
                rhs = np.rint(rhs).astype(np.int64)
                ok = np.array_equal(normal(T[c] * Gvec), normal(rhs))
                check(ok, "coset_identity", {"p": p, "d": d, "c": c})
                sec["identity_checks"] += 1
                if c != 0 and (-c) % p in set(H.tolist()): sec["c_in_minus_H_checks"] += 1
                if c == 0: sec["c_zero_checks"] += 1
            # exact l2 identities over cosets:  sum_j Hhat(xi_j) conj(Hhat(xi_j)) = p - f
            acc = np.zeros(p, dtype=np.int64)
            diff = (H[:, None] - H[None, :]) % p                     # h - h'
            for j in range(d):
                acc += np.bincount(((reps[j] * diff) % p).ravel(), minlength=p)
            check(np.array_equal(normal(acc), normal((p - f) * np.eye(1, p, 0, dtype=np.int64)[0])),
                  "l2_periods_exact", {"p": p, "d": d})
            # twisted: sum_j Hchi(c xi_j) conj(Hchi(c xi_j)) = p - iota f   (take c = 1; multiset is c-independent)
            acc = np.zeros(p, dtype=np.int64)
            wchi = (chiH[:, None] * chiH[None, :]).ravel().astype(float)
            for j in range(d):
                acc += np.rint(np.bincount((reps[j] * diff).ravel() % p, weights=wchi, minlength=p)).astype(np.int64)
            check(np.array_equal(normal(acc), normal((p - iota * f) * np.eye(1, p, 0, dtype=np.int64)[0])),
                  "l2_twisted_exact", {"p": p, "d": d})
            sec["l2_checks"] += 2
            # cyclotomic alternating identity.  e = d (d even) or 2d (d odd); K = <g^e> = H^+.
            e = d if d % 2 == 0 else 2 * d
            K = subgroup(p, g, e)
            for a in range(e):
                ga = pow(g, a, p)
                # (a,k)_e = #{h in K : 1 + g^a h in g^k K}; alternating sum over k = sum_{h in K} chi(1 + g^a h)
                def alt(aa):
                    x = (1 + pow(g, aa, p) * K) % p; x = x[x != 0]
                    return int(((-1) ** (logs[x] % e)).sum())
                if d % 2 == 0:
                    val = alt(a)
                else:
                    val = int(chi[ga]) * (alt((-a) % e) + alt((d - a) % e))
                check(val == T[ga], "cyclotomic_alternating", {"p": p, "d": d, "a": a, "val": val, "T": int(T[ga])})
                sec["cyclotomic_alt_checks"] += 1
            # inversion symmetry for H subset QR: T_H(c) = chi(c) T_H(1/c)
            if iota:
                for c in range(1, p):
                    check(T[c] == chi[c] * T[pow(c, p - 2, p)], "inversion_symmetry", {"p": p, "d": d, "c": c})
                    sec["inversion_symmetry_checks"] += 1
            # exact product relation of Gauss periods (p <= P_PRODUCT):
            # eta_i eta_{i+a} = f [ -1 in g^a H ] + sum_k (a,k)_d eta_{i+k},   (a,k)_d = #{h in H: 1 + g^a h in g^k H}
            if p <= P_PRODUCT:
                Hset = set(H.tolist())
                for a in range(d):
                    ga = pow(g, a, p)
                    x = (1 + ga * H) % p
                    delta = int((x == 0).sum())          # 1 iff -1 in g^a H
                    xk = logs[x[x != 0]] % d
                    cyc = np.bincount(xk, minlength=d)   # (a,k)_d, k in Z/d
                    for i in range(d):
                        lhs = np.bincount(((reps[i] * H)[:, None] + (reps[(i + a) % d] * H)[None, :]).ravel() % p,
                                          minlength=p)
                        rhs = np.zeros(p, dtype=np.int64); rhs[0] += f * delta
                        for k in range(d):
                            if cyc[k]:
                                rhs += cyc[k] * np.bincount((reps[(i + k) % d] * H) % p, minlength=p)
                        check(np.array_equal(normal(lhs), normal(rhs)), "period_product_relation",
                              {"p": p, "d": d, "i": i, "a": a})
                        sec["product_relation_checks"] += 1
    RES["sections"]["1_exact_identities"] = sec
    print("section 1 done", sec, round(time.time() - T0, 1), "s"); sys.stdout.flush()

# ------------------------------------------------------------------ Section 2: Jacobi / Gauss-sum dual (float, 1e-8)
def section2():
    sec = {"jacobi_expansion_checks": 0, "gauss_product_checks": 0, "jacobi_modulus_checks": 0,
           "jacobi_l2_checks": 0, "reflection_checks": 0}
    for p in primes_in(3, P_EXACT):
        chi = legendre(p); g = primitive_root(p); logs = dlog_table(p, g)
        e_p = np.exp(2j * np.pi * np.arange(p) / p)
        Gsum = (chi * e_p).sum()
        check(abs(abs(Gsum) - math.sqrt(p)) < 1e-9, "gauss_modulus", {"p": p})
        xs = np.arange(1, p)
        for d in divisors(p - 1):
            f = (p - 1) // d
            H = subgroup(p, g, d)
            omega = np.exp(2j * np.pi / d)
            psi = omega ** (logs[xs] % d)                           # psi(x), x = 1..p-1, psi(g) = omega
            T = lambda c: int(chi[(H + c) % p].sum())
            # Gauss sums G(psi^i) (G(psi^0) = -1) and G(psi^i chi)
            Gpsi = np.array([(psi ** i * e_p[xs]).sum() for i in range(d)])
            Gpsichi = np.array([(psi ** i * chi[xs] * e_p[xs]).sum() for i in range(d)])
            # Jacobi sums J(psi^i, chi) = sum_{y != 0,1} psi^i(y) chi(1-y)
            ys = np.arange(2, p)
            psiy = omega ** (logs[ys] % d)
            J = np.array([(psiy ** i * chi[(1 - ys) % p]).sum() for i in range(d)])
            for i in range(d):
                if i == 0: ok = abs(J[i] + 1) < 1e-9
                elif d % 2 == 0 and i == d // 2: ok = abs(J[i] + chi[p - 1]) < 1e-9
                else: ok = abs(abs(J[i]) - math.sqrt(p)) < 1e-9
                check(ok, "jacobi_modulus", {"p": p, "d": d, "i": i}); sec["jacobi_modulus_checks"] += 1
            l2 = (np.abs(J) ** 2).sum()
            target = (d - 2) * p + 2 if d % 2 == 0 else (d - 1) * p + 1
            check(abs(l2 - target) < 1e-6 * p * d, "jacobi_l2", {"p": p, "d": d}); sec["jacobi_l2_checks"] += 1
            # reflection J(psi^i, psi^{-i} chi) = psi^i(-1) J(psi^i, chi)  (used to match the two duals)
            # reflection J(psi^i, psi^{-i} chi) = psi^i(-1) J(psi^i, chi), with J(a,b) = sum_{y != 0,1} a(y) b(1-y)
            psi1y = omega ** (logs[(1 - ys) % p] % d)                  # psi(1-y)
            for i in range(d):
                Jr = (psiy ** i * psi1y ** (-i) * chi[(1 - ys) % p]).sum()
                sgn = omega ** ((i * logs[p - 1]) % d)                  # psi^i(-1)
                check(abs(Jr - sgn * J[i]) < 1e-8, "jacobi_reflection", {"p": p, "d": d, "i": i})
                sec["reflection_checks"] += 1
            cs = sorted({1, 2, 3, p - 1, p - 2, (p + 1) // 2} & set(range(1, p)))
            for c in cs:
                lc = logs[c]; lmc = logs[(-c) % p]
                # (i) Jacobi expansion  T_H(c) = (chi(c)/d) sum_i psi^i(-c) J(psi^i, chi)
                val = chi[c] / d * (omega ** ((lmc * np.arange(d)) % d) * J).sum()
                check(abs(val - T(c)) < 1e-8, "jacobi_expansion", {"p": p, "d": d, "c": c, "val": complex(val)})
                sec["jacobi_expansion_checks"] += 1
                # (ii) coset sum = (chi(c)/d) sum_i psi^i(c) G(psi^i) G(psi^{-i} chi)  (Fourier dual on Z/d)
                coset_sum = 0j
                for j in range(d):
                    xi = pow(g, j, p)
                    Hh = e_p[(xi * H) % p].sum(); Hc = (chi[H] * e_p[(c * xi * H) % p]).sum()
                    coset_sum += chi[xi] * Hh * Hc
                dual = chi[c] / d * (omega ** ((lc * np.arange(d)) % d) * Gpsi * Gpsichi[(-np.arange(d)) % d]).sum()
                check(abs(coset_sum - dual) < 1e-7 * p, "gauss_product_dual", {"p": p, "d": d, "c": c})
                check(abs(coset_sum - Gsum * T(c)) < 1e-7 * p, "coset_identity_float", {"p": p, "d": d, "c": c})
                sec["gauss_product_checks"] += 2
    RES["sections"]["2_jacobi_dual"] = sec
    print("section 2 done", sec, round(time.time() - T0, 1), "s"); sys.stdout.flush()

# ------------------------------------------------------------------ Section 3: rho scan
def scan_prime(p, rows, sec):
    chi = legendre(p); g = primitive_root(p)
    for d in divisors(p - 1):
        f = (p - 1) // d; H = subgroup(p, g, d); iota = 1 if d % 2 == 0 else 0
        ind = np.zeros(p); ind[H] = 1.0
        Hhat = np.fft.ifft(ind) * p                         # Hhat(xi) = sum_h e(+xi h/p)
        indchi = np.zeros(p); indchi[H] = chi[H]
        Hchi = np.fft.ifft(indchi) * p
        negH = np.zeros(p); negH[(-H) % p] = 1.0
        T = np.rint(np.fft.ifft(np.fft.fft(chi.astype(float)) * np.fft.fft(negH)).real).astype(np.int64)
        reps = np.array([pow(g, j, p) for j in range(d)], dtype=np.int64)
        for c in sorted({1, 2, p - 1} & set(range(1, p))):
            check(T[c] == chi[(H + c) % p].sum(), "T_fft_direct", {"p": p, "d": d, "c": c})
        u = np.abs(Hhat[reps]) ** 2; v = np.abs(Hchi[reps]) ** 2
        # coset invariance of both moduli families
        check(np.allclose(np.abs(Hhat[(reps[:, None] * H[None, :]) % p]) ** 2, u[:, None], atol=1e-6 * p),
              "coset_invariance_Hhat", {"p": p, "d": d})
        check(np.allclose(np.abs(Hchi[(reps[:, None] * H[None, :]) % p]) ** 2, v[:, None], atol=1e-6 * p),
              "coset_invariance_Hchi", {"p": p, "d": d})
        if iota:
            check(np.allclose(np.sort(u), np.sort(v), atol=1e-6 * p), "same_multiset_when_H_in_QR", {"p": p, "d": d})
        check(abs(u.sum() - (p - f)) < 1e-6 * p, "l2_periods_float", {"p": p, "d": d})
        check(abs(v.sum() - (p - iota * f)) < 1e-6 * p, "l2_twisted_float", {"p": p, "d": d})
        # D^2(a) = sum_j u_j v_{j+a}, all a in Z/d, by cyclic cross-correlation
        D2 = np.fft.ifft(np.conj(np.fft.fft(u)) * np.fft.fft(v)).real
        for a in (0, 1 % d):
            check(abs(D2[a] - (u * np.roll(v, -a)).sum()) < 1e-6 * p * f, "D2_fft_direct", {"p": p, "d": d, "a": a})
        check(abs(D2.sum() - (p - f) * (p * f - iota * f * f) / f) < 1e-6 * p * p, "D2_average_identity", {"p": p, "d": d})
        Tabs = np.abs(T[reps]).astype(float)
        check(int(round((Tabs ** 2).sum() * f)) == p * f - f * f - iota * f * f, "T_second_moment", {"p": p, "d": d})
        rho = np.sqrt(p) * Tabs / np.sqrt(D2)
        a_star = int(np.argmax(Tabs)); a_rho = int(np.argmax(rho))
        # exact signed energy R_H(c) = sum_s r(-cs) rchi(s)
        diff = (H[:, None] - H[None, :]).ravel() % p
        r = np.bincount(diff, minlength=p)
        rchi = np.rint(np.bincount(diff, weights=(chi[H][:, None] * chi[H][None, :]).ravel().astype(float),
                                   minlength=p)).astype(np.int64)
        s = np.arange(p)
        def R_exact(c): return int((r[(-c * s) % p] * rchi[s]).sum())
        def E_cH(c): return int((r[(-c * s) % p] * r[s]).sum())     # E(H, cH) = #{h1 + c h3 = h2 + c h4}
        E_H = int((r * r).sum())                                        # additive energy E(H)
        rec = {}
        for a in sorted({a_star, a_rho, 0}):
            c = int(reps[a]); Rc = R_exact(c)
            D2_exact = Fraction(p * Rc, f) - iota * f ** 3
            check(abs(D2[a] - float(D2_exact)) < 1e-6 * max(1.0, float(D2_exact)), "D2_energy_formula",
                  {"p": p, "d": d, "a": a, "float": float(D2[a]), "exact": str(D2_exact)})
            check(abs(Rc) <= E_cH(c), "signed_energy_bound", {"p": p, "d": d, "a": a})
            check(E_cH(c) <= E_H, "E_cH_le_E_H", {"p": p, "d": d, "a": a})
            rho2_exact = Fraction(p * int(T[c]) ** 2) / D2_exact if D2_exact != 0 else None
            rec[a] = {"c": c, "T": int(T[c]), "R": Rc, "E_cH": E_cH(c), "D2_exact": str(D2_exact),
                      "rho2_exact": str(rho2_exact), "rho": float(rho[a]), "eps": Rc / f ** 2}
        # abstract norm obstruction Phi = (1/sqrt p) min(alpha^2 B/A, beta^2 A/B)
        A = math.sqrt(u.max()); B = math.sqrt(v.max()); alpha2 = u.sum(); beta2 = v.sum()
        Phi = min(alpha2 * B / A, beta2 * A / B) / math.sqrt(p)
        check(A >= math.sqrt(f * (p - f) / (p - 1)) - 1e-9, "pigeonhole_A", {"p": p, "d": d})
        check(B >= math.sqrt((p * f - iota * f * f) / (p - 1)) - 1e-9, "pigeonhole_B", {"p": p, "d": d})
        if iota: check(abs(Phi - (p - f) / math.sqrt(p)) < 1e-6 * p, "Phi_exact_when_H_in_QR", {"p": p, "d": d})
        # sandwich (Cor. 3.2): rho >= |T| sqrt(f / E(H,cH)) for all a (checked at the three witnesses),
        # and, for H in QR and f^2 < p, rho <= (|T|/sqrt f) (1 - f^2/p)^{-1/2}
        for a, w in rec.items():
            if w["E_cH"] > 0:
                check(w["rho"] >= abs(w["T"]) * math.sqrt(f / w["E_cH"]) - 1e-9, "sandwich_lower", {"p": p, "d": d, "a": a})
        if iota and f * f < p:
            check(np.all(rho <= Tabs / math.sqrt(f) / math.sqrt(1 - f * f / p) + 1e-9), "sandwich_upper", {"p": p, "d": d})
        M = int(Tabs.max())
        rows.append({"p": p, "f": f, "d": d, "iota": iota, "M": M, "c_star": rec[a_star]["c"], "T_star": rec[a_star]["T"],
                     "rho_star": rec[a_star]["rho"], "eps_star": rec[a_star]["eps"], "rho2_star_exact": rec[a_star]["rho2_exact"],
                     "rho_max": float(rho.max()), "c_rho": rec[a_rho]["c"], "T_rho": rec[a_rho]["T"], "eps_rho": rec[a_rho]["eps"],
                     "rho2_max_exact": rec[a_rho]["rho2_exact"], "R_star": rec[a_star]["R"], "E_cH_star": rec[a_star]["E_cH"],
                     "E_H": E_H, "rho_median": float(np.median(rho)), "rho_rms": float(math.sqrt((rho ** 2).mean())),
                     "D2max_over_pf": float(D2.max() / (p * f)), "D2min_over_pf": float(D2.min() / (p * f)),
                     "Phi_over_sqrtp": Phi / math.sqrt(p), "A": A, "B": B,
                     "Mnorm": M / math.sqrt(f), "ev": math.sqrt(2 * math.log(max(d, 2)))})
        sec["subgroups"] += 1

def section3():
    sec = {"primes": 0, "subgroups": 0}
    rows = []
    for p in primes_in(3, P_SCAN):
        scan_prime(p, rows, sec); sec["primes"] += 1
    RES["sections"]["3_rho_scan"] = sec
    print("section 3 done", sec, round(time.time() - T0, 1), "s"); sys.stdout.flush()
    return rows

# ------------------------------------------------------------------ Section 4: aggregate + witnesses
def section4(rows):
    W = RES["witnesses"]
    proper = [r for r in rows if r["d"] > 1 and r["f"] > 1]
    # largest rho overall, with exact rational rho^2
    top = sorted(proper, key=lambda r: -r["rho_max"])[:15]
    W["largest_rho"] = [{k: r[k] for k in ("p", "f", "d", "iota", "M", "c_rho", "T_rho", "rho_max", "rho2_max_exact",
                                            "eps_rho", "rho_star", "Mnorm", "ev")} for r in top]
    # ratio rho_star / (M/sqrt f) in the range 2 < |H| < sqrt p
    rat = [(r["rho_star"] / r["Mnorm"], r) for r in proper if 2 < r["f"] < math.sqrt(r["p"])]
    W["ratio_rho_over_Mnorm_range_small_H"] = {
        "n": len(rat), "min": min(rat)[0], "max": max(rat)[0], "mean": sum(x for x, _ in rat) / len(rat),
        "argmin": {k: min(rat)[1][k] for k in ("p", "f", "d", "M", "c_star", "eps_star", "rho_star")},
        "argmax": {k: max(rat)[1][k] for k in ("p", "f", "d", "M", "c_star", "eps_star", "rho_star")}}
    # bands by |H|
    bands = {}
    for r in proper:
        b = min(int(math.log2(r["f"])), 11)
        bands.setdefault(b, []).append(r)
    W["bands_by_log2_H"] = []
    for b in sorted(bands):
        rs = bands[b]
        W["bands_by_log2_H"].append({
            "band": "2^%d<=|H|<2^%d" % (b, b + 1), "n": len(rs),
            "max_rho": max(r["rho_max"] for r in rs),
            "mean_rho_max": sum(r["rho_max"] for r in rs) / len(rs),
            "max_rho_over_ev": max(r["rho_max"] / r["ev"] for r in rs),
            "mean_rho_over_ev": sum(r["rho_max"] / r["ev"] for r in rs) / len(rs),
            "max_Mnorm": max(r["Mnorm"] for r in rs),
            "mean_rho_median": sum(r["rho_median"] for r in rs) / len(rs),
            "mean_rho_rms": sum(r["rho_rms"] for r in rs) / len(rs),
            "max_D2max_over_pf": max(r["D2max_over_pf"] for r in rs),
            "max_eps_star": max(r["eps_star"] for r in rs)})
    # monochromatic shifts (M = |H|): rho there, and whether R = |H|^2
    mono = [r for r in proper if r["M"] == r["f"]]
    W["monochromatic"] = {"n": len(mono), "max_H": max(r["f"] for r in mono),
                          "list_H_ge_8": [{k: r[k] for k in ("p", "f", "d", "c_star", "rho_star", "eps_star", "R_star", "rho2_star_exact")}
                                          for r in sorted(mono, key=lambda r: (-r["f"], r["p"])) if r["f"] >= 8]}
    # Theorem 3.3 witnesses: for each m, smallest p in the scan with a monochromatic shift having R = m^2
    thm = {}
    for r in sorted(mono, key=lambda r: r["p"]):
        m = r["f"]
        if r["R_star"] == m * m and m not in thm:
            thm[m] = {"p": r["p"], "c": r["c_star"], "rho": r["rho_star"], "sqrt_m": math.sqrt(m), "rho2_exact": r["rho2_star_exact"]}
            check(r["rho_star"] >= math.sqrt(m) - 1e-9, "thm33_rho_ge_sqrt_m", {"p": r["p"], "m": m})
    W["theorem33_witnesses"] = thm
    # abstract norm obstruction
    W["Phi_min"] = min(({"Phi_over_sqrtp": r["Phi_over_sqrtp"], "p": r["p"], "f": r["f"], "d": r["d"], "A": r["A"], "B": r["B"]}
                        for r in proper), key=lambda w: w["Phi_over_sqrtp"])
    check(all(r["Phi_over_sqrtp"] >= 0.4 for r in proper), "Phi_ge_0.4_sqrtp")
    # M_H < sqrt p and rho_max <= sqrt(p/(f - f^3/p)) sanity for H in QR
    check(all(r["M"] < math.sqrt(r["p"]) for r in proper), "M_lt_sqrt_p")
    # distribution summaries of typical rho
    W["typical_rho"] = {"median_of_rho_median": float(np.median([r["rho_median"] for r in proper if r["f"] >= 8])),
                        "median_of_rho_rms": float(np.median([r["rho_rms"] for r in proper if r["f"] >= 8]))}
    RES["rows"] = rows

# ------------------------------------------------------------------ main
if __name__ == "__main__":
    sympy_crosscheck()
    section1()
    section2()
    rows = section3()
    section4(rows)
    RES["elapsed_s"] = round(time.time() - T0, 1)
    with open(OUT, "w") as fh:
        json.dump(RES, fh, indent=1, default=str)
    print("checks", RES["checks"], "failures", len(RES["failures"]), "elapsed", RES["elapsed_s"], "s")
    for w in RES["witnesses"]["largest_rho"][:6]:
        print("rho_max %.4f  p=%d |H|=%d d=%d c=%d T=%d eps=%.3f  M/sqrt|H|=%.3f  sqrt(2 log d)=%.3f" %
              (w["rho_max"], w["p"], w["f"], w["d"], w["c_rho"], w["T_rho"], w["eps_rho"], w["Mnorm"], w["ev"]))
    print("ratio rho*/(M/sqrt f), 2<|H|<sqrt p:", RES["witnesses"]["ratio_rho_over_Mnorm_range_small_H"])
    for b in RES["witnesses"]["bands_by_log2_H"]:
        print(b)
    print("Phi_min:", RES["witnesses"]["Phi_min"])
    print("thm 3.3 witnesses:", RES["witnesses"]["theorem33_witnesses"])
    print("typical:", RES["witnesses"]["typical_rho"])
    sys.exit(1 if RES["failures"] else 0)
