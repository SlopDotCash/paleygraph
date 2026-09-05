#!/usr/bin/env python3
"""
sigma / structured  (2026-09-05): exact verifier for the note
research/sigma-structured-2026-09-05.md.

Everything here is exact integer arithmetic (Python ints / numpy int64 with
range checks). Writes results/sigma_structured_2026_09_05.json with the
number of checks performed and every witness.

Notation.  chi = Legendre symbol mod p, chi(0)=0.  For B subset F_p and
r in F_p:   T_B(r) = sum_{b in B} chi(r+b).   S(A,B) = sum_{a in A} T_B(a).
M(B) = max_{r != 0} |T_B(r)|.   M_{2k}(B) = sum_{x in F_p} T_B(x)^{2k}.
"""
import json, os, sys, time, random, math
from math import isqrt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
OUT = os.path.join(ROOT, "results", "sigma_structured_2026_09_05.json")
T0 = time.time()
RNG = random.Random(20260905)

RES = {"checks": 0, "failures": [], "witnesses": {}, "sections": {}}

def check(cond, name, info=None):
    RES["checks"] += 1
    if not cond:
        RES["failures"].append({"name": name, "info": info})
        print("FAIL", name, info)

# ---------------------------------------------------------------- helpers
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

def legendre_table(p):
    chi = np.zeros(p, dtype=np.int64)
    for x in range(1, p):
        v = pow(x, (p - 1) // 2, p)
        chi[x] = 1 if v == 1 else -1
    return chi

def subgroup(p, g, d):
    """the unique subgroup of F_p^* of order d (d | p-1), as a sorted list"""
    assert (p - 1) % d == 0
    e = (p - 1) // d
    return sorted({pow(g, e * i, p) for i in range(d)})

def T_of(chi, B):
    """T_B(x) = sum_b chi[(x+b) % p], as int64 array over x in F_p"""
    p = len(chi); T = np.zeros(p, dtype=np.int64)
    for b in B: T += np.roll(chi, -int(b))     # roll(a,-b)[x] = a[(x+b)%n]
    return T

def S_of(T, A):
    return int(sum(int(T[a]) for a in A))

def product_set(p, B, C):
    return {(b * c) % p for b in B for c in C}

def double_fact_odd(k):   # (2k-1)!!
    r = 1
    for i in range(1, 2 * k, 2): r *= i
    return r

def moment(T, k):
    return sum(int(t) ** (2 * k) for t in T.tolist())

# =====================================================================
# Section 1: identities for subgroups and the equivalence of Task 2
# =====================================================================
sec1 = {"primes": [], "identity_checks": 0, "equivalence_witnesses": []}
for p in [101, 197, 229, 257, 1009, 1093]:
    g = primitive_root(p); chi = legendre_table(p)
    sec1["primes"].append(p)
    xs = np.arange(p)
    for d in divisors(p - 1):
        if d == 1 or d > 400: continue
        H = subgroup(p, g, d); T = T_of(chi, H)
        # (I1)  T_H(c h) = chi(h) T_H(c) for all c, h in H
        ok = all(bool(np.array_equal(T[(xs * h) % p], int(chi[h]) * T)) for h in H)
        check(ok, "I1 coset covariance", (p, d)); sec1["identity_checks"] += 1
        # (I2)  |T_H| constant on cosets => exceptional set is H-invariant
        absT = np.abs(T)
        ok2 = all(bool(np.array_equal(absT[(xs * h) % p], absT)) for h in H)
        check(ok2, "I2 |T_H| H-invariant", (p, d))
        # (I3) second moment identity sum_x T^2 = p|H| - |H|^2 and M_H^2 <= p - |H|
        check(moment(T, 1) == p * d - d * d, "I3 second moment", (p, d))
        M = int(absT[1:].max())
        check(M * M <= p - d, "I3' elementary bound M_H^2 <= p-|H|", (p, d, M))
        # (E) equivalence witness: c attaining M_H, A_c = {c h : chi(h)=1}
        c = int(np.argmax(absT[1:])) + 1
        A_c = [(c * h) % p for h in H if chi[h] == 1]
        S = S_of(T, A_c)
        check(abs(S) == len(A_c) * M, "E converse: S(A_c,H) = |A_c| T_H(c)", (p, d, c))
        # forward direction on a union of few cosets and random A
        for trial in range(2):
            ncos = RNG.randint(1, 3)
            reps = RNG.sample(range(1, p), ncos)
            B = sorted({(r * h) % p for r in reps for h in H})
            if len(B) != ncos * d: continue
            TB = T_of(chi, B)
            A = RNG.sample(range(p), RNG.randint(1, p - 1))
            check(abs(S_of(TB, A)) <= len(A) * ncos * M,
                  "Task2 forward |S(A,B)| <= |A| k M_H", (p, d, ncos))
        if len(sec1["equivalence_witnesses"]) < 12:
            sec1["equivalence_witnesses"].append(
                {"p": p, "|H|": d, "c": c, "M_H": M, "|A_c|": len(A_c), "S(A_c,H)": S})
RES["sections"]["1_subgroup_identities_and_equivalence"] = sec1
print("section 1 done", round(time.time() - T0, 1), "s")

# =====================================================================
# Section 2: top-m reformulation, brute force at tiny p
# =====================================================================
from itertools import combinations
sec2 = {"bruteforce": []}
for p in [13, 17, 19]:
    chi = legendre_table(p)
    for trial in range(3):
        B = RNG.sample(range(p), RNG.randint(2, 6)); T = T_of(chi, B)
        vals = sorted(T.tolist(), reverse=True)
        for m in range(1, 4):
            best = max(S_of(T, A) for A in combinations(range(p), m))
            worst = min(S_of(T, A) for A in combinations(range(p), m))
            check(best == sum(vals[:m]), "top-m max", (p, B, m))
            check(worst == sum(vals[-m:]), "bottom-m min", (p, B, m))
            sec2["bruteforce"].append({"p": p, "|B|": len(B), "m": m, "max_S": best,
                                       "top_m_sum": sum(vals[:m])})

# ---- 2b: |BB| = |B|  iff  B is a coset of a subgroup (brute force, all small subsets)
sec2["coset_iff_K1"] = []
for p in [13, 17]:
    g = primitive_root(p)
    cosets = set()
    for d in divisors(p - 1):
        H = subgroup(p, g, d)
        for c in range(1, p): cosets.add(frozenset((c * h) % p for h in H))
    n_sub = 0; n_k1 = 0
    for m in range(2, 5):
        for B in combinations(range(1, p), m):
            n_sub += 1
            K1 = (len(product_set(p, B, B)) == m)
            is_coset = frozenset(B) in cosets
            check(K1 == is_coset, "|BB|=|B| iff coset", (p, B))
            n_k1 += K1
    sec2["coset_iff_K1"].append({"p": p, "subsets_tested": n_sub, "with_|BB|=|B|": n_k1})
RES["sections"]["2_top_m_reformulation"] = sec2
print("section 2 done", round(time.time() - T0, 1), "s")

# =====================================================================
# Section 3: Weil moment bound, Chebyshev exceptional set, amplification
#            inequality, propagation lemma, symmetric-set bound
# =====================================================================
sec3 = {"moment_checks": 0, "amplification_checks": 0, "propagation_checks": 0,
        "symmetric_bound_checks": 0, "max_moment_ratio": {}, "examples": []}
def sqrt_le(lhs, coef, p):
    """exact test of  lhs <= coef*sqrt(p)  for integers lhs, coef>=0"""
    if lhs <= 0: return True
    return lhs * lhs <= coef * coef * p

for p in [101, 257, 401, 1009, 12289]:
    g = primitive_root(p); chi = legendre_table(p)
    sets = []
    for d in divisors(p - 1):
        if 2 <= d <= 64: sets.append(("subgroup", subgroup(p, g, d)))
    for n in [5, 17, 40]:
        if n < p: sets.append(("random", RNG.sample(range(p), n)))
    for tag, B in sets:
        T = T_of(chi, B); n = len(B)
        for k in [1, 2, 3]:
            Mk = moment(T, k)
            main = double_fact_odd(k) * n**k * p
            ok = sqrt_le(Mk - main, (2 * k - 1) * n**(2 * k), p)
            check(ok, "Weil moment bound", (p, tag, n, k)); sec3["moment_checks"] += 1
            key = f"k={k}"
            ratio = Mk / (main + (2 * k - 1) * n**(2 * k) * math.sqrt(p))
            sec3["max_moment_ratio"][key] = max(sec3["max_moment_ratio"].get(key, 0), ratio)
            # Chebyshev: #{x: |T|>theta n} <= M_2k/(theta n)^{2k}; theta = 1/2
            cnt = int((np.abs(T) * 2 > n).sum())
            check(cnt * (n ** (2 * k)) <= Mk * (2 ** (2 * k)), "Chebyshev exceptional set", (p, tag, n, k))
        # propagation lemma |T_B(rs)| >= |T_B(r)| - |B sym sB|, all r, all s
        Bs = set(B); xs = np.arange(p); absT = np.abs(T)
        for s in range(1, p):
            sB = {(s * b) % p for b in B}
            sym = len(Bs ^ sB)
            ok = bool(np.all(absT[(xs * s) % p] >= absT - sym))
            check(ok, "propagation lemma", (p, tag, n, s)); sec3["propagation_checks"] += 1
        # symmetric-set bound: with Sym_t = {s : |B & sB| >= n - t}, sigma=|Sym_t|,
        #   (M - 2t)_+^2 * sigma <= p n - n^2
        M = int(absT[1:].max())
        for t in [0, max(1, n // 8), n // 4]:
            sigma = sum(1 for s in range(1, p) if len(Bs & {(s * b) % p for b in B}) >= n - t)
            lhs = max(M - 2 * t, 0) ** 2 * sigma
            check(lhs <= p * n - n * n, "symmetric-set bound", (p, tag, n, t, M, sigma))
            sec3["symmetric_bound_checks"] += 1
            if tag == "subgroup" and t == 0:
                check(sigma == n, "Sym_0(H) = H", (p, n))
        if len(sec3["examples"]) < 10:
            sec3["examples"].append({"p": p, "type": tag, "|B|": n, "M(B)": M,
                                     "M2": moment(T, 1), "M4": moment(T, 2)})
    # amplification inequality (exact):  with nu(y) = #{(a,s): a = y s},
    #  |S|*|S(A,B)| <= |A| * sum_s |B sym sB| + sum_y nu(y) |T_B(y)|,
    #  (sum_y nu |T|)^{2k} <= (|A||S|)^{2k-2} * E(A,S) * M_{2k}(B)
    if p <= 1009:
        inv = [0] + [pow(s, p - 2, p) for s in range(1, p)]
        for trial in range(3):
            A = RNG.sample(range(p), RNG.randint(3, 30))
            B = RNG.sample(range(p), RNG.randint(3, 30))
            Sset = RNG.sample(range(1, p), RNG.randint(2, 20))
            T = T_of(chi, B); Bs = set(B)
            nu = np.zeros(p, dtype=np.int64)
            for a in A:
                for s in Sset: nu[(a * inv[s]) % p] += 1
            E = int((nu * nu).sum())
            lhs1 = len(Sset) * abs(S_of(T, A))
            rhs1 = len(A) * sum(len(Bs ^ {(s * b) % p for b in B}) for s in Sset) \
                   + int((nu * np.abs(T)).sum())
            check(lhs1 <= rhs1, "amplification step 1", (p, trial))
            for k in [1, 2, 3]:
                lhs = int((nu * np.abs(T)).sum()) ** (2 * k)
                rhs = (len(A) * len(Sset)) ** (2 * k - 2) * E * moment(T, k)
                check(lhs <= rhs, "amplification Hoelder", (p, trial, k))
                sec3["amplification_checks"] += 1
        # E(A,H) = |A| |H|^2 when A is H-invariant (obstruction to the gain)
        for d in divisors(p - 1):
            if 2 <= d <= 40:
                H = subgroup(p, g, d)
                reps = RNG.sample(range(1, p), 3)
                A = sorted({(r * h) % p for r in reps for h in H})
                nu = np.zeros(p, dtype=np.int64)
                for a in A:
                    for s in H: nu[(a * inv[s]) % p] += 1
                check(int((nu * nu).sum()) == len(A) * d * d, "E(A,H)=|A||H|^2 for H-invariant A", (p, d))
RES["sections"]["3_moments_amplification_propagation"] = sec3
print("section 3 done", round(time.time() - T0, 1), "s")

# =====================================================================
# Section 4: geometric progressions reduce to consecutive-power sums
#   S(aP, P) = sum_{|k|<N} (N-|k|) chi(1 + a g^k)  when chi(g) = 1
# =====================================================================
sec4 = {"checks": 0, "examples": []}
for p in [1009, 12289]:
    chi = legendre_table(p); g0 = primitive_root(p); g = (g0 * g0) % p   # a square of order (p-1)/2
    for N in [8, 30, 100]:
        P = [pow(g, i, p) for i in range(1, N + 1)]
        assert len(set(P)) == N
        TP = T_of(chi, P)
        for a in RNG.sample(range(1, p), 4):
            aP = [(a * x) % p for x in P]
            lhs = S_of(TP, aP)
            rhs = sum((N - abs(k)) * int(chi[(1 + a * pow(g, k % (p - 1), p)) % p]) for k in range(-N + 1, N))
            check(lhs == rhs, "geometric progression identity", (p, N, a)); sec4["checks"] += 1
            if len(sec4["examples"]) < 6:
                sec4["examples"].append({"p": p, "N": N, "a": a, "S(aP,P)": lhs, "|BB|/|B|": len(product_set(p, P, P)) / N})
RES["sections"]["4_geometric_progression_identity"] = sec4
print("section 4 done", round(time.time() - T0, 1), "s")

# =====================================================================
# Section 5: scan of all subgroups, all primes in a range:
#   M_H vs |H|, |H|^2-cliques (c+H monochromatic), elementary bound M_H^2<=p-|H|
# =====================================================================
sec5 = {"prime_range": [101, 4000], "n_primes": 0, "n_subgroups": 0,
        "largest_monochromatic_shift_per_prime": [], "gong_sqrt_p_checks": 0,
        "max_ratio_M_over_H_by_size_band": {}, "index2_witnesses": []}
for p in primes_in(101, 4000):
    g = primitive_root(p); chi = legendre_table(p); sec5["n_primes"] += 1
    data = {}
    for d in divisors(p - 1):
        if d == 1: continue
        H = subgroup(p, g, d); T = T_of(chi, H); absT = np.abs(T)
        M = int(absT[1:].max()); c = int(np.argmax(absT[1:])) + 1
        data[d] = (M, c); sec5["n_subgroups"] += 1
        check(M * M < p, "Gong Thm 2: M_H < sqrt p", (p, d, M)); sec5["gong_sqrt_p_checks"] += 1
        check(M * M <= p - d, "elementary M_H^2 <= p - |H|", (p, d, M))
        band = "small(|H|<=log p)" if d <= math.log(p) else ("mid(log p<|H|<=sqrt p)" if d * d <= p else "large(|H|>sqrt p)")
        sec5["max_ratio_M_over_H_by_size_band"][band] = max(sec5["max_ratio_M_over_H_by_size_band"].get(band, 0), M / d)
        for mult in [1, 2, 3, 4, 6]:
            if d > mult * math.log(p):
                key = f"max M_H/|H| over |H|>{mult}log p"
                sec5["max_ratio_M_over_H_by_size_band"][key] = max(sec5["max_ratio_M_over_H_by_size_band"].get(key, 0), M / d)
    mono = [d for d, (M, c) in data.items() if M == d]
    if mono:
        d = max(mono); sec5["largest_monochromatic_shift_per_prime"].append(
            {"p": p, "|H|": d, "c": data[d][1], "log_p": round(math.log(p), 2), "|H|/log p": round(d / math.log(p), 2)})
    # index-2 (or 3) dense-subset witnesses: H' < H, M_{H'} = |H'|, M_H <= |H|/3
    for d, (M, c) in data.items():
        for k in [2, 3]:
            dd = d // k
            if d % k == 0 and dd in data and data[dd][0] == dd and 3 * M <= d and len(sec5["index2_witnesses"]) < 15:
                sec5["index2_witnesses"].append({"p": p, "|H|": d, "M_H": M, "index": k, "|H'|": dd,
                                                 "M_H'": data[dd][0], "c'": data[dd][1]})
RES["sections"]["5_subgroup_scan"] = sec5
print("section 5 done", round(time.time() - T0, 1), "s", "subgroups:", sec5["n_subgroups"])

# =====================================================================
# Section 6: Task 4 coset-decomposition obstruction, numeric witnesses
# =====================================================================
sec6 = {}
# (i) generic B of size ~p^0.4 vs subgroup of size ~p^0.2: occupancy of cosets
for p, d in [(1000003, 6), (100003, 21), (12289, 16)]:
    g = primitive_root(p); H = subgroup(p, g, d)
    N = int(round(p ** 0.4)); B = RNG.sample(range(1, p), N)
    # coset label: discrete-log-free labelling via canonical rep b*h min
    Hs = H
    reps = {}
    for b in B:
        key = min((b * h) % p for h in Hs); reps[key] = reps.get(key, 0) + 1
    occ = sorted(reps.values(), reverse=True)
    sec6[f"generic_p{p}"] = {"p": p, "|H|": d, "|B|": N, "cosets_met": len(reps), "max_occupancy": occ[0],
                             "expected_collisions_N^2|H|/2p": round(N * N * d / (2 * p), 3)}
    check(len(reps) >= N - N * N * d / p - 5, "occupancy sanity", (p, d))
# (ii) exact witness: the decomposition identity S(A,B) = sum_C chi(c_C) S(A/c_C, B_C)
p = 12289; g = primitive_root(p); chi = legendre_table(p); inv = [0] + [pow(s, p - 2, p) for s in range(1, p)]
H = subgroup(p, g, 64); Hs = set(H)
B = RNG.sample(range(1, p), 150); A = RNG.sample(range(p), 200)
TB = T_of(chi, B); lhs = S_of(TB, A)
cos = {}
for b in B:
    key = min((b * h) % p for h in H); cos.setdefault(key, []).append(b)
rhs = 0
for cC, part in cos.items():
    BC = [(b * inv[cC]) % p for b in part]; assert all(x in Hs for x in BC)
    AC = [(a * inv[cC]) % p for a in A]
    rhs += int(chi[cC]) * S_of(T_of(chi, BC), AC)
check(lhs == rhs, "coset decomposition identity", (p,))
sec6["decomposition_identity"] = {"p": p, "|H|": 64, "|B|": 150, "|A|": 200, "S": lhs, "cosets": len(cos)}
# (iii) dense-subset obstruction, exact: subgroup H, index-2 subgroup H' (dense subset of H),
#       A = c' H' cap {chi=+1 shifted}: |S(A,H')| = |A||H'| while M_H small. Use scan witnesses.
sec6["dense_subset_witnesses"] = sec5["index2_witnesses"][:5]
for w in sec6["dense_subset_witnesses"]:
    pp = w["p"]; gg = primitive_root(pp); cc = legendre_table(pp)
    Hp = subgroup(pp, gg, w["|H'|"]); T = T_of(cc, Hp); c = w["c'"]
    A = [(c * h) % pp for h in Hp if cc[h] == 1]
    check(abs(S_of(T, A)) == len(A) * len(Hp), "dense-subset witness exact", (pp, w["|H'|"]))
    w["|A|"] = len(A); w["|S(A,H')|/(|A||H'|)"] = 1.0

# (iv) smoothing S(A,BH) is unrelated to S(A,B): A={a}, B = {b: chi(a+b)=1}
a0 = 3; B0 = [b for b in range(p) if chi[(a0 + b) % p] == 1]
TB0 = T_of(chi, B0); s_ab = S_of(TB0, [a0])
H = subgroup(p, g, 64)
s_abh = sum(int(chi[(a0 + b * h) % p]) for b in B0 for h in H)   # with multiplicity
check(s_ab == len(B0), "clique-type S(A,B)=|B|", (p,))
sec6["smoothing_witness"] = {"p": p, "a": a0, "|B|": len(B0), "|H|": 64, "S(A,B)": s_ab,
                             "S(A,B*H) with multiplicity": s_abh, "|B||H|": len(B0) * 64}
RES["sections"]["6_task4_obstruction"] = sec6
print("section 6 done", round(time.time() - T0, 1), "s")

# =====================================================================
# Section 7: Task 3 test table at p = 12289 (2^12*3+1): multiplicative doubling K,
#   single-shift maximum M(B)/|B|, top-|B| average, exceptional set, symmetry set
# =====================================================================
p = 12289; g = primitive_root(p); chi = legendre_table(p)
inv = [0] + [pow(s, p - 2, p) for s in range(1, p)]
def profile(tag, B, t_frac=4):
    B = sorted(set(B)); n = len(B); Bs = set(B)
    T = T_of(chi, B); absT = np.abs(T); M = int(absT[1:].max())
    vals = sorted(absT.tolist(), reverse=True)
    K = len(product_set(p, B, B)) / n
    t = n // t_frac
    sigma = sum(1 for s in range(1, p) if len(Bs & {(s * b) % p for b in B}) >= n - t)
    thr = 2 * isqrt(n)
    exc = int((absT > thr).sum())
    return {"type": tag, "|B|": n, "|BB|/|B|": round(K, 3), "M(B)": M, "M/|B|": round(M / n, 4),
            "top_|B|_avg/|B|": round(sum(vals[:n]) / n / n, 4), "#{r:|T|>2sqrt|B|}": exc,
            "|Sym_{|B|/4}(B)|": sigma, "M2": moment(T, 1), "M4": moment(T, 2)}
table = []
H64 = subgroup(p, g, 64); H128 = subgroup(p, g, 128); H256 = subgroup(p, g, 256)
table.append(profile("coset of |H|=64", [(5 * h) % p for h in H64]))
table.append(profile("subgroup |H|=128", H128))
table.append(profile("subgroup |H|=256", H256))
table.append(profile("union of 2 cosets of |H|=64", [(5 * h) % p for h in H64] + [(7 * h) % p for h in H64]))
gp = (g * g) % p
table.append(profile("geometric progression N=128 (ratio square)", [pow(gp, i, p) for i in range(1, 129)]))
table.append(profile("geometric progression N=128 (ratio primitive)", [pow(g, i, p) for i in range(1, 129)]))
table.append(profile("subgroup 128 + 8 random", H128 + RNG.sample(range(1, p), 8)))
table.append(profile("random 128", RNG.sample(range(1, p), 128)))
table.append(profile("random 256", RNG.sample(range(1, p), 256)))
RES["sections"]["7_task3_profiles_p12289"] = table
# exceptional cosets for subgroups at p = 12289: #{cosets: |T_H| > |H|^{3/4}} and Chebyshev bounds
sec7b = []
for d in [16, 32, 48, 64, 96, 128, 192, 256, 384, 512, 768, 1024]:
    if (p - 1) % d: continue
    H = subgroup(p, g, d); T = T_of(chi, H); absT = np.abs(T)
    thr = d ** 0.75
    exc_pts = int((absT[1:] > thr).sum()); exc_cosets = exc_pts // d
    check(exc_pts % d == 0, "exceptional set is a union of cosets", (p, d))
    M4 = moment(T, 2); M6 = moment(T, 3)
    sec7b.append({"|H|": d, "M_H": int(absT[1:].max()), "M_H/|H|": round(int(absT[1:].max()) / d, 4),
                  "exceptional_cosets(|T|>|H|^0.75)": exc_cosets, "total_cosets": (p - 1) // d,
                  "Chebyshev_k2_bound_on_cosets": round(M4 / thr ** 4 / d, 2),
                  "Chebyshev_k3_bound_on_cosets": round(M6 / thr ** 6 / d, 2)})
RES["sections"]["7b_exceptional_cosets_p12289"] = sec7b

print("section 7 done", round(time.time() - T0, 1), "s")

# =====================================================================
RES["runtime_seconds"] = round(time.time() - T0, 1)
RES["n_failures"] = len(RES["failures"])
with open(OUT, "w") as f: json.dump(RES, f, indent=1, default=str)
print("checks:", RES["checks"], "failures:", RES["n_failures"], "runtime:", RES["runtime_seconds"], "s")
print("written", OUT)
