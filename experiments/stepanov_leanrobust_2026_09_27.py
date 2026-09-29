#!/usr/bin/env python3
"""Stepanov wave, worker `leanrobust` (2026-09-27): verifier.

1. Lean: re-elaborates experiments/stepanov_leanrobust_lean/StepanovRobust.lean with the `lean`
   binary directly, LEAN_PATH set to prebuilt package oleans (read-only; no lake command; nothing is
   written outside this repository):
     E1 = Lean v4.30.0-rc2 / Mathlib 5450b53  (~/proximityprize/.lake/packages)
     E2 = Lean v4.29.1     / Mathlib 5e932f9  (~/TheLeaningOfEverything/.lake/packages)
   Records exit code, errors/warnings, `sorry`, and parses every `#print axioms` line (must be
   exactly propext, Classical.choice, Quot.sound).
2. Exact modular checks (integers mod p only) of the formal statements and of the internal claims
   of the formal proof:
   - (*) of `hankel_stepanov_paley` exactly as stated in Lean (sum over all b of
     sum_{i=e_b}^{e} (|A| - delta_b - (i+e)) <= (e+1)((p-1)/2 - e)), exhaustive for p <= 13;
   - the polynomial-level claims for the formal definitions `uPoly`, `hankel`: natDegree
     <= (e+1)(d-e), top coefficient = hankelLead, and (X-b)^{ord_b} | det;
   - `hankelLead_ne_zero` (Lambda != 0 mod p) on all admissible (d, m, e) for small p, plus the
     factorisation `choose_factor` and the triangular witness evaluations used in its proof;
   - `bias_inequality`, `charSum_eq`, `bias_bound`, `subset_bound` and `constant_bias` on random and
     structured sets.
No network access.  Writes results/stepanov_leanrobust_2026_09_27.json.
Flags: --no-lean (skip Lean), --e1-only (skip E2), --timeout=SECONDS per Lean run (default 3600).
On an idle machine the whole run takes about two to three minutes (each Lean run needs roughly
1-2 CPU-minutes; E1 and E2 run in parallel); on a heavily loaded machine the wall time is longer.
"""
import itertools
import json
import math
import os
import random
import re
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEAN_DIR = os.path.join(ROOT, "experiments", "stepanov_leanrobust_lean")
LEAN_FILES = ["StepanovRobust.lean"]
OUT = os.path.join(ROOT, "results", "stepanov_leanrobust_2026_09_27.json")
ENVS = {
    "E1": ("/Users/shawwalters/proximityprize/.lake/packages",
           os.path.expanduser("~/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/bin/lean"),
           "Lean v4.30.0-rc2 / Mathlib 5450b53 (proximityprize packages, read-only)"),
    "E2": ("/Users/shawwalters/TheLeaningOfEverything/.lake/packages",
           os.path.expanduser("~/.elan/toolchains/leanprover--lean4---v4.29.1/bin/lean"),
           "Lean v4.29.1 / Mathlib 5e932f9 (TheLeaningOfEverything packages, read-only)"),
}
EXPECTED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
TIMEOUT = [3600]

results = {"checks": {}, "failures": [], "witnesses": {}, "lean": {}, "notes": []}
random.seed(20260927)


def fail(msg):
    results["failures"].append(msg)


def bump(key, n=1):
    results["checks"][key] = results["checks"].get(key, 0) + n


# ------------------------------------------------------------------------------------ Lean
def run_lean(env_name, fname):
    pk, binary, label = ENVS[env_name]
    dirs = []
    if os.path.isdir(pk):
        for name in sorted(os.listdir(pk)):
            d = os.path.join(pk, name, ".lake", "build", "lib", "lean")
            if os.path.isdir(d):
                dirs.append(d)
    env = dict(os.environ)
    env["LEAN_PATH"] = ":".join(dirs)
    path = os.path.join(LEAN_DIR, fname)
    t0 = time.time()
    rec = {"env": label, "file": os.path.relpath(path, ROOT)}
    if not os.path.exists(binary) or not dirs:
        rec.update({"exit": None, "error": "toolchain or packages missing"})
        return rec
    timeout = TIMEOUT[0]
    try:
        proc = subprocess.Popen([binary, path], cwd=LEAN_DIR, env=env, stdout=subprocess.PIPE,
                                stderr=subprocess.STDOUT, text=True)
        chunks = []
        reader = threading.Thread(target=lambda: chunks.append(proc.stdout.read()))
        reader.start()
        timer = threading.Timer(timeout, proc.kill)
        timer.start()
        _, status, ru = os.wait4(proc.pid, 0)
        timer.cancel()
        reader.join()
        proc.returncode = os.waitstatus_to_exitcode(status)
        out = "".join(chunks)
        code = proc.returncode if proc.returncode >= 0 else None
        rec["cpu_seconds"] = round(ru.ru_utime + ru.ru_stime, 1)
        rec["max_rss_mb"] = round(ru.ru_maxrss / 2 ** 20, 1)  # macOS: bytes
        if code is None:
            out += "\nTIMEOUT/KILLED"
    except OSError as exc:
        out, code = f"OSError {exc}", None
    rec["exit"] = code
    rec["wall_seconds"] = round(time.time() - t0, 1)
    rec["errors"] = [l for l in out.splitlines() if ": error" in l]
    rec["warnings"] = [l for l in out.splitlines() if ": warning" in l]
    rec["mentions_sorry"] = ("sorry" in out)
    axioms = {}
    for m in re.finditer(r"'([^']+)' depends on axioms: \[([^\]]*)\]", out):
        axioms[m.group(1)] = [a.strip() for a in m.group(2).split(",") if a.strip()]
    for m in re.finditer(r"'([^']+)' does not depend on any axioms", out):
        axioms[m.group(1)] = []
    rec["axioms"] = axioms
    n_expected = open(path).read().count("#print axioms")
    rec["print_axioms_lines_in_source"] = n_expected
    rec["axioms_ok"] = (len(axioms) == n_expected and n_expected > 0 and
                        all(set(v) <= EXPECTED_AXIOMS for v in axioms.values()))
    rec["raw_output_tail"] = out[-3000:]
    return rec


def lean_part(second):
    src_ok = True
    for f in LEAN_FILES:
        txt = open(os.path.join(LEAN_DIR, f)).read()
        # the words may appear in comments; forbid the tactic/term in code only
        # (root fix 2026-09-29: strip block and doc comments too; the English word
        # "admit" in a doc comment had triggered a false failure)
        code = re.sub(r"/-.*?-/", " ", txt, flags=re.S)
        code = "\n".join(l.split("--", 1)[0] for l in code.splitlines())
        if re.search(r"\bsorry\b", code) or re.search(r"\badmit\b", code):
            fail(f"{f}: source contains sorry/admit")
            src_ok = False
    bump("lean_source_no_sorry", 1)
    jobs = [("E1", f) for f in LEAN_FILES]
    if second:
        jobs += [("E2", f) for f in LEAN_FILES]
    with ThreadPoolExecutor(max_workers=len(jobs)) as ex:
        recs = list(ex.map(lambda j: (j, run_lean(*j)), jobs))
    for (envn, f), rec in recs:
        results["lean"][f"{envn}:{f}"] = rec
        bump("lean_runs")
        if rec.get("exit") != 0 or rec.get("errors") or rec.get("mentions_sorry"):
            fail(f"Lean {envn} {f}: exit={rec.get('exit')} errors={rec.get('errors')[:3] if rec.get('errors') else None}")
        if not rec.get("axioms_ok"):
            fail(f"Lean {envn} {f}: axioms {rec.get('axioms')}")
    return src_ok


# --------------------------------------------------------------------- arithmetic helpers
def primes_upto(n):
    return [q for q in range(3, n + 1) if all(q % r for r in range(2, int(q ** 0.5) + 1))]


def is_square_table(p):
    sq = [False] * p
    for x in range(p):
        sq[x * x % p] = True
    return sq  # sq[0] = True


def chi_table(p):
    sq = is_square_table(p)
    return [0 if x == 0 else (1 if sq[x] else -1) for x in range(p)]


def weights(A, p):
    """c_a = [x^{|A|-1}] L_a = prod_{a' != a} (a - a')^{-1}."""
    w = {}
    for a in A:
        prod = 1
        for b in A:
            if b != a:
                prod = prod * (a - b) % p
        w[a] = pow(prod, p - 2, p)
    return w


def e_b_list(A, p, sq):
    return [sum(1 for a in A if not sq[(a + b) % p]) for b in range(p)]


def star_lhs(A, p, e, sq):
    m = len(A)
    Aset = set(A)
    tot = 0
    for b in range(p):
        eb = sum(1 for a in A if not sq[(a + b) % p])
        delta = 1 if (-b) % p in Aset else 0
        for i in range(eb, e + 1):
            tot += max(m - delta - (i + e), 0)  # natural-number subtraction as in Lean
    return tot


# polynomials over F_p: list of coefficients, index = degree
def padd(f, g, p):
    n = max(len(f), len(g))
    return trim([((f[i] if i < len(f) else 0) + (g[i] if i < len(g) else 0)) % p for i in range(n)])


def pmul(f, g, p):
    if not f or not g:
        return []
    r = [0] * (len(f) + len(g) - 1)
    for i, a in enumerate(f):
        if a:
            for j, b in enumerate(g):
                r[i + j] = (r[i + j] + a * b) % p
    return trim(r)


def pscale(f, c, p):
    return trim([a * c % p for a in f])


def trim(f):
    f = list(f)
    while f and f[-1] == 0:
        f.pop()
    return f


def shift_pow(a, N, p):
    """(X + a)^N."""
    return trim([math.comb(N, k) * pow(a, N - k, p) % p for k in range(N + 1)])


def taylor_shift(f, b, p):
    """f(X + b)."""
    r = []
    for k, c in enumerate(f):
        if c:
            r = padd(r, pscale(shift_pow(b, k, p), c, p), p)
    return r


def order_at_zero(f):
    for i, c in enumerate(f):
        if c:
            return i
    return None


def u_poly(A, w, d, s, p):
    m = len(A)
    D = d + m - 1
    r = [p - 1] if s == 0 else []
    for a in A:
        r = padd(r, pscale(shift_pow(a, D - s, p), w[a], p), p)
    return r


def det_poly(M, p):
    n = len(M)
    tot = []
    for perm in itertools.permutations(range(n)):
        sign = 1
        seen = list(perm)
        for i in range(n):
            for j in range(i + 1, n):
                if seen[i] > seen[j]:
                    sign = -sign
        term = [1]
        for i in range(n):
            term = pmul(term, M[perm[i]][i], p)
        tot = padd(tot, pscale(term, sign % p, p), p)
    return tot


def det_mod(M, p):
    M = [row[:] for row in M]
    n = len(M)
    det = 1
    for c in range(n):
        piv = next((r for r in range(c, n) if M[r][c] % p), None)
        if piv is None:
            return 0
        if piv != c:
            M[c], M[piv] = M[piv], M[c]
            det = -det
        det = det * M[c][c] % p
        inv = pow(M[c][c], p - 2, p)
        for r in range(c + 1, n):
            f = M[r][c] * inv % p
            if f:
                for k in range(c, n):
                    M[r][k] = (M[r][k] - f * M[c][k]) % p
    return det % p


def lead_matrix(d, m, e):
    D = d + m - 1
    return [[math.comb(D - (i + j), m - 1) for j in range(e + 1)] for i in range(e + 1)]


def desc(n, k):
    r = 1
    for t in range(k):
        r *= (n - t)
    return r


# ------------------------------------------------------------------------------ checks
def check_star_exhaustive():
    """hankel_stepanov_paley, exactly as stated, for every A (p <= 13) and admissible e."""
    worst = (0, None)
    for p in [3, 5, 7, 11, 13]:
        sq = is_square_table(p)
        d = (p - 1) // 2
        for m in range(1, (p + 1) // 2 + 1):
            for A in itertools.combinations(range(p), m):
                for e in range(0, (m - 1) // 2 + 1):
                    lhs = star_lhs(A, p, e, sq)
                    rhs = (e + 1) * (d - e)
                    bump("star_exhaustive")
                    if lhs > rhs:
                        fail(f"(*) fails p={p} A={A} e={e}: {lhs} > {rhs}")
                    if e >= 1 and rhs and Fraction(lhs, rhs) > (worst[0] or 0):
                        worst = (Fraction(lhs, rhs), (p, A, e, lhs, rhs))
    results["witnesses"]["star_tightest_e_ge_1_exhaustive"] = (
        {"ratio": str(worst[0]), "p": worst[1][0], "A": list(worst[1][1]), "e": worst[1][2],
         "lhs": worst[1][3], "rhs": worst[1][4]} if worst[1] else None)


def check_star_random():
    worst = (0, None)
    for p in primes_upto(211):
        if p < 17:
            continue
        sq = is_square_table(p)
        d = (p - 1) // 2
        trials = 40 if p < 100 else 12
        for _ in range(trials):
            m = random.randint(1, (p + 1) // 2)
            kind = random.random()
            if kind < 0.3:
                start = random.randrange(p)
                A = sorted({(start + t) % p for t in range(m)})
            elif kind < 0.5:
                A = sorted({x * x % p for x in range(1, p)} | {0})[:m]
            else:
                A = sorted(random.sample(range(p), m))
            m = len(A)
            for e in sorted({0, 1, 2, (m - 1) // 4, (m - 1) // 2}):
                if 2 * e + 1 > m:
                    continue
                lhs = star_lhs(A, p, e, sq)
                rhs = (e + 1) * (d - e)
                bump("star_random")
                if lhs > rhs:
                    fail(f"(*) fails p={p} A={A} e={e}")
                if e >= 1 and rhs and Fraction(lhs, rhs) > (worst[0] or 0):
                    worst = (Fraction(lhs, rhs), (p, m, e, lhs, rhs))
    results["witnesses"]["star_tightest_e_ge_1_random"] = (
        {"ratio": float(worst[0]), "p": worst[1][0], "m": worst[1][1], "e": worst[1][2],
         "lhs": worst[1][3], "rhs": worst[1][4]} if worst[1] else None)


def check_polynomial_claims():
    """Formal definitions uPoly/hankel: degree, top coefficient, orders at every b."""
    for p in [7, 11, 13, 17, 19, 23]:
        sq = is_square_table(p)
        d = (p - 1) // 2
        for _ in range(6 if p < 20 else 3):
            m = random.randint(3, (p + 1) // 2)
            A = sorted(random.sample(range(p), m))
            w = weights(A, p)
            Aset = set(A)
            for e in range(1, min(2, (m - 1) // 2) + 1):
                M = [[u_poly(A, w, d, i + j, p) for j in range(e + 1)] for i in range(e + 1)]
                H = det_poly(M, p)
                N = (e + 1) * (d - e)
                deg = len(H) - 1
                lead = det_mod(lead_matrix(d, m, e), p)
                bump("poly_degree_and_top_coefficient")
                if deg > N or (H[N] if N < len(H) else 0) != lead:
                    fail(f"degree/top coefficient p={p} A={A} e={e}: deg={deg} N={N}")
                if lead == 0:
                    fail(f"Lambda = 0 mod p for p={p} m={m} e={e}")
                for b in range(p):
                    eb = sum(1 for a in A if not sq[(a + b) % p])
                    delta = 1 if (-b) % p in Aset else 0
                    need = sum(max(m - delta - (i + e), 0) for i in range(eb, e + 1))
                    o = order_at_zero(taylor_shift(H, b, p))
                    bump("poly_order_at_b")
                    if o is None or o < need:
                        fail(f"order p={p} A={A} e={e} b={b}: {o} < {need}")
                # Step 1 internal claim: (X-b)^{m-delta-s} | u_s - 2 rho_s
                for b in range(p):
                    E = [a for a in A if pow((b + a) % p, d, p) == p - 1]
                    delta = 1 if (-b) % p in Aset else 0
                    for s in range(2 * e + 1):
                        rho = []
                        for a in E:
                            rho = padd(rho, pscale(shift_pow(a, d + m - 1 - s, p), w[a], p), p)
                        eps = padd(u_poly(A, w, d, s, p), pscale(rho, p - 2, p), p)
                        o = order_at_zero(taylor_shift(eps, b, p))
                        bump("step1_eps_order")
                        if o is not None and o < m - delta - s:
                            fail(f"eps order p={p} A={A} b={b} s={s}")


def check_lead_nonvanishing():
    """hankelLead_ne_zero: all (d, m, e) with 2e <= d, 2e+1 <= m, d+m-1 < p, p <= 43 (every d),
    and the Paley d = (p-1)/2 for p <= 97; plus the proof's factorisation and witnesses."""
    for p in primes_upto(97):
        ds = range(1, p) if p <= 43 else [(p - 1) // 2]
        for d in ds:
            for m in range(1, p - d + 1):
                if d + m - 1 >= p:
                    continue
                for e in range(0, min(d // 2, (m - 1) // 2) + 1):
                    bump("lambda_nonzero_mod_p")
                    if det_mod(lead_matrix(d, m, e), p) == 0:
                        fail(f"Lambda = 0 mod p: p={p} d={d} m={m} e={e}")
    # necessity of D < p: a witness where Lambda = 0 mod p once D >= p
    wit = None
    for p in primes_upto(40):
        for d in range(1, 2 * p):
            for m in range(3, 2 * p):
                if d + m - 1 < p:
                    continue
                for e in range(1, min(d // 2, (m - 1) // 2) + 1):
                    if det_mod(lead_matrix(d, m, e), p) == 0:
                        wit = {"p": p, "d": d, "m": m, "e": e, "D": d + m - 1}
                        break
                if wit:
                    break
            if wit:
                break
        if wit:
            break
    results["witnesses"]["lambda_vanishes_when_D_ge_p"] = wit
    # choose_factor (integer identity) and the triangular witnesses
    for d in range(1, 26):
        for K in range(0, 26):
            for e in range(0, d // 2 + 1):
                if 2 * e + 1 > K + 1:
                    continue
                for i in range(e + 1):
                    for j in range(e + 1):
                        lhs = math.comb(d + K - (i + j), K) * math.factorial(K) * math.factorial(d - i)
                        rhs = math.factorial(d + K - i - e) * desc(d - i, j) * desc(d + K - i - j, e - j)
                        bump("choose_factor_identity")
                        if lhs != rhs:
                            fail(f"choose_factor d={d} K={K} e={e} i={i} j={j}")
                for k in range(e + 1):
                    for j in range(e + 1):
                        v = desc(k, j) * desc(K + k - j, e - j)
                        # leadPoly evaluated at x = d - k, as a polynomial identity check
                        direct = desc(d - (d - k), j) * desc(d + K - j - (d - k), e - j)
                        bump("witness_evaluation")
                        if v != direct or (j > k and v != 0) or (j == k and v == 0):
                            fail(f"witness d={d} K={K} e={e} k={k} j={j}")


def check_bias_and_constant():
    worst_c25 = None
    worst_sub = None
    for p in primes_upto(163):
        if p < 11:
            continue
        chi = chi_table(p)
        sq = is_square_table(p)
        d = (p - 1) // 2
        trials = 30 if p < 60 else 10
        for _ in range(trials):
            mA = random.randint(1, (p + 1) // 2)
            A = sorted(random.sample(range(p), mA))
            B = sorted(random.sample(range(p), random.randint(0, p)))
            S = sum(chi[(a + b) % p] for a in A for b in B)
            n = len(B)
            r = sum(1 for b in B if (-b) % p in set(A))
            Nm = sum(sum(1 for a in A if not sq[(a + b) % p]) for b in B)
            bump("charSum_eq")
            if S != mA * n - r - 2 * Nm:
                fail(f"charSum_eq p={p}")
            for e in range(0, (mA - 1) // 2 + 1):
                bump("bias_inequality")
                if (mA - 2 * e) * ((e + 1) * n - Nm) > (e + 1) * (d - e + r):
                    fail(f"bias_inequality p={p} A={A} B={B} e={e}")
                bound = Fraction(mA * n - r) - 2 * (e + 1) * (n - Fraction(d - e + r, mA - 2 * e))
                bump("bias_bound")
                if S > bound:
                    fail(f"bias_bound p={p} e={e}")
        # Theorem 2.5 and the per-subset bound on sets with |A||B| >= (1/2 + kappa) p
        for _ in range(trials):
            kappa = random.choice([0.01, 0.05, 0.1, 0.25, 0.5, 1.0, 1.5])
            need = (0.5 + kappa) * p
            mA = random.randint(1, p)
            nmin = math.ceil(need / mA)
            if nmin > p:
                continue
            nB = random.randint(nmin, p)
            A = random.sample(range(p), mA)
            B = random.sample(range(p), nB)
            if random.random() < 0.3:
                # structured: A = {0,1}, B = {b : b, b+1 squares} plus padding (near the threshold)
                A = [0, 1]
                B0 = [b for b in range(p) if sq[b] and sq[(b + 1) % p]]
                rest = [x for x in range(p) if x not in B0]
                random.shuffle(rest)
                B = B0 + rest[:max(0, math.ceil(need / 2) - len(B0))]
            mA, nB = len(A), len(B)
            if mA * nB < need:
                continue
            S = sum(chi[(a + b) % p] for a in A for b in B)
            u = 1 / math.sqrt(1 + 2 * kappa)
            c = 1 - (1 - u) ** 2 + (math.sqrt(need) + 1) / (2 * (p - 1))
            slack = c * mA * nB - abs(S)
            bump("constant_bias")
            if slack < -1e-9:
                fail(f"constant_bias p={p} kappa={kappa} |A|={mA} |B|={nB} S={S}")
            rel = abs(S) / (mA * nB) / c
            if worst_c25 is None or rel > worst_c25["ratio_to_bound"]:
                worst_c25 = {"p": p, "kappa": kappa, "m": mA, "n": nB, "S": S,
                             "bound_factor": c, "ratio_to_bound": rel}
            if 1 <= mA <= (p + 1) // 2:
                cs = 1 - (1 - u) ** 2 + u * (1 - u) * mA / d
                bump("subset_bound")
                if S > cs * mA * nB + 1e-9:
                    fail(f"subset_bound p={p} kappa={kappa}")
                rel2 = S / (mA * nB) / cs
                if worst_sub is None or rel2 > worst_sub["ratio_to_bound"]:
                    worst_sub = {"p": p, "kappa": kappa, "m": mA, "n": nB, "S": S,
                                 "bound_factor": cs, "ratio_to_bound": rel2}
    results["witnesses"]["constant_bias_tightest"] = worst_c25
    results["witnesses"]["subset_bound_tightest"] = worst_sub


def main():
    t0 = time.time()
    args = sys.argv[1:]
    with_lean = "--no-lean" not in args
    for a in args:
        if a.startswith("--timeout="):
            TIMEOUT[0] = int(a.split("=", 1)[1])
    second = "--e1-only" not in args
    lean_thread = None
    if with_lean:
        ex = ThreadPoolExecutor(max_workers=1)
        lean_thread = ex.submit(lean_part, second)
    check_star_exhaustive()
    check_star_random()
    check_lead_nonvanishing()
    check_polynomial_claims()
    check_bias_and_constant()
    if lean_thread is not None:
        lean_thread.result()
    results["total_checks"] = sum(results["checks"].values())
    results["wall_seconds"] = round(time.time() - t0, 1)
    results["ok"] = not results["failures"]
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as f:
        json.dump(results, f, indent=1, default=str)
    print(json.dumps({"ok": results["ok"], "total_checks": results["total_checks"],
                      "checks": results["checks"], "failures": results["failures"][:10],
                      "wall_seconds": results["wall_seconds"],
                      "lean": {k: {"exit": v.get("exit"), "axioms_ok": v.get("axioms_ok"),
                                   "wall": v.get("wall_seconds"), "cpu": v.get("cpu_seconds"),
                                   "warnings": len(v.get("warnings", []))}
                               for k, v in results["lean"].items()}}, indent=1))


if __name__ == "__main__":
    main()
