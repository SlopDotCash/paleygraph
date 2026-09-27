#!/usr/bin/env python3
"""Stepanov wave, worker `leanhp` (2026-09-26): verifier.

1. Re-elaborates every Lean file of this worker (Solutions/Theorems/development file) and
   parses the `#print axioms` output (must be exactly propext, Classical.choice, Quot.sound).
   Environment: if ~/prove2me_workspace exists (Lean v4.30.0, Mathlib c5ea003) it runs
   `lake env lean <file>` there; otherwise it runs the Lean v4.30.0-rc2 binary directly with
   LEAN_PATH set to the prebuilt package oleans of /Users/shawwalters/proximityprize
   (Mathlib 5450b53) -- read-only, no lake command, no write outside this repository.
2. Exact small-p checks of every formal statement (integers mod p only):
   - HP Theorem 1.2 (all proper divisors d of p-1), Paley case d=(p-1)/2,
   - the sharp example |B| = (p+3)/4 and |A||B| = d + r for A={0,1},
   - clique bound |A|(|A|-1) <= (p-1)/2 and the real form |A| <= (sqrt(2p-1)+1)/2,
   - the proof's internal claims: Lagrange weights sum identity, deg F = d exactly,
     root multiplicities >= |A| (resp. |A|-1).
No network access.  Writes results/stepanov_leanhp_2026_09_26.json.
"""
import itertools
import json
import os
import random
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor
from fractions import Fraction
from math import comb, isqrt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LEAN_DIR = os.path.join(ROOT, "experiments", "stepanov_leanhp_lean")
OUT = os.path.join(ROOT, "results", "stepanov_leanhp_2026_09_26.json")
P2M = os.path.expanduser("~/prove2me_workspace")
PP_PACKAGES = "/Users/shawwalters/proximityprize/.lake/packages"
RC2_LEAN = os.path.expanduser("~/.elan/toolchains/leanprover--lean4---v4.30.0-rc2/bin/lean")
EXPECTED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}

results = {"checks": {}, "failures": [], "witnesses": {}, "lean": {}}


def fail(msg):
    results["failures"].append(msg)


def bump(key, n=1):
    results["checks"][key] = results["checks"].get(key, 0) + n


# ---------------------------------------------------------------- Lean part
TLOE_PACKAGES = "/Users/shawwalters/TheLeaningOfEverything/.lake/packages"
V4291_LEAN = os.path.expanduser("~/.elan/toolchains/leanprover--lean4---v4.29.1/bin/lean")


def lean_command(path, second=False):
    if not second and os.path.isdir(P2M) and os.path.exists(os.path.join(P2M, "lakefile.lean")):
        return ["lake", "env", "lean", path], P2M, "prove2me_workspace (c5ea003, v4.30.0)"
    pk, binary, label = ((TLOE_PACKAGES, V4291_LEAN, "TheLeaningOfEverything packages (Mathlib 5e932f9, v4.29.1)")
                         if second else
                         (PP_PACKAGES, RC2_LEAN, "proximityprize packages (Mathlib 5450b53, v4.30.0-rc2)"))
    env_dirs = []
    if os.path.isdir(pk):
        for name in sorted(os.listdir(pk)):
            d = os.path.join(pk, name, ".lake", "build", "lib", "lean")
            if os.path.isdir(d):
                env_dirs.append(d)
    return [binary, path], LEAN_DIR, ("LEAN_PATH=" + ":".join(env_dirs), label)


def run_lean(path, second=False):
    cmd, cwd, env_info = lean_command(path, second)
    env = dict(os.environ)
    label = env_info
    if isinstance(env_info, tuple):
        env["LEAN_PATH"] = env_info[0].split("=", 1)[1]
        label = env_info[1]
    t0 = time.time()
    try:
        proc = subprocess.run(cmd, cwd=cwd, env=env, capture_output=True, text=True, timeout=540)
        out = proc.stdout + proc.stderr
        code = proc.returncode
    except subprocess.TimeoutExpired:
        out, code = "TIMEOUT", -1
    wall = time.time() - t0
    axioms = {}
    # '#print axioms X' prints: 'X' depends on axioms: [propext, Classical.choice, Quot.sound]
    for m in re.finditer(r"'([^']+)' depends on axioms: \[([^\]]*)\]", out):
        axioms[m.group(1)] = sorted(a.strip() for a in m.group(2).split(",") if a.strip())
    for m in re.finditer(r"'([^']+)' does not depend on any axioms", out):
        axioms[m.group(1)] = []
    errors = [ln for ln in out.splitlines() if "error" in ln]
    sorry = re.search(r"declaration uses [`']?sorry", out) is not None
    return {"file": os.path.relpath(path, ROOT), "exit": code, "wall_s": round(wall, 1),
            "env": label, "axioms": axioms, "errors": errors[:20], "uses_sorry": sorry,
            "output_tail": out[-1500:]}


def lean_part(files, expect_sorry=(), second=False):
    with ThreadPoolExecutor(max_workers=min(7, len(files))) as ex:
        reps = list(ex.map(lambda f: run_lean(f, second), files))
    for rep in reps:
        name = rep["file"] + (" [second env]" if second else "")
        results["lean"][name] = rep
        base = os.path.basename(name)
        if rep["exit"] != 0 or rep["errors"]:
            fail(f"lean: {name} exit={rep['exit']} errors={rep['errors'][:3]}")
        if base in expect_sorry:
            if not rep["uses_sorry"]:
                fail(f"lean: statement file {name} expected exactly a sorry placeholder")
        else:
            if rep["uses_sorry"]:
                fail(f"lean: {name} uses sorry")
            if not rep["axioms"]:
                fail(f"lean: {name} printed no axioms report")
            for decl, ax in rep["axioms"].items():
                if not set(ax) <= EXPECTED_AXIOMS:
                    fail(f"lean: {name}: {decl} depends on {ax}")
        bump("lean_files")


# ------------------------------------------------------------- arithmetic part
def primes_upto(n):
    return [q for q in range(3, n + 1) if all(q % r for r in range(2, isqrt(q) + 1))]


def squares0(p):
    """Set of squares in F_p including 0 (IsSquare x <-> x in this set)."""
    return {(x * x) % p for x in range(p)}


def zd0(p, d):
    return {0} | {z for z in range(1, p) if pow(z, d, p) == 1}


def hp_rhs(A, B, d, p):
    r = sum(1 for b in B if (-b) % p in A)
    return d + r


def check_hp_all(p, d, S, exhaustive):
    """Check |A||B| <= d + r for A + B subset of S = Z_d cup {0}."""
    elems = list(range(p))
    worst = None
    if exhaustive:
        As = (frozenset(c) for m in range(0, p + 1) for c in itertools.combinations(elems, m))
    else:
        rng = random.Random(1000 * p + d)
        As = []
        for _ in range(400):
            m = rng.randint(1, min(p, 6))
            As.append(frozenset(rng.sample(elems, m)))
    for A in As:
        Bmax = [b for b in elems if all((a + b) % p in S for a in A)]
        # maximal B is the binding case: each b contributes |A| - [-b in A] >= 0 to lhs - r
        for B in ([Bmax] + ([] if len(Bmax) > 12 else [])):
            lhs = len(A) * len(B)
            rhs = hp_rhs(A, B, d, p)
            bump("hp_bound")
            if lhs > rhs:
                fail(f"HP bound fails p={p} d={d} A={sorted(A)} B={sorted(B)}")
            slack = rhs - lhs
            if len(A) >= 2 and (worst is None or slack < worst[0]):
                worst = (slack, sorted(A), sorted(B))
        # random sub-B checks
        if Bmax:
            rng2 = random.Random(hash((p, d, tuple(sorted(A)))) & 0xFFFF)
            for _ in range(2):
                k = rng2.randint(0, len(Bmax))
                B = rng2.sample(Bmax, k)
                bump("hp_bound_subB")
                if len(A) * len(B) > hp_rhs(A, B, d, p):
                    fail(f"HP bound (sub-B) fails p={p} d={d} A={sorted(A)} B={sorted(B)}")
    return worst


def poly_mul(f, g, p):
    h = [0] * (len(f) + len(g) - 1)
    for i, x in enumerate(f):
        if x:
            for j, y in enumerate(g):
                h[i + j] = (h[i + j] + x * y) % p
    return h


def poly_trim(f):
    while len(f) > 1 and f[-1] == 0:
        f = f[:-1]
    return f


def aux_poly(A, d, p):
    """F(x) = -1 + sum_a c_a (x+a)^D, c_a = 1/prod_{a'!=a}(a-a'), D = d+|A|-1 (coeff list)."""
    M = len(A)
    D = d + M - 1
    F = [0] * (D + 1)
    cs = {}
    for a in A:
        den = 1
        for a2 in A:
            if a2 != a:
                den = den * (a - a2) % p
        c = pow(den, p - 2, p)
        cs[a] = c
        for j in range(D + 1):  # (x+a)^D = sum_j C(D,j) a^{D-j} x^j
            F[j] = (F[j] + c * comb(D, j) * pow(a, D - j, p)) % p
    F[0] = (F[0] - 1) % p
    return poly_trim(F), cs


def root_mult(F, b, p):
    """Multiplicity of root b of F over F_p via repeated synthetic division."""
    m = 0
    f = F[:]
    while len(f) > 1 or f[0] != 0:
        # evaluate / divide by (x - b)
        n = len(f) - 1
        q = [0] * n
        acc = 0
        for i in range(n, -1, -1):
            acc = (acc * b + f[i]) % p
            if i > 0:
                q[i - 1] = acc
        if acc != 0 or n == 0:
            break
        m += 1
        f = q if q else [0]
    return m


def check_construction(p, d, S, trials):
    rng = random.Random(7 * p + d)
    elems = list(range(p))
    for _ in range(trials):
        M = rng.randint(1, min(d + 1, p - d, 6))
        A = rng.sample(elems, M)
        # Lagrange weight identity sum_a c_a (t+a)^n = [n = M-1], n < M
        _, cs = aux_poly(A, d, p)
        for t in range(p):
            for n in range(M):
                s = sum(cs[a] * pow(t + a, n, p) for a in A) % p
                bump("lagrange_identity")
                if s != (1 if n == M - 1 else 0):
                    fail(f"Lagrange identity fails p={p} A={A} t={t} n={n}")
        F, _ = aux_poly(A, d, p)
        bump("aux_degree")
        if len(F) - 1 != d:
            fail(f"deg F != d: p={p} d={d} A={A} deg={len(F)-1}")
        if F[-1] != comb(d + M - 1, d) % p:
            fail(f"leading coeff != C(D,d): p={p} d={d} A={A}")
        for b in elems:
            if all((a + b) % p in S for a in A):
                need = M - 1 if (-b) % p in A else M
                mult = root_mult(F, b, p)
                bump("aux_root_mult")
                if mult < need:
                    fail(f"root multiplicity {mult} < {need}: p={p} d={d} A={A} b={b}")


def taylor_coeffs(F, b, p):
    """Coefficients of F(x + b) over F_p."""
    n = len(F)
    out = [0] * n
    for i, c in enumerate(F):
        if c:
            for j in range(i + 1):
                out[j] = (out[j] + c * comb(i, j) * pow(b, i - j, p)) % p
    return out


def check_bad_partner_formula(p, trials):
    """Brief's formula (Paley case d=(p-1)/2): for b with b+a != 0 for all a in A,
    [x^j] F(x+b) = -2 C(D,j) sum_{a in E(b)} c_a (b+a)^{M-1-j} for 0 <= j <= M-1,
    E(b) = {a : chi(b+a) = -1}; equivalently F^{(j)}(b) = -2 (D)_j sum_E ... ."""
    d = (p - 1) // 2
    Q = {(x * x) % p for x in range(1, p)}
    rng = random.Random(99 * p)
    for _ in range(trials):
        M = rng.randint(1, min(d + 1, 6))
        A = rng.sample(range(p), M)
        F, cs = aux_poly(A, d, p)
        D = d + M - 1
        for b in range(p):
            if any((b + a) % p == 0 for a in A):
                continue
            T = taylor_coeffs(F, b, p)
            E = [a for a in A if (b + a) % p not in Q]
            for j in range(M):
                rhs = (-2 * comb(D, j) * sum(cs[a] * pow(b + a, M - 1 - j, p) for a in E)) % p
                bump("bad_partner_formula")
                if (T[j] if j < len(T) else 0) != rhs:
                    fail(f"bad-partner formula fails p={p} A={A} b={b} j={j}")


def check_sharp(p):
    Q0 = squares0(p)
    B = [b for b in range(p) if b in Q0 and (b + 1) % p in Q0]
    A = {0, 1}
    ok_card = len(B) == (p + 3) // 4
    ok_sum = all((a + b) % p in Q0 for a in A for b in B)
    r = sum(1 for b in B if (-b) % p in A)
    ok_eq = len(A) * len(B) == (p - 1) // 2 + r
    bump("sharp_example")
    if not (ok_card and ok_sum and ok_eq):
        fail(f"sharp example fails p={p}: |B|={len(B)} r={r}")
    return {"p": p, "B_card": len(B), "r": r, "AB": 2 * len(B), "d_plus_r": (p - 1) // 2 + r}


def max_clique(p, Q):
    """Exact clique number of the Paley graph (Bron-Kerbosch with pivoting)."""
    adj = {v: {u for u in range(p) if u != v and (u - v) % p in Q} for v in range(p)}
    best = [0]

    def bk(R, P, X):
        if not P and not X:
            best[0] = max(best[0], len(R))
            return
        if len(R) + len(P) <= best[0]:
            return
        u = max(P | X, key=lambda w: len(adj[w] & P))
        for v in list(P - adj[u]):
            bk(R | {v}, P & adj[v], X & adj[v])
            P = P - {v}
            X = X | {v}

    bk(set(), set(range(p)), set())
    return best[0]


def check_clique(p):
    Q = {(x * x) % p for x in range(1, p)}
    w = max_clique(p, Q)
    ok_nat = w * (w - 1) <= (p - 1) // 2
    # real form: w <= (sqrt(2p-1)+1)/2  <=>  (2w-1)^2 <= 2p-1  (w >= 1)
    ok_real = (2 * w - 1) ** 2 <= 2 * p - 1
    bump("clique_bound")
    if not (ok_nat and ok_real):
        fail(f"clique bound fails p={p} omega={w}")
    return {"p": p, "omega": w, "w(w-1)": w * (w - 1), "(p-1)/2": (p - 1) // 2}


def arithmetic_part():
    t0 = time.time()
    worst = {}
    for p in primes_upto(47):
        divs = [d for d in range(1, p - 1) if (p - 1) % d == 0]
        for d in divs:
            S = zd0(p, d)
            w = check_hp_all(p, d, S, exhaustive=(p <= 13))
            if w is not None:
                worst[f"p={p},d={d}"] = w
            check_construction(p, d, S, trials=6 if p <= 23 else 2)
        check_bad_partner_formula(p, trials=4 if p <= 23 else 1)
        # Paley statement uses IsSquare, i.e. S = squares incl. 0 = Z_{(p-1)/2} cup {0}
        bump("paley_set_equals_Zd0")
        if squares0(p) != zd0(p, (p - 1) // 2):
            fail(f"squares0 != Z_d0 at p={p}")
    results["witnesses"]["hp_min_slack_|A|>=2"] = {k: v for k, v in worst.items() if v[0] == 0}
    results["witnesses"]["sharp_example"] = [check_sharp(p) for p in primes_upto(200) if p % 4 == 1]
    results["witnesses"]["clique"] = [check_clique(p) for p in primes_upto(109) if p % 4 == 1]
    # p = 2 counterexample showing hp : p != 2 is needed in the Paley statement
    bump("p2_counterexample")
    A = B = {0, 1}
    results["witnesses"]["p2_counterexample"] = {
        "A": [0, 1], "B": [0, 1], "AB": 4,
        "rhs": (2 - 1) // 2 + sum(1 for b in B if (-b) % 2 in A)}
    if not 4 > results["witnesses"]["p2_counterexample"]["rhs"]:
        fail("p=2 counterexample check wrong")
    results["arith_wall_s"] = round(time.time() - t0, 1)


def norm_stmt(text):
    return " ".join(text.split())


def statement_match(sol_dir, thm_dir):
    """Solution type == target statement, textually (binders and type after the theorem name)."""
    for f in sorted(os.listdir(thm_dir)):
        if not f.startswith("Thm_") or not f.endswith(".lean"):
            continue
        name = f[len("Thm_"):-len(".lean")]
        thm = open(os.path.join(thm_dir, f)).read()
        m = re.search(r"theorem " + re.escape(name) + r"\s(.*?):= by sorry", thm, re.S)
        solp = os.path.join(sol_dir, "Sol_" + name + ".lean")
        bump("statement_match")
        if not m or not os.path.exists(solp):
            fail(f"statement match: missing statement or solution for {name}")
            continue
        sol = open(solp).read()
        m2 = re.search(r"theorem solution\s(.*?):=\n", sol, re.S)
        if not m2 or norm_stmt(m.group(1)) != norm_stmt(m2.group(1)):
            fail(f"statement match: solution type differs from target for {name}")
        # the solution file must not import platform theorems and must not contain sorry
        if re.search(r"^import Theorems", sol, re.M) or "sorry" in sol:
            fail(f"solution file {name} imports Theorems or contains sorry")


def main():
    t0 = time.time()
    files = []
    sol_dir = os.path.join(LEAN_DIR, "Solutions")
    thm_dir = os.path.join(LEAN_DIR, "Theorems")
    stmt_files = set()
    for sub in (sol_dir, thm_dir):
        if os.path.isdir(sub):
            for f in sorted(os.listdir(sub)):
                if f.endswith(".lean"):
                    files.append(os.path.join(sub, f))
                    if sub == thm_dir:
                        stmt_files.add(f)
    if "--with-dev" in sys.argv:
        for dev in ("StepanovHP.lean", "StepanovSharp.lean"):
            if os.path.exists(os.path.join(LEAN_DIR, dev)):
                files.append(os.path.join(LEAN_DIR, dev))
    if os.path.isdir(sol_dir) and os.path.isdir(thm_dir):
        statement_match(sol_dir, thm_dir)
    if "--no-lean" not in sys.argv:
        lean_part(files, expect_sorry=stmt_files)
    if "--second-env" in sys.argv:
        lean_part([f for f in files if "/Solutions/" in f], second=True)
    arithmetic_part()
    results["total_checks"] = sum(results["checks"].values())
    results["n_failures"] = len(results["failures"])
    results["wall_s"] = round(time.time() - t0, 1)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w") as fh:
        json.dump(results, fh, indent=1, default=str)
    print(json.dumps({"total_checks": results["total_checks"], "failures": results["failures"][:10],
                      "wall_s": results["wall_s"],
                      "lean": {k: (v["exit"], v["axioms"], v["uses_sorry"], v["wall_s"])
                               for k, v in results["lean"].items()}}, indent=1))


if __name__ == "__main__":
    main()
