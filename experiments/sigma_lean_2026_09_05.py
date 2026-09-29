#!/usr/bin/env python3
"""sigma pass, direction `lean` (2026-09-05): verifier.

Re-runs `lake env lean` on every Solutions/Sol_paley_*.lean and Theorems/Thm_paley_*.lean in the
prove2me workspace, prints the axioms of every `solution`, checks each solution's `theorem solution`
against its theorem statement (normalized text + a Lean `rfl` defeq check), checks the proposal JSON
against the Lean files, and re-tests every formalized identity/inequality by exact integer
arithmetic at small primes.  Downloads nothing.  Standard library only (numpy/scipy/sympy unused).

Writes results/sigma_lean_2026_09_05.json.  Wall time on this machine: ~2 minutes.

Environment variables:  PROVE2ME_WORKSPACE (default ~/prove2me_workspace),
                        SIGMA_LEAN_SKIP_LEAN=1 to skip the Lean runs (integer checks only).
"""
import json, os, pathlib, random, re, subprocess, sys, time
from itertools import combinations

HERE = pathlib.Path(__file__).resolve().parent
PALEY = HERE.parent
WS = pathlib.Path(os.environ.get("PROVE2ME_WORKSPACE", pathlib.Path.home() / "prove2me_workspace"))
SCRATCH = PALEY / "tmp" / "sigma_lean"          # brief: never use /tmp
OUT = PALEY / "results" / "sigma_lean_2026_09_05.json"
ALLOWED_AXIOMS = {"propext", "Classical.choice", "Quot.sound"}
SOLVED = ["shift_orthogonality", "second_moment", "chung_bound", "interval_nonresidue"]
OPEN = ["two_set_conjecture", "karatsuba_amplification",
        "hanson_petridis_difference_bound", "hanson_petridis_clique_number"]
T0 = time.time()
R = {"date": "2026-09-05", "direction": "lean", "workspace": str(WS), "checks": 0, "failures": [],
     "lean": {}, "integer_checks": {}, "witnesses": {}}

def note(msg):
    print(f"[{time.time()-T0:7.1f}s] {msg}", flush=True)

def check(cond, what, witness=None):
    R["checks"] += 1
    if not cond:
        R["failures"].append({"what": what, "witness": witness})
        note(f"FAIL: {what} {witness if witness is not None else ''}")
    return cond

def lean(path: pathlib.Path, timeout=600):
    t = time.time()
    p = subprocess.run(["lake", "env", "lean", str(path)], cwd=WS, capture_output=True, text=True, timeout=timeout)
    return {"file": str(path.relative_to(WS)) if path.is_relative_to(WS) else str(path),
            "returncode": p.returncode, "seconds": round(time.time() - t, 2),
            "stdout": p.stdout[-4000:], "stderr": p.stderr[-4000:]}

# ----------------------------------------------------------------------------- Lean runs
def lean_runs():
    SCRATCH.mkdir(parents=True, exist_ok=True)
    r = lean(WS / "Solutions" / "SmokeTest.lean")
    R["lean"]["smoke"] = r
    check(r["returncode"] == 0, "SmokeTest.lean compiles", r)
    for n in SOLVED:
        sol = WS / "Solutions" / f"Sol_paley_{n}.lean"
        thm = WS / "Theorems" / f"Thm_paley_{n}.lean"
        stxt, ttxt = sol.read_text(), thm.read_text()
        check("sorry" not in stxt, f"no sorry in {sol.name}")
        check(re.search(r"^theorem solution\b", stxt, re.M) is not None, f"top-level `theorem solution` in {sol.name}")
        check(not re.search(r"^\s*import\s+Theorems\.", stxt, re.M), f"{sol.name} does not import Theorems.*")
        r = lean(sol); R["lean"][f"Sol_{n}"] = r
        check(r["returncode"] == 0 and r["stdout"].strip() == "" and r["stderr"].strip() == "",
              f"Sol_paley_{n}.lean compiles with empty output", r)
        ax = SCRATCH / f"Axioms_{n}.lean"
        ax.write_text(stxt.rstrip() + "\n\n#print axioms solution\n")
        r = lean(ax); R["lean"][f"Axioms_{n}"] = r
        m = re.search(r"'solution' depends on axioms: \[([^\]]*)\]", r["stdout"])
        axioms = sorted(a.strip() for a in m.group(1).split(",")) if m else None
        R["lean"][f"Axioms_{n}"]["axioms"] = axioms
        check(r["returncode"] == 0 and axioms is not None and set(axioms) <= ALLOWED_AXIOMS,
              f"axioms of solution in Sol_paley_{n} within {sorted(ALLOWED_AXIOMS)}", axioms)
        # statement match: normalized text
        tm = re.search(r"^theorem\s+(\S+)(.*?):=\s*by\s+sorry\s*$", ttxt, re.S | re.M)
        sm = re.search(r"^theorem\s+solution(.*?):=\s*by\s*$", stxt, re.S | re.M)
        norm = lambda s: re.sub(r"\s+", " ", s).strip()
        check(tm is not None and sm is not None and norm(tm.group(2)) == norm(sm.group(1)),
              f"binders+conclusion of Thm_paley_{n} equal those of solution (normalized text)",
              {"thm": norm(tm.group(2)) if tm else None, "sol": norm(sm.group(1)) if sm else None})
        # statement match: Lean defeq (`example : @thm = @solution := rfl` succeeds iff the types are defeq)
        mt = SCRATCH / f"Match_{n}.lean"
        mt.write_text(stxt.rstrip() + "\n\n" + tm.group(0) + f"\n\nexample : @{tm.group(1)} = @solution := rfl\n")
        r = lean(mt); R["lean"][f"Match_{n}"] = r
        check(r["returncode"] == 0 and "error" not in r["stdout"], f"Lean defeq match Thm_paley_{n} vs solution", r)
    for n in SOLVED + OPEN:
        thm = WS / "Theorems" / f"Thm_paley_{n}.lean"
        ttxt = thm.read_text()
        check(ttxt.rstrip().endswith(":= by sorry"), f"Thm_paley_{n}.lean ends with ':= by sorry'")
        r = lean(thm); R["lean"][f"Thm_{n}"] = r
        warn = r["stdout"].count("declaration uses `sorry`")
        check(r["returncode"] == 0 and warn == 1 and "error" not in r["stdout"],
              f"Thm_paley_{n}.lean elaborates with exactly one sorry warning", r)

# ----------------------------------------------------------------------------- proposal JSON
def proposal_checks():
    pj = WS / "proposals" / "paley-mission-proposal.json"
    if not check(pj.exists(), "proposal JSON exists"):
        return
    d = json.loads(pj.read_text())
    names = [i["theorem_name"] for i in d["items"]]
    check(d["main_item_theorem_name"] in names, "goal item present in items")
    check(d["item_order_theorem_names"] == names, "item_order matches items")
    check(all(m["theorem_name"] in names and m["theorem_name"] != d["main_item_theorem_name"] for m in d["milestones"]),
          "every milestone is a non-goal item")
    for it in d["items"]:
        ttxt = (WS / it["local_lean_file"]).read_text()
        tm = re.search(r"^theorem\s+(\S+).*?:=\s*by\s+sorry\s*$", ttxt, re.S | re.M)
        check(tm and tm.group(1) == it["theorem_name"] and tm.group(0).rstrip() == it["formal_statement"].rstrip(),
              f"proposal formal_statement of {it['theorem_name']} equals the Lean file", it["local_lean_file"])
        check(it["formal_statement"].rstrip().endswith(":= by sorry"), f"{it['theorem_name']} formal_statement ends in := by sorry")
        for line in it["preamble"].splitlines():
            check(line in ttxt, f"preamble line of {it['theorem_name']} present in Lean file", line)
        for key in ("theorem_title", "natural_language_statement", "source"):
            check(bool(it[key].strip()), f"{it['theorem_name']} has non-empty {key}")
        check(len(it["theorem_title"]) <= 200, f"{it['theorem_name']} title <= 200 chars")
    for m in d["milestones"]:
        check(len(m["milestone_title"]) <= 200 and m["milestone_description"].strip() != "", f"milestone {m['theorem_name']} well-formed")
    R["proposal"] = {"items": len(d["items"]), "milestones": len(d["milestones"]), "env": d["proposal"]["env"],
                     "description_words": len(d["proposal"]["description"].split())}

# ----------------------------------------------------------------------------- exact integer checks
def chi_table(p):
    sq = {(x * x) % p for x in range(1, p)}
    return [0] + [1 if x in sq else -1 for x in range(1, p)]

def integer_checks():
    rng = random.Random(20260905)
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53]
    wit = R["witnesses"]
    # (i) shift orthogonality, all p, all c != 0
    cnt = 0; wit["shift"] = []
    for p in primes:
        chi = chi_table(p)
        for c in range(1, p):
            s = sum(chi[x] * chi[(x + c) % p] for x in range(p))
            check(s == -1, "shift orthogonality", (p, c, s)); cnt += 1
            if c == 1: wit["shift"].append([p, c, s])
    R["integer_checks"]["shift_orthogonality"] = cnt
    # (ii) second moment: exhaustive over all B for p <= 13, random B otherwise
    cnt = 0; wit["second_moment"] = []
    for p in primes:
        chi = chi_table(p)
        subsets = ([[b for b in range(p) if mask >> b & 1] for mask in range(1 << p)] if p <= 13
                   else [[b for b in range(p) if rng.random() < rng.random()] for _ in range(300)])
        for B in subsets:
            m2 = sum(sum(chi[(x - b) % p] for b in B) ** 2 for x in range(p))
            check(m2 == len(B) * (p - len(B)), "second moment", (p, B, m2)); cnt += 1
        wit["second_moment"].append([p, len(subsets), "exhaustive" if p <= 13 else "random"])
    R["integer_checks"]["second_moment"] = cnt
    # (iii) Chung bound: exhaustive pairs for p <= 7, random pairs otherwise; record max ratio
    cnt = 0; wit["chung"] = []
    for p in primes:
        chi = chi_table(p); best = (0, 1, None)
        if p <= 7:
            pairs = [([a for a in range(p) if ma >> a & 1], [b for b in range(p) if mb >> b & 1])
                     for ma in range(1 << p) for mb in range(1 << p)]
        else:
            pairs = [([a for a in range(p) if rng.random() < rng.random()], [b for b in range(p) if rng.random() < rng.random()])
                     for _ in range(400)]
        for A, B in pairs:
            S = sum(chi[(a - b) % p] for a in A for b in B)
            rhs = len(A) * len(B) * (p - len(B))
            check(S * S <= rhs, "Chung bound", (p, A, B, S)); cnt += 1
            if rhs and S * S * best[1] > best[0] * rhs: best = (S * S, rhs, (A, B, S))
        wit["chung"].append([p, len(pairs), "max S^2/(|A||B|(p-|B|)) =", f"{best[0]}/{best[1]}", best[2]])
    R["integer_checks"]["chung_bound"] = cnt
    # (iv) interval => nonresidue, all N with 2N < p (and the sum identity when no nonresidue exists)
    cnt = 0; wit["interval"] = []
    for p in primes:
        chi = chi_table(p)
        for N in range(1, (p - 1) // 2 + 1):
            S = sum(chi[(a + b) % p] for a in range(1, N + 1) for b in range(1, N + 1))
            has_nr = any(chi[n % p] == -1 for n in range(2, 2 * N + 1))
            check((not S < N * N) or has_nr, "interval => nonresidue", (p, N, S)); cnt += 1
            check(has_nr or S == N * N, "no nonresidue => sum = N^2", (p, N, S)); cnt += 1
        wit["interval"].append([p, "least nonresidue", next(n for n in range(2, p) if chi[n] == -1)])
    R["integer_checks"]["interval_nonresidue"] = cnt
    # (v) Hanson-Petridis difference bound: exhaustive over subsets A for p <= 13 and every proper divisor d of p-1
    cnt = 0; wit["hp_difference"] = []
    for p in [p for p in primes if p <= 13]:
        for d in [d for d in range(1, p - 1) if (p - 1) % d == 0]:
            Zd = {z for z in range(1, p) if pow(z, d, p) == 1}
            mx = 0
            for mask in range(1 << p):
                A = [a for a in range(p) if mask >> a & 1]
                if all(((a - b) % p) in Zd for a, b in combinations(A, 2) for (a, b) in [(a, b), (b, a)]):
                    check(len(A) * (len(A) - 1) <= d, "HP difference bound", (p, d, A)); cnt += 1
                    mx = max(mx, len(A))
            wit["hp_difference"].append([p, d, "max |A| with A-A in Z_d u {0}", mx])
    R["integer_checks"]["hp_difference_bound"] = cnt
    # (vi) HP clique bound: exact clique number by brute force for p = 1 mod 4, p <= 53
    cnt = 0; wit["hp_clique"] = []
    for p in [p for p in primes if p % 4 == 1]:
        chi = chi_table(p)
        adj = [{y for y in range(p) if y != x and chi[(x - y) % p] == 1} for x in range(p)]
        best = 0
        def bk(Rset, P, X):
            nonlocal best
            if not P and not X: best = max(best, len(Rset)); return
            for v in list(P):
                bk(Rset | {v}, P & adj[v], X & adj[v]); P = P - {v}; X = X | {v}
        bk(set(), set(range(p)), set())
        bound = ((2 * p - 1) ** 0.5 + 1) / 2
        check(best <= bound, "HP clique bound", (p, best, bound)); cnt += 1
        wit["hp_clique"].append([p, "omega(G_p) =", best, "bound", round(bound, 3)])
    R["integer_checks"]["hp_clique_number"] = cnt
    R["integer_checks"]["not_numerically_testable"] = ["paley_two_set_conjecture (asymptotic, open)",
                                                       "paley_karatsuba_amplification (asymptotic with unspecified constant)"]

def main():
    if os.environ.get("SIGMA_LEAN_SKIP_LEAN") != "1":
        note("Lean runs"); lean_runs(); note("proposal checks"); proposal_checks()
    note("integer checks"); integer_checks()
    R["seconds"] = round(time.time() - T0, 1)
    R["status"] = "ALL CHECKS PASSED" if not R["failures"] else f"{len(R['failures'])} FAILURES"
    OUT.parent.mkdir(exist_ok=True)
    OUT.write_text(json.dumps(R, indent=2) + "\n")
    note(f"{R['status']}: {R['checks']} checks, wrote {OUT}")
    sys.exit(0 if not R["failures"] else 1)

if __name__ == "__main__":
    main()
