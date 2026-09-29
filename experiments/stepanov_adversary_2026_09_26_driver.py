#!/usr/bin/env python3
"""Driver for the Stepanov-wave `adversary` searches (2026-09-26).

Runs the C++ helper experiments/stepanov_adversary_2026_09_26.cpp in at most
six concurrent processes, appends every result line to a JSONL log in the
scratch directory (so a killed run loses nothing), and merges the log into
results/stepanov_adversary_2026_09_26_search.json.  The search results are
NOT trusted: experiments/stepanov_adversary_2026_09_26.py re-checks every
stored witness by exact arithmetic.

Usage:  python3 stepanov_adversary_2026_09_26_driver.py <phase> [<phase> ...]
phases: exh exhq rhp rhpq rect rectq merge   (merge is always run at the end)
"""
import json
import math
import os
import random
import subprocess
import sys
import threading
import time
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
SRC = ROOT / "experiments" / "stepanov_adversary_2026_09_26.cpp"
SCR = Path(os.environ.get(
    "ADV_SCRATCH",
    "/private/tmp/claude-501/-Users-shawwalters-Desktop-paleygraph/"
    "1b733267-6e28-4860-bf0c-a4c3dcfdbfb5/scratchpad/adversary"))
BIN = SCR / "stepanov_adversary"
LOG = SCR / "runs.jsonl"
OUT = ROOT / "results" / "stepanov_adversary_2026_09_26_search.json"
MAXPROC = int(os.environ.get("ADV_MAXPROC", "5"))
LOCK = threading.Lock()


def compile_bin():
    if not BIN.exists() or BIN.stat().st_mtime < SRC.stat().st_mtime:
        subprocess.run(["clang++", "-O2", "-std=c++17", "-o", str(BIN), str(SRC)], check=True)


def done_jobs():
    s = set()
    if LOG.exists():
        for line in LOG.open():
            try:
                s.add(json.loads(line)["job"])
            except Exception:
                pass
    return s


def run_job(job, args, timeout):
    t0 = time.time()
    try:
        out = subprocess.run([str(BIN)] + [str(a) for a in args], capture_output=True, text=True,
                             timeout=timeout).stdout
    except subprocess.TimeoutExpired:
        out = ""
    recs = []
    for line in out.splitlines():
        line = line.strip()
        if line.startswith("{"):
            r = json.loads(line)
            r["job"] = job
            recs.append(r)
    with LOCK:
        with LOG.open("a") as f:
            for r in recs:
                f.write(json.dumps(r) + "\n")
    return job, len(recs), time.time() - t0


def run_all(jobs, label):
    have = done_jobs()
    todo = [j for j in jobs if j[0] not in have]
    est = sum(j[2] for j in todo)
    print(f"[{label}] {len(todo)} jobs to run ({len(jobs) - len(todo)} done), "
          f"nominal {est:.0f}s cpu", flush=True)
    t0 = time.time()
    with ThreadPoolExecutor(MAXPROC if label != "exh" else 1) as ex:
        futs = [ex.submit(run_job, j[0], j[1], j[3]) for j in todo]
        for i, f in enumerate(futs):
            job, n, dt = f.result()
            if i % 25 == 0 or i == len(futs) - 1:
                print(f"[{label}] {i + 1}/{len(futs)} {job} recs={n} {dt:.1f}s "
                      f"elapsed={time.time() - t0:.0f}s", flush=True)


# ------------------------------------------------------------------ job lists

def primes_upto(n):
    s = [True] * (n + 1)
    s[0] = s[1] = False
    for i in range(2, int(n ** 0.5) + 1):
        if s[i]:
            s[i * i::i] = [False] * len(s[i * i::i])
    return [i for i in range(n + 1) if s[i]]


def comb(n, k):
    return math.comb(n, k) if 0 <= k <= n else 0


def exh_mmax(p):
    if p <= 23:
        return min(p, 31)
    # largest mmax with 2 * sum_{k <= mmax-2} C(p-2, k) <= 3e8
    m = 2
    while m < 31 and 2 * sum(comb(p - 2, k) for k in range(m - 1)) <= 1.5e8:
        m += 1
    return m - 1


def jobs_exh():
    jobs = []
    for p in primes_upto(61):
        if p < 3:
            continue
        mm = exh_mmax(p)
        jobs.append((f"exh:{p}:{mm}", ["exh", p, mm, MAXPROC], 30, 3600))
    return jobs


def jobs_exhq():
    return [(f"exhq:{p}:{mm}", ["exhq", p, mm, MAXPROC], 30, 3600) for p, mm in [(3, 9), (5, 12), (7, 8)]]


RHP_PRIMES = [503, 1009, 2003, 4001, 8009, 8011, 16001]
RHP_M = [4, 5, 6, 7, 8, 10, 12, 16, 20, 24, 32, 48, 64, 96, 128]
FAMS_BASE = ["random", "rect", "gp", "coset", "clique", "nbhd"]
FAMS_EXTRA = ["b01", "m4", "interval", "squares", "b01core"]


def crux_seeds():
    """Sign-corrected A and B sides of the near-extremal rectangles of sigma-crux."""
    path = ROOT / "results" / "sigma_crux_2026_09_05_search.json"
    seeds = {}
    if not path.exists():
        return seeds
    d = json.load(path.open())
    for key in ("exhaustive", "anneal", "kstar"):
        for ent in d[key]:
            p = ent["p"]
            if p not in RHP_PRIMES:
                continue
            ws = ent.get("witnesses") or [v for v in ent.get("per_k", {}).values()]
            for w in ws:
                nu = next(x for x in range(2, p) if pow(x, (p - 1) // 2, p) == p - 1)
                sgn = 1 if w["S"] > 0 else -1
                for side in ("A", "B"):
                    S = [x % p for x in w[side]]
                    if sgn < 0:
                        S = [(nu * x) % p for x in S]
                    seeds.setdefault(p, []).append(sorted(set(S)))
    return seeds


def rhp_secs(p, m):
    """CPU seconds per search job."""
    return round(0.3 + p / 40000.0 + m / 400.0, 3)


def rhp_es(m):
    return [e for e in sorted({0, 1, 2, m // 16, m // 8, m // 4}) if e < m / 2]


def jobs_rhp():
    jobs = []
    cs = crux_seeds()
    for p in RHP_PRIMES:
        for m in RHP_M:
            if m * m > p:
                continue
            for e in rhp_es(m):
                if e == 0:
                    fams = ["random", "rect"]
                elif m < 8:
                    fams = ["random", "rect", "gp", "clique"]
                elif m >= 32:
                    fams = ["random", "rect", "gp", "clique"]
                else:
                    fams = list(FAMS_BASE)
                if m in (8, 16, 32) and e >= 1:
                    fams += FAMS_EXTRA
                secs = rhp_secs(p, m)
                for fam in fams:
                    jobs.append((f"rhp:fp:{p}:{m}:{e}:{fam}:1",
                                 ["rhp", "fp", p, m, e, fam, 1, secs], secs, 300 * secs + 600))
                # crux seeds (only at primes present in sigma-crux data)
                if p in cs and e >= 1:
                    good = [s for s in cs[p] if abs(len(s) - m) <= m // 2]
                    for i, s in enumerate(good[:2]):
                        jobs.append((f"rhp:fp:{p}:{m}:{e}:crux{i}:1",
                                     ["rhp", "fp", p, m, e, "given", 1, secs, "A=" + ",".join(map(str, s))],
                                     secs, 300 * secs + 600))
    prio = {f: i for i, f in enumerate(["random", "rect", "gp", "clique", "coset", "nbhd"] + FAMS_EXTRA)}
    jobs.sort(key=lambda j: (prio.get(j[0].split(":")[5], 20), int(j[0].split(":")[2]), int(j[0].split(":")[3])))
    return jobs


def jobs_rhp2():
    """Supplement: the common-neighbourhood family for large m (the unbalanced-biclique construction)."""
    jobs = []
    for p in RHP_PRIMES:
        for m in RHP_M:
            if m * m > p or m < 32:
                continue
            for e in rhp_es(m):
                secs = rhp_secs(p, m)
                for seed in (1, 2):
                    jobs.append((f"rhp:fp:{p}:{m}:{e}:nbhd:{seed}", ["rhp", "fp", p, m, e, "nbhd", seed, secs], secs,
                                 300 * secs + 600))
    return jobs


def jobs_rhpcal():
    """Search-quality calibration against the exact maxima of the exh phase."""
    jobs = []
    for p in (53, 59, 61):
        for m in (7, 8):
            for e in (1, 2, 3):
                for fam in ("random", "rect", "gp"):
                    jobs.append((f"rhp:fp:{p}:{m}:{e}:{fam}:cal", ["rhp", "fp", p, m, e, fam, 7, 0.3],
                                 0.3, 600))
    return jobs


def jobs_rhpq():
    """Calibration over F_{p^2}: does the same searcher find the subfield?"""
    jobs = []
    for p in [5, 7, 11, 13, 17, 19]:
        for m in sorted({p, (p + 1) // 2}):
            for e in (0, 1, 2):
                if e >= m / 2:
                    continue
                secs = round(0.3 + p * p / 1000.0, 2)
                jobs.append((f"rhp:fq:{p}:{m}:{e}:random:1",
                             ["rhp", "fq", p, m, e, "random", 1, secs], secs, 300 * secs + 600))
    return jobs


RECT_PRIMES = [1009, 2003, 4001, 8009, 16001]
RECT_RATIOS = [0.25, 0.5, 0.55, 0.6, 0.75, 1.0, 1.5]


def jobs_rect():
    jobs = []
    for p in RECT_PRIMES:
        for rho in RECT_RATIOS:
            shapes = []
            k = round(math.sqrt(rho * p))
            shapes.append(("bal", k, k))
            m = max(2, round(math.sqrt(p) / 2))
            shapes.append(("1:4", m, max(1, round(rho * p / m))))
            for tag, m, n in shapes:
                secs = round(0.3 + p / 12000.0, 2)
                for fam, seed in (("random", 1), ("random", 2), ("gp", 1)):
                    jobs.append((f"rect:fp:{p}:{rho}:{tag}:{m}:{n}:{fam}:{seed}",
                                 ["rect", "fp", p, m, n, fam, seed, secs], secs, 300 * secs + 600))
    return jobs


def jobs_rectq():
    jobs = []
    for p in [11, 13, 17, 19, 23, 29, 31]:
        q = p * p
        for rho in (0.5, 1.0):
            k = round(math.sqrt(rho * q))
            secs = round(0.3 + q / 2000.0, 2)
            for seed in (1, 2):
                jobs.append((f"rect:fq:{p}:{rho}:bal:{k}:{k}:random:{seed}",
                             ["rect", "fq", p, k, k, "random", seed, secs], secs, 300 * secs + 600))
    return jobs


TREND_PRIMES = [32003, 64007, 128021]
TREND_ME = [(4, 1), (6, 1), (7, 1), (8, 1), (8, 2), (12, 2), (16, 2), (16, 4)]


def jobs_trend():
    """Larger primes at small (m, e): does max R_e approach the Weil limit g(m,e)?"""
    jobs = []
    for p in TREND_PRIMES:
        for m, e in TREND_ME:
            for fam in ("random", "rect"):
                secs = round(0.8 + p / 64000.0, 2)
                jobs.append((f"rhp:fp:{p}:{m}:{e}:{fam}:1", ["rhp", "fp", p, m, e, fam, 1, secs], secs,
                             300 * secs + 600))
    return jobs


SHARE_PRIMES = [1009, 2003, 4001, 8009, 16001]


def share_thrs(m):
    """thr = 1 (positive part of the second moment) and thr = m - 2e - 1 for e = floor(m/8), floor(m/4)."""
    return sorted({1, m - 2 * (m // 8) - 1, m - 2 * (m // 4) - 1})


def jobs_share():
    jobs = []
    for p in SHARE_PRIMES:
        for m in sorted({8, 16, 32, math.isqrt(p)}):
            for thr in share_thrs(m):
                secs = round(0.3 + p / 20000.0, 2)
                for fam in ("random", "rect", "clique", "nbhd"):
                    jobs.append((f"share:fp:{p}:{m}:{thr}:{fam}:1", ["share", "fp", p, m, thr, fam, 1, secs], secs,
                                 300 * secs + 600))
    for p in NULL_PRIMES_SHARE:
        fk = "fr:" + str(null_path(p, 1))
        for m in sorted({8, 16, 32, math.isqrt(p)}):
            for thr in share_thrs(m):
                secs = round(0.3 + p / 20000.0, 2)
                for fam in ("random", "rect"):
                    jobs.append((f"share:fr1:{p}:{m}:{thr}:{fam}:1", ["share", fk, p, m, thr, fam, 1, secs], secs,
                                 300 * secs + 600))
    for p in (11, 13, 17, 19, 23):
        for thr in share_thrs(p):
            secs = round(0.3 + p * p / 2000.0, 2)
            jobs.append((f"share:fq:{p}:{p}:{thr}:random:1", ["share", "fq", p, p, thr, "random", 1, secs], secs,
                         300 * secs + 600))
    return jobs


NULL_PRIMES_SHARE = [2003, 8009]


# ------------------------------------------------------------------ null model (random +-1 function)

def null_table(p, seed):
    """Balanced random +-1 function on F_p^*, value 0 at 0; regenerated identically by the verifier."""
    rng = random.Random(1000003 * p + seed)
    d = (p - 1) // 2
    signs = [1] * d + [-1] * d
    rng.shuffle(signs)
    return [0] + signs


def null_path(p, seed):
    path = SCR / f"fr_{p}_{seed}.txt"
    if not path.exists():
        path.write_text(" ".join(map(str, null_table(p, seed))) + "\n")
    return path


NULL_PRIMES = [2003, 8009]


def jobs_null():
    jobs = []
    for p in NULL_PRIMES:
        fk = "fr:" + str(null_path(p, 1))
        for m in [7, 8, 12, 16, 24, 32, 64]:
            if m * m > p:
                continue
            for e in rhp_es(m):
                fams = ["random", "rect"] if e == 0 else (["random", "rect", "gp", "clique"] if m >= 32 or m < 8
                                                           else list(FAMS_BASE))
                secs = rhp_secs(p, m)
                for fam in fams:
                    jobs.append((f"rhp:fr1:{p}:{m}:{e}:{fam}:1", ["rhp", fk, p, m, e, fam, 1, secs], secs,
                                 300 * secs + 600))
        for rho in RECT_RATIOS:
            k = round(math.sqrt(rho * p))
            secs = round(0.3 + p / 12000.0, 2)
            for seed in (1, 2, 3):
                jobs.append((f"rect:fr1:{p}:{rho}:bal:{k}:{k}:random:{seed}", ["rect", fk, p, k, k, "random", seed, secs],
                             secs, 300 * secs + 600))
    return jobs


# ------------------------------------------------------------------ merge

def merge():
    recs = []
    if LOG.exists():
        for line in LOG.open():
            try:
                recs.append(json.loads(line))
            except Exception:
                pass
    out = {
        "description": ("Stepanov wave, worker adversary (2026-09-26). Witness data produced by "
                        "experiments/stepanov_adversary_2026_09_26.cpp via the driver. obj = m|B_e(A)| - r_e(A); "
                        "R_e = obj/((q-1)/2). F_{p^2} elements are indexed u + p*v for u + v*x, x^2 = nu, "
                        "nu the least non-residue mod p. exh/exhq records are exact maxima over all A of size m "
                        "(normalised to contain {0,1} or {0,nonsq}); rhp/rect records are search lower bounds. "
                        "Every witness is re-checked by experiments/stepanov_adversary_2026_09_26.py."),
        "exh": [r for r in recs if r.get("mode") == "exh"],
        "exh_done": [r for r in recs if r.get("mode") == "exh_done"],
        "exhq": [r for r in recs if r.get("mode") == "exhq"],
        "exhq_done": [r for r in recs if r.get("mode") == "exhq_done"],
        "rhp": [r for r in recs if r.get("mode") == "rhp" and r.get("field") == "fp" and not r["job"].endswith(":cal")],
        "rhpcal": [r for r in recs if r.get("mode") == "rhp" and r["job"].endswith(":cal")],
        "rhpq": [r for r in recs if r.get("mode") == "rhp" and r.get("field") == "fq"],
        "rect": [r for r in recs if r.get("mode") == "rect" and r.get("field") == "fp"],
        "rhpnull": [r for r in recs if r.get("mode") == "rhp" and str(r.get("field", "")).startswith("fr:")],
        "rectnull": [r for r in recs if r.get("mode") == "rect" and str(r.get("field", "")).startswith("fr:")],
        "share": [r for r in recs if r.get("mode") == "share"],
        "rectq": [r for r in recs if r.get("mode") == "rect" and r.get("field") == "fq"],
    }
    OUT.write_text(json.dumps(out, separators=(",", ":")))
    print("merged", {k: len(v) for k, v in out.items() if isinstance(v, list)}, "->", OUT)


def main():
    SCR.mkdir(parents=True, exist_ok=True)
    compile_bin()
    phases = sys.argv[1:] or ["merge"]
    table = {"exh": jobs_exh, "exhq": jobs_exhq, "rhp": jobs_rhp, "rhpq": jobs_rhpq, "rhpcal": jobs_rhpcal,
             "rect": jobs_rect, "rectq": jobs_rectq, "null": jobs_null,
             "trend": jobs_trend, "rhp2": jobs_rhp2, "share": jobs_share}
    for ph in phases:
        if ph in table:
            run_all(table[ph](), ph)
            merge()
    merge()


if __name__ == "__main__":
    main()
