#!/usr/bin/env python3
"""sigma / p2m verifier (2026-09-05): re-read, by id, what the sigma-p2m pass created on prove2me.

Given a prove2me API key (``~/prove2me_workspace/credentials.json`` = ``{"api_key": "..."}``, or the
path in ``$P2M_CREDENTIALS``), this script

  1. exchanges the key for a 1-hour bearer token at POST /api/v1/agent/refresh (the ONLY place the key
     is sent; the token is kept in memory and sent ONLY to https://prove2.me/api/v1);
  2. re-reads the four proof submissions by ``submission_id`` (GET /verify?submission_id=...) and records
     each verdict verbatim;
  3. re-reads the eight private theorems by ``theorem_id`` and the private mission proposal by id
     (status, visibility, items, milestones), and checks they are unchanged: private / Draft / never
     launched;
  4. independently re-checks, by exact integer arithmetic at small primes, the four identities whose
     proofs were ACCEPTED (so the verdicts are not the only evidence);
  5. writes results/sigma_p2m_2026_09_05.json with the counts of checks performed and every witness.

It never prints or writes a secret. Standard library only. Read-only: it performs no POST except the
token exchange. Runs in well under a minute.
"""
import json, os, sys, time, random, urllib.request, urllib.error, urllib.parse

BASE = "https://prove2.me/api/v1"
HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, "results", "sigma_p2m_2026_09_05.json")
CRED = os.environ.get("P2M_CREDENTIALS", os.path.expanduser("~/prove2me_workspace/credentials.json"))

ENV = "c5ea00351c28e24afc9f0f84379aa41082b1188f"
PROPOSAL_ID = "b5e7121a-eb3c-48f1-a495-29fffd24e71d"
THEOREMS = {  # theorem_name -> (theorem_id, expected status after this pass)
    "paley_two_set_conjecture": ("973dac7a-fac5-4e8e-961a-936f59818faa", "Open"),
    "paley_shift_orthogonality": ("088f1fc8-5c05-4c9a-88ba-06444c4e3bef", "Proved"),
    "paley_second_moment": ("3b819822-953a-4aec-8409-5f1197fa26a8", "Proved"),
    "paley_chung_bound": ("63bb5a7e-af3a-405c-ab99-53cb381ac121", "Proved"),
    "paley_interval_nonresidue": ("3168650d-eaa9-4a62-b26e-2b51a554f39c", "Proved"),
    "paley_karatsuba_amplification": ("ca0607b1-a20c-4460-a5ab-61d9de0e8ad3", "Open"),
    "paley_hanson_petridis_difference_bound": ("58b3552d-bb5c-42f6-8595-ecf402103ee1", "Open"),
    "paley_hanson_petridis_clique_number": ("cb279d48-6789-4966-9b30-fb95f7de7bc6", "Open"),
}
SUBMISSIONS = {  # theorem_name -> submission_id (all expected ACCEPTED)
    "paley_shift_orthogonality": "965f3276-7f70-4916-9865-b8ca6c5ada62",
    "paley_second_moment": "6a5656c4-7ede-40bc-a230-2b366fba29aa",
    "paley_chung_bound": "35802501-9a0a-4c58-84cc-c3bc2d6a8907",
    "paley_interval_nonresidue": "f31c31b7-907c-477a-b77f-6c00a2e986ad",
}
MILESTONE_NAMES = ["paley_shift_orthogonality", "paley_second_moment", "paley_chung_bound",
                   "paley_interval_nonresidue", "paley_karatsuba_amplification",
                   "paley_hanson_petridis_difference_bound", "paley_hanson_petridis_clique_number"]

res = {"date": "2026-09-05", "direction": "p2m", "checks": 0, "failures": 0, "witnesses": [],
       "requests": [], "verdicts": {}, "theorems": {}, "proposal": {}, "exact_arithmetic": {}}

def check(name, ok, witness=None):
    res["checks"] += 1
    if not ok:
        res["failures"] += 1
        res["witnesses"].append({"check": name, "witness": witness})

# ----------------------------------------------------------------------------- HTTP (secret-safe)
def _raw(method, url, headers, body=None):
    req = urllib.request.Request(url, data=body, method=method, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=120) as r:
            return r.status, r.read()
    except urllib.error.HTTPError as e:
        return e.code, e.read()

def get_token():
    key = json.load(open(CRED))["api_key"]
    st, raw = _raw("POST", BASE + "/agent/refresh", {"Content-Type": "application/json"},
                   json.dumps({"api_key": key}).encode())
    res["requests"].append({"method": "POST", "path": "/agent/refresh", "status": st})
    if st != 200:
        raise SystemExit("token exchange failed: HTTP %d" % st)
    d = json.loads(raw)
    res["platform_version"] = d.get("version")
    return d["access_token"]

TOKEN = None
def get(path, params=None):
    global TOKEN
    if TOKEN is None:
        TOKEN = get_token()
    url = BASE + path + (("?" + urllib.parse.urlencode(params)) if params else "")
    st, raw = _raw("GET", url, {"Authorization": "Bearer " + TOKEN})
    try:
        d = json.loads(raw.decode("utf-8", "replace") or "{}")
    except json.JSONDecodeError:
        d = {"_raw": raw[:300].decode("utf-8", "replace")}
    res["requests"].append({"method": "GET", "path": path + (("?" + urllib.parse.urlencode(params)) if params else ""),
                            "status": st})
    return st, d

# ----------------------------------------------------------------------------- server re-reads
def server_checks():
    t0 = time.time()
    # (a) verdicts by submission id
    for name, sid in SUBMISSIONS.items():
        st, d = get("/verify", {"submission_id": sid})
        v = {"http": st, "status": d.get("status"), "error_message": d.get("error_message"),
             "theorem_id": d.get("theorem_id"), "submission_id": sid}
        res["verdicts"][name] = v
        check("verdict %s is ACCEPTED" % name, st == 200 and v["status"] == "ACCEPTED", v)
        check("verdict %s targets the recorded theorem" % name, v["theorem_id"] == THEOREMS[name][0], v)
    # (b) theorems by id
    for name, (tid, expected) in THEOREMS.items():
        st, d = get("/theorems/%s" % tid)
        t = {"http": st, "theorem_name": d.get("theorem_name"), "status": d.get("status"),
             "mathlib_rev": d.get("mathlib_rev"), "deprecated_at": d.get("deprecated_at"),
             "visibility": d.get("visibility"), "is_private": d.get("is_private"), "private": d.get("private"),
             "created_by_username": d.get("created_by_username")}
        res["theorems"][name] = t
        check("theorem %s readable by id" % name, st == 200, t)
        check("theorem %s has the recorded name" % name, t["theorem_name"] == name, t)
        check("theorem %s status == %s" % (name, expected), t["status"] == expected, t)
        check("theorem %s in pinned env" % name, t["mathlib_rev"] == ENV, t)
        check("theorem %s not deprecated" % name, t["deprecated_at"] in (None, ""), t)
        vis = t["visibility"] if t["visibility"] is not None else (
            "private" if (t["is_private"] or t["private"]) else None)
        if vis is not None:
            check("theorem %s is private" % name, vis == "private", t)
        else:
            res.setdefault("notes", []).append("GET /theorems/:id exposes no visibility field for %s; privacy was "
                                               "confirmed at publish time (publish job visibility=private)" % name)
    # (c) proposal by id
    st, d = get("/mission-proposals/%s" % PROPOSAL_ID)
    items = d.get("items", [])
    p = {"http": st, "name": d.get("name"), "status": d.get("status"), "visibility": d.get("visibility"),
         "env": d.get("env"), "mission_type": d.get("mission_type"), "main_item_id": d.get("main_item_id"),
         "n_items": len(items), "fields": [f.get("slug") for f in d.get("fields", [])],
         "items": [{"id": i.get("id"), "kind": i.get("kind"), "theorem_id": i.get("theorem_id"),
                    "theorem_name": i.get("theorem_name")} for i in items]}
    res["proposal"] = p
    check("proposal readable", st == 200, p)
    check("proposal status is Draft (never launched)", p["status"] == "Draft", p)
    check("proposal visibility is private", p["visibility"] == "private", p)
    check("proposal env is pinned", p["env"] == ENV, p)
    check("proposal has 8 items", p["n_items"] == 8, p)
    by_tid = {i["theorem_id"]: i for i in p["items"]}
    for name, (tid, _) in THEOREMS.items():
        check("proposal references %s" % name, tid in by_tid and by_tid[tid]["kind"] == "reference", p["items"])
    goal = next((i for i in p["items"] if i["id"] == p["main_item_id"]), None)
    check("proposal goal is paley_two_set_conjecture", goal is not None and goal["theorem_id"] == THEOREMS["paley_two_set_conjecture"][0], goal)
    st, d = get("/mission-proposals/%s/milestones" % PROPOSAL_ID)
    ms = d.get("milestones", d if isinstance(d, list) else [])
    res["proposal"]["milestones"] = [{"sort_order": m.get("sort_order"), "theorem_id": m.get("theorem_id"),
                                      "item_id": m.get("item_id"), "milestone_title": m.get("milestone_title")} for m in ms]
    check("proposal has 7 milestones", st == 200 and len(ms) == 7, res["proposal"]["milestones"])
    ms_tids = [m.get("theorem_id") for m in sorted(ms, key=lambda m: m.get("sort_order", 0))]
    check("milestones are the 7 non-goal items in order", ms_tids == [THEOREMS[n][0] for n in MILESTONE_NAMES], ms_tids)
    # (d) the mission catalogue must NOT contain this mission (never launched)
    st, d = get("/missions", {"limit": 200, "offset": 0})
    names = [m.get("name") for m in d.get("missions", [])]
    check("no live mission named like the proposal", st == 200 and not any("paley" in (n or "").lower() for n in names), names[:5])
    res["server_wall_seconds"] = round(time.time() - t0, 2)

# ----------------------------------------------------------------------------- exact arithmetic
def chi_table(p):
    sq = {(x * x) % p for x in range(1, p)}
    return [0] + [1 if x in sq else -1 for x in range(1, p)]

def exact_checks():
    rng = random.Random(20260905)
    primes = [3, 5, 7, 11, 13, 17, 19, 23, 29, 31]
    counts = {"shift_orthogonality": 0, "second_moment": 0, "chung_bound": 0, "interval_nonresidue": 0}
    for p in primes:
        chi = chi_table(p)
        for c in range(1, p):                                    # (i) sum_x chi(x)chi(x+c) = -1
            s = sum(chi[x] * chi[(x + c) % p] for x in range(p))
            check("shift p=%d c=%d" % (p, c), s == -1, {"p": p, "c": c, "sum": s}); counts["shift_orthogonality"] += 1
        subsets = [[x for x in range(p) if (m >> x) & 1] for m in range(1 << p)] if p <= 11 else \
                  [sorted(rng.sample(range(p), rng.randint(0, p))) for _ in range(120)]
        for B in subsets:                                        # (ii) sum_x F_B(x)^2 = |B|(p-|B|)
            F = [sum(chi[(x - b) % p] for b in B) for x in range(p)]
            m2 = sum(f * f for f in F)
            check("second moment p=%d" % p, m2 == len(B) * (p - len(B)), {"p": p, "B": B, "m2": m2}); counts["second_moment"] += 1
            A = sorted(rng.sample(range(p), rng.randint(0, p)))  # (iii) S^2 <= |A||B|(p-|B|)
            S = sum(F[a] for a in A)
            check("chung p=%d" % p, S * S <= len(A) * len(B) * (p - len(B)), {"p": p, "A": A, "B": B, "S": S}); counts["chung_bound"] += 1
        for N in range(1, (p - 1) // 2 + 1):                     # (iv) interval => non-residue in [2,2N]
            s = sum(chi[(a + b) % p] for a in range(1, N + 1) for b in range(1, N + 1))
            has_nr = any(chi[n] == -1 for n in range(2, 2 * N + 1))
            check("interval p=%d N=%d" % (p, N), (s >= N * N) or has_nr, {"p": p, "N": N, "sum": s, "has_nonresidue": has_nr})
            check("interval converse p=%d N=%d" % (p, N), has_nr or s == N * N, {"p": p, "N": N, "sum": s}); counts["interval_nonresidue"] += 2
    res["exact_arithmetic"] = {"primes": primes, "counts": counts}

def main():
    t0 = time.time()
    exact_checks()
    try:
        server_checks()
    except Exception as e:  # network failure must still leave a results file
        res["server_error"] = repr(e)
        check("server re-read completed", False, repr(e))
    res["wall_seconds"] = round(time.time() - t0, 2)
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(res, open(OUT, "w"), indent=1, ensure_ascii=False)
    print("checks=%d failures=%d wall=%.1fs -> %s" % (res["checks"], res["failures"], res["wall_seconds"], OUT))
    for name, v in res["verdicts"].items():
        print("  %s: %s (submission %s)" % (name, v["status"], v["submission_id"]))
    print("  proposal %s: status=%s visibility=%s items=%s milestones=%d" % (
        PROPOSAL_ID, res["proposal"].get("status"), res["proposal"].get("visibility"),
        res["proposal"].get("n_items"), len(res["proposal"].get("milestones", []))))
    for w in res["witnesses"]:
        print("  WITNESS:", json.dumps(w)[:300])

if __name__ == "__main__":
    main()
