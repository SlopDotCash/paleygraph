#!/usr/bin/env python3
"""Data sweep for research/stepanov-balanced-2026-09-29.md (task 2).

P_K(p) := max{|A||B| : A + B subset of Q u {0} in F_p, min(|A|,|B|) >= K},  p = 1 mod 4.

Method (exact where flagged):
  * M_m(p) = max_{|A|=m} |B(A)| from results/sigma_biclique_2026_09_05_search.json (exact values,
    computed by the `biclique` worker) for m <= 8 (p <= 769) or m <= 6 (p <= 997), and from our own
    C helper experiments/stepanov_balanced_2026_09_29.c (mode `profile`) otherwise.  For p <= 997 the
    helper re-derives M_3..M_6 and they are compared with the table.
  * T_K = max over K <= m <= mk (M_m >= K) of m*M_m is attained (witness A, B = B(A)).
  * Any biclique whose smaller side has size m <= mk has product <= m*M_m <= T_K.  Bicliques whose
    smaller side has size >= mk+1 are searched exhaustively by the helper (mode `cert`), with the
    smaller side at most bmax: bmax = b(p) (exact balanced number, table, p <= 509) or the
    Hanson--Petridis bound floor((1+sqrt(2p-1))/2) otherwise.  If the search finishes, P_K is exact.
Writes results/stepanov_balanced_2026_09_29_data.json.  Run time: about an hour with 3 processes.
"""
import json, os, subprocess, sys, math, time
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.join(HERE, '..')
SRC = os.path.join(HERE, 'stepanov_balanced_2026_09_29.c')
OUT = os.path.join(ROOT, 'results', 'stepanov_balanced_2026_09_29_data.json')
TABLE = os.path.join(ROOT, 'results', 'sigma_biclique_2026_09_05_search.json')
import tempfile
BIN = os.environ.get('SB_BIN', os.path.join(os.environ.get('SB_TMP') or tempfile.mkdtemp(), 'sb_balanced'))
NPROC = int(os.environ.get('SB_NPROC', '3'))
TLIM = float(os.environ.get('SB_TLIM', '120'))


def compile_helper():
    subprocess.check_call(['clang', '-O3', '-o', BIN, SRC])


def run(args, timeout=None):
    out = subprocess.run([BIN] + [str(a) for a in args], capture_output=True, text=True, timeout=timeout)
    return json.loads(out.stdout.strip().splitlines()[-1])


def isprime(n):
    if n < 2:
        return False
    i = 2
    while i * i <= n:
        if n % i == 0:
            return False
        i += 1
    return True


def hp_bmax(p):
    # Hanson--Petridis: a k x k biclique has k^2 <= (p-1)/2 + k, i.e. k <= (1+sqrt(2p-1))/2.
    k = 1
    while (k + 1) * (k + 1) <= (p - 1) // 2 + (k + 1):
        k += 1
    return k


def one_prime(p, table, extra_profiles):
    e = table.get(str(p))
    rec = {'p': p, 'M': {}, 'M_src': {}, 'witness': {}}
    if e is not None:
        for k, v in e['M'].items():
            if v['exact']:
                rec['M'][int(k)] = v['M']; rec['M_src'][int(k)] = 'table'; rec['witness'][int(k)] = v['A']
    # our own profiles (independent re-derivation for small k, extension for large p)
    for k in extra_profiles:
        r = run([p, 'profile', k])
        if k in rec['M'] and rec['M'][k] != r['M']:
            rec.setdefault('mismatch', []).append((k, rec['M'][k], r['M']))
        if k not in rec['M']:
            rec['M'][k] = r['M']; rec['M_src'][k] = 'helper'; rec['witness'][k] = r['A']
        else:
            rec['M_src'][k] += '+helper'
    mk = max(rec['M'])
    b_exact = None
    if e is not None and e.get('balanced', {}).get('exact'):
        b_exact = e['balanced']['b']
    bmax = b_exact if b_exact is not None else hp_bmax(p)
    rec['mk'] = mk; rec['bmax'] = bmax; rec['bmax_src'] = 'exact b(p)' if b_exact is not None else 'HP bound'
    caps = [rec['M'].get(m, p) for m in range(1, mk + 1)]
    caps[0] = (p + 1) // 2
    rec['P'] = {}
    for K in (3, 4, 5, 6):
        cands = [(m * rec['M'][m], m) for m in range(K, mk + 1) if m in rec['M'] and rec['M'][m] >= K]
        if not cands:
            continue
        T, am = max(cands)
        entry = {'T_profile': T, 'argm': am, 'A': rec['witness'][am], 'n': rec['M'][am]}
        if bmax <= mk:
            entry.update({'P': T, 'exact': True, 'cert': 'bmax <= mk'})
        elif p > 509 and K >= 5:
            entry.update({'P': T, 'exact': False, 'cert': 'not attempted (lower bound)'})
        else:
            r = run([p, 'cert', mk + 1, T, bmax, TLIM if p <= 509 else TLIM / 4] + caps)
            entry.update({'P': r['P'], 'exact': r['exact'], 'cert_nodes': r['nodes'], 'cert_time': r['time']})
            if r['improved']:
                entry['A_large'] = r['A']
        entry['ratio'] = entry['P'] / p
        rec['P'][K] = entry
    return rec


def main():
    compile_helper()
    table = json.load(open(TABLE))['table']
    jobs = []
    for p in range(13, 1000):
        if isprime(p) and p % 4 == 1:
            jobs.append((p, [3, 4, 5, 6]))
    big = [p for p in range(1000, 3001) if isprime(p) and p % 4 == 1]
    for i, p in enumerate(big):
        if i % 5 == 0 or p == big[-1]:
            jobs.append((p, [3, 4, 5]))
    res = {}
    t0 = time.time()
    with ThreadPoolExecutor(NPROC) as ex:
        futs = {ex.submit(one_prime, p, table, ks): p for p, ks in jobs}
        for f in futs:
            p = futs[f]
            try:
                res[p] = f.result()
            except Exception as exc:  # record, do not hide
                res[p] = {'p': p, 'error': repr(exc)}
            json.dump({'description': __doc__, 'elapsed': time.time() - t0,
                       'rows': {str(k): v for k, v in sorted(res.items())}}, open(OUT, 'w'), indent=0)
    print('done', len(res), time.time() - t0)


if __name__ == '__main__':
    main()
