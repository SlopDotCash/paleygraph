#!/usr/bin/env python3
"""Complete finite quartic-window D6 census; no asymptotic extrapolation."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from datetime import datetime, timezone
from fractions import Fraction
from hashlib import sha256
import importlib.util
import json
from math import isqrt
from pathlib import Path
import time

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'six_distinct', ROOT/'experiments/parallel27_six_distinct_2026_09_05.py')
M = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(M)

def sieve(limit):
    flags = bytearray(b'\1') * (limit + 1)
    flags[:2] = b'\0\0'
    for d in range(2, isqrt(limit) + 1):
        if flags[d]:
            flags[d*d::d] = b'\0' * ((limit-d*d)//d+1)
    return flags

def run_case(args):
    p,n = args
    row = M.sparse(p,n)
    for k in ('seconds','balanced_samples','repeated_patterns'):
        del row[k]
    return row

def main():
    start = time.monotonic()
    flags = sieve(64**4)
    orders = [4,8,16,32,64]
    tasks = [(p,n) for n in orders for p in range(n**4//4+1,n**4+1,n)
             if flags[p]]
    old = json.loads((ROOT/'results/sigma_subgroup_2026_09_05.json').read_text())
    old_rows = old['compact_records']['by_order']
    old_lookup = {(r[0],int(n)):r for n,rows in old_rows.items() for r in rows}
    # The old energy window covers the full quartic window at these orders.
    for n in [16,32,64]:
        assert {p for p,k in tasks if k==n} == {
            r[0] for r in old_rows[str(n)] if n**4//4 <= r[0] <= n**4}
    rows=[]
    with ProcessPoolExecutor(max_workers=8) as pool:
        for row in pool.map(run_case,tasks,chunksize=64):
            prior=old_lookup.get((row['p'],row['n']))
            if prior:
                assert row['E2']==prior[1] and row['E3']==prior[2]
            rows.append(row)
            if len(rows)%4000==0:
                print(json.dumps({'completed':len(rows),'total':len(tasks),
                    'seconds':round(time.monotonic()-start,2)}),flush=True)
    summary={}
    for n in orders:
        part=[r for r in rows if r['n']==n]
        top=max(r['distinct_unbalanced_R6'] for r in part)
        summary[str(n)]={
            'quartic_primes':len(part),
            'max_D6':top,'max_D6_over_n_cubed':str(Fraction(top,n**3)),
            'maximizers':[r['p'] for r in part if r['distinct_unbalanced_R6']==top],
            'D6_orbit_histogram':dict(sorted(Counter(r['D6_scaling_orbits'] for r in part).items())),
            'D6_le_n_cubed_failures':[r['p'] for r in part if r['distinct_unbalanced_R6']>n**3],
            'E3_le_15_n_cubed_failures':[r['p'] for r in part if r['E3']>15*n**3],
            'nonintrinsic_E2_primes':[r['p'] for r in part if r['E2']>3*n*n-3*n],
        }
    files=['experiments/parallel28_complete_remainder_2026_09_05.py',
           'experiments/parallel27_six_distinct_2026_09_05.py',
           'results/sigma_subgroup_2026_09_05.json']
    result={'status':'passed','checked_at_utc':datetime.now(timezone.utc).isoformat(),
        'previous_turn_classification':'progress','orders':orders,'total_cases':len(rows),
        'prime_enumeration':'Independent Eratosthenes sieve through 64^4; exact trial-division primality rechecked in sparse counter.',
        'prior_energy_comparisons':sum((r['p'],r['n']) in old_lookup for r in rows),
        'per_order':summary,'cases':rows,'seconds':time.monotonic()-start,
        'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in files},
        'scope':'Complete finite census at five fixed orders. Does not imply a uniform bound at growing orders or prove Paley.'}
    (ROOT/'results/parallel28_complete_remainder_2026_09_05.json').write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'status':'passed','total_cases':len(rows),'per_order':summary,
                      'seconds':time.monotonic()-start}),flush=True)

if __name__=='__main__': main()
