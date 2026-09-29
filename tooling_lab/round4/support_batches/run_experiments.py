#!/usr/bin/env python3
"""Exhaustive family operations, census equality and structured/generic cost audit."""
from collections import Counter
from itertools import combinations
from pathlib import Path
from hashlib import sha256
import json
import random
import sys
import time

sys.dont_write_bytecode = True
from symbolic_batches import *
import compressed_support as baseline

HERE = Path(__file__).resolve().parent
LAB = HERE.parents[1]


def family_audit():
    checks = 0
    for n in range(1, 9):
        for k in range(n + 1):
            z = FamilyZDD()
            root = z.choose(range(n), k)
            full = set(combinations(range(n), k))
            assert z.materialize(root) == full
            for allowed in range(1 << n):
                after = z.subtract_subsets(root, allowed)
                exact = {B for B in full if any(not (allowed & (1 << x)) for x in B)}
                assert z.materialize(after) == exact
                assert z.counts[after] == len(exact)
                checks += 1
                if exact:
                    assert z.first(after) == min(exact)
            rng = random.Random(n * 123 + k)
            current, exact = root, full
            for _ in range(20):
                mask = rng.randrange(1 << n)
                current = z.subtract_subsets(current, mask)
                exact = {B for B in exact if any(not (mask & (1 << x)) for x in B)}
                assert z.materialize(current) == exact
                checks += 1
    return checks


def compare(F, dom, k, s, u0, u1, name, materialize_limit=1000):
    old = baseline.run_census(F, dom, k, s, u0, u1, materialize_limit=materialize_limit)
    new = run_census(F, dom, k, s, u0, u1, materialize_limit=materialize_limit)
    assert node_set(old) == node_set(new)
    for key in ['finite_bad_scalar_count', 'whole_field_correlated', 'bad_scalars',
                'node_ledger_materialized', 'finite_track_certificates', 'all_field_track_certificates']:
        assert old[key] == new[key], (name, key)
    verify = verify_batches(F, new)
    out = {'name': name, 'full_node_and_track_certificate_equality': True,
           'baseline_seconds': old['elapsed_seconds'], 'batch_seconds': new['elapsed_seconds'],
           'baseline_interpolations': old['stats']['bases_processed'],
           'batch_interpolations': new['stats']['tracks_interpolated'], 'verification': verify,
           'result': new}
    print(name, 'bases', out['baseline_interpolations'], 'tracks', out['batch_interpolations'],
          'seconds', out['baseline_seconds'], out['batch_seconds'], flush=True)
    return out


def piecewise(F, dom, k, s, layout='aligned', parts=3):
    n = len(dom)
    plan = cover_plan(n, k, s)
    labels = [i % parts for i in range(n)]
    if layout == 'aligned':
        for j, block in enumerate(plan['blocks']):
            for i in block:
                labels[i] = j % parts
    rng = random.Random(491823 + n + k)
    polynomials = [[rng.randrange(F.q) for _ in range(k)] for _ in range(parts)]
    destination = [rng.randrange(F.q) for _ in range(k)]
    evaluate = lambda cs, x: F.sum(F.mul(c, F.power(x, j)) for j, c in enumerate(cs))
    u0 = [evaluate(polynomials[labels[i]], x) for i, x in enumerate(dom)]
    target = [evaluate(destination, x) for x in dom]
    u1 = [F.sub(c, a) for a, c in zip(u0, target)]
    return u0, u1, {'layout': layout, 'parts': parts, 'labels': labels,
                   'intercept_polynomials': polynomials, 'common_scalar': 1,
                   'target_polynomial': destination, 'target_codeword': target}


def main():
    start = time.perf_counter()
    output = {'schema': 'certified-support-batches/v1', 'family_subtraction_checks': family_audit(),
              'round1_and_extension_equalities': [], 'cost_experiments': [], 'large_structured': []}
    old = json.loads((LAB / 'proximity' / 'stack_results.json').read_text())
    for r in old['records']:
        p = r['p']
        dom = baseline.domain(r['n'], p)
        name = str(r['configuration']) + '/' + str(r['sample'])
        output['round1_and_extension_equalities'].append(compare(PrimeField(p), dom, r['k'], r['s'], r['u0'], r['u1'], name))
    extension = json.loads((LAB / 'proximity' / 'extension_results.json').read_text())
    for r in extension['experiments']:
        F = Field(r['p'], r['modulus'])
        output['round1_and_extension_equalities'].append(compare(F, r['domain'], r['k'], r['s'], r['u0'], r['u1'], r['name']))
    # Force the symbolic and materialized branches on the same extension case.
    F = Field(3, [1, 0, 1]); dom = list(range(9))
    for limit in [0, 1000]:
        output['round1_and_extension_equalities'].append(compare(F, dom, 2, 5,
            [F.add(1, x) for x in dom], [F.mul(2, x) for x in dom],
            'GF9-whole-field-limit' + str(limit), materialize_limit=limit))
    F = PrimeField(65537)
    for n, k, s, layout in [(32, 4, 12, 'generic'), (48, 6, 28, 'aligned'),
                            (48, 6, 28, 'interleaved'), (128, 8, 80, 'aligned')]:
        dom = list(range(n))
        if layout == 'generic':
            rng = random.Random(7843)
            u0 = [rng.randrange(F.q) for _ in dom]
            u1 = [rng.randrange(F.q) for _ in dom]
            planting = {'layout': 'generic', 'seed': 7843}
        else:
            u0, u1, planting = piecewise(F, dom, k, s, layout)
        r = compare(F, dom, k, s, u0, u1, f'{layout}-n{n}-k{k}-s{s}')
        r['planting'] = planting
        output['cost_experiments'].append(r)
    # These selected-base counts are too large for the baseline. Every skipped
    # basis is instead accounted for by exact disjoint same-track certificates.
    for n, k, s in [(256, 16, 160), (512, 32, 320), (1024, 64, 410)]:
        dom = baseline.domain(n, F.q)
        u0, u1, planting = piecewise(F, dom, k, s)
        result = run_census(F, dom, k, s, u0, u1)
        verification = verify_batches(F, result)
        assert any(x['scalar'] == 1 and x['codeword'] == planting['target_codeword'] for x in result['nodes'])
        output['large_structured'].append({'planting': planting, 'verification': verification,
            'baseline_not_run': True, 'result': result})
        print('large', n, k, s, 'covered', result['stats']['bases_certified'],
              'tracks', result['stats']['tracks_interpolated'], 'seconds', result['elapsed_seconds'], flush=True)
    output['elapsed_seconds'] = round(time.perf_counter() - start, 6)
    # Round1 inputs remain archived in their original ledger. Keep validation
    # summaries here; full certificates remain returned by run_census and are
    # exported for the new structured/generic experiments below.
    for r in output['round1_and_extension_equalities']:
        full = r.pop('result')
        r['result_summary'] = {key: full[key] for key in ['field', 'n', 'k', 's', 'stats',
            'whole_field_correlated', 'finite_bad_scalar_count', 'node_ledger_materialized']}
        stable = {key: full[key] for key in ['batch_certificates', 'nodes',
            'cover_certificate', 'all_field_track_certificates', 'finite_track_certificates']}
        r['deterministic_certificate_sha256'] = sha256(json.dumps(stable, sort_keys=True).encode()).hexdigest()
    output['source_sha256'] = {name: sha256((HERE / name).read_bytes()).hexdigest()
                               for name in ['symbolic_batches.py', 'run_experiments.py']}
    (HERE / 'results.json').write_text(json.dumps(output, indent=2) + '\n')
    print('saved results.json', output['elapsed_seconds'], 'seconds', flush=True)


if __name__ == '__main__':
    main()
