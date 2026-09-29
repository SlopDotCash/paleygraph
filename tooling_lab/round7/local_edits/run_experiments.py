#!/usr/bin/env python3
"""Direct tiny oracles, information-loss search, then streamed critical-size inputs."""
from collections import defaultdict
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import json
from math import comb
from pathlib import Path
import random
from time import perf_counter
import numpy as np
import local_edits as tool

HERE = Path(__file__).resolve().parent


def direct_signs(q):
    residues = {x*x % q for x in range(1, q)}
    return [[0 if x == y else 1 if (x-y) % q in residues else -1 for y in range(q)] for x in range(q)]


def direct_target(matrix, selected, d):
    total = 0
    for row in matrix:
        coefficients = [1] + [0]*d
        for column in selected:
            for j in range(d, 0, -1):
                coefficients[j] += row[column]*coefficients[j-1]
        total += coefficients[d]
    return total


def direct_neighbours(matrix, selected, d, cache=None):
    selected = tuple(selected)
    q = len(matrix)
    outside = set(range(q))-set(selected)
    values = []
    for a in selected:
        for b in outside:
            neighbour = tuple(sorted((set(selected)-{a}) | {b}))
            values.append(cache[neighbour] if cache is not None else direct_target(matrix, neighbour, d))
    mean = F(sum(values), len(values))
    second = F(sum(x*x for x in values), len(values))
    return {'mean': mean, 'second': second, 'variance': second-mean*mean,
            'minimum': min(values), 'maximum': max(values), 'count': len(values)}


def check(record, oracle):
    for saved, direct in (('neighbour_mean', 'mean'), ('neighbour_second_moment', 'second'), ('neighbour_variance', 'variance')):
        assert F(*record[saved]) == oracle[direct], (saved, record, oracle)


def tiny_suite():
    rng = random.Random(73017)
    cases = targets = 0
    for q in (5, 13, 17, 29):
        matrix = direct_signs(q)
        for n in sorted({1, min(4,q-1), min(6,q-1), min(8,q-1), q-1}):
            for d in sorted({0, 1, min(3,n), min(6,n)}):
                for _ in range(2):
                    selected = sorted(rng.sample(range(q), n))
                    record = tool.from_matrix(matrix, selected, d, chunk_size=3)
                    oracle = direct_neighbours(matrix, selected, d)
                    check(record, oracle)
                    assert record['target'] == direct_target(matrix, selected, d)
                    prime = tool.from_prime(q, selected, d, chunk_size=7)
                    for key in ('target','internal_contraction','deletion_records','neighbour_mean','neighbour_variance'):
                        assert prime[key] == record[key]
                    cases += 1
                    targets += oracle['count']
    return {'cases': cases, 'direct_neighbour_targets': targets, 'matrix_prime_backend_equalities': cases}


def collision_search(q=17, n=6):
    matrix = direct_signs(q)
    cache = {c: direct_target(matrix, c, 6) for c in combinations(range(q), n)}
    fibres = defaultdict(list)
    examined = 0
    for tail in combinations(range(2, q), n-2):
        selected = (0, 1) + tail
        record = tool.from_prime(q, list(selected), 6)
        oracle = direct_neighbours(matrix, selected, 6, cache)
        check(record, oracle)
        signature = tuple(tuple(row) for row in record['signed_row_count_histogram'])
        # The full signed/zero row-count histogram determines every scalar row-sum
        # moment, T_d for every d, and the exact swap mean of each T_d.
        fibres[signature].append((record, oracle))
        examined += 1
    ambiguous = []
    for signature, rows in fibres.items():
        low = min(rows, key=lambda x: F(*x[0]['neighbour_variance']))
        high = max(rows, key=lambda x: F(*x[0]['neighbour_variance']))
        if low[0]['neighbour_variance'] != high[0]['neighbour_variance']:
            assert low[0]['target'] == high[0]['target'] and low[0]['neighbour_mean'] == high[0]['neighbour_mean']
            ambiguous.append((F(*high[0]['neighbour_variance'])-F(*low[0]['neighbour_variance']), low, high))
    ambiguous.sort(key=lambda row: row[0], reverse=True)
    witness = None
    if ambiguous:
        gap, low, high = ambiguous[0]
        witness = {'lower': low[0], 'upper': high[0], 'variance_gap': tool.frac(gap),
                   'direct_neighbour_ranges': [[row[1]['minimum'],row[1]['maximum']] for row in (low, high)]}
    return {'q':q, 'n':n, 'degree':6, 'normalization':'All sets containing 0 and 1, not all affine equivalence classes.',
            'all_set_targets_cached':len(cache), 'normalized_inputs':examined, 'histogram_fibres':len(fibres),
            'ambiguous_variance_fibres':len(ambiguous), 'direct_neighbour_targets':examined*n*(q-n), 'witness':witness}


def scale_cases():
    rng = random.Random(73018)
    rows = []
    for q,n in ((1297,6),(65537,16),(1000033,31),(6700417,50)):
        families = [('arithmetic_progression',list(range(n))),('seeded_uniform',sorted(rng.sample(range(q),n)))]
        for label, selected in families:
            row = tool.from_prime(q, selected)
            row['family'] = label
            rows.append(row)
            print(json.dumps({'q':q,'n':n,'family':label,'target':row['target'],
                              'variance':float(F(*row['neighbour_variance'])),'seconds':row['work']['seconds']}), flush=True)
    return rows


def main():
    start = perf_counter()
    tiny = tiny_suite()
    print('Tiny direct oracles: '+str(tiny), flush=True)
    collision = collision_search()
    print('Actual histogram ambiguity: '+str({k:v for k,v in collision.items() if k!='witness'}), flush=True)
    # Preserve the completed toy evidence before starting the larger streaming runs.
    bindings = {name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('local_edits.py','run_experiments.py')}
    toy = {'tiny':tiny,'collision_search':collision,'source_sha256':bindings}
    (HERE/'toy_results.json').write_text(json.dumps(toy,indent=2)+'\n')
    rows = scale_cases()
    out = {'scope':'Exact local one-swap moments of actual selected sets; neither a uniform bound nor a new general concentration theorem.',
           'tiny':tiny,'collision_search':collision,'scale_cases':rows,'source_sha256':bindings,'elapsed_seconds':perf_counter()-start}
    (HERE/'results.json').write_text(json.dumps(out,indent=2)+'\n')
    print('Completed local-edit experiments.',flush=True)


if __name__ == '__main__':
    main()
