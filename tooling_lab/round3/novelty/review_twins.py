#!/usr/bin/env python3
"""Independent GF49, pair-type, and translation-reduced full histogram audit."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
from math import comb
from pathlib import Path
import json
import subprocess

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent/'pair_type_twins'


def mul(x, y):
    a, b, c, d = x % 7, x // 7, y % 7, y // 7
    return (a*c-b*d) % 7 + 7*((a*d+b*c) % 7)


def sub(x, y):
    return (x % 7-y % 7) % 7 + 7*((x//7-y//7) % 7)


def power(x, k):
    out = 1
    for _ in range(k):
        out = mul(out, x)
    return out


def direct_kernel(matrix, subset):
    out = 0
    for row in matrix:
        value = 1
        for i in subset:
            value *= row[i]
        out += value
    return out


def main():
    saved = json.loads((SOURCE/'results.json').read_text())
    pair_saved = json.loads((SOURCE/'pair_type_certificates.json').read_text())
    assert all((a*a+1) % 7 for a in range(7)), 'X^2+1 must be irreducible'
    for x in range(1, 49):
        a, b = x % 7, x // 7
        den_inv = pow((a*a+b*b) % 7, -1, 7)
        inverse = a*den_inv % 7 + 7*((-b*den_inv) % 7)
        assert mul(x, inverse) == 1
    g = saved['construction']['primitive_element']
    assert g == 9 and len({power(g, k) for k in range(48)}) == 48
    assert power(g, 48) == 1 and power(g, 16) != 1 and power(g, 24) != 1
    classes = [sorted(power(g, k) for k in range(48) if k % 4 == i) for i in range(4)]
    assert classes == saved['construction']['cyclotomic_classes']
    conn = [set(classes[0]+classes[2]), set(classes[0]+classes[1])]
    assert conn[0] == {mul(x, x) for x in range(1, 49)}
    matrices = [[[0 if x == y else (1 if sub(x, y) in connection else -1)
                  for y in range(49)] for x in range(49)] for connection in conn]
    assert matrices == saved['sign_matrices']
    patterns = []
    checked = 0
    for matrix, name in zip(matrices, ('paley', 'peisert')):
        counts = Counter()
        for x, row in enumerate(matrix):
            assert Counter(row) == {0: 1, -1: 24, 1: 24}
            for y in range(49):
                assert matrix[x][y] == matrix[y][x]
                histogram = Counter(zip(row, matrix[y]))
                entry = pair_saved[name][x*49+y]
                assert entry['rows'] == [x, y]
                assert entry['types'] == [[a, b, n] for (a, b), n in sorted(histogram.items())]
                inner = sum(a*b*n for (a, b), n in histogram.items())
                assert inner == entry['inner_product'] == (48 if x == y else -1)
                counts[tuple(sorted(histogram.items()))] += 1
                checked += 1
        patterns.append(counts)
    assert patterns[0] == patterns[1]
    removed = [[a, b] for a, b in combinations(range(49), 2) if matrices[0][a][b] == 1 and matrices[1][a][b] == -1]
    added = [[a, b] for a, b in combinations(range(49), 2) if matrices[0][a][b] == -1 and matrices[1][a][b] == 1]
    assert removed == saved['construction']['removed_edges']
    assert added == saved['construction']['added_edges']
    binary = HERE/'twin_origin_census'
    subprocess.run(['clang++', '-O3', '-std=c++17', str(HERE/'twin_origin_census.cpp'), '-o', str(binary)], check=True)
    payload = '\n'.join(' '.join(map(str, row)) for matrix in matrices for row in matrix)+'\n'
    output = subprocess.check_output([str(binary)], input=payload, text=True)
    hist = {'P': {}, 'Q': {}}
    joint = {}
    origin_count = None
    for line in output.splitlines():
        kind, *entries = line.split()
        nums = list(map(int, entries))
        if kind == 'COUNT':
            origin_count = nums[0]
            continue
        count = nums[-1]
        assert count*49 % 6 == 0
        full_count = count*49//6
        if kind in ('P', 'Q'):
            hist[kind][nums[0]] = full_count
        else:
            assert kind == 'J'
            joint[tuple(nums[:2])] = full_count
    assert origin_count == comb(48, 5)
    census = saved['census']
    assert sum(hist['P'].values()) == sum(hist['Q'].values()) == census['six_sets'] == comb(49, 6)
    witness_checks = 0
    for symbol, name, matrix in zip(('P', 'Q'), ('paley', 'peisert'), matrices):
        assert hist[symbol] == {int(k): v for k, v in census[name].items()}
        for value, subset in census['witnesses'][name].items():
            assert direct_kernel(matrix, subset) == int(value)
            witness_checks += 1
        for k, stored in census['moments'][name].items():
            expected = Fraction(sum(count*value**int(k) for value, count in hist[symbol].items()), comb(49, 6))
            assert expected == Fraction(*stored)
    assert joint == {(r['paley_value'], r['peisert_value']): r['count'] for r in census['joint']}
    w = census['largest_pointwise_difference']
    values = [direct_kernel(matrix, w['C']) for matrix in matrices]
    assert values == [w['paley_value'], w['peisert_value']]
    assert abs(values[0]-values[1]) == max(abs(a-b) for a, b in joint)
    assert hist['P'] != hist['Q']
    for name, h in saved['source_sha256'].items():
        assert sha256((SOURCE/name).read_bytes()).hexdigest() == h
    result = {'date': '2026-09-05', 'passed': True, 'field': 'F7[X]/(X^2+1)',
              'primitive_element': g, 'nonzero_inverses': 48, 'ordered_row_pairs_checked': checked,
              'matching_patterns_as_multiset': len(patterns[0]), 'removed_edges': len(removed), 'added_edges': len(added),
              'origin_subsets_enumerated': origin_count, 'full_subsets_certified_by_translation': comb(49, 6),
              'histogram_and_joint_bins_exactly_match': True, 'saved_single_witnesses_checked': witness_checks,
              'largest_pointwise_difference': {'C': w['C'], 'values': values},
              'nonisomorphism_certificate': 'Different distributions of an isomorphism-invariant six-subset kernel',
              'candidate_source_sha256': saved['source_sha256'],
              'review_sources_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                                        for name in ('review_twins.py', 'twin_origin_census.cpp')},
              'results_sha256': sha256((SOURCE/'results.json').read_bytes()).hexdigest()}
    (HERE/'twins_review.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
