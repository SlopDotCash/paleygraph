#!/usr/bin/env python3
"""Reproduce the marked-moment compiler's scale, null, and degree checks."""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import json
from pathlib import Path
import sys
import time

from conditional_moments import (conditional_moments, elementary,
                                 matrix_counts, paley_prime_matrix)

HERE = Path(__file__).resolve().parent
COUNTS = HERE.parent/'marked_counts'
sys.path.insert(0, str(HERE.parents[1]/'round2'/'observables'))
from exchange_grouped import grouped_spectrum


def empty_counts(q):
    return {'q': q, 'marks': [], 'cells': [{'pattern': [], 'size': q}],
            'relations': [{'left_pattern': [], 'right_pattern': [], 'sign': s,
                           'count': q if s == 0 else q*(q-1)//2} for s in (-1, 0, 1)]}


def exact_target(matrix, subset, degree):
    return sum(elementary((int(row[i]) for i in subset), degree)[degree] for row in matrix)


def main():
    started = time.perf_counter()
    scale = []
    for q, n in ((13, 6), (17, 8), (101, 6), (101, 8), (1297, 6), (1297, 8),
                 (65537, 16), (1000033, 31)):
        path = COUNTS/f'counts_p{q}_marks_0_1_2.json'
        data = json.loads(path.read_text())
        result = conditional_moments(data, n, audit_sign_dependence=True)
        baseline = conditional_moments(empty_counts(q), n)
        # This is an independent earlier overlap-kernel implementation,
        # with different polynomial variables and a Johnson spectral solve.
        if n <= q//2:
            previous = grouped_spectrum(q, n)
            assert baseline['mean'] == [previous['mean'].numerator, previous['mean'].denominator]
            assert Fraction(*baseline['second_moment']) == previous['distance_correlations'][0]
        result['global_second_moment'] = baseline['second_moment']
        relative = Fraction(*result['second_moment'])/Fraction(*baseline['second_moment'])-1
        result['relative_second_moment_change'] = [relative.numerator, relative.denominator]
        result['relative_second_moment_change_float'] = float(relative)
        result['source_counts'] = str(path.relative_to(HERE.parents[1]))
        result['source_counts_sha256'] = sha256(path.read_bytes()).hexdigest()
        scale.append(result)
    nulls = []
    for path in sorted(COUNTS.glob('counts_p101_marks_*.json')):
        data = json.loads(path.read_text())
        result = conditional_moments(data, 8, audit_sign_dependence=True)
        if len(data['marks']) <= 2:
            baseline = conditional_moments(empty_counts(101), 8)
            for field in ('mean', 'second_moment', 'variance'):
                assert result[field] == baseline[field]
        nulls.append(result)
    lookup = {tuple(x['marks']): x for x in nulls}
    for marks in ((1,2,3), (0,4,8), (0,2,4)):
        for field in ('mean', 'second_moment', 'variance'):
            assert lookup[marks][field] == lookup[0,1,2][field]
    # Degree-general checks avoid silently building an even-degree-only API.
    matrix = paley_prime_matrix(13)
    direct_checks = []
    for degree, n, marks in ((0,0,()), (1,3,(0,)), (2,4,(0,1)),
                             (3,3,(0,1)), (3,5,(0,1)), (4,5,(0,1,2)),
                             (5,6,(0,1,2)), (6,6,tuple(range(6))),
                             (6,8,tuple(range(8))), (6,13,(0,1,2))):
        outside = [x for x in range(13) if x not in marks]
        values = [exact_target(matrix, (*marks, *tail), degree)
                  for tail in combinations(outside, n-len(marks))]
        expected_mean = Fraction(sum(values), len(values))
        expected_second = Fraction(sum(x*x for x in values), len(values))
        actual = conditional_moments(matrix_counts(matrix, list(marks)), n, degree)
        assert Fraction(*actual['mean']) == expected_mean
        assert Fraction(*actual['second_moment']) == expected_second
        direct_checks.append({'degree': degree, 'n': n, 'marks': marks, 'sets': len(values),
                              'mean': actual['mean'], 'second_moment': actual['second_moment']})
    out = {'scope': 'Exact conditional finite averages; no worst-case or historical novelty claim',
           'scale_cases': scale, 'affine_and_low_mark_controls': nulls,
           'direct_degree_and_boundary_checks': direct_checks,
           'elapsed_seconds': time.perf_counter()-started,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ('conditional_moments.py', 'run_experiments.py')}}
    (HERE/'results.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps({'scale_cases': len(scale), 'nulls': len(nulls),
                      'degree_boundary_checks': len(direct_checks),
                      'scale_relative_second_moment_changes':
                      [[x['q'],x['n'],x['relative_second_moment_change_float']] for x in scale]}, indent=2))


if __name__ == '__main__':
    main()
