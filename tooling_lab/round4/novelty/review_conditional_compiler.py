#!/usr/bin/env python3
"""Independent direct-enumeration comparisons for the marked-moment compiler."""
import sys
sys.dont_write_bytecode = True
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from importlib.util import module_from_spec, spec_from_file_location
from itertools import combinations
from math import comb, prod
from pathlib import Path
import json
import time

HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'marked_moments' / 'conditional_moments.py'


def load(path):
    spec = spec_from_file_location(path.stem, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def matrices():
    out = {}
    for p in (5, 13, 17):
        chi = [0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1)
               for x in range(p)]
        out[f'Paley{p}'] = [[chi[(y-x) % p] for x in range(p)] for y in range(p)]
    oracle = load(HERE / 'direct_conditional_oracle.py')
    powers = [1]
    for _ in range(47):
        powers.append(oracle.mul(powers[-1], 9))
    assert len(set(powers)) == 48 and oracle.mul(powers[-1], 9) == 1
    for name, classes in (('Paley49', (0, 2)), ('Peisert49', (0, 1))):
        connections = {powers[j] for j in range(48) if j % 4 in classes}
        out[name] = [[0 if x == y else (1 if oracle.sub(x, y) in connections else -1)
                      for y in range(49)] for x in range(49)]
    return out


def direct_generic(S, marks, n, degree):
    rest = sorted(set(range(len(S))) - set(marks))
    values = []
    for outside in combinations(rest, n-len(marks)):
        C = tuple(marks)+outside
        values.append(sum(prod(row[x] for x in D) for row in S for D in combinations(C, degree)))
    mean = Fraction(sum(values), len(values))
    second = Fraction(sum(x*x for x in values), len(values))
    return {'mean': [mean.numerator, mean.denominator],
            'second_moment': [second.numerator, second.denominator],
            'variance': [(second-mean*mean).numerator, (second-mean*mean).denominator]}


def main():
    started = time.perf_counter()
    source_hash = sha256(SOURCE.read_bytes()).hexdigest()
    compiler = load(SOURCE)
    all_matrices = matrices()
    oracle_path = HERE / 'direct_conditional_oracle.json'
    oracle = json.loads(oracle_path.read_text())
    for filename, digest in oracle['source_sha256'].items():
        assert sha256((HERE / filename).read_bytes()).hexdigest() == digest
    cache = {}
    cases = []
    for case in oracle['cases']:
        S = all_matrices[case['graph']]
        compared = 0
        full_records = case['records'] + ([case['full_mark_boundary']] if 'full_mark_boundary' in case else [])
        for expected in full_records:
            marks = expected['marks']
            key = (case['graph'], tuple(marks))
            if key not in cache:
                cache[key] = compiler.matrix_counts(S, marks, validate_matrix=True)
            actual = compiler.conditional_moments(cache[key], case['n'], 6)
            assert expected['count'] == comb(case['q']-len(marks), case['n']-len(marks))
            for field in ('mean', 'second_moment', 'variance'):
                assert actual[field] == expected[field], (case['graph'], case['n'], marks, field, actual[field], expected[field])
            compared += 1
        cases.append({'graph': case['graph'], 'n': case['n'], 'compared': compared})
        print(json.dumps(cases[-1]), flush=True)
    generic = 0
    S = all_matrices['Paley5']
    for n in range(6):
        marksets = {M for m in range(min(n, 2)+1) for M in combinations(range(5), m)}
        marksets.add(tuple(range(n)))
        for marks in sorted(marksets):
            counts = compiler.matrix_counts(S, list(marks))
            for degree in range(n+1):
                actual = compiler.conditional_moments(counts, n, degree)
                expected = direct_generic(S, marks, n, degree)
                assert all(actual[key] == expected[key] for key in expected), (n, marks, degree, actual, expected)
                generic += 1
    reversed_marks = 0
    for name, S in all_matrices.items():
        if len(S) < 13:
            continue
        for n in (6, 7):
            a = compiler.conditional_moments(compiler.matrix_counts(S, [0, 1, 2]), n)
            b = compiler.conditional_moments(compiler.matrix_counts(S, [2, 0, 1]), n)
            assert all(a[k] == b[k] for k in ('mean', 'second_moment', 'variance'))
            reversed_marks += 1
    assert sha256(SOURCE.read_bytes()).hexdigest() == source_hash, 'Candidate source changed during review'
    result = {'status': 'passed', 'date': '2026-09-05', 'cases': cases,
              'saved_oracle_records': sum(len(c['records']) for c in oracle['cases']),
              'full_mark_boundary_records': sum('full_mark_boundary' in c for c in oracle['cases']),
              'additional_generic_degree_checks': generic, 'mark_permutation_checks': reversed_marks,
              'compiler_sha256': source_hash, 'oracle_sha256': sha256(oracle_path.read_bytes()).hexdigest(),
              'reviewer_sha256': sha256(Path(__file__).read_bytes()).hexdigest(),
              'elapsed_seconds': time.perf_counter()-started}
    (HERE / 'conditional_compiler_review.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
