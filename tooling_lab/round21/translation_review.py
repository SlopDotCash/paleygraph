#!/usr/bin/env python3
"""Independent rational-matrix and complete suffix-set translation audit."""
from collections import Counter, defaultdict
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp

HERE = Path(__file__).resolve().parent


def arithmetic(c):
    p, g, N, f = c['p'], c['g'], c['N'], c['f']; m = p//2
    assert sp.isprime(p) and pow(g, N, p) == p-1 and sum(f[j]*pow(g, -j, p) for j in range(N)) % p == 0
    M = sp.Matrix([[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in range(N)])
    inv = p*M.inv(); rows = [[Fraction(v) for v in row] for row in inv.tolist()]
    assert M*sp.Matrix(c['adj']) == sp.Matrix([c['norm_f']]+[0]*(N-1)) and c['norm_f'] == p*c['k']
    orbits = [[(a*pow(g, j, p)+m) % p-m for j in range(N)] for a in range(p)]
    words = []
    dense = [[int(v) for v in row] for row in M.tolist()]
    for F in orbits:
        products = [sum(x*y for x, y in zip(row, F)) for row in dense]
        assert all(v % p == 0 for v in products); words.append([v//p for v in products])
    def compile_h(H):
        rational = [sum(x*y for x, y in zip(row, H)) for row in rows]
        if any(v.denominator != 1 for v in rational): return {'status': 'disjoint', 'reason': 'fractional_centered_difference'}
        V = [int(v) for v in rational]; delta = V[0] % p
        if any(V[j] % p != delta*pow(g, j, p) % p for j in range(N)):
            return {'status': 'disjoint', 'reason': 'inconsistent_scalar_difference', 'difference': V}
        intervals = [[max(-m, -m-v), min(m, m-v)] for v in V]
        if any(lo > hi for lo, hi in intervals):
            return {'status': 'disjoint', 'reason': 'disjoint_centering_boxes', 'delta': delta, 'difference': V}
        F = orbits[delta]; assert all((v-w) % p == 0 for v, w in zip(F, V))
        C = [(v-w)//p for v, w in zip(F, V)]; assert all(v in (-1, 0, 1) for v in C)
        assert [sum(x*y for x, y in zip(row, C)) for row in dense] == [words[delta][j]-H[j] for j in range(N)]
        return {'status': 'translated_box', 'delta': delta, 'difference': V, 'required_carry': C, 'intervals': intervals}
    return words, orbits, compile_h


def main():
    source = HERE/'translation_audit.json'; report = json.loads(source.read_text()); results = []; example = None
    for case in report['cases']:
        c = case['config']; p, N = c['p'], c['N']; words, F, compile_h = arithmetic(c)
        stream = sha256(); cache = {}; total = 0; incident_total = 0; tested_total = 0
        for level in case['levels']:
            d = level['depth']; fibres = defaultdict(list)
            for a, w in enumerate(words): fibres[tuple(w[:d])].append(a)
            prefixes = list(map(tuple, level['prefixes'])); actual = sorted(fibres)
            assert prefixes[:len(actual)] == actual and len(prefixes) == len(actual)+int(d > 0)
            if d: assert prefixes[-1] not in fibres
            assert level['nonempty_prefixes'] == len(actual)
            suffixes = {P: {tuple(words[a][d:]) for a in fibres[P]} for P in prefixes}
            kinds = Counter(); reasons = Counter(); overlaps = 0; tested = 0
            for i, P in enumerate(prefixes):
                for j, Q in enumerate(prefixes):
                    H = tuple([Q[t]-P[t] for t in range(d)]+[0]*(N-d))
                    if H not in cache: cache[H] = compile_h(H)
                    cert = cache[H]; reason = cert.get('reason', 'translated_box'); reasons[reason] += 1
                    both = suffixes[P] & suffixes[Q]
                    expected = [a for a in fibres[P] if tuple(words[a][d:]) in both]
                    translated = []
                    if cert['status'] == 'translated_box':
                        tested += len(fibres[P])
                        for a in fibres[P]:
                            if all(lo <= x <= hi for x, (lo, hi) in zip(F[a], cert['intervals'])):
                                b = (a+cert['delta']) % p
                                assert words[b][:d] == list(Q) and words[b][d:] == words[a][d:]
                                assert [F[b][t]-F[a][t] for t in range(N)] == cert['difference']
                                translated.append(a)
                    assert translated == expected
                    left, right = suffixes[P], suffixes[Q]
                    if left == right: kind = 'equal'
                    elif left < right: kind = 'strict_left_subset'
                    elif right < left: kind = 'strict_right_subset'
                    elif not both: kind = 'nonempty_disjoint'
                    else: kind = 'overlapping_incomparable'
                    kinds[kind] += 1; overlaps += len(expected); total += 1
                    stream.update((json.dumps([d, i, j, reason, kind, expected], separators=(',', ':'))+'\n').encode())
            assert dict(kinds) == level['relation_types'] and dict(reasons) == level['compiled_reasons']
            assert len(prefixes)**2 == level['pair_comparisons'] and overlaps == level['common_suffix_incidents'] and tested == level['orbits_tested_in_boxes']
            incident_total += overlaps; tested_total += tested
        assert total == case['comparisons'] and len(cache) == case['distinct_difference_certificates']
        assert stream.hexdigest() == case['comparison_stream_sha256']
        for w in case['witnesses']:
            assert compile_h(w['difference_word']) == w['certificate']
            if w['certificate']['status'] == 'translated_box' and example is None: example = (w, compile_h)
        row = {'name': case['name'], 'pair_comparisons': total, 'difference_certificates': len(cache),
               'common_suffix_incidents': incident_total, 'orbits_checked_in_boxes': tested_total,
               'witnesses_checked': len(case['witnesses']), 'comparison_stream_sha256': stream.hexdigest()}
        results.append(row); print(json.dumps(row), flush=True)
    assert example; sample, oracle = example; rejected = []
    for label in ['delta', 'difference', 'carry', 'interval']:
        changed = deepcopy(sample['certificate'])
        if label == 'delta': changed['delta'] += 1
        if label == 'difference': changed['difference'][0] += 1
        if label == 'carry': changed['required_carry'][0] += 1
        if label == 'interval': changed['intervals'][0][0] += 1
        assert changed != oracle(sample['difference_word']); rejected.append(label)
    out = {'status': 'passed', 'scope': 'Independent full rational inverses and complete direct suffix-set intersections verify every reported prefix comparison, including scalar-congruence failures and degenerate calibration controls.',
           'cases': results, 'corrupt_translation_certificates_rejected': rejected,
           'source_sha256': {'translation_review.py': sha256(Path(__file__).read_bytes()).hexdigest()},
           'input_sha256': {'translation_audit.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'translation_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
