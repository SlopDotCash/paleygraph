#!/usr/bin/env python3
"""Exhaust every pair of nonempty residuals plus an empty representative."""
from collections import Counter, defaultdict
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
from translation import prepare, compile_difference, centered_orbit
from lattice_erasure import encode

HERE = Path(__file__).resolve().parent


def classify(left, right, intersection):
    if intersection == left == right: return 'equal'
    if intersection == left: return 'strict_left_subset'
    if intersection == right: return 'strict_right_subset'
    if intersection == 0: return 'nonempty_disjoint'
    return 'overlapping_incomparable'


def main():
    source = HERE.parent/'round17/realizability.json'; saved = json.loads(source.read_text())['cases'][0]['certificate']
    configs = [('p257_saved', saved['p'], saved['g'], saved['f']),
               ('p257_canonical', 257, 2, [2]+[0]*6+[1]),
               ('p17_squared', 17, 2, [4, 0, -1, 4]),
               ('p17_scalar_p', 17, 2, [17, 0, 0, 0]),
               ('p17_scalar_p_squared', 17, 2, [289, 0, 0, 0])]
    cases = []
    for name, p, g, f in configs:
        c = prepare(p, g, f); N = len(f); B = c['bound']
        words = [encode(p, g, f, a)[0] for a in range(p)]
        orbits = [centered_orbit(p, g, a, N) for a in range(p)]
        levels = []; witnesses = {}; stream = sha256(); cache = {}; calls = 0
        for d in range(N+1):
            fibres = defaultdict(list)
            for a, w in enumerate(words): fibres[tuple(w[:d])].append(a)
            prefixes = sorted(fibres)
            if d:
                # Scan only until the first absent bounded pattern. If every
                # bounded pattern occurs, use one outside the digit interval.
                missing = next((P for P in product(range(-B, B+1), repeat=d) if P not in fibres), None)
                if missing is None: missing = tuple([B+1]+[0]*(d-1))
                prefixes.append(missing); fibres[missing] = []
            suffixes = {P: {tuple(words[a][d:]): a for a in fibres[P]} for P in prefixes}
            kinds = Counter(); reasons = Counter(); overlap_total = 0; tested_orbits = 0
            for i, P in enumerate(prefixes):
                for j, Q in enumerate(prefixes):
                    H = tuple([y-x for x, y in zip(P, Q)]+[0]*(N-d))
                    if H not in cache: cache[H] = compile_difference(c, list(H))
                    cert = cache[H]; reason = cert.get('reason', 'translated_box'); reasons[reason] += 1
                    common = suffixes[P].keys() & suffixes[Q].keys()
                    expected = sorted(suffixes[P][S] for S in common)
                    actual = []
                    if cert['status'] == 'translated_box':
                        for a in fibres[P]:
                            tested_orbits += 1
                            if all(lo <= v <= hi for v, (lo, hi) in zip(orbits[a], cert['intervals'])):
                                actual.append(a)
                                b = (a+cert['delta']) % p
                                assert tuple(words[b][:d]) == Q and words[b][d:] == words[a][d:]
                                assert [orbits[b][j]-orbits[a][j] for j in range(N)] == cert['difference']
                    assert actual == expected
                    kind = classify(len(fibres[P]), len(fibres[Q]), len(actual)); kinds[kind] += 1
                    overlap_total += len(actual); calls += 1
                    item = [d, i, j, reason, kind, actual]
                    stream.update((json.dumps(item, separators=(',', ':'))+'\n').encode())
                    key = reason+'|'+kind
                    if key not in witnesses:
                        no = next((a for a in fibres[P] if tuple(words[a][d:]) not in common), None)
                        witnesses[key] = {'depth': d, 'P': P, 'Q': Q, 'difference_word': H,
                                          'certificate': cert, 'intersection_scalars': actual,
                                          'left_not_right_scalar': no}
            levels.append({'depth': d, 'prefixes': [list(P) for P in prefixes],
                           'nonempty_prefixes': sum(bool(fibres[P]) for P in prefixes),
                           'pair_comparisons': len(prefixes)**2, 'relation_types': dict(kinds),
                           'compiled_reasons': dict(reasons), 'common_suffix_incidents': overlap_total,
                           'orbits_tested_in_boxes': tested_orbits})
        out = {'name': name, 'config': c, 'levels': levels, 'comparisons': calls,
               'distinct_difference_certificates': len(cache), 'comparison_stream_sha256': stream.hexdigest(),
               'witnesses': list(witnesses.values())}
        cases.append(out); print(json.dumps({'name': name, 'comparisons': calls, 'difference_certificates': len(cache),
                                             'witness_types': sorted(witnesses)}), flush=True)
    report = {'status': 'produced', 'scope': 'Every pair of nonempty prefixes at every depth, plus one explicit empty prefix at each nonzero depth, in five complete scalar codebooks. Translated centering boxes give exact shared-suffix sets, including zero-calibration controls.',
              'cases': cases, 'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['translation_audit.py', 'translation.py', '../round17/realizability.py', '../round17/norm_carry.py']},
              'input_sha256': {'../round17/realizability.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'translation_audit.json').write_text(json.dumps(report, separators=(',', ':'))+'\n')


if __name__ == '__main__': main()
