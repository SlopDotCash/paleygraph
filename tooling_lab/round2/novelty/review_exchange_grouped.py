#!/usr/bin/env python3
"""Independent character-type and marked-containment checks for grouped energies."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations
import importlib.util
import json
from math import comb
from pathlib import Path
import sys

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OBS = HERE.parent / 'observables'
sys.path.insert(0, str(OBS))
import exchange_grouped as candidate


def coefficient(rows, r, d):
    lim = (r, d-r, d-r)
    coefficients = {(0, 0, 0): 1}
    for a, b in rows:
        new = dict(coefficients)
        for key, value in coefficients.items():
            for axis, factor in enumerate((a*b, a, b)):
                nxt = list(key)
                nxt[axis] += 1
                if nxt[axis] <= lim[axis]:
                    nxt = tuple(nxt)
                    new[nxt] = new.get(nxt, 0) + value * factor
        coefficients = {key: value for key, value in new.items() if value}
    return coefficients.get(lim, 0)


def main():
    type_checks, coefficient_checks, row_pair_checks = 0, 0, 0
    for p in (5, 13, 17, 29, 37, 41, 61, 73, 89, 97, 101):
        chi = [0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1) for x in range(p)]
        K, groups = candidate.grouped_kernel_sums(p, min(6, p//2))
        diag = [(chi[-x % p], chi[-x % p]) for x in range(p)]
        off = [(chi[-x % p], chi[(1-x) % p]) for x in range(p)]
        for rows, name in ((diag, 'diagonal'), (off, 'distinct')):
            assert Counter(rows) == Counter({(a, b): m for a, b, m in groups[name] if m})
            type_checks += 1
            for d in range(1, min(6, p//2)+1):
                for r in range(d+1):
                    assert coefficient(rows, r, d) == candidate.grouped_coefficient(groups[name], r, d)
                    coefficient_checks += 1
        if p <= 17:
            base_hist = Counter(off)
            for a in range(p):
                for b in range(p):
                    if a == b:
                        continue
                    eps = chi[(b-a) % p]
                    actual = Counter((chi[(a-x) % p], chi[(b-x) % p]) for x in range(p))
                    assert actual == Counter({(eps*u, eps*v): count for (u, v), count in base_hist.items()})
                    row_pair_checks += 1
    containment_checks = 0
    input_pairs_examined = 0
    for p in (5, 7, 8, 9):
        for n in range(1, p//2+1):
            sets = [sum(1 << x for x in s) for s in combinations(range(p), n)]
            pairs = [(A, B, n-(A&B).bit_count()) for A in sets for B in sets]
            input_pairs_examined += len(pairs)
            denominators = [comb(p, n)*comb(n, j)*comb(p-n, j) for j in range(n+1)]
            for d in range(1, n+1):
                for r in range(d+1):
                    S = (1 << d)-1
                    T = (1 << r)-1 | sum(1 << x for x in range(d, 2*d-r))
                    counts = [0]*(n+1)
                    for A, B, j in pairs:
                        if A&S == S and B&T == T:
                            counts[j] += 1
                    for j, count in enumerate(counts):
                        expected = Fraction(count, denominators[j])
                        assert expected == candidate.containment_falling(p, n, j, r, d)
                        containment_checks += 1
    odd = []
    for p, n, d in ((13, 6, 3), (13, 6, 5), (17, 6, 3)):
        new = candidate.grouped_spectrum(p, n, d)
        old = candidate.exact_spectrum(p, n, d)
        assert new['energies'] == old['energies']
        assert new['energies'][2] > 0
        odd.append({'p': p, 'n': n, 'd': d, 'level_two_energy': str(new['energies'][2])})
    result = {'date': '2026-09-05', 'passed': True, 'row_type_histograms': type_checks,
              'single_coordinate_coefficient_checks': coefficient_checks,
              'all_distinct_row_pair_histograms': row_pair_checks,
              'exact_containment_probabilities': containment_checks,
              'ordered_input_pairs_enumerated': input_pairs_examined,
              'odd_degree_cases': odd,
              'source_sha256': sha256((OBS/'exchange_grouped.py').read_bytes()).hexdigest(),
              'review_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE/'exchange_grouped_review.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(result, indent=2))


if __name__ == '__main__':
    main()
