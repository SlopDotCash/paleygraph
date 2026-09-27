#!/usr/bin/env python3
"""Bounded exact checks for the pass41 character smoothing identities."""
from collections import Counter
from fractions import Fraction
from itertools import product
import json
from pathlib import Path


def check_case(p, V, A, t):
    chi = [0] + [1 if pow(x, (p - 1) // 2, p) == 1 else -1
                 for x in range(1, p)]
    F = [sum(chi[(x - v) % p] for v in V) for x in range(p)]
    M = {j: sum(x ** j for x in F) for j in range(1, 7)}
    n, m = len(V), len(A)
    assert M[1] == 0 and M[2] == p * n - n * n
    mixed = sum(c * c for c in Counter((a - v) % p
                                     for a in A for v in V).values())
    variance = mixed - Fraction(m * m * n * n, p)
    S = sum(F[a] for a in A)
    totals = dict(h6=0, err6=0, err2=0, paired_h=0,
                  paired_h2=0, paired_err=0, paired_err2=0)
    for shifts in product(range(p), repeat=t):
        Hnum = [sum(F[(x - s) % p] for s in shifts) for x in range(p)]
        Dnum = [t * F[x] - Hnum[x] for x in range(p)]
        totals['h6'] += sum(x ** 6 for x in Hnum)
        totals['err6'] += sum(x ** 6 for x in Dnum)
        totals['err2'] += sum(x * x for x in Dnum)
        paired_h = sum(Hnum[a] for a in A)
        paired_err = t * S - paired_h
        totals['paired_h'] += paired_h
        totals['paired_h2'] += paired_h * paired_h
        totals['paired_err'] += paired_err
        totals['paired_err2'] += paired_err * paired_err
    cases = p ** t
    avg_h6 = Fraction(totals['h6'], cases * t ** 6)
    predicted = (Fraction(M[6], t ** 5)
                 + Fraction(15 * (t - 1) * M[4] * M[2], p * t ** 5)
                 + Fraction(10 * (t - 1) * M[3] ** 2, p * t ** 5)
                 + Fraction(15 * (t - 1) * (t - 2) * M[2] ** 3,
                            p * p * t ** 5))
    assert avg_h6 == predicted
    assert Fraction(totals['err2'], cases * t * t) == Fraction(t + 1, t) * M[2]
    assert totals['paired_h'] == 0
    assert Fraction(totals['paired_h2'], cases * t * t) == variance / t
    assert Fraction(totals['paired_err'], cases * t) == S
    assert Fraction(totals['paired_err2'], cases * t * t) == S * S + variance / t
    assert Fraction(totals['err6'], cases * t ** 6) >= M[6]
    if n ** 4 <= p:
        bound = Fraction(p * n ** 3 * (15 * (t * t + 7 * t - 7) + 5 * n),
                         t ** 5)
        assert avg_h6 <= bound
        if t >= 2 and t ** 5 >= n:
            assert bound <= Fraction(325, 32) * p * n ** 3
    return dict(p=p, V=V, A=A, t=t, shift_tuples=cases,
                M2=M[2], M3=M[3], M4=M[4], M6=M[6],
                expected_smoothed_M6=str(avg_h6),
                expected_smoothed_ratio=str(avg_h6 / (p * n ** 3)),
                original_pairing=S, mixed_additive_energy=mixed,
                identities_passed=True)


def main():
    specs = [(17, [0, 1], [0, 1, 2, 4, 8], t) for t in (1, 2, 3)]
    specs += [(29, [0, 1], [0, 2, 3, 5, 7, 11], t) for t in (1, 2, 3)]
    specs += [(97, [0, 1, 3], [0, 1, 2, 4, 8, 16, 32], 2)]
    cases = [check_case(*spec) for spec in specs]
    result = dict(scope='bounded exact verification; not an original-moment upper bound',
                  cases=cases, total_shift_tuples=sum(c['shift_tuples'] for c in cases),
                  all_passed=all(c['identities_passed'] for c in cases))
    dest = Path(__file__).resolve().parents[1] / 'results' / 'parallel41_classical_smoothing_2026_09_06.json'
    dest.write_text(json.dumps(result, indent=2) + '\n')
    print(json.dumps(dict(all_passed=result['all_passed'],
                          total_shift_tuples=result['total_shift_tuples'],
                          cases=len(cases), result=str(dest))))


if __name__ == '__main__':
    main()
