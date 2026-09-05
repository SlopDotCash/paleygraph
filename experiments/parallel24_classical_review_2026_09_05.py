#!/usr/bin/env python3
"""Root's independent actual-field checks of the local sixth-moment drift."""
from collections import Counter
from hashlib import sha256
from math import comb
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]
counts = Counter()
reports = []
for p, C in ((17, (0, 1, 3, 7, 9, 14)),
             (29, (1, 3, 8, 9, 13, 21, 25)),
             (97, (0, 2, 5, 11, 36, 80))):
    n = len(C)
    chi = [0]+[1 if pow(x, (p-1)//2, p) == 1 else -1 for x in range(1, p)]
    inside, outside = list(C), [x for x in range(p) if x not in C]
    M, total_actual, total_drift = 0, 0, 0
    for x in range(p):
        U = [chi[(x-a) % p] for a in inside]
        V = [chi[(x-b) % p] for b in outside]
        F, N = sum(U), sum(u*u for u in U)
        den = n*(p-n)
        S = n*(p-1)+(p-2*n)*N
        P = N*(p-1-N)
        mu = [-p*F, S+2*F*F, -(4*p-3)*F,
              S+8*F*F+6*P, -(16*p-15)*F, S+32*F*F+30*P]
        for j in range(1, 7):
            assert sum((v-u)**j for u in U for v in V) == mu[j-1]
            counts['independently_enumerated_row_increment_moments'] += 1
        expanded = den*F**6+sum(comb(6, j)*F**(6-j)*mu[j-1] for j in range(1, 7))
        direct = sum((F+v-u)**6 for u in U for v in V)
        drift = ((den-6*(p-5))*F**6+(15*S-80*p+180)*F**4
                 +(15*S+90*P-96*p+122)*F*F+S+30*P)
        assert expanded == direct == drift
        counts['independent_sixth_row_expansions'] += 1
        M += F**6
        total_actual += direct
        total_drift += drift
    assert total_actual == total_drift
    counts['complete_swap_averages'] += 1
    reports.append({'p': p, 'C': C, 'M6': M, 'swap_count': n*(p-n),
                    'total_swap_M6': total_actual,
                    'scope': 'Outside the quartic slice; exact identities only.'})
inputs = ['research/parallel24-classical-exceptions-2026-09-05.md',
          'experiments/parallel24_classical_exceptions_2026_09_05.py',
          'results/parallel24_classical_exceptions_2026_09_05.json',
          'experiments/parallel24_classical_review_2026_09_05.py']
report = {'status': 'passed', 'counts': dict(counts), 'fields': reports,
          'input_sha256': {f: sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs},
          'scope': 'Independent root review checks. No exceptional-input upper bound, no formal verification.'}
out = ROOT/'results/parallel24_classical_review_2026_09_05.json'
out.write_text(json.dumps(report, indent=2)+'\n')
print(json.dumps(report, indent=2))
