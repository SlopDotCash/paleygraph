#!/usr/bin/env python3
"""Independent literal-column derivation of the three-mark degree-six coefficient.

No import from the generating compiler, coefficient lane, or its polynomial
arithmetic. Symbolic factorization uses SymPy after independent integer DP.
"""
from collections import defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import product
from math import factorial
from pathlib import Path
import json
import sys
import time
import sympy as sp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'coefficient_structure'


def multiply(a, b):
    out = defaultdict(Fraction)
    for (i, j, k), x in a.items():
        for (ii, jj, kk), y in b.items():
            if i + ii <= 6 and j + jj <= 6 and k + kk <= 12:
                out[i + ii, j + jj, k + kk] += x * y
    return {index: value for index, value in out.items() if value}


def column(a, b, outside):
    power = int(outside)
    return {index: value for index, value in
            [((0, 0, 0), 1), ((1, 0, power), a),
             ((0, 1, power), b), ((1, 1, power), a * b)] if value}


def literal_pair(q, u, w, s):
    counts = {(0, s): 1, (s, 0): 1,
              (1, 1): (q - 3 - 2*s)//4, (-1, -1): (q - 3 + 2*s)//4,
              (1, -1): (q - 1)//4, (-1, 1): (q - 1)//4}
    polynomial = {(0, 0, 0): Fraction(1)}
    for a, b in zip(u, w):
        counts[a, b] -= 1
        polynomial = multiply(polynomial, column(a, b, False))
    assert min(counts.values()) >= 0 and sum(counts.values()) == q - 3
    for (a, b), multiplicity in counts.items():
        # Deliberately multiply literal columns, not grouped multinomials.
        for _ in range(multiplicity):
            polynomial = multiply(polynomial, column(a, b, True))
    return polynomial


def literal_walsh(q):
    out = defaultdict(Fraction)
    for u in product((-1, 1), repeat=3):
        for w in product((-1, 1), repeat=3):
            mark_character = u[0] * u[1] * w[0] * w[1] * w[2]
            for sign in (-1, 1):
                for index, coefficient in literal_pair(q, u, w, sign).items():
                    out[index] += Fraction(sign * mark_character, 64) * coefficient
    return {index: value for index, value in out.items() if value}


def falling(x, count):
    out = 1
    for j in range(count):
        out *= x - j
    return out


def rational(x):
    return sp.Rational(x.numerator, x.denominator)


def main():
    start = time.perf_counter()
    saved = json.loads((SOURCE / 'symbolic_formula.json').read_text())
    fixed = json.loads((SOURCE / 'fixed_q_polynomials.json').read_text())
    q, r, t, w, v = sp.symbols('q r t w v')
    B = sp.expand(sp.prod(1 + v*(a*t + b*w + a*b*t*w)
                          for a, b in product((-1, 1), repeat=2)))
    explicit_B = 1 - 2*v**2*(t**2 + w**2 + t**2*w**2) + 8*v**3*t**2*w**2 \
                 + v**4*((t**2*w**2-t**2-w**2)**2 - 4*t**2*w**2)
    assert sp.expand(B - explicit_B) == 0
    bdict = {(int(i), int(j), int(k)): Fraction(int(value))
             for (i, j, k), value in sp.Poly(B, t, w, v).terms()}
    assert bdict.pop((0, 0, 0)) == 1
    assert min(k for _, _, k in bdict) == 2
    W = literal_walsh(17)
    T, current = [], W
    for power in range(7):
        T.append([current.get((6, 6, k), Fraction()) for k in range(13)])
        current = multiply(current, bdict)
    assert T == [[Fraction(*value) for value in row] for row in saved['binomial_basis_T']]
    assert all(not value for row in T[3:] for value in row)
    print('independent W17 and seven binomial rows match', flush=True)
    h = (q - 17)/4
    G = [sp.expand(sum(rational(T[j][k]) * falling(h, j)/factorial(j)
                       for j in range(7))) for k in range(13)]
    assert all(sp.expand(g - sp.sympify(export, locals={'q': q})) == 0
               for g, export in zip(G, saved['union_coefficients_G_as_polynomials_in_q']))
    assert all(g == 0 or sp.degree(g, q) <= 2 for g in G)
    denominator = falling(q - 3, 12)
    numerator = sp.expand(sum(G[k] * falling(r, k) * falling(q - 3 - k, 12 - k)
                              for k in range(13)))
    endpoints = falling(r, 4) * falling(q - r - 3, 4)
    quotient, remainder = sp.div(numerator, endpoints, r)
    assert sp.expand(remainder) == 0 and sp.degree(quotient, r) == 4
    formula = numerator / denominator
    exported = sp.sympify(saved['formula_in_q_r'], locals={'q': q, 'r': r})
    assert sp.cancel(formula - exported) == 0
    assert sp.expand(sp.Poly(quotient, r).LC() - 2*(q - 21)*(q - 19)) == 0
    numeric_replays = []
    # A second literal-column reconstruction verifies the q-step prediction
    # without any call to the original compiler or canonicalization routine.
    literal29 = literal_walsh(29)
    direct29 = [literal29.get((6, 6, k), Fraction()) for k in range(13)]
    assert all(rational(value) == g.subs(q, 29) for value, g in zip(direct29, G))
    for qvalue, coefficients in [(17, T[0]), (29, direct29)]:
        for n in [6, 7, 8, qvalue - 3, qvalue - 2, qvalue - 1, qvalue]:
            outside = n - 3
            direct = sum((coefficient * Fraction(falling(outside, k), falling(qvalue - 3, k))
                          for k, coefficient in enumerate(coefficients)), Fraction())
            assert rational(direct) == formula.subs({q: qvalue, r: outside})
            numeric_replays.append({'q': qvalue, 'n': n, 'coefficient': [direct.numerator, direct.denominator]})
    assert formula.subs({q: 29, r: 4}) == -sp.Rational(2, 7475)
    fixed_replays = []
    for record in fixed['cases']:
        qvalue = record['q']
        coefs = [Fraction(*x) for x in record['union_coefficients']]
        assert all(rational(a) == g.subs(q, qvalue) for a, g in zip(coefs, G))
        original = sp.sympify(record['factored_in_r'], locals={'r': r})
        assert sp.cancel(formula.subs(q, qvalue) - original) == 0
        fixed_replays.append(qvalue)
    exceptional = sp.degree(numerator.subs(q, 21), r)
    assert exceptional == 11
    result = {'status': 'passed', 'proof_method': 'independent literal-column expansion plus exact finite binomial identity; no fitted coefficients',
        'literal_column_walsh_reconstructions': [17, 29], 'signed_marked_pair_polynomials_per_q': 128,
        'base_W_nonzero_terms': len(W), 'all_seven_binomial_rows_equal': True,
        'binomial_rows_3_through_6_vanish_exactly': True,
        'common_factor_B_expanded': str(B),
        'union_coefficients_degree_in_q_at_most': 2,
        'degree_bound_in_r': 12, 'degree_over_Qq': 12,
        'eight_endpoint_factor_product': str(sp.expand(endpoints)),
        'quartic_after_division': str(sp.factor(quotient)),
        'quartic_leading_coefficient': str(sp.factor(sp.Poly(quotient, r).LC())),
        'exceptional_formal_specialization': {'q': 21, 'degree_in_r': int(exceptional)},
        'numeric_coefficients_from_literal_columns': numeric_replays,
        'fixed_q_export_factorizations_rechecked': fixed_replays,
        'coefficient_formula_exactly_equal': True,
        'validity_domain': 'q>=17,q=1mod4,3<=n<=q; degree6; graph interpretation additionally requires a realized conference matrix',
        'elapsed_seconds': round(time.perf_counter() - start, 6),
        'source_sha256': {name: sha256((SOURCE / name).read_bytes()).hexdigest()
                         for name in ['symbolic_formula.py', 'symbolic_formula.json', 'coefficient_polynomial.py', 'fixed_q_polynomials.json']},
        'review_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'review_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print('passed factorization, q29 literal replay, numeric coefficients; seconds', result['elapsed_seconds'], flush=True)


if __name__ == '__main__':
    main()
