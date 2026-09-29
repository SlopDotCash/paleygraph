#!/usr/bin/env python3
"""Independent critical-scaling and conference-norm check for c(q,n,6)Q."""
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from pathlib import Path
import json
import sys
import sympy as sp

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
SOURCE = HERE.parent / 'coefficient_structure'


def matrix(p):
    def character(x):
        x %= p
        return 0 if x == 0 else (1 if pow(x, (p-1)//2, p) == 1 else -1)
    S = [[character(x-y) for y in range(p)] for x in range(p)]
    assert all(sum(row) == 0 for row in S)
    assert all(S[i][j] == S[j][i] for i in range(p) for j in range(p))
    assert all(sum(a*b for a, b in zip(S[i], S[j])) == (p-1 if i == j else -1)
               for i in range(p) for j in range(p))
    return S


def main():
    formula_record = json.loads((SOURCE / 'symbolic_formula.json').read_text())
    review = json.loads((HERE / 'review_results.json').read_text())
    assert review['coefficient_formula_exactly_equal']
    assert review['source_sha256']['symbolic_formula.json'] == sha256((SOURCE / 'symbolic_formula.json').read_bytes()).hexdigest()
    q, r, t, alpha = sp.symbols('q r t alpha')
    numerator = sp.sympify(formula_record['numerator_factored'], locals={'q': q, 'r': r})
    denominator = sp.sympify(formula_record['common_denominator'], locals={'q': q})
    substitution = {q: t**4, r: alpha*t - 3}
    top = sp.Poly(numerator.subs(substitution), t)
    bottom = sp.Poly(denominator.subs(substitution), t)
    assert top.degree() == 40 and bottom.degree() == 48
    assert sp.simplify(top.LC()/bottom.LC() + alpha**4) == 0
    quotient = sp.sympify(review['quartic_after_division'], locals={'q': q, 'r': r})
    P = sp.Poly(sp.expand(-2*quotient), q, r)
    weighted_terms = []
    for (i, j), coefficient in P.terms():
        assert i + sp.Rational(j, 2) <= 5
        if j == 0 and i == 5:
            assert coefficient == 2
        weighted_terms.append({'q_degree': i, 'r_degree': j, 'coefficient': str(coefficient)})
    # With r=o(sqrt(q)), every term except 2q^5 is o(q^5): terms
    # at weighted degree5 have a positive power of r/sqrt(q).
    assert [(i, j) for (i, j), _ in P.terms() if i == 5 and j == 0] == [(5, 0)]
    for signs in product((-1, 1), repeat=3):
        h2 = sum(signs[i]*signs[j] for i, j in combinations(range(3), 2))
        h3 = signs[0]*signs[1]*signs[2]
        assert h2*h2 == 3 + 2*h2 and h3*h3 == 1
    Lvalues = {a*b + a*c + b*c for a, b, c in product((-1, 1), repeat=3)}
    assert Lvalues == {-1, 3}
    numeric = []
    for p in [17, 29]:
        S = matrix(p); checked = 0; largest_ratio = Fraction()
        for marks in combinations(range(p), 3):
            markset = set(marks)
            h2, h3 = [0]*p, [0]*p
            for x in range(p):
                if x not in markset:
                    a, b, c = [S[x][m] for m in marks]
                    h2[x], h3[x] = a*b + a*c + b*c, a*b*c
            a, b, c = S[marks[0]][marks[1]], S[marks[0]][marks[2]], S[marks[1]][marks[2]]
            L = a*b + a*c + b*c
            assert sum(h2) == -3-L
            assert sum(x*x for x in h2) == 3*p - 15 - 2*L
            assert sum(x*x for x in h3) == p - 3
            Q = sum(h2[x]*sum(S[x][y]*h3[y] for y in range(p)) for x in range(p))
            bound2 = p*(p-3)*(3*p-13)
            assert Q*Q <= bound2
            largest_ratio = max(largest_ratio, Fraction(Q*Q, bound2))
            checked += 1
        numeric.append({'q': p, 'all_marked_triples_checked': checked,
            'largest_observed_squared_ratio_to_uniform_bound': [largest_ratio.numerator, largest_ratio.denominator]})
    exact_examples = []
    expression = numerator/denominator
    for field_order, size in [(65537, 16), (1000033, 31), (6700417, 50)]:
        coeff = sp.cancel(expression.subs({q: field_order, r: size-3}))
        squared = sp.cancel(coeff**2 * field_order*(field_order-3)*(3*field_order-13))
        exact_examples.append({'q': field_order, 'n': size,
            'c': [int(sp.numer(coeff)), int(sp.denom(coeff))],
            'uniform_squared_absolute_correction_bound': [int(sp.numer(squared)), int(sp.denom(squared))]})
    output = {'status': 'passed',
        'critical_assumptions': 'fixed alpha>0; admissible q tends to infinity; n/q^(1/4) tends to alpha',
        'polynomial_substitution': 'q=t^4,r=alpha*t-3', 'numerator_degree_in_t': 40,
        'denominator_degree_in_t': 48, 'leading_coefficient_ratio': '-alpha^4',
        'coefficient_asymptotic': 'c(q,n,6) ~ -alpha^4 q^(-2)',
        'weighted_terms_proving_P_asymptotic': weighted_terms,
        'triangle_L_values': sorted(Lvalues), 'sum_h2': '-3-L',
        'h2_squared_norm': '3q-15-2L', 'h3_squared_norm': 'q-3',
        'uniform_contraction_bound_squared': 'q(q-3)(3q-13)',
        'uniform_critical_correction_bound': '|cQ| <= (sqrt(3)*alpha^4+o(1))*q^(-1/2)',
        'norm_proof': 'S symmetric and S^2=qI-J imply operator norm sqrt(q); apply Cauchy-Schwarz and the exact vector norms',
        'scope': 'uniform over realized conference matrices and marked triples; bounds only the additive conditional-second-moment correction, not a relative baseline error or individual T6 values',
        'all_triple_numeric_checks': numeric, 'finite_absolute_correction_bounds': exact_examples,
        'source_sha256': {'symbolic_formula.json': sha256((SOURCE / 'symbolic_formula.json').read_bytes()).hexdigest(),
                          'review_results.json': sha256((HERE / 'review_results.json').read_bytes()).hexdigest()},
        'review_sha256': sha256(Path(__file__).read_bytes()).hexdigest()}
    (HERE / 'critical_scaling_review.json').write_text(json.dumps(output, indent=2) + '\n')
    print('Critical scaling and all4334 marked-triple norm checks passed.')


if __name__ == '__main__':
    main()
