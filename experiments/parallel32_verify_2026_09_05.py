#!/usr/bin/env python3
"""Exact finite checks of the dyadic trinomial inverse bounds.

The uniform Fourier argument is in the companion note; finite tests do
not prove it. No Lean or hosted proof job is launched.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from pathlib import Path
import importlib.util
import json
import time

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'pass31', ROOT / 'experiments/parallel31_short_multiples_2026_09_05.py')
MOD = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(MOD)
CHECKS = Counter()


def check(value, label):
    assert value, label
    CHECKS[label] += 1


def sparse_action(vector, terms):
    """Signed cyclic row action, separate from the norm-descent multiplication."""
    d = len(vector)
    out = [0] * d
    for exponent, coefficient in terms:
        for i, value in enumerate(vector):
            turns, j = divmod(i + exponent, d)
            out[j] += coefficient * value * (-1 if turns % 2 else 1)
    return out


def base_inverse(d):
    # Coefficients solve (1+X+X^2)b=1 modulo X^d+1.
    pattern = (1, 0, -1) if d % 3 == 1 else (0, -1, 1)
    return [pattern[i % 3] for i in range(d)]


def encode(value):
    return {'numerator': value.numerator, 'denominator': value.denominator}


def main():
    started = time.monotonic()
    rows = []
    largest = {}
    for n in [4, 8, 16, 32, 64, 128]:
        d = n // 2
        base = base_inverse(d)
        epsilon = 1 if n % 3 == 1 else -1
        t = (n - epsilon) // 3
        check(sparse_action(base, [(0, 1), (1, 1), (2, 1)]) == [1] + [0] * (d-1),
              'explicit_cyclotomic_unit_inverse')
        check(sum(x*x for x in base) == sum(map(abs, base)) == t,
              'exact_base_inverse_norms')
        check(81 * t <= 32 * n, 'uniform_constant_four_thirds')
        maximum = (Fraction(0), None)
        for b in range(n):
            f = sparse_action([1] + [0]*(d-1), [(0, 1), (1, 1), (b, 1)])
            norm, adjugate = MOD.norm_adjugate(f)
            check(norm > 0, 'nonzero_trinomial_norm')
            check(sparse_action(adjugate, [(0, 1), (1, 1), (b, 1)]) == [norm]+[0]*(d-1),
                  'independent_inverse_row_identity')
            l1 = sum(map(abs, adjugate))
            l2sq = sum(x*x for x in adjugate)
            check(l2sq <= 9*t*norm*norm, 'fourier_dominance_l2_consequence')
            check(l1*l1 <= d*l2sq, 'coefficient_cauchy')
            check(l1*l1 <= 9*d*t*norm*norm, 'sharper_inverse_l1_bound')
            check(3*l1 <= 4*n*norm, 'uniform_inverse_l1_bound')
            # Check a six-term signed target, allowing reduction and repeats.
            terms = [(0, 1), (1, -1), (b, 1), (2*b+1, -1), (3*b+2, 1), (n-1, 1)]
            target = sparse_action([1]+[0]*(d-1), terms)
            numerator = sparse_action(adjugate, terms)
            target_l1 = sum(map(abs, target))
            check(sparse_action(numerator, [(0, 1), (1, 1), (b, 1)]) == [norm*x for x in target],
                  'six_term_rational_quotient_identity')
            check(sum(map(abs, numerator)) <= target_l1*l1,
                  'signed_convolution_bound')
            ratio = Fraction(l1, n*norm)
            if ratio > maximum[0]:
                maximum = ratio, b
            rows.append({'n': n, 'b': b, 'norm': norm,
                         'inverse_l1': encode(Fraction(l1, norm)),
                         'inverse_l2_squared': encode(Fraction(l2sq, norm*norm))})
        largest[str(n)] = {'b': maximum[1], 'inverse_l1_over_n': encode(maximum[0])}

    # Lifting an inverse from a smaller dyadic conductor preserves coefficient
    # norms exactly. No assumption of integrality for all field relations.
    lifts = []
    for m, b, n in [(4, 2, 32), (8, 3, 64), (16, 5, 128), (32, 7, 128), (64, 19, 256)]:
        small_d = m//2
        f = sparse_action([1]+[0]*(small_d-1), [(0, 1), (1, 1), (b, 1)])
        norm, adjugate = MOD.norm_adjugate(f)
        c = n//m
        lifted = [0]*(n//2)
        lifted[::c] = adjugate
        check(sparse_action(lifted, [(0, 1), (c, 1), (c*b, 1)]) == [norm]+[0]*(n//2-1),
              'conductor_lift_inverse')
        check(sum(map(abs, lifted)) == sum(map(abs, adjugate)), 'conductor_lift_norm_preserved')
        check(3*sum(map(abs, lifted)) <= 4*m*norm, 'bound_uses_conductor_not_ambient_order')
        lifts.append({'conductor': m, 'ambient_order': n, 'b': b})

    certificate_path = 'results/parallel31_short_multiples_2026_09_05.json'
    known = json.loads((ROOT/certificate_path).read_text())
    inverse_l1 = Fraction(sum(map(abs, known['adjugate'])), known['norm'])
    for row in known['orbits']:
        length = row['quotient_l1']
        check(length <= 6*inverse_l1, 'all_119_minima_obey_actual_inverse_bound')
        check(length <= 8*known['n'], 'all_119_minima_obey_uniform_bound')
        if row['kind'] == 'primitive':
            check(row['cancellation_graph']['cycle_rank'] <= 4*known['n']-2,
                  'primitive_cycle_upper_bound')
    inputs = ['experiments/parallel32_verify_2026_09_05.py',
              'experiments/parallel31_short_multiples_2026_09_05.py', certificate_path]
    result = {
        'status': 'exact finite inverse-bound checks passed; uniform argument is ordinary mathematics; full goal unproved',
        'check_counts': dict(CHECKS), 'trinomial_rows': rows,
        'largest_inverse_l1_over_order': largest, 'conductor_lifts': lifts,
        'fixed_field': {'p': known['p'], 'n': known['n'],
                        'inverse_l1': encode(inverse_l1),
                        'six_term_bound': encode(6*inverse_l1),
                        'integer_six_term_bound': (6*inverse_l1).__floor__(),
                        'uniform_six_term_bound': 8*known['n']},
        'elapsed_seconds': time.monotonic()-started,
        'input_sha256': {p: sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
        'scope': 'No bound on the number of six-term outputs, no stronger Paley exponent, no Lean verification or separate-author review.'}
    (ROOT/'results/parallel32_verification_2026_09_05.json').write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({k: result[k] for k in ['status', 'check_counts', 'largest_inverse_l1_over_order',
                                           'fixed_field', 'elapsed_seconds']}, indent=2))


if __name__ == '__main__':
    main()
