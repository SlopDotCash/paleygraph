#!/usr/bin/env python3
"""Bound the one-contraction contribution using the exact conference norm.

This is the standard centered spectral Cauchy-Schwarz inequality applied to
the newly isolated diagnostic. It bounds a conditional average correction,
not individual sets, eigenvalue outliers, or either prize.
"""
from fractions import Fraction
from hashlib import sha256
import json
from math import isqrt, prod
from pathlib import Path

from contraction_reduction import fraction_pair, reduce_three_marks

HERE = Path(__file__).resolve().parent


def rational_sqrt_upper(value, digits=20):
    assert value >= 0
    scale = 10**digits
    numerator = value.numerator*scale*scale
    root = isqrt(numerator//value.denominator)
    if root*root*value.denominator < numerator:
        root += 1
    out = Fraction(root, scale)
    assert out*out >= value
    return out


def envelope(record, n, degree=6):
    reduced = reduce_three_marks(record, n, degree)
    q = reduced['q']
    norm2 = norm3 = sum2 = sum3 = 0
    for cell in record['cells']:
        u, size = cell['pattern'], cell['size']
        if 0 in u:
            continue
        h2 = u[0]*u[1]+u[0]*u[2]+u[1]*u[2]
        h3 = prod(u)
        norm2 += size*h2*h2
        norm3 += size*h3*h3
        sum2 += size*h2
        sum3 += size*h3
    # On 1-perp, S has norm sqrt(q), since S^2=qI-J.
    q_bound_squared = Fraction((q*norm2-sum2*sum2)*(q*norm3-sum3*sum3), q)
    assert reduced['Q']**2 <= q_bound_squared
    c = Fraction(*reduced['Q_coefficient'])
    radius_squared = c*c*q_bound_squared
    base = Fraction(*reduced['cell_only_term'])
    radius_upper = rational_sqrt_upper(radius_squared)
    actual = Fraction(*reduced['second_moment'])
    assert (actual-base)**2 <= radius_squared
    relative_squared = radius_squared/(base*base) if base else None
    return {**reduced,
            'vector_statistics': {'h2_norm_squared': norm2, 'h3_norm_squared': norm3,
                                  'h2_sum': sum2, 'h3_sum': sum3},
            'Q_bound_squared': fraction_pair(q_bound_squared),
            'cell_only_envelope_radius_squared': fraction_pair(radius_squared),
            'envelope_radius_rational_upper': fraction_pair(radius_upper),
            'envelope_radius_upper_float': float(radius_upper),
            'relative_radius_rational_upper': fraction_pair(rational_sqrt_upper(relative_squared))
                                              if relative_squared is not None else None,
            'relative_radius_upper_float': float(rational_sqrt_upper(relative_squared))
                                           if relative_squared is not None else None,
            'scope': 'All symmetric balanced conference matrices realizing the same marked cell data satisfy this interval; an outer bound, not a realization claim.'}


def main():
    root = HERE.parent
    witness = json.loads((root/'marked_counts_ablation'/'witness.json').read_text())
    rows = [envelope(item['record'], 7) for item in witness['pair']]
    for p, n in ((1297,6), (1297,8), (65537,16), (1000033,31)):
        path = root/'marked_counts'/f'counts_p{p}_marks_0_1_2.json'
        row = envelope(json.loads(path.read_text()), n)
        row['counts_sha256'] = sha256(path.read_bytes()).hexdigest()
        rows.append(row)
    out = {'status': 'exact squared-radius certificates with rational outward rounding',
           'cases': rows,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ('conditional_moments.py','contraction_reduction.py','certified_envelope.py')}}
    (HERE/'contraction_envelopes.json').write_text(json.dumps(out, indent=2)+'\n')
    print(json.dumps([{'q':x['q'], 'n':x['n'], 'Q':x['Q'],
                       'relative_radius_upper':x['relative_radius_upper_float']} for x in rows], indent=2))


if __name__ == '__main__':
    main()
