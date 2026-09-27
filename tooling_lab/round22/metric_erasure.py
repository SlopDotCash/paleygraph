#!/usr/bin/env python3
"""Propose an erasure basis reduced in inverse-image geometry, check exactly."""
from copy import deepcopy
from fractions import Fraction as F
from math import prod
from fpylll import IntegerMatrix, LLL
import sympy as sp
from conditioned_erasure import pair
from pullback_erasure import improve as pullback


def change_basis(base):
    c = deepcopy(base); d = len(c['erased']); N = c['N']; adj = c['adj']
    if not d: return c
    # Columns are the inverse-adjugate images of erased unit vectors.
    inverse_image = sp.Matrix([[adj[i-j] if i >= j else -adj[N+i-j] for i in range(N)] for j in c['erased']])
    old_basis = sp.Matrix(c['basis']); lifted = old_basis*inverse_image
    matrix = IntegerMatrix.from_matrix([[int(v) for v in row] for row in lifted.tolist()])
    transform = IntegerMatrix.identity(d); LLL.reduction(matrix, U=transform)
    S = sp.Matrix([[int(transform[i,j]) for j in range(d)] for i in range(d)])
    assert abs(S.det()) == 1
    new_basis = S*old_basis; new_lift = new_basis*inverse_image
    assert new_lift == sp.Matrix([[int(matrix[i,j]) for j in range(N)] for i in range(d)])
    T = new_basis.inv(); w = [c['projection_weights'][i] for i in c['erased']]
    steps = []
    for row in new_basis.tolist():
        numerator = sum(v*x for v,x in zip(row, w)); assert numerator % c['k'] == 0
        steps.append(int(numerator//c['k'] % c['p']))
    radii = [F(c['digit_bound'])*sum(abs(F(T[i,j])) for i in range(d)) for j in range(d)]
    c['basis'] = [[int(v) for v in row] for row in new_basis.tolist()]
    c['inverse'] = [[pair(T[i,j]) for j in range(d)] for i in range(d)]
    c['coordinate_radii'] = [pair(v) for v in radii]
    c['scalar_steps'] = steps
    c['visible_directions'] = [j for j,s in enumerate(steps) if s]
    c['invisible_directions'] = [j for j,s in enumerate(steps) if not s]
    c['universal_candidate_box_cap'] = prod(int(2*radii[j])+1 for j in c['visible_directions'])
    c['universal_unique_completion'] = all(2*radii[j] < 1 for j in c['visible_directions'])
    c['metric_change'] = {'rule': 'LLL_on_integer_adjugate_inverse_images', 'old_basis': base['basis'],
                          'unimodular_transform': [[int(v) for v in row] for row in S.tolist()]}
    return pullback(c)
