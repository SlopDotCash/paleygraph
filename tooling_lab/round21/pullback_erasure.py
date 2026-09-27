#!/usr/bin/env python3
"""Retain cancellations by bounding dual coordinates in centered orbit space.

For x=E(a)_A and inverse-basis column v, x.v = (M_A^T v).F_a/p.
The centered cube gives a second exact radius, often much smaller than the
sum of independent digit bounds. This changes no candidate-coverage logic.
"""
from copy import deepcopy
from fractions import Fraction
from math import prod


def improve(base):
    c = deepcopy(base); N, p, f, A = c['N'], c['p'], c['relation'], c['erased']; d = len(A)
    inv = [[Fraction(*x) for x in row] for row in c['inverse']]
    forms = []; centered_radii = []; final_radii = []
    def pair(x): return [x.numerator, x.denominator]
    for j in range(d):
        form = [sum((f[i-t] if i >= t else -f[N+i-t])*inv[a][j] for a, i in enumerate(A)) for t in range(N)]
        radius = Fraction(p//2, p)*sum(map(abs, form))
        forms.append([pair(v) for v in form]); centered_radii.append(pair(radius))
        final_radii.append(min(Fraction(*c['coordinate_radii'][j]), radius))
    c['digit_coordinate_radii'] = c['coordinate_radii']; c['coordinate_radii'] = [pair(v) for v in final_radii]
    c['centered_coordinate_forms'] = forms; c['centered_coordinate_radii'] = centered_radii
    c['radius_rule'] = 'minimum_of_digit_cube_and_centered_orbit_pullback'
    c['universal_candidate_box_cap'] = prod(int(2*final_radii[j])+1 for j in c['visible_directions'])
    c['universal_unique_completion'] = all(2*final_radii[j] < 1 for j in c['visible_directions'])
    return c
