#!/usr/bin/env python3
"""Exact known-digit cuts, with floating LP used only for proposal discovery."""
from copy import deepcopy
from fractions import Fraction as F
from itertools import product
from math import prod, ceil, floor
from pathlib import Path
import sys
import numpy as np
import sympy as sp
from scipy.optimize import linprog

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round21'))
from lattice_erasure import encode


def pair(x):
    x = F(x); return [x.numerator, x.denominator]


def discover(q, rows):
    """Find lambda reducing ||q-M_K^T lambda||1; optionally certify optimality.

    Normalize before the floating solve. Reconstruct a proposed zero-residual
    subsystem over Q. A dual vector is accepted only after exact feasibility
    and primal/dual equality checks; LP status alone certifies nothing.
    """
    N, h = len(q), len(rows)
    A = sp.Matrix(rows).T if h else sp.zeros(N, 0)
    zero = [F(0)]*h
    scale = max(map(abs, q), default=F(0))
    if h == 0 or scale == 0:
        residual = list(q); lam = zero; method = 'trivial'
    else:
        a = np.array(A.tolist(), dtype=float)
        b = np.array([float(x/scale) for x in q])
        lp = linprog(np.r_[np.zeros(h), np.ones(N)],
                     A_ub=np.vstack([np.c_[-a, -np.eye(N)], np.c_[a, -np.eye(N)]]),
                     b_ub=np.r_[-b, b], bounds=[(None, None)]*h+[(0, None)]*N,
                     method='highs', options={'time_limit': 20})
        candidates = [('zero', zero)]
        if lp.success:
            proposed = lp.x[:h]
            candidates.append(('rational_approximation', [F(float(v)).limit_denominator(10**8)*scale for v in proposed]))
            active = [i for i, x in enumerate(b-a@proposed) if abs(x) < 1e-7]
            if active:
                pivots = A[active, :].T.rref()[1]
                chosen = [active[i] for i in pivots]
                if len(chosen) == h:
                    exact = A[chosen, :].inv()*sp.Matrix([q[i] for i in chosen])
                    candidates.append(('exact_active_subsystem', [F(v) for v in exact]))
        def evaluate(lam):
            return [q[i]-sum(F(A[i,j])*lam[j] for j in range(h)) for i in range(N)]
        method, lam = min(candidates, key=lambda c: sum(map(abs, evaluate(c[1]))))
        residual = evaluate(lam)
    # Derive a feasible dual from the exact residual's zero/nonzero split.
    Z = [i for i, x in enumerate(residual) if x == 0]
    nz = [i for i, x in enumerate(residual) if x != 0]
    dual = [F((x > 0)-(x < 0)) for x in residual]
    optimal = False
    try:
        if Z and h:
            rhs = -A[nz, :].T*sp.Matrix([dual[i] for i in nz]) if nz else sp.zeros(h, 1)
            solution, params = A[Z, :].T.gauss_jordan_solve(rhs)
            solution = solution.subs({t: 0 for t in params})
            for i, value in zip(Z, solution): dual[i] = F(value)
        optimal = (all(abs(v) <= 1 for v in dual)
                   and all(sum(F(A[i,j])*dual[i] for i in range(N)) == 0 for j in range(h))
                   and sum(q[i]*dual[i] for i in range(N)) == sum(map(abs, residual)))
    except ValueError:
        pass
    return {'lambda': [pair(v) for v in lam], 'residual': [pair(v) for v in residual],
            'l1': pair(sum(map(abs, residual))), 'discovery': method,
            'optimality_certified': optimal, 'dual': [pair(v) for v in dual] if optimal else None}


def improve(base):
    c = deepcopy(base); N, f, U = c['N'], c['relation'], c['erased']
    K = [i for i in range(N) if i not in U]
    rows = [[f[i-t] if i >= t else -f[N+i-t] for t in range(N)] for i in K]
    cuts = []
    for j, encoded in enumerate(c['centered_coordinate_forms']):
        cut = discover([F(*v) for v in encoded], rows)
        cut['radius'] = pair(F(c['p']//2, c['p'])*F(*cut['l1']))
        cuts.append(cut)
    c['conditioning'] = {'known_coordinates': K, 'cuts': cuts,
                         'rule': 'intersect_distinct_centered_integer_intervals'}
    c['unconditioned_universal_cap'] = c['universal_candidate_box_cap']
    c['unconditioned_unique_completion'] = c['universal_unique_completion']
    radii = [min(F(*old), F(*cut['radius'])) for old, cut in zip(c['coordinate_radii'], cuts)]
    c['universal_candidate_box_cap'] = prod(floor(2*radii[j])+1 for j in c['visible_directions'])
    c['universal_unique_completion'] = all(2*radii[j] < 1 for j in c['visible_directions'])
    return c


def intervals(c, known):
    U, k = c['erased'], c['k']; w = c['projection_weights']
    value = sum(w[i]*y for i, y in known.items()); rhs = -value % k
    h = c['kernel_gcd']
    if rhs % h: return None, None
    x0 = [rhs//h*v for v in c['bezout_lift']]
    boxes = []
    for j, cut in enumerate(c['conditioning']['cuts']):
        oldcenter = -sum(F(*c['inverse'][i][j])*x0[i] for i in range(len(U)))
        oldradius = F(*c['coordinate_radii'][j])
        center = oldcenter+sum(F(*v)*known[i] for i, v in zip(c['conditioning']['known_coordinates'], cut['lambda']))
        radius = F(*cut['radius'])
        boxes.append([max(ceil(oldcenter-oldradius), ceil(center-radius)),
                      min(floor(oldcenter+oldradius), floor(center+radius))])
    numerator = value+sum(w[i]*x for i, x in zip(U, x0)); assert numerator % k == 0
    return boxes, numerator//k % c['p']


def decode(c, known, candidate_budget=50000):
    if type(candidate_budget) is not int or candidate_budget < 1: raise ValueError('invalid candidate budget')
    if set(known) != set(c['conditioning']['known_coordinates']) or any(type(i) is not int or type(y) is not int for i, y in known.items()):
        raise ValueError('known coordinates must be exactly the complement of erasures')
    def empty(reason, boxes=None):
        r = {'status': 'complete', 'reason': reason, 'count': 0, 'completions': [],
             'candidate_box_size': 0, 'scalar_candidates_checked': 0}
        if boxes is not None: r['intervals'] = boxes
        return r
    if any(abs(y) > c['digit_bound'] for y in known.values()): return empty('known_digit_height')
    boxes, a0 = intervals(c, known)
    if boxes is None: return empty('syndrome_divisibility')
    if any(lo > hi for lo, hi in boxes): return empty('empty_integer_coordinate_interval', boxes)
    visible = c['visible_directions']; size = prod(boxes[j][1]-boxes[j][0]+1 for j in visible)
    if size > candidate_budget:
        return {'status': 'budget_exceeded', 'reason': 'candidate_box_too_large',
                'candidate_box_size': size, 'candidate_budget': candidate_budget, 'intervals': boxes}
    seen = set(); actual = []
    for values in product(*(range(boxes[j][0], boxes[j][1]+1) for j in visible)):
        a = (a0+sum(c['scalar_steps'][j]*z for j, z in zip(visible, values))) % c['p']
        if a in seen: continue
        seen.add(a); word, _ = encode(c['p'], c['g'], c['relation'], a)
        if all(word[i] == y for i, y in known.items()): actual.append({'scalar': a, 'digits': word})
    return {'status': 'complete', 'reason': 'exhausted_projected_candidate_box',
            'count': len(actual), 'completions': sorted(actual, key=lambda x: x['scalar']),
            'candidate_box_size': size, 'scalar_candidates_checked': len(seen), 'intervals': boxes}
