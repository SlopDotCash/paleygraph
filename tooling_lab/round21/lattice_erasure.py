#!/usr/bin/env python3
"""Certified bounded erasure decoding with a nested lattice scalar projection.

LLL chooses coordinates only. Exact inverse identities and scalar steps carry
the proof. Coordinates that do not change the scalar are never enumerated.
The resulting finite candidate box is a superset, filtered by re-encoding.
"""
from fractions import Fraction
from itertools import product
from math import gcd, prod
import sys
from pathlib import Path
import sympy as sp
from fpylll import IntegerMatrix, LLL

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from realizability import prepare as inverse_prepare
sys.path.insert(0, str(HERE.parent/'round19'))
from projection_codec import encode


def egcd(a, b):
    old_r, r, old_s, s, old_t, t = a, b, 1, 0, 0, 1
    while r:
        q = old_r//r
        old_r, r, old_s, s, old_t, t = r, old_r-q*r, s, old_s-q*s, t, old_t-q*t
    if old_r < 0: return -old_r, -old_s, -old_t
    return old_r, old_s, old_t


def modular_kernel(weights, modulus):
    """Basis of {x: weights.x == 0 mod modulus}, plus a Bezout lift."""
    d = len(weights); column = [w % modulus for w in weights]+[modulus]
    U = [[int(i == j) for j in range(d+1)] for i in range(d+1)]
    for i in range(1, d+1):
        a, b = column[0], column[i]
        if a == b == 0: continue
        g, s, t = egcd(a, b); row0, rowi = U[0], U[i]
        U[0] = [s*x+t*y for x, y in zip(row0, rowi)]
        U[i] = [(-b//g)*x+(a//g)*y for x, y in zip(row0, rowi)]
        column[0], column[i] = g, 0
    return column[0], U[0][:d], [row[:d] for row in U[1:]]


def pair(x):
    x = Fraction(x); return [x.numerator, x.denominator]


def prepare(p, g, f, erased):
    base = inverse_prepare(p, g, f); N = len(f); A = sorted(erased)
    if len(set(A)) != len(A) or any(type(i) is not int or i not in range(N) for i in A):
        raise ValueError('invalid erasure set')
    adj, k = base['adj'], base['k']; weights = [adj[0]]+[-adj[N-j] for j in range(1, N)]
    d = len(A); wA = [weights[j] for j in A]
    if d:
        divisor, lift, raw = modular_kernel(wA, k)
        mat = IntegerMatrix.from_matrix(raw); LLL.reduction(mat)
        B = sp.Matrix([[int(mat[i,j]) for j in range(d)] for i in range(d)])
        inv = B.inv(); steps = []
        for row in B.tolist():
            value = sum(int(x)*w for x, w in zip(row, wA)); assert value % k == 0
            steps.append(value//k % p)
        radii = [Fraction(base['bound'])*sum(Fraction(abs(inv[i,j])) for i in range(d)) for j in range(d)]
        inverse = [[pair(inv[i,j]) for j in range(d)] for i in range(d)]
        basis = [[int(v) for v in row] for row in B.tolist()]
    else:
        divisor, lift, basis, inverse, steps, radii = k, [], [], [], [], []
    visible = [i for i, step in enumerate(steps) if step]
    caps = [int(2*radii[j])+1 for j in visible]
    return {'p': p, 'g': g, 'N': N, 'relation': f, 'norm_f': base['norm_f'], 'k': k,
            'adj': adj, 'digit_bound': base['bound'], 'erased': A,
            'projection_weights': weights, 'kernel_gcd': divisor, 'kernel_index': k//divisor,
            'bezout_lift': lift, 'basis': basis, 'inverse': inverse,
            'coordinate_radii': [pair(v) for v in radii], 'scalar_steps': steps,
            'visible_directions': visible, 'invisible_directions': [i for i in range(d) if not steps[i]],
            'universal_candidate_box_cap': prod(caps),
            'universal_unique_completion': all(2*radii[j] < 1 for j in visible)}


def decode(c, known, candidate_budget=100000):
    N, p, k, B = c['N'], c['p'], c['k'], c['digit_bound']; A = c['erased']; d = len(A)
    if type(candidate_budget) is not int or candidate_budget < 1: raise ValueError('invalid candidate budget')
    if set(known) != set(range(N))-set(A) or any(type(j) is not int or type(v) is not int for j, v in known.items()):
        raise ValueError('known coordinates must be exactly the complement of the erasures')
    if any(abs(v) > B for v in known.values()):
        return {'status': 'complete', 'reason': 'known_digit_height', 'count': 0, 'completions': [], 'candidate_box_size': 0, 'scalar_candidates_checked': 0}
    weights = c['projection_weights']; known_sum = sum(weights[j]*v for j, v in known.items())
    rhs = -known_sum % k; divisor = c['kernel_gcd']
    if rhs % divisor:
        return {'status': 'complete', 'reason': 'syndrome_divisibility', 'count': 0, 'completions': [], 'candidate_box_size': 0, 'scalar_candidates_checked': 0}
    x0 = [(rhs//divisor)*v for v in c['bezout_lift']]
    intervals = []
    for j in range(d):
        center = -sum(Fraction(*c['inverse'][i][j])*x0[i] for i in range(d))
        radius = Fraction(*c['coordinate_radii'][j])
        lower, upper = center-radius, center+radius
        lo = -((-lower.numerator)//lower.denominator); hi = upper.numerator//upper.denominator
        intervals.append([lo, hi])
    if any(lo > hi for lo, hi in intervals):
        return {'status': 'complete', 'reason': 'empty_integer_coordinate_interval', 'count': 0, 'completions': [],
                'candidate_box_size': 0, 'scalar_candidates_checked': 0, 'intervals': intervals}
    visible = c['visible_directions']; size = prod(intervals[j][1]-intervals[j][0]+1 for j in visible)
    if size > candidate_budget:
        return {'status': 'budget_exceeded', 'reason': 'candidate_box_too_large', 'candidate_box_size': size,
                'candidate_budget': candidate_budget, 'intervals': intervals}
    numerator = known_sum+sum(weights[j]*x for j, x in zip(A, x0)); assert numerator % k == 0
    base_scalar = numerator//k % p; seen = set(); actual = []
    for values in product(*(range(intervals[j][0], intervals[j][1]+1) for j in visible)):
        a = (base_scalar+sum(c['scalar_steps'][j]*z for j, z in zip(visible, values))) % p
        if a in seen: continue
        seen.add(a); word, _ = encode(p, c['g'], c['relation'], a)
        if all(word[j] == v for j, v in known.items()): actual.append({'scalar': a, 'digits': word})
    actual.sort(key=lambda r: r['scalar'])
    return {'status': 'complete', 'reason': 'exhausted_projected_candidate_box', 'count': len(actual),
            'completions': actual, 'candidate_box_size': size, 'scalar_candidates_checked': len(seen),
            'intervals': intervals}


def decode_cyclic(c, start, known, candidate_budget=100000):
    """Transport a consecutive cyclic erasure block to the certified first block."""
    N, p, g = c['N'], c['p'], c['g']; d = len(c['erased'])
    if c['erased'] != list(range(d)) or type(start) is not int or start not in range(N):
        raise ValueError('cyclic decoder requires an initial consecutive block and a valid start')
    missing = {(start+j) % N for j in range(d)}
    if set(known) != set(range(N))-missing: raise ValueError('wrong cyclic known-coordinate set')
    rotated = {j: known[(j+start) % N]*(1 if j+start < N else -1) for j in range(d, N)}
    out = decode(c, rotated, candidate_budget)
    if out['status'] != 'complete': return out
    result = {**out, 'completions': []}
    for item in out['completions']:
        a = item['scalar']*pow(g, -start, p) % p; original = [0]*N
        for j, v in enumerate(item['digits']): original[(j+start) % N] = v*(1 if j+start < N else -1)
        assert original == encode(p, g, c['relation'], a)[0]
        result['completions'].append({'scalar': a, 'digits': original})
    result['completions'].sort(key=lambda r: r['scalar'])
    return result
