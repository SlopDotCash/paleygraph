#!/usr/bin/env python3
"""Exact residual-intersection reduction, including degenerate calibrations."""
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from realizability import prepare
from norm_carry import multiply


def centered_orbit(p, g, a, N):
    return [(a*pow(g, j, p)+p//2) % p-p//2 for j in range(N)]


def compile_difference(c, H):
    """For E(b)-E(a)=H, determine b-a and the exact box allowed for F_a."""
    p, N, k = c['p'], c['N'], c['k']; m = p//2
    if len(H) != N or any(type(v) is not int for v in H): raise ValueError('invalid difference word')
    numerator = multiply(c['adj'], H)
    if any(v % k for v in numerator): return {'status': 'disjoint', 'reason': 'fractional_centered_difference'}
    V = [v//k for v in numerator]; delta = V[0] % p
    if any(V[j] % p != delta*pow(c['g'], j, p) % p for j in range(N)):
        return {'status': 'disjoint', 'reason': 'inconsistent_scalar_difference', 'difference': V}
    intervals = [[max(-m, -m-v), min(m, m-v)] for v in V]
    if any(lo > hi for lo, hi in intervals):
        return {'status': 'disjoint', 'reason': 'disjoint_centering_boxes', 'delta': delta, 'difference': V}
    F_delta = centered_orbit(p, c['g'], delta, N)
    assert all((x-y) % p == 0 for x, y in zip(F_delta, V))
    carry = [(x-y)//p for x, y in zip(F_delta, V)]; assert set(carry) <= {-1, 0, 1}
    return {'status': 'translated_box', 'delta': delta, 'difference': V, 'required_carry': carry, 'intervals': intervals}


def contains_orbit(c, certificate, a):
    if certificate['status'] == 'disjoint': return False
    F = centered_orbit(c['p'], c['g'], a, c['N'])
    return all(lo <= x <= hi for x, (lo, hi) in zip(F, certificate['intervals']))
