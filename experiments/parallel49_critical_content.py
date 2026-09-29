#!/usr/bin/env python3
"""Exact critical-content, p-adic and nonsplit-field checks.

These check the accompanying reductions, not a uniform Paley bound.
"""
import argparse
from collections import Counter
from hashlib import sha256
from math import factorial, gcd
from pathlib import Path
import json
import sys

from parallel2_subgroup_2026_09_04 import primitive_polynomial, derivative
from quadruple_orbits import Ring, find_generator
from dyadic_energy_descent import QuadraticTower, order_mod

ROOT = Path(__file__).resolve().parents[1]


def pencil(coefficients, flint):
    f = flint.fmpz_poly(coefficients)
    a, d = f.derivative(), f.degree()
    b = a.derivative()
    values = [int(f.resultant(a + t*b)) for t in range(d+1)]
    row, differences = values[:], []
    while row:
        differences.append(row[0])
        row = [row[j+1]-row[j] for j in range(len(row)-1)]
    denominator = factorial(d)
    out, falling = [0]*(d+1), [1]
    for k in range(d+1):
        for j, c in enumerate(falling):
            out[j] += differences[k]*(denominator//factorial(k))*c
        following = [0]*(len(falling)+1)
        for j, c in enumerate(falling):
            following[j] -= k*c
            following[j+1] += c
        falling = following
    assert all(c % denominator == 0 for c in out)
    out = [c//denominator for c in out]
    assert out[0] == int(f.resultant(a))
    assert out[-1] == int(f.resultant(b))
    assert sum(c*(d+1)**j for j, c in enumerate(out)) == int(f.resultant(a+(d+1)*b))
    content, coarse = gcd(*out), gcd(out[0], out[-1])
    assert content and coarse % content == 0
    return out, content, coarse, len(values)+3


def valuation(value, p):
    assert value
    result = 0
    while value % p == 0:
        result += 1
        value //= p
    return result


def truncated_valuation(value, p, precision):
    return precision if value == 0 else valuation(value, p)


def evaluate(f, x, modulus):
    out = 0
    for c in reversed(f):
        out = (out*x+c) % modulus
    return out


def padic_check(f, content, n, case):
    p, initial = case['p'], case['generator']
    a, b = derivative(f), derivative(derivative(f))
    expected = valuation(content, p)
    for precision in [1, 2, 4, 8, 16, 32, 64]:
        modulus = p**precision
        g = pow(initial, p**(precision-1), modulus)
        assert pow(g, n, modulus) == 1 and pow(g, n//2, modulus) == modulus-1
        roots = [pow(1+pow(g, j, modulus), n, modulus) for j in range(1, n//2, 2)]
        mins = [min(truncated_valuation(evaluate(a, x, modulus), p, precision),
                    truncated_valuation(evaluate(b, x, modulus), p, precision)) for x in roots]
        if max(mins) >= precision or len(set(roots)) != len(roots):
            continue
        assert sum(mins) == expected
        pair_sums, nearest = [], []
        for i, x in enumerate(roots):
            distances = [valuation((x-y) % modulus, p) for j, y in enumerate(roots) if j != i]
            pair_sums.append(sum(distances))
            nearest.append(max(distances))
        baseline = sum(x-y for x, y in zip(pair_sums, nearest))
        levels = []
        for r in range(1, max(nearest)+2):
            counts = Counter(x % (p**r) for x in roots)
            Y = sum(c*(c-2) for c in counts.values() if c >= 3)
            levels.append(Y)
        assert levels[-1] == 0 and sum(levels) == baseline
        assert levels[0] == case['triple_mass_Y']
        corrections = [v-A+M for v, A, M in zip(mins, pair_sums, nearest)]
        assert min(corrections) >= 0 and sum(corrections) == expected-baseline
        return dict(n=n, p=p, precision=precision, content_valuation=expected,
            derivative_min_sum=sum(mins), triple_mass_by_precision=levels,
            cancellation_correction=sum(corrections),
            derivative_min_histogram=dict(sorted(Counter(mins).items())), all_passed=True)
    raise AssertionError('Increase the verified precision before claiming a valuation.')


def build_field(p, n):
    degree = order_mod(p, n)
    if degree <= 2:
        nonresidue = None if degree == 1 else next(a for a in range(2, p) if pow(a,(p-1)//2,p)==p-1)
        ring = Ring(p, 1, nonresidue)
        return ring, find_generator(ring, n)
    base, g = build_field(p, n//2)
    ring = QuadraticTower(base, g, n//2)
    return ring, (base.zero, base.one)


def ring_evaluate(f, x, ring):
    out = ring.zero
    for c in reversed(f):
        out = ring.add(ring.mul(out, x), ring.scalar(c))
    return out


def nonsplit_check(p, n, f, content):
    assert p % n != 1
    ring, g = build_field(p, n)
    assert ring.power(g, n) == ring.one and ring.power(g, n//2) == ring.scalar(-1)
    H = [ring.power(g, j) for j in range(n)]
    K, L = H[::2], H[1::2]
    assert len(set(H)) == n and set(K).isdisjoint(L)
    pairs = Counter(ring.add(a, b) for a in K for b in L)
    assert ring.zero not in pairs
    roots = [ring.power(ring.add(ring.one, H[j]), n) for j in range(1, n//2, 2)]
    counts = Counter(roots)
    assert max(pairs.values()) == max(counts.values())
    assert sum(c*c for c in pairs.values()) == n*sum(c*c for c in counts.values())
    a, b = derivative(f), derivative(derivative(f))
    for x in roots:
        assert ring_evaluate(f, x, ring) == ring.zero
        first = ring_evaluate(a, x, ring)
        assert (first == ring.zero) == (counts[x] >= 2)
        assert first != ring.zero or ring_evaluate(b, x, ring) != ring.zero
    assert content % p
    doubled = order_mod(p, n) == 2*order_mod(p, n//2)
    norm_checks = 0
    if doubled:
        assert max(pairs.values()) == 1
        branch = 'field degree doubles'
    else:
        assert ring.q == p*p
        epsilon = 1 if p % n == n-1 else -1
        assert epsilon == 1 or p % n == n//2-1
        assert all(ring.power(a,p+1) == ring.one for a in K)
        assert all(ring.power(b,p+1) == ring.scalar(epsilon) for b in L)
        norms = {x:ring.power(x,p+1) for x in pairs}
        for a0 in K:
            for b0 in L:
                x = ring.add(a0,b0)
                norm_x = norms[x]
                constant = ring.add(norm_x,ring.scalar(1-epsilon))
                # Norm(x)*a^2-x*(Norm(x)+1-epsilon)*a+x^2=0.
                expression = ring.add(ring.sub(ring.mul(norm_x,ring.mul(a0,a0)),
                    ring.mul(ring.mul(x,constant),a0)),ring.mul(x,x))
                assert expression == ring.zero
                norm_checks += 1
        assert max(pairs.values()) <= (1 if epsilon == 1 else 2)
        branch = 'both norms one' if epsilon == 1 else 'opposite norms'
    return dict(p=p,n=n,field_degree=order_mod(p,n),field_size=ring.q,
        branch=branch,maximum_fiber=max(pairs.values()),
        root_multiplicity_histogram=dict(sorted(Counter(counts.values()).items())),
        balanced_pair_checks=(n//2)**2,norm_quadratic_checks=norm_checks,
        primitive_content_nonzero_mod_p=True,all_passed=True)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--flint-path',default=str(ROOT/'tmp/arithmetic/python-flint-0.9.0-cp313'))
    args=parser.parse_args()
    sys.path.insert(0,args.flint_path)
    import flint
    assert flint.__version__ == '0.9.0'
    prior=json.loads((ROOT/'results/parallel48_collision_eliminants_2026_09_06.json').read_text())
    orders,polynomials,contents,padic=[],{}, {},[]
    for previous in prior['orders']:
        n=previous['n'];f=primitive_polynomial(n)
        coefficients,content,coarse,calls=pencil(f,flint)
        assert coarse == int(previous['gcd_hex'],16)
        polynomials[n],contents[n]=f,content
        for case in previous['fields']:
            padic.append(padic_check(f,content,n,case))
        orders.append(dict(n=n,content_hex=hex(content),coarse_gcd_hex=hex(coarse),
            content_equals_coarse_gcd=content==coarse,
            pencil_coefficient_hex=[hex(c) for c in coefficients],resultant_evaluations=calls))
        print(json.dumps({'n':n,'content_equals_coarse':content==coarse,'padic_cases':len(previous['fields'])}),flush=True)
    polynomials[4],contents[4]=[4,1],1
    examples=[]
    for p,n in [(3,4),(3,8),(7,8),(7,16),(11,4),(11,8),(23,8),(23,16),
                (31,32),(31,64),(47,16),(47,32),(127,128),(191,128),
                (5,8),(13,8),(17,32),(3,16),(5,16),(7,32),(3,32)]:
        examples.append(nonsplit_check(p,n,polynomials[n],contents[n]))
    toy1,content1,coarse1,_=pencil([0,-20,29,-10,1],flint)
    assert content1==16 and coarse1==400 and content1%5 and coarse1%5==0
    toy2,content2,coarse2,_=pencil([0,175,-40,1],flint)
    assert content2==2500 and valuation(content2,5)==4
    toy_levels=[]
    for precision in [1,2]:
        mult=Counter(x%5**precision for x in [0,5,35])
        toy_levels.append(sum(c*(c-2) for c in mult.values() if c>=3))
    assert toy_levels==[3,0]
    out=dict(scope='Exact critical-content valuation and nonsplit prime-support theorem; no improved uniform estimate at the splitting target primes.',
        orders=orders,padic_checks=padic,nonsplit_checks=examples,
        toy_refinement=dict(polynomial=[0,-20,29,-10,1],p=5,content=content1,coarse_gcd=coarse1,
            pencil_coefficients=toy1,actual_subgroup=False),
        toy_cancellation=dict(polynomial=[0,175,-40,1],p=5,content=content2,content_valuation=4,
            triple_mass_by_precision=toy_levels,pencil_coefficients=toy2,actual_subgroup=False),
        total_padic_cases=len(padic),nonsplit_field_cases=len(examples),
        prior_certificate_sha256=sha256((ROOT/'results/parallel48_collision_eliminants_2026_09_06.json').read_bytes()).hexdigest(),
        python_flint_version=flint.__version__,all_passed=True)
    path=ROOT/'results/parallel49_critical_content_2026_09_06.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'all_passed':True,'output':str(path)}))


if __name__ == '__main__':
    main()
