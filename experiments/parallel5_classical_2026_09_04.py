#!/usr/bin/env python3
"""Exact bounded checks: one-size transfers and B_h forced-row witnesses."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb, isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def prime(p):
    return p >= 2 and all(p % d for d in range(2, isqrt(p) + 1))


def character(p):
    assert prime(p)
    values = [-1] * p
    values[0] = 0
    for a in range(1, p):
        values[a*a % p] = 1
    return values


def rows(p, chi, B):
    return [sum(chi[(a-b) % p] for b in B) for a in range(p)]


def subset_transfer():
    records = []
    for p in (5, 7, 11, 13):
        chi = character(p)
        sets = [B for n in range(1, min(p, 4)+1)
                for B in combinations(range(p), n)]
        cache = {B: {r: sum(v**(2*r) for v in rows(p, chi, B))
                     for r in (2, 3)} for B in sets}
        count = 0
        for B in sets:
            n = len(B)
            for k in range(1, n+1):
                for r in (2, 3):
                    total = sum(cache[C][r] for C in combinations(B, k))
                    assert k**(2*r)*comb(n, k)*cache[B][r] <= n**(2*r)*total
                    count += 1
        records.append(dict(p=p, sets=len(sets), exact_average_checks=count))
    return records


def rectangle_transfer():
    records = []
    for p in (5, 7):
        chi = character(p)
        sets = [A for n in range(1, p+1) for A in combinations(range(p), n)]
        row_cache = {B: rows(p, chi, B) for B in sets}
        maxima = {(k, r): max(sum(v**(2*r) for v in row_cache[C])
                              for C in sets if len(C) == k)
                  for k in range(1, 4) for r in (2, 3)}
        count = 0
        for A in sets:
            m = len(A)
            for B in sets:
                n = len(B)
                total = sum(row_cache[B][a] for a in A)
                transposed = sum(row_cache[A][b] for b in B)
                assert transposed == chi[p-1]*total
                for k in range(1, min(3, n)+1):
                    for r in (2, 3):
                        lhs = abs(total)**(2*r)*k**(2*r)
                        rhs = m**(2*r-1)*n**(2*r)*maxima[k, r]
                        assert lhs <= rhs
                        count += 1
        records.append(dict(p=p, rectangles=len(sets)**2,
                            exact_transfer_checks=count,
                            moment_maxima=[dict(k=k, r=r, maximum=L)
                                           for (k, r), L in maxima.items()]))
    return records


def energies(p, B, h):
    output = []
    for j in range(1, h+1):
        sums = Counter()
        permutations = Counter()
        for t in product(B, repeat=j):
            sums[sum(t) % p] += 1
            permutations[tuple(sorted(t))] += 1
        actual = sum(v*v for v in sums.values())
        minimum = sum(v*v for v in permutations.values())
        assert actual == minimum
        if j == 2:
            assert actual == 2*len(B)**2-len(B)
        if j == 3:
            assert actual == 6*len(B)**3-9*len(B)**2+4*len(B)
        output.append(dict(j=j, actual=actual, permutation_minimum=minimum,
                           ordered_tuples=len(B)**j))
    return output


def structured_witness(h, q, p):
    assert prime(q) and q > h and prime(p)
    base = h*q
    Q = tuple(sum(pow(t, j, q)*base**(j-1) for j in range(1, h+1))
              for t in range(q))
    assert len(set(Q)) == q and h*max(Q) < p
    # For finite examples, this exact no-wrap condition suffices. The note's
    # more conservative p >= 2 h^(h+1) q^h guarantees it uniformly.
    container_energies = energies(p, Q, h)
    chi = character(p)
    target = max(h+1, q//4)
    current = Q
    U = []
    stages = []
    for step in range(8):
        candidates = [(sum(chi[(a-b) % p] == 1 for b in current), a)
                      for a in range(p) if a not in U]
        best, a = max(candidates)
        if best < target:
            break
        expected_numerator = len(current)*((p-1)//2-step)
        assert best*(p-step) >= expected_numerator
        U.append(a)
        current = tuple(b for b in current if chi[(a-b) % p] == 1)
        stages.append(dict(row=a, surviving_columns=len(current)))
    assert U and len(current) >= target
    B = current[:target]
    assert set(U).isdisjoint(B)
    F = rows(p, chi, B)
    assert all(F[a] == len(B) for a in U)
    moments = []
    for r in (h+1, h+2, h+3):
        actual = sum(v**(2*r) for v in F)
        lower = len(U)*len(B)**(2*r)
        assert actual >= lower
        ratio = Fraction(actual, p*len(B)**r)
        moments.append(dict(r=r, actual=actual, forced_row_lower_bound=lower,
                            gaussian_ratio_numerator=ratio.numerator,
                            gaussian_ratio_denominator=ratio.denominator))
    k = len(U)
    average = Fraction(q*comb((p-1)//2, k), comb(p, k))
    return dict(p=p, auxiliary_prime=q, h=h, base=base, container=Q,
                no_wrap_margin=p-h*max(Q), container_energies=container_energies,
                row_stages=stages, B=B, U=U,
                witness_energies=energies(p, B, h), moments=moments,
                container_row_average=[average.numerator, average.denominator])


def exponent_audit():
    records = []
    for r in range(2, 101):
        theta = Fraction(1, r+1)
        beta = Fraction(r-1, 2*(r+1))
        assert theta+beta == Fraction(1, 2)
        assert r*theta < 1
        for h in range(2, r):
            alpha = (Fraction(1, r)+Fraction(1, h))/2
            assert alpha*r > 1 and alpha < Fraction(1, h)
        records.append(dict(r=r, theta=str(theta), termwise_beta=str(beta)))
    examples = []
    for eps in (Fraction(3, 4), Fraction(1, 2), Fraction(1, 3), Fraction(1, 4),
                Fraction(1, 10), Fraction(1, 20), Fraction(1, 50), Fraction(1, 100)):
        r = (4*eps.denominator+eps.numerator-1)//eps.numerator
        theta = beta = Fraction(1, r+1)
        delta = eps/(8*r)
        assert theta+beta < eps/2
        assert 2*r*delta < eps-theta-beta
        assert Fraction(1, r) < eps and Fraction(1, r) < eps/2
        examples.append(dict(epsilon=str(eps), r=r, SS_beta=str(beta),
                             CS_alpha=str(Fraction(1, r)), delta=str(delta)))
    return dict(termwise_ranks=records, conditional_examples=examples)


def main():
    transfer = subset_transfer()
    rectangles = rectangle_transfer()
    witnesses = [structured_witness(*case) for case in
                 ((2, 5, 101), (2, 7, 127), (2, 11, 509),
                  (2, 31, 5003), (2, 61, 20011),
                  (3, 5, 3001), (3, 7, 10009), (3, 11, 131071))]
    algebra = exponent_audit()
    paths = ('research/parallel5-classical-2026-09-04.md',
             'experiments/parallel5_classical_2026_09_04.py')
    record = dict(status='Exact reductions and structured obstruction proved; SS, CS, LM and Paley remain unproved.',
                  arithmetic='Integer character sums and counts; rational exponents and ratios.',
                  subset_transfers=transfer, rectangle_transfers=rectangles,
                  structured_witnesses=witnesses, exponent_audit=algebra,
                  source_sha256={f: sha256((ROOT/f).read_bytes()).hexdigest() for f in paths})
    destination = ROOT/'results/parallel5_classical_2026_09_04.json'
    destination.write_text(json.dumps(record, indent=2)+'\n')
    print(json.dumps(dict(subset_average_checks=sum(x['exact_average_checks'] for x in transfer),
                          rectangle_transfer_checks=sum(x['exact_transfer_checks'] for x in rectangles),
                          structured_witnesses=len(witnesses),
                          energy_checks=sum(len(w['container_energies'])+len(w['witness_energies']) for w in witnesses),
                          moment_lower_checks=sum(len(w['moments']) for w in witnesses),
                          exponent_ranks=len(algebra['termwise_ranks']),
                          result=str(destination))))


if __name__ == '__main__':
    main()
