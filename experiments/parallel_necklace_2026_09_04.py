#!/usr/bin/env python3
"""Integer checks of projective label transfer and the three-label wheel identity."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from paley_exact import character
from localized_necklace_identities import convolution_moments
from planar_necklace_reductions import graph_sum

ROOT = Path(__file__).resolve().parents[1]
PI = str.maketrans('ABC', 'ACB')


def trace(a):
    return sum(int(x) for x in a.diagonal())


def matrices(p):
    chi = character(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    a = np.array(chi, dtype=object)
    b = np.array([chi[(x-1) % p] for x in range(p)], dtype=object)
    diag = dict(A=a, B=b, C=a*b)
    return chi, diag, {label: d[:, None]*s for label, d in diag.items()}


def direct_trace(word, cs, restrict=False):
    size = len(cs['A']) - int(restrict)
    result = np.eye(size, dtype=object)
    for label in word:
        a = cs[label][1:, 1:] if restrict else cs[label]
        result = result @ a
    return trace(result)


def literal(p, chi, diag, word, restrict=False):
    total = 0
    domain = range(int(restrict), p)
    for xs in product(domain, repeat=len(word)):
        value = 1
        for i, label in enumerate(word):
            value *= int(diag[label][xs[i]]) * chi[(xs[i]-xs[(i+1) % len(xs)]) % p]
        total += value
    return total


def all_traces(cs, max_k, restrict=False):
    size = len(cs['A']) - int(restrict)
    values = {}

    def visit(word, matrix):
        if word:
            values[word] = trace(matrix)
        if len(word) == max_k:
            return
        for label in 'ABC':
            a = cs[label][1:, 1:] if restrict else cs[label]
            visit(word+label, matrix @ a)
    visit('', np.eye(size, dtype=object))
    return values


def field_case(p):
    chi, diag, cs = matrices(p)
    # Independently check the conjugations on F_p*, where inversion is defined.
    s = np.array([[chi[(x-y) % p] for y in range(1,p)] for x in range(1,p)], dtype=object)
    inversion = np.zeros((p-1,p-1),dtype=object)
    for x in range(1,p):
        inversion[x-1,pow(x,-1,p)-1]=1
    quadratic = np.diag(diag['A'][1:])
    assert np.array_equal(inversion @ s @ inversion, quadratic @ s @ quadratic)
    for label in 'ABC':
        other=label.translate(PI)
        assert np.array_equal(inversion @ np.diag(diag[label][1:]) @ inversion,
                              np.diag(diag[other][1:]))
        assert np.array_equal(inversion @ cs[label][1:,1:] @ inversion,
                              quadratic @ cs[other][1:,1:] @ quadratic)
    t = convolution_moments(p, chi, 14)
    max_k = 6 if p <= 17 else 4
    full = all_traces(cs, max_k)
    restricted = all_traces(cs, max_k, restrict=True)
    for w, value in full.items():
        k = len(w)
        assert restricted[w] == restricted[w.translate(PI)]
        assert (value-restricted[w])**2 <= w.count('B')**2*p**k
        assert (value-full[w.translate(PI)])**2 <= (w.count('B')+w.count('C'))**2*p**k
        assert value == full[w.translate(str.maketrans('ABC','BAC'))]
    pair_values = {}
    for letters in product('ABC', repeat=6):
        w = ''.join(letters)
        if len(set(w)) <= 2:
            value = full[w] if w in full else direct_trace(w, cs)
            assert value**2 <= 42**2*p**7
            pair_values[w] = value
    assert len(pair_values) == 189
    blocks = 0
    for a, b in (('A','B'),('A','C'),('B','C')):
        for r in range(1,8):
            for s in range(1,8):
                k=r+s
                value = direct_trace(a*r+b*s, cs)
                assert value**2 <= (k*k+k)**2*p**k
                if b == 'C':
                    # Check the sharper transfer error against the exact binary formula.
                    reference = p*t[r]*t[s]-t[k]
                    assert (value-reference)**2 <= s*s*p**k
                blocks += 1
    triples = {}
    for letters in product('ABC', repeat=3):
        if len(set(letters)) == 3:
            w=''.join(letters)
            value=full[w]
            assert value == t[4]
            triples[w]=value
    literals = 0
    if p <= 13:
        for k in (1,2,3):
            for letters in product('ABC',repeat=k):
                w=''.join(letters)
                assert literal(p,chi,diag,w) == full[w]
                assert literal(p,chi,diag,w,restrict=True) == restricted[w]
                literals += 2
        # H is summed from its edges, not via either claimed formula.
        # a,b,x,y,z = 0,1,2,3,4.
        edges=[(2,3),(3,4),(2,4),(0,2),(0,4),(1,3),(1,4),(0,1)]
        wheel=graph_sum(p,chi,range(5),edges)
        assert wheel == p*(p-1)*full['ABC'] == p*(p-1)*t[4]
    else:
        wheel=None
    return dict(p=p, operator_conjugation_checks=7,
                transfer_checks=len(full), restricted_inversion_checks=len(full),
                full_restricted_rank_one_checks=len(full), reflection_checks=len(full),
                two_label_length_six_checks=len(pair_values), two_block_checks=blocks,
                triple_checks=len(triples), literal_coordinate_sums=literals,
                wheel_partition=wheel, triples=triples,
                two_label_length_six_values=pair_values, t0_through_t14=t)


def main():
    cases=[]
    for p in (5,13,17,29,41):
        cases.append(field_case(p))
        print(json.dumps({'p':p,'status':'all exact checks passed'}),flush=True)
    inputs=['research/parallel-necklace-2026-09-04.md',
            'experiments/parallel_necklace_2026_09_04.py',
            'experiments/planar_necklace_reductions.py',
            'experiments/localized_necklace_identities.py',
            'sources/kunisky-2303.16475v1.html',
            'sources/lu-zheng-zheng-1305.3405v3.html']
    result=dict(status='Uniform proofs in accompanying note; finite exact checks supplement them.',
                arithmetic='Python integers; NumPy object arrays; no floating-point comparisons.',
                cases=cases, input_sha256={path:sha256((ROOT/path).read_bytes()).hexdigest() for path in inputs})
    path=ROOT/'results/parallel_necklace_2026_09_04.json'
    path.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'written':str(path)}),flush=True)


if __name__=='__main__':
    main()
