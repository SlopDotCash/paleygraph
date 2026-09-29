#!/usr/bin/env python3
"""Exact pole-inclusive Schur-kernel and AABBCC checks; no open-class bounds."""
from collections import Counter
from hashlib import sha256
from itertools import permutations, product
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, direct_trace, literal, trace
from parallel2_necklace_2026_09_04 import positive_principal_minors
from planar_necklace_reductions import graph_sum

ROOT = Path(__file__).resolve().parents[1]


def dihedral(word):
    return {w[i:]+w[:i] for w in (word, word[::-1]) for i in range(len(w))}


def orbit(word):
    return {v for labels in permutations('ABC')
            for v in dihedral(word.translate(str.maketrans('ABC', ''.join(labels))))}


def catalog():
    words = {''.join(w) for w in product('ABC', repeat=6)}
    balanced = {w for w in words if all(w.count(a) == 2 for a in 'ABC')}
    classes = []
    remaining = set(balanced)
    while remaining:
        representative = min(remaining)
        members = orbit(representative)
        assert members <= remaining
        classes.append(dict(representative=representative, count=len(members),
                            status='proved' if representative == 'AABBCC' else 'open',
                            words=sorted(members)))
        remaining -= members
    assert [(c['representative'], c['count']) for c in classes] == [
        ('AABBCC', 12), ('AABCBC', 36), ('AABCCB', 18),
        ('ABACBC', 18), ('ABCABC', 6)]
    blocks = orbit('AABBCC')
    assert blocks == dihedral('AABBCC')
    multiplicities = Counter(tuple(sorted(Counter(w).values(), reverse=True)) for w in words)
    old = {w for w in words if len(set(w)) <= 2 or min(Counter(w).values()) == 1}
    # Prior two-label, singleton-pair and three-gap results give exactly these words.
    assert len(old) == 639 and words-old == balanced
    assert multiplicities[(4, 1, 1)] == 90 and multiplicities[(3, 2, 1)] == 360
    assert sum(n for m, n in multiplicities.items() if len(m) <= 2) == 189
    assert len(old | blocks) == 651 and len(words-old-blocks) == 78
    return dict(total_words=len(words), prior_coverage=len(old), new_coverage=len(blocks),
                total_coverage=len(old | blocks), remaining_count=len(words-old-blocks),
                balanced_classes=classes, remaining_words=sorted(words-old-blocks))


def necklace_edges(word):
    edges = [(i, (i+1) % 6) for i in range(6)]
    edges += [(i, 6) for i, a in enumerate(word) if a in 'AC']
    edges += [(i, 7) for i, a in enumerate(word) if a in 'BC']
    assert len(edges) == 14
    return edges


def literal_path(p, chi, q, pattern):
    """Four free vertices after fixing the final two C vertices to 1 and 0."""
    pairs = [tuple(i for i, a in enumerate(pattern) if a == letter) for letter in 'AB']
    total = 0
    for xs in product(range(p), repeat=4):
        value = chi[xs[0]]*chi[(xs[3]-1) % p]
        for i in range(3):
            value *= chi[(xs[i]-xs[i+1]) % p]
        for i, j in pairs:
            value *= int(q[xs[i], xs[j]])
        total += value
    return total


def field_case(p, block_words):
    chi, diag, cs = matrices(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    identity = np.eye(p, dtype=object)
    assert np.array_equal(s@s, p*identity-np.ones((p, p), dtype=object))
    da, db, dc = (np.diag(diag[a]) for a in 'ABC')
    k2 = s@da@s
    q = s@dc@s
    wa, wc = s*k2, s*q
    ones = np.ones(p-1, dtype=object)
    assert np.array_equal(wa[0, 1:], -ones)
    assert np.array_equal(wa[1:, 1:]@ones, 2*ones)
    assert np.array_equal(q[0], np.array([p*int(y == 1)-1-chi[y] for y in range(p)], dtype=object))
    augmented = np.zeros((p+1, p+1), dtype=object)
    augmented[:p, :p] = wa
    border = np.array([p*int(y == 0)-1 for y in range(p)], dtype=object)
    augmented[p, :p] = augmented[:p, p] = border
    e0 = np.zeros(p+1, dtype=object); e0[0] = 1
    einf = np.zeros(p+1, dtype=object); einf[p] = 1
    v = np.zeros(p+1, dtype=object); v[1:p] = 1
    for vec, value in ((einf-e0, -(p-1)), (einf+e0+v, 0),
                       ((p-1)*(einf+e0)-2*v, p+1)):
        assert np.array_equal(augmented@vec, value*vec)
    # A basis of the remaining invariant subspace, with both border coordinates zero.
    for j in range(2, p):
        vector = np.zeros(p+1, dtype=object); vector[j] = 1; vector[1] = -1
        out = augmented@vector
        assert out[0] == out[p] == sum(out[1:p]) == 0
        assert np.array_equal(out[1:p], wa[1:, 1:]@vector[1:p])
    index = [p if x == 0 else (pow(x, -1, p)-1) % p for x in range(p)]
    assert set(index) == set(range(p+1))-{p-1}
    assert np.array_equal(wc+s, augmented[np.ix_(index, index)])
    # Exact positive definite certificates, not numerical eigenvalue comparisons.
    aug_cert = positive_principal_minors(4*p*p*np.eye(p+1, dtype=object)-augmented@augmented)
    wc_cert = positive_principal_minors(25*p*p*identity-4*wc@wc)
    # Check the general-anchor covariance in square and nonsquare scales.
    nonsquare = next(x for x in range(1, p) if chi[x] == -1)
    anchor_pairs = sorted({(0, 1), (1, 0), (2 % p, (2+nonsquare) % p),
                           (3 % p, (3+nonsquare) % p)})
    for a, b in anchor_pairs:
        delta = (b-a) % p
        relabel = [(a+delta*x) % p for x in range(p)]
        dab = np.diag([chi[(x-a) % p]*chi[(x-b) % p] for x in range(p)])
        qab = s@dab@s
        wab = s*qab
        assert np.array_equal(qab[np.ix_(relabel, relabel)], q)
        assert np.array_equal((wab+s)[np.ix_(relabel, relabel)], chi[delta]*(wc+s))
    f0, f1 = s[:, 0], s[:, 1]
    ca4 = np.linalg.matrix_power(cs['A'], 4)
    assert trace(ca4) % (p-1) == 0
    t4 = trace(ca4)//(p-1)
    assert abs(t4) <= p*p
    pzero = identity.copy(); pzero[0, 0] = 0
    nzero = trace(ca4@pzero@s@pzero@s)
    cstar, sstar = cs['A'][1:, 1:], s[1:, 1:]
    u = diag['A'][1:]
    assert np.array_equal(cstar@ones, -ones) and np.array_equal(cstar@u, -u)
    assert np.array_equal(sstar@sstar, p*np.eye(p-1, dtype=object)-np.outer(ones, ones)-np.outer(u, u))
    assert nzero == trace(np.linalg.matrix_power(cstar, 4)@sstar@sstar)
    assert nzero == (p-1)*(p*t4-2)
    contractions = {
        'AABB': int(f0@wc@s@wc@f1),
        'ABAB': sum(int(z) for z in (s*(s@da@q)*(q@db@s)).flat),
        'ABBA': sum(int(z) for z in (np.outer(f0, f1)*q*(s@wc@s)).flat),
    }
    path_checks = {}
    for pattern, value in contractions.items():
        actual = direct_trace(pattern+'CC', cs)
        assert (p-1)*actual == (p-1)*value-nzero
        if p <= 13:
            path_checks[pattern] = literal_path(p, chi, q, pattern)
            assert path_checks[pattern] == value
    n = direct_trace('AABBCC', cs)
    assert n == contractions['AABB']-p*t4+2
    assert n*n <= 49*p**7
    block_values = {w: direct_trace(w, cs) for w in block_words}
    assert all(value == n and value*value <= 49*p**7 for value in block_values.values())
    graph_checks, original_literals = {}, {}
    if p == 5:
        for pattern, value in contractions.items():
            word = pattern+'CC'
            edges = necklace_edges(word)
            z = graph_sum(p, chi, range(8), edges)
            fixed_pair = graph_sum(p, chi, range(8), edges, fixed={4: 1, 5: 0})
            distinct_anchors = graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 1})
            same_anchors = graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 0})
            assert fixed_pair == value and same_anchors == nzero
            assert distinct_anchors == direct_trace(word, cs)
            assert z == p*(p-1)*distinct_anchors+p*same_anchors == p*(p-1)*fixed_pair
            graph_checks[word] = dict(full_graph=z, fixed_C_pair=fixed_pair,
                                     distinct_anchors=distinct_anchors, coincident_anchors=same_anchors)
        for word in ('AABBCC', 'AABCBC', 'AABCCB', 'ABACBC', 'ABCABC'):
            original_literals[word] = literal(p, chi, diag, word)
            assert original_literals[word] == direct_trace(word, cs)
    return dict(p=p, full_projective_compression_entries=p*p, augmented_eigenvectors=3,
                zero_sum_basis_vectors=p-2, pole_row_entries=p,
                anchor_covariance_pairs=len(anchor_pairs), norm_certificates=2,
                augmented_norm_leading_minors=aug_cert, schur_norm_leading_minors=wc_cert,
                t4=t4, coincident_anchor_term=nzero, N_AABBCC=n,
                block_orbit_values=block_values, exact_contractions=contractions,
                unresolved_contractions_are_bounds=False, literal_four_vertex_paths=path_checks,
                graph_partition_checks=graph_checks, literal_original_necklaces=original_literals)


def main():
    coverage = catalog()
    cases = []
    for p in (5, 13, 17, 29, 41):
        cases.append(field_case(p, sorted(dihedral('AABBCC'))))
        print(json.dumps({'p': p, 'status': 'all exact operator, identity and bounded-family checks passed'}), flush=True)
    files = ['research/parallel6-necklace-2026-09-04.md',
             'experiments/parallel6_necklace_2026_09_04.py',
             'research/parallel2-necklace-2026-09-04.md',
             'research/parallel4-necklace-2026-09-04.md',
             'research/parallel5-necklace-2026-09-04.md',
             'experiments/parallel_necklace_2026_09_04.py',
             'experiments/parallel2_necklace_2026_09_04.py',
             'experiments/planar_necklace_reductions.py',
             'sources/katz-finite-field-mellin.pdf']
    result = dict(status='Uniform Schur-kernel bound and the twelve-word AABBCC class proved; four balanced classes remain open.',
                  arithmetic='Python integers, NumPy object arrays, exact Bareiss positivity certificates.',
                  spectral_scope='No localized adjacency spectral-edge or signed aggregate claim.',
                  primary_source=dict(url='https://web.math.princeton.edu/~nmk/mellin186.pdf',
                      archive='sources/katz-finite-field-mellin.pdf',
                      locations=['Theorem 15.1, printed p. 52', 'Corollary 4.2, printed p. 22']),
                  finite_word_coverage=coverage, cases=cases,
                  input_sha256={f: sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    output = ROOT/'results/parallel6_necklace_2026_09_04.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'written': str(output)}), flush=True)


if __name__ == '__main__':
    main()
