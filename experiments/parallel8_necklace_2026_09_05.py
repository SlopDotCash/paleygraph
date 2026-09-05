#!/usr/bin/env python3
"""Exact nested/direct contractions and rank bookkeeping; uniform proof is prose."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, direct_trace, literal, trace
from parallel6_necklace_2026_09_04 import orbit, dihedral, necklace_edges, literal_path
from parallel7_necklace_2026_09_04 import catalog as previous_catalog
from planar_necklace_reductions import graph_sum

ROOT = Path(__file__).resolve().parents[1]


def surd_sign(a, b, p):
    """Sign of a+b*sqrt(p), using integers only; p is a positive nonsquare."""
    a, b = int(a), int(b)
    if a == 0:
        return (b > 0)-(b < 0)
    if b == 0:
        return (a > 0)-(a < 0)
    if a*b > 0:
        return (a > 0)-(a < 0)
    delta = a*a-p*b*b
    return ((delta > 0)-(delta < 0))*((a > 0)-(a < 0))


def rank(rep):
    return sum(size for char, size in rep)


def inv(rep):
    return sum(char == 0 for char, size in rep)


def twist(rep):
    return sorted((1-char, size) for char, size in rep)


def quotient(rep):
    return sorted((char, size-(char == 0)) for char, size in rep
                  if char != 0 or size > 1)


def reconstruct(quot, total):
    """Recover tame monodromy from its quotient by invariants and total rank."""
    rep = [(char, size+(char == 0)) for char, size in quot]
    assert rank(rep) <= total
    rep += [(0, 1)]*(total-rank(rep))
    return sorted(rep)


def mc(local):
    finite = (0, 1, 'y')
    codims = {str(x): rank(local[x])-inv(local[x]) for x in finite}
    ambient = sum(codims.values())
    output_rank = ambient-inv(twist(local['infinity']))
    output = {x: reconstruct(twist(quotient(local[x])), output_rank) for x in finite}
    aux = reconstruct(local['infinity'], ambient)
    output['infinity'] = quotient(twist(aux))
    assert rank(output['infinity']) == output_rank
    return output, dict(finite_codimensions=codims, auxiliary_rank=ambient,
                       twisted_infinity_invariants=inv(twist(local['infinity'])),
                       output_rank=output_rank, local_types=output)


def monodromy_checks():
    cases = []
    for a, b in product((0, 1), repeat=2):
        g = {b: [(1, 2)], 1-b: [(1, 1), (1, 1)],
             'y': [(0, 2)], 'infinity': [(1, 2)]}
        middle, first = mc(g)
        assert first['output_rank'] == 4
        assert middle[b] == [(0, 1), (0, 3)]
        assert middle[1-b] == [(0, 2), (0, 2)]
        assert middle['y'] == [(0, 1)]*3+[(1, 1)]
        assert middle['infinity'] == [(0, 1)]+[(1, 1)]*3
        second_input = dict(middle)
        second_input[a] = twist(middle[a])
        second_input['infinity'] = twist(middle['infinity'])
        assert inv(second_input[a]) == 0  # Zero mask is the actual new stalk.
        final, second = mc(second_input)
        assert second['output_rank'] == 6
        assert second['finite_codimensions'] == {str(a): 4, str(1-a): 2, 'y': 1}
        assert second['twisted_infinity_invariants'] == 1
        assert final['infinity'] == [(1, 2)]*3
        assert final['y'] == [(0, 1)]*4+[(0, 2)]
        cases.append(dict(a=a, b=b, first=first, second=second))
    return cases


def catalog():
    previous = previous_catalog()
    reps = ('ABBACC', 'ABCBAC', 'ABCABC')
    classes = {w: sorted(orbit(w)) for w in reps}
    assert [len(classes[w]) for w in reps] == [18, 18, 6]
    assert orbit('ABBACC') == orbit('AABCCB')
    assert orbit('ABCBAC') == orbit('ABACBC')
    assert orbit('ABCABC') == dihedral('ABCABC')
    new = set().union(*(set(c) for c in classes.values()))
    old_remaining = set().union(*(set(c) for c in previous['remaining_classes'].values()))
    assert len(new) == 42 and new == old_remaining
    assert new <= {''.join(w) for w in product('ABC', repeat=6)}
    return dict(previous_coverage=previous['total_coverage'], new_count=42,
                total_coverage=729, remaining_count=0, new_classes=classes)


def field_case(p, classes):
    chi, d, cs = matrices(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    eye = np.eye(p, dtype=object)
    da = {a: np.diag([chi[(x-a) % p] for x in range(p)]) for a in (0, 1)}
    dc = np.diag(d['C'])
    q = s@dc@s
    assert np.array_equal(s@s, p*eye-np.ones((p, p), dtype=object))
    kernels = {}
    finite_stalks = pointwise = compact_vectors = core_bounds = formal_coeffs = 0
    fiber_bounds = exceptional = 0
    fiber_sums = {}
    for a, b in product((0, 1), repeat=2):
        kb = s@da[b]@s
        v = q@da[b]@s
        ell = s@da[a]@v
        kernels[a, b] = ell
        for y in range(2, p):
            e = -q[:, y]-1
            assert (e[0], e[1], e[y]) == (chi[y], chi[(1-y) % p], chi[y*(y-1) % p])
            assert all(int(t)**2 <= 4*p for t in e)
            g = -d['C']*kb[:, y]
            assert g[0] == g[1] == 0
            finite_stalks += 5
            assert np.array_equal(-s@g, v[:, y])
            spike = np.array([p*int(x == a)-1 for x in range(p)], dtype=object)
            assert np.array_equal(s@da[a]@np.ones(p, dtype=object), spike)
            assert np.array_equal(-s@da[a]@v[:, y], -ell[:, y])
            compact_vectors += 3
            for x in range(p):
                # Every complex tau,kappa with these modulus bounds is allowed.
                assert surd_sign(abs(ell[x, y])+p, abs(spike[x])-6*p, p) <= 0
                pointwise += 1
            for c in (0, 1, y):
                h = np.array([chi[(x-c) % p] for x in range(p)], dtype=object)
                assert sum(h) == 0
                actual = int(h@(q[:, y]*ell[:, y]))
                open_x = [x for x in range(p) if x not in (0, 1, y)]
                core_zero = sum(int(-h[x]*e[x]*ell[x, y]) for x in open_x)
                core_slope = sum(int(h[x]*e[x]) for x in open_x)
                # Supremum of |core_zero-(tau+kappa)*core_slope|.
                assert surd_sign(abs(core_zero)+p*abs(core_slope),
                                 abs(core_slope)-24*p*p, p) <= 0
                core_bounds += 1
                # Coefficients of 1,tau,kappa in the full expansion (9).
                he = int(h@e)
                coeff_const = int(h@(e*(-ell[:, y]))+h@(-ell[:, y]))
                coeff_tau = int(h@(e*spike)+he+h@spike-p*h[a]*(e[a]+1))
                coeff_kappa = -he+he-int(sum(h))
                assert (coeff_const, coeff_tau, coeff_kappa) == (actual, 0, 0)
                formal_coeffs += 3
                assert surd_sign(abs(actual)-26*p*p, -(32*p*p+2*p), p) <= 0
                fiber_bounds += 1
                fiber_sums[a, b, y, c] = actual
        for y in (0, 1):
            for c in (0, 1, y):
                value = sum(int(chi[(x-c) % p]*q[x, y]*ell[x, y]) for x in range(p))
                assert abs(value) <= p**3
                fiber_sums[a, b, y, c] = value
                exceptional += 1
    f0, f1 = s[:, 0], s[:, 1]
    wc = s*q
    nested_original = int(sum((np.outer(f0, f1)*q*(s@wc@s)).flat))
    nested_rearranged = sum(fiber_sums[0, 1, y, y] for y in range(p))
    assert nested_original == nested_rearranged
    ca4 = np.linalg.matrix_power(cs['A'], 4)
    assert trace(ca4) % (p-1) == 0
    t4 = trace(ca4)//(p-1)
    assert abs(t4) <= p*p
    pzero = eye.copy()
    pzero[0, 0] = 0
    nzero = trace(ca4@pzero@s@pzero@s)
    assert nzero == (p-1)*(p*t4-2)
    direct = {
        'ABBACC': nested_rearranged-p*t4+2,
        'ABCBAC': sum(chi[y]*fiber_sums[1, 1, y, 0] for y in range(p)),
        'ABCABC': sum(chi[(y-1) % p]*fiber_sums[1, 0, y, 0] for y in range(p)),
    }
    constants = {'ABBACC': 46, 'ABCBAC': 45, 'ABCABC': 45}
    orbit_constants = {'ABBACC': 48, 'ABCBAC': 47, 'ABCABC': 45}
    orbits = {}
    for representative, value in direct.items():
        assert value == direct_trace(representative, cs)
        assert value*value <= constants[representative]**2*p**7
        vals = {w: direct_trace(w, cs) for w in classes[representative]}
        assert all(abs(v-value) <= 4*p**3 for v in vals.values())
        assert all(v*v <= orbit_constants[representative]**2*p**7 for v in vals.values())
        if representative == 'ABCABC':
            assert set(vals.values()) == {value}
        orbits[representative] = vals
    graphs = {}
    literals = {}
    if p == 5:
        edges = necklace_edges('ABBACC')
        graphs = {
            'full': graph_sum(p, chi, range(8), edges),
            'distinct_anchors': graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 1}),
            'coincident_anchors': graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 0}),
            'fixed_C_pair': graph_sum(p, chi, range(8), edges, fixed={4: 1, 5: 0}),
        }
        assert graphs['full'] == p*(p-1)*direct['ABBACC']+p*nzero
        assert graphs['full'] == p*(p-1)*nested_original
        assert graphs['distinct_anchors'] == direct['ABBACC']
        assert graphs['coincident_anchors'] == nzero
        assert graphs['fixed_C_pair'] == literal_path(p, chi, q, 'ABBA') == nested_original
        for w in direct:
            literals[w] = literal(p, chi, d, w)
            assert literals[w] == direct[w]
    # Normalized constant comparisons attain their largest error at p=5.
    assert surd_sign(28*p**3, 2*p*p-13*p**3, p) < 0
    assert surd_sign(29*p**3+2, 2*p*p-14*p**3, p) < 0
    return dict(p=p, finite_stalk_checks=finite_stalks,
                compact_vector_checks=compact_vectors, interval_pointwise_checks=pointwise,
                interval_core_bounds=core_bounds, formal_coefficient_checks=formal_coeffs,
                generic_fiber_bounds=fiber_bounds, exceptional_fiber_bounds=exceptional,
                direct_identities=direct, nested_contraction=nested_original,
                t4=t4, coincident_anchor_term=nzero, orbit_values=orbits,
                literal_graphs=graphs, literal_necklaces=literals)


def main():
    ranks = monodromy_checks()
    coverage = catalog()
    cases = []
    for p in (5, 13, 17, 29, 41):
        cases.append(field_case(p, coverage['new_classes']))
        print(json.dumps(dict(p=p, status='all exact checks passed')), flush=True)
    paths = [
        'research/parallel8-necklace-2026-09-05.md',
        'experiments/parallel8_necklace_2026_09_05.py',
        'research/parallel7-necklace-2026-09-04.md',
        'research/parallel6-necklace-2026-09-04.md',
        'research/parallel-necklace-2026-09-04.md',
        'experiments/parallel7_necklace_2026_09_04.py',
        'experiments/parallel6_necklace_2026_09_04.py',
        'experiments/parallel_necklace_2026_09_04.py',
        'experiments/planar_necklace_reductions.py',
        'sources/katz-rigid-local-systems.pdf',
        'sources/katz-gauss-kloosterman-monodromy.pdf',
    ]
    counts = {k: sum(c[k] for c in cases) for k in (
        'finite_stalk_checks', 'compact_vector_checks', 'interval_pointwise_checks',
        'interval_core_bounds', 'formal_coefficient_checks',
        'generic_fiber_bounds', 'exceptional_fiber_bounds')}
    counts.update(monodromy_cases=len(ranks), direct_trace_identities=3*len(cases),
                  new_orbit_values=sum(len(v) for c in cases for v in c['orbit_values'].values()),
                  original_literal_necklaces=sum(len(c['literal_necklaces']) for c in cases))
    result = dict(
        status='All exact checks passed. Uniform proof is in the note; full Paley goal remains open.',
        arithmetic='Integers and exact signs in Q(sqrt(p)); no floating-point acceptance.',
        scope='All fixed length-six degree-two words; no growing-depth aggregate or Paley proof.',
        infinity_correction_policy='No guessed kappa. Finite bounds allow all complex tau,kappa with |tau|<=sqrt(p), |kappa|<=p.',
        monodromy=ranks, coverage=coverage, counts=counts, cases=cases,
        input_sha256={p: sha256((ROOT/p).read_bytes()).hexdigest() for p in paths})
    out = ROOT/'results/parallel8_necklace_2026_09_05.json'
    out.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps(dict(written=str(out), counts=counts)), flush=True)


if __name__ == '__main__':
    main()
