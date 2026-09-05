#!/usr/bin/env python3
"""Exact crossing checks, including raw traces, all boundaries, and F_(5^n)."""
from hashlib import sha256
from itertools import product
from pathlib import Path
import json

import numpy as np

from parallel_necklace_2026_09_04 import matrices, direct_trace, literal, trace
from parallel6_necklace_2026_09_04 import orbit, dihedral, necklace_edges
from planar_necklace_reductions import graph_sum

ROOT = Path(__file__).resolve().parents[1]


def catalog():
    words = {''.join(w) for w in product('ABC', repeat=6)}
    balanced = {w for w in words if all(w.count(a) == 2 for a in 'ABC')}
    old = (words-balanced) | orbit('AABBCC')
    new = orbit('ABABCC')
    assert new == orbit('AABCBC') and len(new) == 36 and not old & new
    remaining = words-old-new
    assert len(old) == 651 and len(old | new) == 687 and len(remaining) == 42
    classes = {w: sorted(orbit(w)) for w in ('AABCCB', 'ABACBC', 'ABCABC')}
    assert set().union(*(set(c) for c in classes.values())) == remaining
    return dict(total_words=729, prior_coverage=651, newly_covered=36, total_coverage=687,
                new_words=sorted(new), remaining_count=42, remaining_classes=classes)


def field_case(p, new_words):
    chi, d, cs = matrices(p)
    s = np.array([[chi[(x-y) % p] for y in range(p)] for x in range(p)], dtype=object)
    eye = np.eye(p, dtype=object)
    assert np.array_equal(s@s, p*eye-np.ones((p, p), dtype=object))
    q = s@np.diag(d['C'])@s
    kb = s@np.diag(d['B'])@s
    u_raw = s@np.diag(d['A'])@q
    v_raw = q@np.diag(d['B'])@s
    assert np.array_equal(v_raw, s@np.diag(d['C'])@kb)
    elliptic = [sum(chi[z]*chi[(z-1) % p]*chi[(z-t) % p] for z in range(p)) for t in range(p)]
    rows = {}; cross_ratio_checks = 0; input_stalk_checks = 0
    for y in range(2, p):
        for x in range(1, p):
            t = y*(1-x)*pow(x*(1-y) % p, -1, p) % p
            assert q[x, y]+1 == chi[x*(1-y) % p]*elliptic[t]
            cross_ratio_checks += 1
        for x in range(p):
            t = (x-1)*pow(y-1, -1, p) % p
            assert kb[x, y] == chi[y-1]*elliptic[t]
        f = -d['A']*(q[:, y]+1)
        g = -d['C']*kb[:, y]
        assert f[0] == g[0] == g[1] == 0
        input_stalk_checks += 3
        u = u_raw[:, y].copy(); u[0] += p; u -= 1
        v = v_raw[:, y]-1
        # Compact-convolution identities are checked from the original input traces.
        assert np.array_equal(-s@f, u)
        assert np.array_equal(-s@g, v+1)
        assert all(abs(int(z)) <= 4*p for z in u)
        assert all(abs(int(z)) <= 4*p for z in v)
        a = s[:, y]
        assert sum(a) == 0
        actual = sum(int(a[x]*u_raw[x, y]*v_raw[x, y]) for x in range(p))
        core = sum(int(a[x]*u[x]*v[x]) for x in range(p) if x not in (0, 1, y))
        restored = sum(int(a[x]*u[x]*v[x]) for x in (0, 1))
        linear_u = int(a@u); linear_v = int(a@v)
        pole = -p*chi[(-y) % p]*(int(v[0])+1)
        assert actual == core+restored+linear_u+linear_v+pole
        assert core*core <= 32**2*p**5
        assert abs(restored) <= 32*p*p
        assert abs(linear_u) <= 4*p*p and abs(linear_v) <= 4*p*p
        assert abs(pole) <= 4*p*p+p
        residual = max(abs(actual)-40*p*p, 0)
        assert residual*residual <= p*(36*p*p+p)**2
        rows[str(y)] = dict(inner=actual, pure_open_core=core, restored_x0_x1=restored,
                           linear_u=linear_u, linear_v=linear_v, pole=pole,
                           infinity_correction=1,
                           max_abs_u=max(abs(int(z)) for z in u),
                           max_abs_v=max(abs(int(z)) for z in v))
    exceptional = {}
    for y in (0, 1):
        value = sum(int(s[x, y]*u_raw[x, y]*v_raw[x, y]) for x in range(p))
        assert abs(value) <= p**3
        exceptional[str(y)] = value
    contraction = sum(row['inner'] for row in rows.values())+sum(exceptional.values())
    ca4 = np.linalg.matrix_power(cs['A'], 4)
    assert trace(ca4) % (p-1) == 0
    t4 = trace(ca4)//(p-1)
    p0 = eye.copy(); p0[0, 0] = 0
    nzero = trace(ca4@p0@s@p0@s)
    assert nzero == (p-1)*(p*t4-2)
    value = direct_trace('ABABCC', cs)
    assert value == contraction-p*t4+2 and value*value <= 56**2*p**7
    orbit_values = {w: direct_trace(w, cs) for w in new_words}
    # Each S3 orbit element uses at most one inversion exchange after exact reflections.
    reference_words = dihedral('ABABCC') | dihedral('ABABCC'.translate(str.maketrans('AB', 'BA')))
    reference_values = {direct_trace(w, cs) for w in reference_words}
    assert reference_values == {value}
    assert all(abs(v-value) <= 4*p**3 for v in orbit_values.values())
    assert all(v*v <= 58**2*p**7 for v in orbit_values.values())
    graphs = {}; literals = {}
    if p == 5:
        edges = necklace_edges('ABABCC')
        graphs = dict(full=graph_sum(p, chi, range(8), edges),
                      distinct_anchors=graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 1}),
                      coincident_anchors=graph_sum(p, chi, range(8), edges, fixed={6: 0, 7: 0}),
                      adjacent_C_pair=graph_sum(p, chi, range(8), edges, fixed={4: 1, 5: 0}))
        assert graphs['distinct_anchors'] == value
        assert graphs['coincident_anchors'] == nzero
        assert graphs['adjacent_C_pair'] == contraction
        assert graphs['full'] == p*(p-1)*value+p*nzero == p*(p-1)*contraction
        for w in ('ABABCC', 'AABCBC'):
            literals[w] = literal(p, chi, d, w)
            assert literals[w] == orbit_values[w]
    return dict(p=p, cross_ratio_entries=cross_ratio_checks, affine_Legendre_entries=p*(p-2),
                actual_zero_stalk_checks=input_stalk_checks,
                compact_trace_vector_checks=2*(p-2), corrected_pointwise_traces=2*p*(p-2),
                generic_inner_expansions_and_bounds=len(rows), open_core_bounds=len(rows),
                exceptional_y=exceptional, generic_rows=rows,
                t4=t4, coincident_anchor_term=nzero, crossing_contraction=contraction,
                N_ABABCC=value, orbit_values=orbit_values,
                literal_graphs=graphs, literal_necklaces=literals)


def extension_case(n):
    """Construct F_(5^n), verifying that every nonzero element is a unit."""
    p = 5; q = p**n
    modulus = {1: [0], 2: [2, 0], 3: [1, 1, 0], 4: [2, 0, 0, 0]}[n]
    digits = [[(x//p**i) % p for i in range(n)] for x in range(q)]

    def multiply(x, y):
        a = [0]*(2*n-1)
        for i in range(n):
            for j in range(n):
                a[i+j] += digits[x][i]*digits[y][j]
        for i in range(2*n-2, n-1, -1):
            for j in range(n):
                a[i-n+j] -= a[i]*modulus[j]
        return sum((a[i] % p)*p**i for i in range(n))

    def power(x, k):
        out = 1
        while k:
            if k & 1:
                out = multiply(out, x)
            x = multiply(x, x); k //= 2
        return out

    def subtract(x, y):
        return sum(((digits[x][i]-digits[y][i]) % p)*p**i for i in range(n))

    chi = [0]
    for x in range(1, q):
        assert power(x, q-1) == 1
        euler = power(x, (q-1)//2)
        assert euler in (1, p-1)
        chi.append(1 if euler == 1 else -1)
    chi = np.array(chi, dtype=np.int64)
    s = np.array([[chi[subtract(x, y)] for y in range(q)] for x in range(q)], dtype=np.int64)
    b = np.array([chi[subtract(x, 1)] for x in range(q)], dtype=np.int64)
    traces = {}
    for y in (2, 3, 4):
        cy = np.array([chi[subtract(x, y)] for x in range(q)], dtype=np.int64)
        qy = s@(chi*b*cy)
        kby = s@(b*cy)
        raw_u = s@(chi*qy); raw_v = s@(chi*b*kby)
        for x in (2, 3, 4):
            if x != y:
                traces[f'{x},{y}'] = dict(u=int(raw_u[x]-1), v=int(raw_v[x]-1))
    return dict(p=p, degree=n, q=q, modulus_low_coefficients=modulus,
                verified_nonzero_units=q-1, generic_points=traces)


def frobenius_polynomials(extension_cases):
    output = {}
    p = 5
    for point in extension_cases[0]['generic_points']:
        output[point] = {}
        for name in ('u', 'v'):
            powers = [0]+[case['generic_points'][point][name] for case in extension_cases]
            elementary = [1]
            for k in range(1, 5):
                numerator = sum((-1)**(i-1)*elementary[k-i]*powers[i] for i in range(1, k+1))
                assert numerator % k == 0
                elementary.append(numerator//k)
            # Coefficients of X^4+c1 X^3+c2 X^2+c3 X+c4.
            coefficients = [(-1)**i*v for i, v in enumerate(elementary)]
            _, a, b, c, d = coefficients
            assert d in (p**4, -p**4)
            if d == p**4:
                assert c == p*p*a
                discriminant = a*a-4*(b-2*p*p)
                assert discriminant >= 0 and 4*p >= abs(a)
                assert discriminant <= (4*p-abs(a))**2
                certificate = 'Both real roots of Z²+aZ+b−2p² lie in [−2p,2p], with Z=X+p²/X.'
            else:
                assert b == 0 and c == -p*p*a and abs(a) <= 2*p
                certificate = '(X²−p²)(X²+aX+p²), with |a|≤2p.'
            output[point][name] = dict(extension_traces=powers[1:], coefficients=coefficients,
                                      exact_modulus_certificate=certificate, root_modulus=p)
    return output


def main():
    coverage = catalog()
    cases = []
    for p in (5, 13, 17, 29, 41):
        cases.append(field_case(p, coverage['new_words']))
        print(json.dumps({'p': p, 'status': 'exact crossing, boundary and orbit checks passed'}), flush=True)
    extension_cases = []
    for n in range(1, 5):
        extension_cases.append(extension_case(n))
        print(json.dumps({'q': 5**n, 'status': 'original extension-field trace sums complete'}), flush=True)
    polynomials = frobenius_polynomials(extension_cases)
    files = ['research/parallel7-necklace-2026-09-04.md',
             'experiments/parallel7_necklace_2026_09_04.py',
             'research/parallel6-necklace-2026-09-04.md',
             'research/parallel2-necklace-2026-09-04.md',
             'research/parallel-necklace-2026-09-04.md',
             'experiments/parallel6_necklace_2026_09_04.py',
             'experiments/parallel_necklace_2026_09_04.py',
             'experiments/planar_necklace_reductions.py',
             'sources/katz-finite-field-mellin.pdf']
    source_archive = ROOT/'sources/katz-rigid-local-systems.pdf'
    if source_archive.exists():
        assert sha256(source_archive.read_bytes()).hexdigest() == 'ca6eb8d5d9e21076dda4e154e83dfa2821f586d6ccbe8c8f1b371b4352f753d1'
        files.append(str(source_archive.relative_to(ROOT)))
    result = dict(status='All exact checks passed; the note proves the 36-word crossing class, while 42 balanced words remain open.',
                  arithmetic='Exact integers throughout. Extension sums use int64 below their checked maximum possible magnitude 625².',
                  scope='Fixed length-six crossing class; no full signed aggregate or Paley conclusion.',
                  primary_source=dict(url='https://web.math.princeton.edu/~nmk/wholebookRLScorr.pdf',
                      bytes=1115958, sha256='ca6eb8d5d9e21076dda4e154e83dfa2821f586d6ccbe8c8f1b371b4352f753d1',
                      sections=['2.9.4', '2.9.7', '3.3.3', '3.3.6', '3.3.7', '5.5.5.10']),
                  finite_word_coverage=coverage, cases=cases, extension_cases=extension_cases,
                  frobenius_polynomials=polynomials,
                  input_sha256={f: sha256((ROOT/f).read_bytes()).hexdigest() for f in files})
    output = ROOT/'results/parallel7_necklace_2026_09_04.json'
    output.write_text(json.dumps(result, indent=2)+'\n')
    print(json.dumps({'written': str(output)}), flush=True)


if __name__ == '__main__':
    main()
