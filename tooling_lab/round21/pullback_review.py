#!/usr/bin/env python3
"""Independent rational matrix identities for the sharper erasure radii."""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from math import prod
import json
from pathlib import Path
import sympy as sp
from erasure_review import certify, replay, rational

HERE = Path(__file__).resolve().parent


def certify_pullback(c):
    assert c['radius_rule'] == 'minimum_of_digit_cube_and_centered_orbit_pullback'
    base = deepcopy(c); base['coordinate_radii'] = c['digit_coordinate_radii']
    digit_radii = [rational(v) for v in base['coordinate_radii']]; visible = c['visible_directions']
    base['universal_candidate_box_cap'] = prod(int(2*digit_radii[j])+1 for j in visible)
    base['universal_unique_completion'] = all(2*digit_radii[j] < 1 for j in visible)
    prepared = certify(base); inv, _, encoder, limits = prepared
    p, N, f, A = c['p'], c['N'], c['relation'], c['erased']; d = len(A)
    M = sp.Matrix([[f[i-j] if i >= j else -f[N+i-j] for j in range(N)] for i in A])
    V = sp.Matrix([[sp.Rational(x.numerator, x.denominator) for x in row] for row in inv])
    expected = M.T*V
    assert len(c['centered_coordinate_forms']) == d
    radii = []
    for j in range(d):
        given = [rational(v) for v in c['centered_coordinate_forms'][j]]
        assert given == [Fraction(expected[i,j]) for i in range(N)]
        radius = Fraction(p//2, p)*sum(map(abs, given))
        assert radius == rational(c['centered_coordinate_radii'][j])
        radii.append(min(digit_radii[j], radius))
    assert radii == [rational(v) for v in c['coordinate_radii']]
    assert c['universal_candidate_box_cap'] == prod(int(2*radii[j])+1 for j in visible)
    assert c['universal_unique_completion'] == all(2*radii[j] < 1 for j in visible)
    return inv, radii, encoder, limits


def main():
    path = HERE/'pullback_summary.json'; summary = json.loads(path.read_text()); inputs = {path.name: sha256(path.read_bytes()).hexdigest()}
    heldout = HERE/'heldout_summary.json'; extra = json.loads(heldout.read_text()); inputs[heldout.name] = sha256(heldout.read_bytes()).hexdigest()
    results = []; sample = None
    for row in summary['cases']+extra['cases']:
        path = HERE/row['artifact']; r = json.loads(path.read_text()); inputs[path.name] = sha256(path.read_bytes()).hexdigest()
        c = r['certificate']; prepared = certify_pullback(c); checks = 0
        for q in r['records']:
            out = replay(c, dict(q['known']), r['candidate_budget'], prepared); assert out == q['output']
            checks += out.get('scalar_candidates_checked', 0)
        for q in r.get('cyclic_controls', []):
            N, s, d = c['N'], q['start'], len(c['erased']); known = dict(q['known'])
            rotated = {j: known[(j+s) % N]*(1 if j+s < N else -1) for j in range(d, N)}
            out = replay(c, rotated, r['candidate_budget'], prepared)
            scalars = sorted(x['scalar']*pow(c['g'], -s, c['p']) % c['p'] for x in out['completions'])
            words = prepared[2](scalars); out['completions'] = [{'scalar': a, 'digits': w} for a, w in zip(scalars, words)]
            assert out == q['output'] and q['anchor_scalar'] in scalars
            if c['universal_unique_completion']: assert scalars == [q['anchor_scalar']]
            assert all(all(w[j] == v for j, v in known.items()) for w in words)
        result = {'name': r['name'], 'dual_functional_identities_checked': len(c['erased']),
                  'queries_replayed': len(r['records']), 'scalar_candidates_checked': checks,
                  'cyclic_controls_checked': len(r.get('cyclic_controls', [])),
                  'universal_candidate_box_cap': c['universal_candidate_box_cap'],
                  'universal_unique_completion': c['universal_unique_completion']}
        results.append(result); print(json.dumps(result), flush=True)
        sample = c
    rejected = []
    for label in ['false_form', 'understated_physical_radius', 'wrong_minimum_radius', 'false_uniqueness']:
        changed = deepcopy(sample); j = changed['visible_directions'][0]
        if label == 'false_form': changed['centered_coordinate_forms'][j][0][0] += 1
        if label == 'understated_physical_radius': changed['centered_coordinate_radii'][j] = [0, 1]
        if label == 'wrong_minimum_radius': changed['coordinate_radii'][j] = [0, 1]
        if label == 'false_uniqueness': changed['universal_unique_completion'] = not changed['universal_unique_completion']
        try: certify_pullback(changed)
        except AssertionError: rejected.append(label)
        else: raise AssertionError('corrupt pullback accepted')
    out = {'status': 'passed', 'scope': 'Every centered-orbit dual identity and radius verified with separate rational matrix multiplication; all scalar candidate sets and cyclic transports independently replayed.',
           'cases': results, 'corrupt_pullback_certificates_rejected': rejected,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['pullback_review.py', 'erasure_review.py']}, 'input_sha256': inputs}
    (HERE/'pullback_review.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
