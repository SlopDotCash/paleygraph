#!/usr/bin/env python3
"""Verify metric transformations and replay queries independently of producer."""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import sympy as sp
from conditioned_review import certify, replay

HERE = Path(__file__).resolve().parent


def check_transform(c):
    change = c['metric_change']; assert change['rule'] == 'LLL_on_integer_adjugate_inverse_images'
    S = sp.Matrix(change['unimodular_transform']); old = sp.Matrix(change['old_basis'])
    assert S.shape == old.shape == (len(c['erased']), len(c['erased']))
    assert all(v.is_Integer for v in S) and abs(S.det()) == 1
    assert S*old == sp.Matrix(c['basis'])


def main():
    path = HERE/'metric_summary.json'; summary = json.loads(path.read_text())
    inputs = {path.name: sha256(path.read_bytes()).hexdigest()}; cases = []; sample = None
    for row in summary['cases']:
        path = HERE/row['artifact']; data = json.loads(path.read_text()); inputs[path.name] = sha256(path.read_bytes()).hexdigest()
        c = data['certificate']; check_transform(c); prepared = certify(c); checks = 0
        source = HERE.parent/'round21'/('heldout_'+data['name'].removeprefix('metric_')+'.json')
        old = json.loads(source.read_text()); inputs['../round21/'+source.name] = sha256(source.read_bytes()).hexdigest()
        assert c['metric_change']['old_basis'] == old['certificate']['basis']
        assert len(data['records']) == len(old['records'])
        for r, original in zip(data['records'], old['records']):
            assert r['known'] == original['known'] and r['kind'] == original['kind'] and r['anchor_scalar'] == original['anchor_scalar']
            out = replay(c, dict(r['known']), data['candidate_budget'], prepared)
            assert out == r['output']; checks += out.get('scalar_candidates_checked', 0)
            if out['status'] == 'complete' and r['kind'] == 'actual': assert r['anchor_scalar'] in [v['scalar'] for v in out['completions']]
        result = {'name': data['name'], 'unimodular_basis_change_verified': True, 'cuts_verified': len(c['erased']),
                  'optimality_certificates': sum(x['optimality_certified'] for x in c['conditioning']['cuts']),
                  'queries_replayed': len(data['records']), 'scalar_candidates_checked': checks}
        cases.append(result); print(json.dumps(result), flush=True); sample = c
    rejected = []
    for label in ['non_unimodular', 'wrong_original_basis']:
        c = deepcopy(sample)
        if label == 'non_unimodular': c['metric_change']['unimodular_transform'][0] = [0]*len(c['erased'])
        else: c['metric_change']['old_basis'][0][0] += 1
        try: check_transform(c)
        except AssertionError: rejected.append(label)
        else: raise AssertionError('bad basis transformation accepted')
    (HERE/'metric_review.json').write_text(json.dumps({'status': 'passed', 'cases': cases,
        'corrupt_transformations_rejected': rejected, 'input_sha256': inputs,
        'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['metric_review.py', 'conditioned_review.py']}}, indent=2)+'\n')


if __name__ == '__main__': main()
