#!/usr/bin/env python3
"""Read back track coverage and reconstructed node exports; reject tampering.

This verifier shares interpolation and field helpers with the generator. It does
not use ZDD operations. Independent census equality is a separate test lane.
"""
from copy import deepcopy
from hashlib import sha256
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
from symbolic_batches import PrimeField, Field, run_census, verify_batches, node_set

HERE = Path(__file__).resolve().parent


def verify_export(F, record):
    n, k, s = record['n'], record['k'], record['s']
    assert record['field'] == {'p': F.p, 'e': F.e, 'q': F.q}
    assert (record['cover_certificate']['n'], record['cover_certificate']['k'],
            record['cover_certificate']['s']) == (n, k, s)
    assert len(record['domain']) == len(record['u0']) == len(record['u1']) == n
    assert all(isinstance(x, int) and 0 <= x < F.q
               for key in ['domain', 'u0', 'u1'] for x in record[key])
    assert len(node_set(record)) == len(record['nodes'])
    coverage = verify_batches(F, record)
    nodes, finite, whole = {}, [], []
    for batch in record['batch_certificates']:
        common, buckets = [], {}
        A, B = batch['intercept'], batch['slope']
        for i in range(n):
            intercept = F.sub(record['u0'][i], A[i])
            slope = F.sub(record['u1'][i], B[i])
            if slope:
                scalar = F.mul(F.sub(0, intercept), F.inv(slope))
                buckets.setdefault(scalar, []).append(i)
            elif not intercept:
                common.append(i)
        assert common == batch['common_coordinates']
        if len(common) >= s:
            assert batch['all_field'] is True
            whole.append(batch)
            scalars = range(F.q) if record['node_ledger_materialized'] else []
        else:
            scalars = [z for z in buckets if len(common) + len(buckets[z]) >= s]
            if scalars:
                assert batch['all_field'] is False
                assert batch['qualifying_scalar_buckets'] == {str(z): buckets[z] for z in sorted(scalars)}
                finite.append(batch)
        for scalar in scalars:
            cw = tuple(F.add(a, F.mul(scalar, b)) for a, b in zip(A, B))
            support = tuple(sorted(common + buckets.get(scalar, [])))
            nodes[(scalar, cw)] = (scalar, cw, support)
    assert set(nodes.values()) == node_set(record)
    assert finite == record['finite_track_certificates']
    assert whole == record['all_field_track_certificates']
    assert bool(whole) == record['whole_field_correlated']
    scalar_set = sorted({z for z, _ in nodes})
    assert record['finite_bad_scalar_count'] == (F.q if whole else len(scalar_set))
    expected = (list(range(F.q)) if record['node_ledger_materialized'] else None) if whole else scalar_set
    assert record['bad_scalars'] == expected
    if not whole:
        assert record['node_ledger_materialized']
    return {**coverage, 'node_exports_checked': len(nodes), 'whole_field_tracks_checked': len(whole)}


def main():
    results = json.loads((HERE / 'results.json').read_text())
    records = [x['result'] for x in results['cost_experiments'] + results['large_structured']]
    F = Field(3, [1, 0, 1]); dom = list(range(9))
    records += [run_census(F, dom, 2, 5, [F.add(1, x) for x in dom],
                          [F.mul(2, x) for x in dom], materialize_limit=limit)
                for limit in [0, 1000]]
    checks = []
    for r in records:
        field = PrimeField(r['field']['p']) if r['field']['e'] == 1 else F
        checks.append(verify_export(field, r))
    original = records[1]
    def mutation(name, function):
        record = deepcopy(original)
        function(record)
        try:
            verify_export(PrimeField(65537), record)
        except AssertionError:
            return name
        raise AssertionError('accepted tampering: ' + name)
    rejected = [
        mutation('basis_count', lambda r: r['batch_certificates'][0].__setitem__('base_multiplicity', 1)),
        mutation('missing_track', lambda r: r['batch_certificates'].pop()),
        mutation('duplicate_track', lambda r: r['batch_certificates'].append(deepcopy(r['batch_certificates'][0]))),
        mutation('false_common_coordinate', lambda r: r['batch_certificates'][0]['common_coordinates'].append(47)),
        mutation('missing_node', lambda r: r['nodes'].clear()),
        mutation('mismatched_cover_problem', lambda r: r.__setitem__('s', r['s'] - 1)),
        mutation('noncanonical_field_coordinate', lambda r: r['domain'].__setitem__(0, r['domain'][0] + 65537)),
        mutation('duplicate_node', lambda r: r['nodes'].append(deepcopy(r['nodes'][0]))),
    ]
    output = {'status': 'passed', 'exports_checked': checks, 'tampering_rejected': rejected,
              'scope': 'shared-helper readback, no ZDD replay; not an independently implemented end-to-end decoder',
              'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ['symbolic_batches.py', 'verify_exports.py', 'results.json']}}
    (HERE / 'verification_results.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
