#!/usr/bin/env python3
"""Full constant-code/scalar oracle in GF(3^6), a tiny degree-six field check."""
from hashlib import sha256
from pathlib import Path
import json
import sys

sys.dont_write_bytecode = True
from scalar_fibers import Field, certify_pencil, materialize_scalar, verify_certificate

HERE = Path(__file__).resolve().parent


def main():
    F = Field(3, [1, 0, 0, 0, 1, 1, 1])
    dom = [0, 1, 3, 9]
    direction = [1, 100, 1, 100]
    zstar, anchor_value = 642, 517
    u0 = [F.sub(anchor_value, F.mul(zstar, d)) for d in direction]
    pieces = [{'coefficients': [1], 'coordinates': [0, 2]},
              {'coefficients': [100], 'coordinates': [1, 3]}]
    certificate = certify_pencil(F, dom, 1, 2, u0, direction, pieces)
    verify_certificate(F, certificate)
    nodes_checked = 0
    for z in range(F.q):
        word = [F.add(a, F.mul(z, b)) for a, b in zip(u0, direction)]
        exact = {(z, (c, c, c, c), tuple(i for i, value in enumerate(word) if value == c))
                 for c in range(F.q) if sum(value == c for value in word) >= 2}
        transported = {(x['scalar'], tuple(x['codeword']), tuple(x['agreement_support']))
                       for x in materialize_scalar(F, certificate, z)}
        assert exact == transported
        nodes_checked += len(exact)
    assert certificate['anchor']['scalar'] == zstar
    assert certificate['bad_scalar_count'] == 729
    assert nodes_checked == certificate['complete_symbolic_node_count'] == 1457
    result = {'status': 'passed', 'field': 'GF(3^6)', 'full_scalar_constant_codeword_pairs': 729**2,
              'exact_nodes_compared': nodes_checked, 'certificate': certificate,
              'scope': 'n4,k1,s2 exact toy; field degree matches the profile, field size does not',
              'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest()
                                for name in ['scalar_fibers.py', 'degree_six_control.py']}}
    (HERE / 'degree_six_results.json').write_text(json.dumps(result, indent=2) + '\n')
    print('Passed GF729: 531441 scalar/constant pairs,1457 exact nodes,anchor642.')


if __name__ == '__main__':
    main()
