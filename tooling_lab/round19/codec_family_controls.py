#!/usr/bin/env python3
"""Exhaust bounded small relations and digit cubes against scalar codebooks."""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import sys
import sympy as sp
from projection_codec import prepare, decode

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent/'round17'))
from carry_review import product as polynomial_product


def main():
    cases = []
    for p, radius in ((17, 2), (41, 2), (97, 3)):
        N = 4; primitive = int(sp.primitive_root(p)); g = pow(primitive, (p-1)//8, p)
        assert pow(g, N, p) == p-1; u = pow(g, -1, p)
        relations = []; modes = Counter(); total = 0
        smaller_count = sum(any(f) and sum(v*pow(u, j, p) for j, v in enumerate(f)) % p == 0
                            for f in product(range(-2, 3), repeat=N))
        if p == 97: assert smaller_count == 0
        for f in product(range(-radius, radius+1), repeat=N):
            if not any(f) or sum(v*pow(u, j, p) for j, v in enumerate(f)) % p: continue
            c = prepare(p, g, list(f)); modes[c['mode']] += 1
            # Independent exact codebook: rational field inverses are not used.
            oracle = {}
            for a in range(p):
                F = []
                for j in range(N):
                    v = a*pow(g, j, p) % p; F.append(v if v <= p//2 else v-p)
                raw = polynomial_product(list(f), F, N); assert all(v % p == 0 for v in raw)
                digits = tuple(v//p for v in raw); assert digits not in oracle
                oracle[digits] = a
            bound = (p-1)*sum(map(abs, f))//(2*p)
            assert c['digit_bound'] == bound
            digest = sha256(); count = 0; accepted = 0
            for digits in product(range(-bound, bound+1), repeat=N):
                result = decode(c, digits); yes = result['status'].startswith('accepted_')
                assert yes == (digits in oracle)
                if yes: assert result['extracted_scalar'] == oracle[digits]
                accepted += yes; count += 1; digest.update(bytes([yes]))
            assert accepted == p
            relations.append({'relation': list(f), 'mode': c['mode'], 'digit_bound': bound,
                              'cube_words_checked': count, 'actual_scalars': accepted,
                              'verdict_stream_digest': digest.hexdigest()})
            total += count
        assert relations
        cases.append({'p': p, 'N': N, 'g': g, 'relation_coefficient_interval': [-radius, radius],
                      'radius_two_nonzero_relations': smaller_count,
                      'relation_vectors_examined': (2*radius+1)**N,
                      'relations_checked': len(relations), 'modes': dict(modes),
                      'digit_words_checked': total, 'records': relations})
        print(json.dumps({key: value for key, value in cases[-1].items() if key != 'records'}), flush=True)
    out = {'status': 'passed', 'scope': 'Every nonzero relation in the stated four-dimensional coefficient box vanishing at the selected inverse root; every word in that relation-specific digit cube checked against a complete direct scalar codebook. p97 has no relation in[-2,2]^4, so its nonempty family uses[-3,3]^4.',
           'cases': cases,
           'source_sha256': {name: sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ['codec_family_controls.py', 'projection_codec.py', '../round17/realizability.py',
                                          '../round17/norm_carry.py', '../round17/carry_review.py']}}
    (HERE/'codec_family_controls.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
