#!/usr/bin/env python3
"""Replay the same large queries with exact cancellation retained."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from lattice_erasure import decode, decode_cyclic, encode
from pullback_erasure import improve

HERE = Path(__file__).resolve().parent


def main():
    rows = []; bindings = {}
    for d in [4, 8, 13, 14, 15, 16, 20, 24, 32]:
        source = HERE/f'p2013265921_erase{d}.json'; original = json.loads(source.read_text()); c = improve(original['certificate'])
        bindings[source.name] = sha256(source.read_bytes()).hexdigest(); records = []
        for r in original['records']:
            out = decode(c, dict(r['known']), 50000)
            assert out['status'] == 'complete'
            if r['output']['status'] == 'complete': assert out['completions'] == r['output']['completions']
            if r['kind'] == 'actual': assert r['anchor_scalar'] in [x['scalar'] for x in out['completions']]
            records.append({'kind': r['kind'], 'anchor_scalar': r['anchor_scalar'], 'known': r['known'], 'output': out})
        extras = {}
        if d == 20:
            rotated = []
            anchors = [r['anchor_scalar'] for r in original['records'] if r['kind'] == 'actual']
            for start in range(c['N']):
                for a in anchors:
                    word = encode(c['p'], c['g'], c['relation'], a)[0]; erased = {(start+j) % c['N'] for j in range(d)}
                    known = {j: v for j, v in enumerate(word) if j not in erased}
                    out = decode_cyclic(c, start, known, 50000)
                    assert out['status'] == 'complete' and [x['scalar'] for x in out['completions']] == [a]
                    rotated.append({'start': start, 'anchor_scalar': a, 'known': list(map(list, known.items())), 'output': out})
            extras['cyclic_controls'] = rotated
        name = f'pullback_erase{d}'; path = HERE/(name+'.json')
        path.write_text(json.dumps({'name': name, 'certificate': c, 'candidate_budget': 50000, 'records': records,
                                    'baseline_artifact': source.name, **extras}, separators=(',', ':'))+'\n')
        row = {'name': name, 'artifact': path.name, 'erasures': d,
               'baseline_universal_cap': original['certificate']['universal_candidate_box_cap'],
               'universal_candidate_box_cap': c['universal_candidate_box_cap'],
               'universal_unique_completion': c['universal_unique_completion'], 'queries': len(records),
               'complete_counts': dict(Counter(str(r['output']['count']) for r in records)),
               'maximum_observed_box': max(r['output']['candidate_box_size'] for r in records),
               'scalar_candidates_checked': sum(r['output']['scalar_candidates_checked'] for r in records),
               'cyclic_controls': len(extras.get('cyclic_controls', []))}
        rows.append(row); print(json.dumps(row), flush=True)
    out = {'status': 'produced', 'scope': 'Same large bounded-erasure queries with exact centered-orbit pullback radii. Uniform caps and uniqueness follow from rational dual identities, independent of the finite query sample.',
           'cases': rows,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['pullback_experiments.py', 'pullback_erasure.py', 'lattice_erasure.py']},
           'input_sha256': bindings}
    (HERE/'pullback_summary.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
