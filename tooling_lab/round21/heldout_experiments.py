#!/usr/bin/env python3
"""Held-out scalar and erasure-geometry controls for the improved decoder."""
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
import random
from lattice_erasure import prepare, decode, decode_cyclic, encode
from pullback_erasure import improve

HERE = Path(__file__).resolve().parent


def main():
    source = HERE/'pullback_erase32.json'; base = json.loads(source.read_text())['certificate']
    p, N, g, f = base['p'], base['N'], base['g'], base['relation']; rng = random.Random(20260907)
    scalars = [rng.randrange(p) for _ in range(16)]; assert len(set(scalars)) == 16
    random_set = sorted(rng.sample(range(N), 32)); rows = []
    for label, A in [('consecutive', list(range(32))), ('alternating', list(range(0, 64, 2))), ('random', random_set)]:
        c = base if label == 'consecutive' else improve(prepare(p, g, f, A)); records = []; cyclic = []
        for a in scalars:
            word = encode(p, g, f, a)[0]; known = {j: v for j, v in enumerate(word) if j not in A}; j = min(known)
            for kind, data in [('actual', known), ('edited', {**known, j: 0 if known[j] else 1})]:
                out = decode(c, data, 50000)
                if kind == 'actual' and out['status'] == 'complete': assert a in [r['scalar'] for r in out['completions']]
                records.append({'kind': kind, 'anchor_scalar': a, 'known': list(map(list, data.items())), 'output': out})
        if label == 'consecutive':
            a = scalars[0]; word = encode(p, g, f, a)[0]
            for start in range(N):
                erased = {(start+j) % N for j in range(32)}; known = {j: v for j, v in enumerate(word) if j not in erased}
                out = decode_cyclic(c, start, known, 50000)
                assert out['status'] == 'complete' and a in [x['scalar'] for x in out['completions']]
                cyclic.append({'start': start, 'anchor_scalar': a, 'known': list(map(list, known.items())), 'output': out})
        name = 'heldout_'+label+'32'; path = HERE/(name+'.json')
        path.write_text(json.dumps({'name': name, 'seed': 20260907, 'certificate': c, 'candidate_budget': 50000,
                                    'records': records, 'cyclic_controls': cyclic}, separators=(',', ':'))+'\n')
        row = {'name': name, 'artifact': path.name, 'erased': A, 'universal_candidate_box_cap': c['universal_candidate_box_cap'],
               'universal_unique_completion': c['universal_unique_completion'], 'queries': len(records),
               'statuses': dict(Counter(r['output']['status'] for r in records)),
               'complete_counts': dict(Counter(str(r['output'].get('count')) for r in records)),
               'maximum_observed_box': max(r['output']['candidate_box_size'] for r in records),
               'cyclic_controls': len(cyclic)}
        rows.append(row); print(json.dumps(row), flush=True)
    out = {'status': 'produced', 'scope': 'Same16 held-out scalar anchors under consecutive, alternating and seeded random32-coordinate erasures, each with a bounded edit; incomplete budgets explicitly retained. Consecutive cyclic transport checked at all64 starts.',
           'cases': rows, 'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['heldout_experiments.py', 'lattice_erasure.py', 'pullback_erasure.py']},
           'input_sha256': {'pullback_erase32.json': sha256(source.read_bytes()).hexdigest()}}
    (HERE/'heldout_summary.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
