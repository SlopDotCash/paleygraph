#!/usr/bin/env python3
"""Discover, exercise and record bounded erasure certificates and failures."""
from collections import Counter
from hashlib import sha256
from itertools import product
import json
from pathlib import Path
import random
import fpylll
from lattice_erasure import prepare, decode, decode_cyclic, encode

HERE = Path(__file__).resolve().parent


def save(name, c, records, scope, extra=None):
    target = HERE/(name+'.json')
    out = {'name': name, 'scope': scope, 'certificate': c, 'records': records, **(extra or {})}
    target.write_text(json.dumps(out, separators=(',', ':'))+'\n')
    stats = {'name': name, 'artifact': target.name, 'p': c['p'], 'N': c['N'], 'erasures': len(c['erased']),
             'visible_directions': len(c['visible_directions']), 'invisible_directions': len(c['invisible_directions']),
             'universal_unique_completion': c['universal_unique_completion'],
             'universal_candidate_box_cap': c['universal_candidate_box_cap'], 'queries': len(records),
             'statuses': dict(Counter(r['output']['status'] for r in records)),
             'complete_counts': dict(Counter(str(r['output'].get('count')) for r in records)),
             'maximum_observed_box': max((r['output']['candidate_box_size'] for r in records), default=0)}
    print(json.dumps(stats), flush=True); return stats


def main():
    inputs = {}; specs = [(17, 2, [2, 0, 0, 1], 'p17_basic')]
    for p in [41, 97]:
        path = HERE.parent/f'round20/p{p}_general.json'; r = json.loads(path.read_text())['config']
        inputs['../round20/'+path.name] = sha256(path.read_bytes()).hexdigest()
        specs.append((r['p'], r['g'], r['relation'], f'p{p}_general'))
    summaries = []
    for p, g, f, label in specs:
        words = [encode(p, g, f, a)[0] for a in range(p)]
        for d in range(5):
            c = prepare(p, g, f, range(d)); records = []; B = c['digit_bound']
            for values in product(range(-B, B+1), repeat=4-d):
                known = dict(zip(range(d, 4), values)); out = decode(c, known, 100000)
                expected = [a for a, w in enumerate(words) if all(w[j] == v for j, v in known.items())]
                assert out['status'] == 'complete' and [r['scalar'] for r in out['completions']] == expected
                records.append({'known': list(map(list, known.items())), 'output': out})
            summaries.append(save(f'{label}_erase{d}', c, records, 'Every bounded assignment on every unerased coordinate; expected completions from the full scalar codebook.'))
    # Degenerate calibration and non-unit syndrome controls do not assume beta.
    p, g, N = 17, 2, 4
    for label, f in [('squared_relation', [4, 0, -1, 4]), ('scalar_p', [17, 0, 0, 0]), ('scalar_p_squared', [289, 0, 0, 0])]:
        for A in [[0], [1, 2, 3]]:
            c = prepare(p, g, f, A); words = [encode(p, g, f, a)[0] for a in range(p)]
            records = []; seen = set()
            for w in words:
                known = {j: v for j, v in enumerate(w) if j not in A}
                candidates = [known, {**known, min(known): 0 if known[min(known)] else 1}]
                for known in candidates:
                    key = tuple(sorted(known.items()))
                    if key in seen: continue
                    seen.add(key); out = decode(c, known, 100000)
                    expected = [a for a, word in enumerate(words) if all(word[j] == v for j, v in known.items())]
                    assert out['status'] == 'complete' and [r['scalar'] for r in out['completions']] == expected
                    records.append({'known': list(map(list, known.items())), 'output': out})
            summaries.append(save(label+'_erase'+str(len(A)), c, records, 'All scalar-origin assignments and one bounded edit per distinct assignment; checked against the full17-word codebook.'))
    source = HERE.parent/'round17/norm_compression.json'; raw = json.loads(source.read_text())
    inputs['../round17/norm_compression.json'] = sha256(source.read_bytes()).hexdigest(); old = raw['cases'][-1]
    p, g, f = old['p'], old['g'], old['relation']; N = len(f); rng = random.Random(20260906)
    scalars = list(dict.fromkeys([0, 1, p//2, (p+1)//2, 1234567]+[rng.randrange(p) for _ in range(59)]))
    for d in [4, 8, 13, 14, 15, 16, 20, 24, 32]:
        c = prepare(p, g, f, range(d)); records = []; first_budget = []
        budget = 100000 if d == 24 else 20000
        count = 64 if d <= 16 else (8 if d == 20 else 3)
        for a in scalars[:count]:
            w = encode(p, g, f, a)[0]; known = {j: w[j] for j in range(d, N)}
            for kind, data in [('actual', known), ('edited', {**known, d: 0 if known[d] else 1})]:
                if d == 24: first_budget.append({'known': list(map(list, data.items())), 'output': decode(c, data, 20000)})
                out = decode(c, data, budget)
                if kind == 'actual' and out['status'] == 'complete': assert a in [r['scalar'] for r in out['completions']]
                records.append({'kind': kind, 'anchor_scalar': a, 'known': list(map(list, data.items())), 'output': out})
        extras = {'candidate_budget': budget}
        if first_budget: extras['initial_budget_controls'] = first_budget
        if d == 15:
            rotated = []
            for start in range(N):
                for a in scalars[:8]:
                    w = encode(p, g, f, a)[0]; erased = {(start+j) % N for j in range(d)}
                    known = {j: v for j, v in enumerate(w) if j not in erased}
                    out = decode_cyclic(c, start, known)
                    assert out['status'] == 'complete' and [v['scalar'] for v in out['completions']] == [a]
                    rotated.append({'start': start, 'anchor_scalar': a, 'known': list(map(list, known.items())), 'output': out})
            extras['cyclic_controls'] = rotated
            # A bounded lattice alias shares all known digits with zero. The
            # scalar quotient identifies it with zero; re-encoding removes it.
            extras['bounded_syndrome_alias'] = {'digits': f, 'known': [[j, 0] for j in range(d, N)],
                                                'lifted_first_coordinate': p,
                                                'decoded': decode(c, {j: 0 for j in range(d, N)})}
        summaries.append(save(f'p{p}_erase{d}', c, records, 'Complete candidate-box decoding or explicitly incomplete budget status; actual/edited bounded queries. Universal uniqueness/candidate caps come from the exact certificate, not these samples.', extras))
    out = {'status': 'produced', 'scope': 'Exact lattice-coordinate erasure decoding with scalar-invisible directions removed. No k-state table, full p-codebook at large p, asymptotic runtime claim or prize result.',
           'cases': summaries, 'fpylll_version': fpylll.__version__,
           'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['erasure_experiments.py', 'lattice_erasure.py', '../round17/realizability.py', '../round19/projection_codec.py']},
           'input_sha256': inputs}
    (HERE/'erasure_summary.json').write_text(json.dumps(out, indent=2)+'\n')


if __name__ == '__main__': main()
