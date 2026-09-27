#!/usr/bin/env python3
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from conditioned_erasure import improve, decode

HERE = Path(__file__).resolve().parent


def main():
    cases = []; inputs = {}
    for name in ['consecutive32', 'alternating32', 'random32']:
        path = HERE.parent/'round21'/('heldout_'+name+'.json')
        inputs[str(path.relative_to(HERE.parent)).replace('round21/', '../round21/')] = sha256(path.read_bytes()).hexdigest()
        old = json.loads(path.read_text()); c = improve(old['certificate'])
        records = []
        for r in old['records']:
            out = decode(c, dict(r['known']), old['candidate_budget'])
            if out['status'] == 'complete' and r['kind'] == 'actual':
                assert r['anchor_scalar'] in [v['scalar'] for v in out['completions']]
            if r['output']['status'] == 'complete': assert out['completions'] == r['output']['completions']
            records.append({**r, 'baseline_output': r['output'], 'output': out})
        data = {'name': name, 'certificate': c, 'candidate_budget': old['candidate_budget'], 'records': records}
        output = HERE/('conditioned_'+name+'.json'); output.write_text(json.dumps(data, separators=(',', ':'))+'\n')
        row = {'name': name, 'artifact': output.name, 'baseline_cap': c['unconditioned_universal_cap'],
               'universal_cap': c['universal_candidate_box_cap'], 'unique_completion': c['universal_unique_completion'],
               'optimality_certificates': sum(cut['optimality_certified'] for cut in c['conditioning']['cuts']),
               'cuts': len(c['conditioning']['cuts']), 'statuses': dict(Counter(r['output']['status'] for r in records)),
               'complete_counts': dict(Counter(str(r['output']['count']) for r in records if r['output']['status'] == 'complete')),
               'maximum_observed_box': max(r['output']['candidate_box_size'] for r in records),
               'scalar_candidates_checked': sum(r['output'].get('scalar_candidates_checked', 0) for r in records),
               'radius_ratios': [float(F(*cut['radius'])/F(*v)) if F(*v) else None for cut, v in zip(c['conditioning']['cuts'], c['coordinate_radii'])]}
        cases.append(row); print(json.dumps(row), flush=True)
    report = {'status': 'produced', 'cases': cases, 'input_sha256': inputs,
              'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['conditioned_experiments.py', 'conditioned_erasure.py']}}
    (HERE/'conditioned_summary.json').write_text(json.dumps(report, indent=2)+'\n')


if __name__ == '__main__': main()
