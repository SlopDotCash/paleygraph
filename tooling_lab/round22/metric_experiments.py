#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from copy import deepcopy
from metric_erasure import change_basis
from conditioned_erasure import improve, decode

HERE = Path(__file__).resolve().parent


def main():
    cases = []; inputs = {}
    for name in ['consecutive32', 'alternating32', 'random32']:
        path = HERE.parent/'round21'/('heldout_'+name+'.json')
        old = json.loads(path.read_text()); inputs['../round21/'+path.name] = sha256(path.read_bytes()).hexdigest()
        base = deepcopy(old['certificate']); base['coordinate_radii'] = base.pop('digit_coordinate_radii')
        # change_basis rewrites every basis-dependent field before pullback.
        c = improve(change_basis(base)); records = []
        for r in old['records']:
            out = decode(c, dict(r['known']), 50000)
            if out['status'] == 'complete' and r['kind'] == 'actual': assert r['anchor_scalar'] in [x['scalar'] for x in out['completions']]
            if r['output']['status'] == 'complete' and out['status'] == 'complete': assert out['completions'] == r['output']['completions']
            records.append({**r, 'output': out})
        path = HERE/('metric_'+name+'.json'); path.write_text(json.dumps({'name': 'metric_'+name, 'certificate': c, 'records': records, 'candidate_budget': 50000}, separators=(',', ':'))+'\n')
        row = {'name': 'metric_'+name, 'artifact': path.name, 'universal_cap': c['universal_candidate_box_cap'],
               'unconditioned_cap': c['unconditioned_universal_cap'], 'visible_directions': len(c['visible_directions']),
               'unique_completion': c['universal_unique_completion'], 'optimality_certificates': sum(cut['optimality_certified'] for cut in c['conditioning']['cuts']),
               'statuses': dict(Counter(r['output']['status'] for r in records)),
               'complete_counts': dict(Counter(str(r['output']['count']) for r in records if r['output']['status'] == 'complete')),
               'maximum_observed_box': max(r['output']['candidate_box_size'] for r in records),
               'scalar_candidates_checked': sum(r['output'].get('scalar_candidates_checked', 0) for r in records)}
        cases.append(row); print(json.dumps(row), flush=True)
    (HERE/'metric_summary.json').write_text(json.dumps({'status': 'produced', 'cases': cases, 'input_sha256': inputs,
        'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['metric_experiments.py', 'metric_erasure.py', 'conditioned_erasure.py']}}, indent=2)+'\n')


if __name__ == '__main__': main()
