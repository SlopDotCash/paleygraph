#!/usr/bin/env python3
from collections import Counter
from hashlib import sha256
import json
from pathlib import Path
from conditioned_erasure import improve, decode
from pullback_erasure import improve as pullback

HERE = Path(__file__).resolve().parent


def main():
    path = HERE.parent/'round21'/'erasure_summary.json'; summary = json.loads(path.read_text())
    inputs = {'../round21/erasure_summary.json': sha256(path.read_bytes()).hexdigest()}; cases = []
    for row in summary['cases']:
        if row['p'] >= 100: continue
        path = HERE.parent/'round21'/row['artifact']; original = json.loads(path.read_text())
        inputs['../round21/'+path.name] = sha256(path.read_bytes()).hexdigest()
        c = improve(pullback(original['certificate'])); records = []
        for r in original['records']:
            out = decode(c, dict(r['known']), 50000)
            assert out['completions'] == r['output']['completions']
            records.append({**r, 'output': out})
        filename = 'small_'+row['name']+'.json'
        (HERE/filename).write_text(json.dumps({'name': row['name'], 'certificate': c, 'candidate_budget': 50000, 'records': records}, separators=(',', ':'))+'\n')
        result = {'name': row['name'], 'artifact': filename, 'queries': len(records),
                  'cuts': len(c['erased']), 'optimality_certificates': sum(x['optimality_certified'] for x in c['conditioning']['cuts']),
                  'complete_counts': dict(Counter(str(r['output']['count']) for r in records))}
        cases.append(result); print(json.dumps(result), flush=True)
    (HERE/'small_summary.json').write_text(json.dumps({'status': 'produced', 'cases': cases, 'input_sha256': inputs,
        'source_sha256': {n: sha256((HERE/n).read_bytes()).hexdigest() for n in ['small_experiments.py', 'conditioned_erasure.py', '../round21/pullback_erasure.py']}}, indent=2)+'\n')


if __name__ == '__main__': main()
