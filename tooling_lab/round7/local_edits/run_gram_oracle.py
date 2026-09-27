#!/usr/bin/env python3
"""Run the separate compiled Gram method on fixed inputs; save every contraction."""
from hashlib import sha256
import json
from pathlib import Path
import random
import subprocess
from time import perf_counter

HERE = Path(__file__).resolve().parent


def run(q, selected, degree=6):
    text = f'{q} {len(selected)} {degree}\n' + ' '.join(map(str, selected)) + '\n'
    start = perf_counter()
    result = subprocess.run([str(HERE/'gram_oracle')], input=text, text=True, capture_output=True, check=True)
    row = json.loads(result.stdout)
    row['selected'] = sorted(selected)
    row['seconds'] = perf_counter()-start
    return row


def main():
    rng = random.Random(73018)
    records = []
    for q,n in ((1297,6),(65537,16),(1000033,31),(6700417,50)):
        for label, selected in [('arithmetic_progression',list(range(n))),('seeded_uniform',sorted(rng.sample(range(q),n)))]:
            record = run(q, selected)
            record['family'] = label
            records.append(record)
            (HERE/'gram_scale_partial.json').write_text(json.dumps(records,indent=2)+'\n')
            print(json.dumps({'q':q,'n':n,'family':label,'seconds':record['seconds']}),flush=True)
    witnesses = json.loads((HERE/'toy_results.json').read_text())['collision_search']['witness']
    for key in ('lower','upper'):
        source = witnesses[key]
        row = run(source['q'],source['selected'])
        assert row['target'] == source['target'] and row['internal_contraction'] == source['internal_contraction']
        assert row['insertion_sums'] == [r['insertion_sum'] for r in source['deletion_records']]
        assert row['insertion_norms_squared'] == [r['insertion_norm_squared'] for r in source['deletion_records']]
    out = {'status':'passed', 'scope':'Separate synthetic-division and symmetric weighted-Gram arithmetic; full local contraction data, not a third direct neighbour census.',
           'scale_cases':records, 'toy_witnesses_checked':2,
           'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('gram_oracle.cpp','gram_oracle','run_gram_oracle.py')},
           'input_sha256':{'toy_results.json':sha256((HERE/'toy_results.json').read_bytes()).hexdigest()}}
    (HERE/'gram_results.json').write_text(json.dumps(out,indent=2)+'\n')
    (HERE/'gram_scale_partial.json').unlink()
    print('Completed independent Gram backend.',flush=True)


if __name__=='__main__':
    main()
