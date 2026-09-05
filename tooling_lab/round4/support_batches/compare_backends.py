#!/usr/bin/env python3
"""Measured iteration: compare the archived and cardinality-pruned ZDD engines."""
from hashlib import sha256
from pathlib import Path
from statistics import median
import json
import sys
import warnings

sys.dont_write_bytecode = True
warnings.filterwarnings('ignore', category=SyntaxWarning)
import symbolic_batches_initial as initial
import symbolic_batches as final

HERE = Path(__file__).resolve().parent


def main():
    saved = json.loads((HERE / 'results.json').read_text())
    comparisons = []
    for case in saved['cost_experiments']:
        r = case['result']
        args = (final.PrimeField(r['field']['p']), r['domain'], r['k'], r['s'], r['u0'], r['u1'])
        timings = {'initial': [], 'pruned': []}
        runs = {}
        for repeat in range(3):
            # Alternate execution order to reduce systematic first-run bias.
            names = ['initial', 'pruned'] if repeat % 2 == 0 else ['pruned', 'initial']
            for name in names:
                engine = initial if name == 'initial' else final
                runs[name] = engine.run_census(*args)
                timings[name].append(runs[name]['elapsed_seconds'])
        for key in ['nodes', 'batch_certificates', 'finite_track_certificates',
                    'all_field_track_certificates', 'bad_scalars']:
            assert runs['initial'][key] == runs['pruned'][key]
        comparisons.append({'case': case['name'], 'all_census_and_batch_certificates_equal': True,
            'initial_node_visits': runs['initial']['stats']['zdd_subtraction_node_visits'],
            'pruned_node_visits': runs['pruned']['stats']['zdd_subtraction_node_visits'],
            'median_seconds': {name: median(values) for name, values in timings.items()},
            'timing_samples': timings})
    output = {'status': 'passed', 'iteration': 'minimum remaining cardinality pruning',
        'comparisons': comparisons,
        'limits': 'Exact operation counts are reproducible; timings depend on concurrent workloads.',
        'source_sha256': {name: sha256((HERE / name).read_bytes()).hexdigest()
                         for name in ['symbolic_batches_initial.py', 'symbolic_batches.py', 'compare_backends.py']}}
    (HERE / 'iteration_results.json').write_text(json.dumps(output, indent=2) + '\n')
    print(json.dumps(output, indent=2))


if __name__ == '__main__':
    main()
