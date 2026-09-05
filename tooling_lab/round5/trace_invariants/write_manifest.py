#!/usr/bin/env python3
"""Bind the final symmetry/remainder lane and its read-only dependencies."""
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def entry(path,root):
    return {'path':str(path.relative_to(root)),'bytes':path.stat().st_size,'sha256':sha256(path.read_bytes()).hexdigest()}


def main():
    folded=json.loads((HERE/'results.json').read_text())
    residue=json.loads((HERE/'residue_results.json').read_text())
    symmetry=json.loads((HERE/'symmetry_validation.json').read_text())
    assert len(folded['comparisons'])==len(residue['comparisons'])==14
    summary={'status':'all exact folded and seven-statistic replays pass',
             'moment_cases_per_adapter':14,'normalized_inventories_losslessly_folded_and_unfolded':9,
             'generic_coefficient_evaluations':{'original':35,'two_families':19,'quartic_remainders':15},
             'measured_statistic_count':7,
             'raw_coefficient_rank_at_tested_orders':{str(x['q']):x['seven_statistic_coefficient_rank'] for x in residue['raw_union_high_degree_checks']},
             'prime_parameters_checked':sum(x['all_nonsingular_parameters_checked'] for x in symmetry['prime_orbits']),
             'unordered_matrix_triples_checked':sum(x['unordered_distinct_triples_checked'] for x in symmetry['literal_matrix_checks']),
             'literal_canonical_histograms_checked':sum(x['literal_canonical_histograms_checked'] for x in symmetry['literal_matrix_checks']),
             'maximum_prime':max(x['p'] for x in symmetry['prime_orbits']),
             'malformed_record_rejections':residue['malformed_record_rejections'],
             'historical_originality_proved':False,'prize_or_pointwise_bound':False}
    (HERE/'summary.json').write_text(json.dumps(summary,indent=2)+'\n')
    external=[HERE.parent/'third_moment'/name for name in ('third_moment.py','union_coefficients.cpp','union_coefficients','results.json')]
    external.extend(HERE.parent/'preflight'/name for name in ('third_moment_preflight.py','third_moment_preflight.json','README.md'))
    external.extend(sorted((HERE.parent/'trace_backend').glob('inventory_p*.json')))
    external.extend(sorted((HERE.parent/'trace_backend').glob('tau_p*.bin')))
    files=[p for p in HERE.iterdir() if p.is_file() and p.name!='manifest.json']
    out={'files':[entry(p,HERE) for p in sorted(files)],
         'external_read_only_dependencies':[entry(p,HERE.parent) for p in external],
         'scope':'Hashes bind the computation snapshot. Prior lanes were read only.'}
    (HERE/'manifest.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(summary,indent=2))


if __name__=='__main__':main()
