#!/usr/bin/env python3
"""Validate the finite shared-insertion witnesses, then extend to frozen scale inputs."""
from fractions import Fraction as F
from hashlib import sha256
import json
from pathlib import Path
from covariance import query

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]


def main():
    toy=json.loads((HERE/'results.json').read_text());tiny=[]
    for pair in toy['pairs']:
        for member in pair['members']:
            r,g=query(member['q'],member['selected'])
            for key in ('target','derivative_gram','insertion_cross_second_sum','covariance_trace_squared'):
                assert r[key]==member[key],key
            for i,row in enumerate(member['insertion_covariance']):
                for j,value in enumerate(row):assert F(*value)==F(r['covariance_numerator'][i][j],r['covariance_denominator'])
            tiny.append(r)
    old_path=LAB/'round7/local_edits/results.json';old=json.loads(old_path.read_text());rows=[]
    for prev in old['scale_cases']:
        row,g=query(prev['q'],prev['selected']);row['family']=prev['family']
        assert row['internal_contraction']==prev['internal_contraction']
        assert g['insertion_sums']==[x['insertion_sum'] for x in prev['deletion_records']]
        assert g['insertion_norms_squared']==[x['insertion_norm_squared'] for x in prev['deletion_records']]
        assert row['marginal_neighbour_mean']==prev['neighbour_mean'] and row['marginal_neighbour_variance']==prev['neighbour_variance']
        rows.append(row)
        (HERE/'scale_partial.json').write_text(json.dumps(rows,indent=2)+'\n')
        print(json.dumps({'q':row['q'],'n':row['n'],'family':row['family'],'seconds':row['seconds'],
                          'off_diagonal_energy_fraction':float(F(*row['off_diagonal_squared_energy_fraction']))}),flush=True)
    names=('covariance_backend.cpp','covariance_backend','covariance.py','run_scale.py')
    out={'status':'passed','scope':'Exact shared-insertion covariance; complete tiny direct checks and inherited marginal comparisons. New large off-diagonal Gram entries receive a separate review.',
         'tiny_witnesses':tiny,'scale_cases':rows,
         'source_sha256':{f'round8/shared_insertions/{name}':sha256((HERE/name).read_bytes()).hexdigest() for name in names},
         'input_sha256':{'round8/shared_insertions/results.json':sha256((HERE/'results.json').read_bytes()).hexdigest(),
                         'round7/local_edits/results.json':sha256(old_path.read_bytes()).hexdigest()}}
    (HERE/'scale_results.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'scale_partial.json').unlink()
    print('Completed covariance scale run.',flush=True)


if __name__=='__main__':main()
