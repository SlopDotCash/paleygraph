#!/usr/bin/env python3
"""Cross-lane scale and information-limit experiments for exact third moments."""
from fractions import Fraction
from hashlib import sha256
import json
import math
from pathlib import Path
import sys

from third_moment import encode,global_third_moment

HERE=Path(__file__).resolve().parent
LAB=HERE.parents[1]
sys.path.insert(0,str(LAB/'round4/marked_moments'))
from conditional_moments import conditional_moments


def first_two(q,n):
    counts={'q':q,'marks':[],'cells':[{'pattern':[],'size':q}],
            'relations':[{'left_pattern':[],'right_pattern':[],'sign':s,
                          'count':q if s==0 else q*(q-1)//2} for s in (-1,0,1)]}
    return conditional_moments(counts,n)


def enrich(row):
    earlier=first_two(row['q'],row['n'])
    mean,second,variance=[Fraction(*earlier[key]) for key in ('mean','second_moment','variance')]
    central=Fraction(*row['third_moment'])-3*mean*second+2*mean**3
    normalized_squared=central**2/variance**3 if variance else None
    row.update({'mean':encode(mean),'second_moment':encode(second),'variance':encode(variance),
                'central_third_moment':encode(central),
                'standardized_third_squared':encode(normalized_squared) if normalized_squared is not None else None,
                'standardized_third_float':math.copysign(math.sqrt(float(normalized_squared)),central)
                                           if normalized_squared is not None and central else 0.0})
    return row


def main():
    scale=[]
    for q,n in ((13,6),(17,8),(29,6),(101,8),(1297,6),(1297,8),(65537,16),(1000033,31)):
        path=HERE.parent/'trace_backend'/f'inventory_p{q}.json'
        inventory=json.loads(path.read_text())
        row=enrich(global_third_moment(inventory,n))
        row.update({'source_inventory':str(path.relative_to(LAB)),
                    'source_inventory_sha256':sha256(path.read_bytes()).hexdigest(),
                    'direct_inventory_coefficient_evaluations':3+len(inventory['records'])})
        scale.append(row)
        print(f'q={q}, n={n}: third={row["third_moment_float"]:.9g}, standardized={row["standardized_third_float"]:.6g}, evaluations={row["coefficient_evaluations"]}',flush=True)
    comparisons=[]
    for q,n in ((101,8),(1297,8)):
        inventory=json.loads((HERE.parent/'trace_backend'/f'inventory_p{q}.json').read_text())
        direct=global_third_moment(inventory,n,method='direct_inventory')
        accelerated=next(x for x in scale if (x['q'],x['n'])==(q,n))
        assert direct['third_moment']==accelerated['third_moment']
        comparisons.append({'q':q,'n':n,'exact_equality':True,
                            'direct_evaluations':direct['coefficient_evaluations'],
                            'interpolation_evaluations':accelerated['coefficient_evaluations'],
                            'direct_elapsed_seconds':direct['elapsed_seconds'],
                            'interpolation_elapsed_seconds':accelerated['elapsed_seconds']})
    preflight=json.loads((HERE.parent/'preflight/third_moment_preflight.json').read_text())['cases']
    twins=[]
    for n in (6,7,8):
        reference=HERE.parent/'review'/f'twins_n{n}_third_oracle.json'
        oracle=json.loads(reference.read_text())['records']
        for graph in ('Paley49','Peisert49'):
            inventory=next(x for x in preflight if x['graph']==graph)
            row=enrich(global_third_moment(inventory,n))
            target=next(x for x in oracle if x['graph']==graph)
            assert [row['mean'],row['second_moment'],row['third_moment']]==target['moments']
            row.update({'graph':graph,'minimum':target['minimum'],'maximum':target['maximum'],
                        'subset_count':target['subset_count'],'oracle_sha256':sha256(reference.read_bytes()).hexdigest()})
            twins.append(row)
    n7=[x for x in twins if x['n']==7]
    n8=[x for x in twins if x['n']==8]
    assert n7[1]['third_moment_float']>n7[0]['third_moment_float'] and n7[1]['maximum']<n7[0]['maximum']
    assert n8[1]['third_moment_float']>n8[0]['third_moment_float'] and n8[1]['maximum']>n8[0]['maximum']
    output={'scope':'exact global moments and finite information-limit examples; no pointwise or prize bound',
            'scale_cases':scale,'direct_interpolation_comparisons':comparisons,'twins':twins,
            'finite_countercheck':'Larger raw third moment accompanies a smaller maximum at n7 and a larger maximum at n8 in these actual twins; no monotone upper-extremum ordering follows.',
            'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest()
                             for name in ('third_moment.py','union_coefficients.cpp','union_coefficients','run_experiments.py')},
            'previous_second_moment_compiler_sha256':sha256((LAB/'round4/marked_moments/conditional_moments.py').read_bytes()).hexdigest()}
    (HERE/'results.json').write_text(json.dumps(output,indent=2)+'\n')


if __name__=='__main__':main()
