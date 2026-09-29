#!/usr/bin/env python3
"""Replay folded powers against the frozen general compiler's exact outputs."""
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path
import sys
sys.dont_write_bytecode=True
from folded_moment import fold_inventory,unfold_inventory,folded_global_third_moment,compiler,encode

HERE=Path(__file__).resolve().parent


def main():
    baseline_path=HERE.parent/'third_moment'/'results.json'
    baseline=json.loads(baseline_path.read_text())
    preflight=json.loads((HERE.parent/'preflight'/'third_moment_preflight.json').read_text())
    inputs={}
    for p in (13,17,29,101,1297,65537,1000033):
        inputs[f'Paley{p}']=json.loads((HERE.parent/'trace_backend'/f'inventory_p{p}.json').read_text())
    for row in preflight['cases']:
        if row['q']==49:inputs[row['graph']]=row
    for name,record in inputs.items():
        folded=fold_inventory(record)
        original=record.get('records',record.get('joint_edges_tau_histogram'))
        assert unfold_inventory(folded)==sorted(original,key=lambda row:(row['edges'],row['tau']))
        (HERE/f'inventory_{name}.json').write_text(json.dumps(folded,indent=2)+'\n')
    outputs=[]
    for expected in baseline['scale_cases']+baseline['twins']:
        name=expected.get('graph',f"Paley{expected['q']}")
        folded=fold_inventory(inputs[name])
        # Deliberately remove the inventory: the consumer really uses powers.
        del folded['records']
        actual=folded_global_third_moment(folded,expected['n'],expected['degree'])
        assert actual['third_moment']==expected['third_moment'],(name,actual,expected)
        actual.update(graph=name,matched_general_compiler=True,general_coefficient_evaluations=expected['coefficient_evaluations'])
        outputs.append(actual)
        print(name,actual['n'],actual['third_moment'],actual['coefficient_evaluations'],flush=True)
    api_checks=[]
    for degree in (0,2,4):
        expected=compiler.global_third_moment(inputs['Paley13'],7,degree)
        actual=folded_global_third_moment(inputs['Paley13'],7,degree)
        assert actual['third_moment']==expected['third_moment']
        api_checks.append({'q':13,'n':7,'degree':degree,'third_moment':actual['third_moment']})
    try:folded_global_third_moment(inputs['Paley13'],7,3)
    except AssertionError:pass
    else:raise AssertionError('odd degree must be rejected')
    # On q101 both families have >=7 nodes, so their degree6 polynomials
    # are uniquely defined, not merely their values on a tiny finite domain.
    comparisons=[]
    for n in (6,7,8):
        result=next((v for v in outputs if v['q']==101 and v['n']==n),None)
        if result is None:result=folded_global_third_moment(inputs['Paley101'],n)
        polys={v['monochromatic']:[Fraction(*x) for x in v['polynomial_in_theta']] for v in result['family_plans']}
        difference=[a-b for a,b in zip(polys[True],polys[False])]
        comparisons.append({'q':101,'n':n,'mono_minus_mixed_polynomial':[encode(c) for c in difference],
                            'coincide':all(c==0 for c in difference),
                            'polynomial_extension_comparison':'The two actual theta domains are disjoint mod8.'})
    assert comparisons[0]['coincide']
    result={'scope':'Exact symmetry compression for conditional third-moment coefficients; no pointwise or prize bound.',
            'all_nine_folded_inventories_unfold_exactly':True,
            'comparisons':outputs,'even_degree_api_checks':api_checks,'odd_degree_rejected':True,
            'family_polynomial_comparison':comparisons,
            'input_sha256':{str(p.relative_to(HERE.parent)):sha256(p.read_bytes()).hexdigest() for p in
                             (baseline_path,HERE.parent/'preflight'/'third_moment_preflight.json')},
            'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'folded_moment.py')}}
    (HERE/'results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
