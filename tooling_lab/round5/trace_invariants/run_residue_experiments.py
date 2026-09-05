#!/usr/bin/env python3
from fractions import Fraction
from hashlib import sha256
import json
from copy import deepcopy
from pathlib import Path
import sys
sys.dont_write_bytecode=True
from folded_moment import compiler,polynomial_from_values,encode
from residue_moment import seven_statistics,seven_stat_global_third_moment,universal_sixth_union

HERE=Path(__file__).resolve().parent


def exact_rank(rows):
    matrix=[[Fraction(x) for x in row] for row in rows];rank=0
    for col in range(len(matrix[0])):
        pivot=next((i for i in range(rank,len(matrix)) if matrix[i][col]),None)
        if pivot is None:continue
        matrix[rank],matrix[pivot]=matrix[pivot],matrix[rank]
        scale=matrix[rank][col];matrix[rank]=[x/scale for x in matrix[rank]]
        for i in range(len(matrix)):
            if i!=rank and matrix[i][col]:
                scale=matrix[i][col];matrix[i]=[a-scale*b for a,b in zip(matrix[i],matrix[rank])]
        rank+=1
    return rank


def main():
    baseline=json.loads((HERE/'results.json').read_text())
    outputs=[]
    for expected in baseline['comparisons']:
        name=expected['graph']
        record=json.loads((HERE/f'inventory_{name}.json').read_text())
        summary=seven_statistics(record)
        (HERE/f'seven_statistics_{name}.json').write_text(json.dumps(summary,indent=2)+'\n')
        actual=seven_stat_global_third_moment(summary,expected['n'])
        assert actual['third_moment']==expected['third_moment']
        actual['graph']=name;actual['matched_four_and_two_family_compilers']=True
        outputs.append(actual);print(name,actual['n'],actual['coefficient_evaluations'],flush=True)
    # Check raw union coefficients before choosing n or applying inclusion.
    audit=[]
    universal=universal_sixth_union()
    for q in (49,101,1297):
        hist=[];plans=[]
        for mono in (False,True):
            edges=(1,1,1) if mono else (1,1,-1)
            taus=compiler.admissible_tau_nodes(q,edges)
            assert len(taus)==8
            plans.append((mono,len(hist),[t if mono else -t for t in taus]))
            hist.extend(compiler.distinct_histogram(q,edges,t) for t in taus)
        raw,elapsed=compiler.union_coefficients_batch(hist,6,max_union=18)
        polynomials={}
        for mono,index,nodes in plans:
            polynomials[mono]=[]
            for k in range(19):
                values=[row['coefficients'][k] for row in raw[index:index+len(nodes)]]
                poly=polynomial_from_values(nodes[:7],values[:7])
                assert sum((c*nodes[-1]**j for j,c in enumerate(poly)),Fraction())==values[-1]
                assert poly[5]==0 and poly[6]==universal[k]
                polynomials[mono].append(poly)
        statistic_coefficients=[]
        for k in range(19):
            a,b=polynomials[True][k],polynomials[False][k]
            statistic_coefficients.append([a[1]+3*b[1],a[2]-b[2],a[3],b[3],a[4],b[4],universal[k]])
        audit.append({'q':q,'families':2,'union_degrees_checked':list(range(19)),
                      'theta5_zero':True,'theta6_equals_universal_polynomial':True,'extra_nodes_pass':True,
                      'seven_statistic_coefficient_rank':exact_rank(statistic_coefficients),
                      'rank_scope':'Exact formal linear dependence at this fixed q after the stated conference identities; not graph-realizability or historical minimality.',
                      'statistic_coefficient_rows_by_union_degree':[[encode(c) for c in row] for row in statistic_coefficients],
                      'backend_seconds':elapsed})
        print('raw union',q,flush=True)
    result={'scope':'Proved sufficient seven-statistic identity for degree6; no minimality or graph-realizability claim.',
            'comparisons':outputs,'raw_union_high_degree_checks':audit,
            'universal_theta6_union_coefficients':[encode(c) for c in universal],
            'source_sha256':{p.name:sha256(p.read_bytes()).hexdigest() for p in
                (Path(__file__),HERE/'residue_moment.py',HERE/'folded_moment.py')},
            'baseline_sha256':sha256((HERE/'results.json').read_bytes()).hexdigest()}
    valid=json.loads((HERE/'inventory_Paley101.json').read_text())
    malformed=[]
    changed=deepcopy(valid);changed['family_power_sums'][0]['powers'][3]+=1;malformed.append(changed)
    changed=deepcopy(valid);changed['family_power_sums'].append(changed['family_power_sums'][0]);malformed.append(changed)
    changed=deepcopy(valid);changed['family_power_sums'][0]['monochromatic']=0;malformed.append(changed)
    for record in malformed:
        try:seven_statistics(record)
        except AssertionError:pass
        else:raise AssertionError('malformed record accepted')
    result['malformed_record_rejections']=len(malformed)
    (HERE/'residue_results.json').write_text(json.dumps(result,indent=2)+'\n')


if __name__=='__main__':main()
