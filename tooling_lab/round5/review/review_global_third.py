#!/usr/bin/env python3
"""Compare the accelerated compiler to independent complete-set oracles."""
import sys
sys.dont_write_bytecode=True
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations
from pathlib import Path
import importlib.util
import json
import time

HERE=Path(__file__).resolve().parent


def load(name,path):
    spec=importlib.util.spec_from_file_location(name,path)
    module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)
    return module


def matrix9():
    def sub(a,b):return (a%3-b%3)%3+3*((a//3-b//3)%3)
    def mul(a,b):return (a%3*(b%3)-a//3*(b//3))%3+3*((a%3*(b//3)+a//3*(b%3))%3)
    squares={mul(x,x) for x in range(1,9)}
    return [[0 if a==b else 1 if sub(a,b) in squares else -1 for b in range(9)] for a in range(9)]


def normalized_ordered_check(S,record):
    """Check the aggregate normalization directly, including all six orders."""
    q=len(S);actual=Counter()
    for a,b,c in combinations(range(q),3):
        tau=sum(S[a][x]*S[b][x]*S[c][x] for x in range(q))
        for u,v,w in permutations((a,b,c)):
            sign=S[u][v]
            actual[(1,sign*S[u][w],sign*S[v][w]),sign*tau]+=1
    wanted=Counter({(tuple(r['edges']),r['tau']):r['count']*q*(q-1) for r in record['records']})
    assert actual==wanted
    return sum(actual.values())


def main():
    started=time.perf_counter()
    source=HERE.parent/'third_moment'/'third_moment.py'
    source_hash=sha256(source.read_bytes()).hexdigest()
    root=load('reviewed_third_moment',source)
    oracle=load('independent_direct_oracle',HERE/'direct_third_oracle.py')
    base=json.loads((HERE/'direct_third_oracle.json').read_text())
    inventories={r['graph']:r for r in base['normalized_joint_records']}
    matrices={'Paley13':oracle.prime(13),'Paley17':oracle.prime(17),**oracle.twins(),'Paley9':matrix9()}
    inventories['Paley9']=oracle.normalized_record('Paley9',matrices['Paley9'])
    ordered_checks=sum(normalized_ordered_check(S,inventories[name]) for name,S in matrices.items())
    records=list(base['records'])
    for n in (6,7,8):records+=json.loads((HERE/f'twins_n{n}_third_oracle.json').read_text())['records']
    reports=[];odd_rejected=0
    for row in records:
        degree=row['degree']
        if degree not in (0,2,4,6):
            if degree%2:
                assert row['moments'][2]==[0,1]
                try:root.global_third_moment(inventories[row['graph']],row['n'],degree)
                except AssertionError:odd_rejected+=1
                else:raise AssertionError('odd normalized degree accepted')
            continue
        pair=[]
        for method in ('direct_inventory','interpolate'):
            out=root.global_third_moment(inventories[row['graph']],row['n'],degree,method)
            assert out['third_moment']==row['moments'][2],(row,out)
            pair.append({'method':method,'coefficient_evaluations':out['coefficient_evaluations'],
                         'interpolated_classes':sum(x['method']=='proved-degree interpolation' for x in out['class_plans']),
                         'extra_node_checks':sum(x['extra_node_checked'] for x in out['class_plans'])})
        reports.append({key:row[key] for key in ('graph','q','n','degree')}|{'third_moment':row['moments'][2],'comparisons':pair})
    # Degree-six finite differences on raw union coefficients at synthetic,
    # nonnegative integral inventories. This regression check uses an
    # independent finite-difference criterion, not the root Lagrange function.
    synthetic=[]
    for b in (-1,1):
        for c in (-1,1):
            edges=(1,b,c)
            valid=[]
            for tau in range(-98,99):
                try:hist=oracle.triple_histogram(101,edges,tau)
                except AssertionError:continue
                valid.append((tau,hist))
            assert len(valid)>=8
            start=(len(valid)-8)//2
            chosen=valid[start:start+8]
            assert all(chosen[i+1][0]-chosen[i][0]==8 for i in range(7))
            raw,_=root.union_coefficients_batch([h for _,h in chosen],degree=6,max_union=8)
            tested=0
            for k in range(9):
                differences=[r['coefficients'][k] for r in raw]
                for _ in range(7):differences=[b-a for a,b in zip(differences,differences[1:])]
                assert differences==[0]
                tested+=1
            synthetic.append({'edges':list(edges),'nodes':[t for t,_ in chosen],
                              'seventh_finite_differences_zero':tested,'graph_realizability_claimed':False})
    # Optional cached moments must be complete and unique if supplied.
    valid=json.loads(json.dumps(inventories['Paley13']))
    valid['edge_power_sums']=[]
    for b in (-1,1):
        for c in (-1,1):
            edges=[1,b,c]
            valid['edge_power_sums'].append({'edges':edges,'powers':[
                sum(r['count']*r['tau']**j for r in valid['records'] if r['edges']==edges)
                for j in range(7)]})
    root.validate_inventory(valid)
    metadata_rejections=0
    for kind in ('missing','duplicate','wrong_value','empty'):
        bad=json.loads(json.dumps(valid))
        if kind=='missing':bad['edge_power_sums'].pop()
        elif kind=='duplicate':bad['edge_power_sums'][-1]=bad['edge_power_sums'][0]
        elif kind=='wrong_value':bad['edge_power_sums'][0]['powers'][6]+=1
        else:bad['edge_power_sums']=[]
        try:root.validate_inventory(bad)
        except AssertionError:metadata_rejections+=1
        else:raise AssertionError(('corrupt metadata accepted',kind))
    assert sha256(source.read_bytes()).hexdigest()==source_hash,'source changed during review; rerun required'
    paths=['direct_third_oracle.json','twins_n6_third_oracle.json','twins_n7_third_oracle.json','twins_n8_third_oracle.json']
    out={'status':'all accelerated third moments match independent exhaustive subset oracles',
         'finite_cases':len(reports),'compiler_comparisons':2*len(reports),'odd_degree_cases_rejected':odd_rejected,
         'ordered_row_triples_normalized_independently':ordered_checks,'records':reports,
         'malformed_optional_metadata_rejections':metadata_rejections,
         'synthetic_degree_bound_regressions':synthetic,
         'source_sha256':{'../third_moment/third_moment.py':source_hash,
                          '../third_moment/union_coefficients.cpp':sha256((source.parent/'union_coefficients.cpp').read_bytes()).hexdigest(),
                          Path(__file__).name:sha256(Path(__file__).read_bytes()).hexdigest()},
         'oracle_artifact_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in paths},
         'elapsed_seconds':time.perf_counter()-started}
    (HERE/'review_global_third.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:v for k,v in out.items() if k not in ('records','oracle_artifact_sha256')},indent=2))


if __name__=='__main__':main()
