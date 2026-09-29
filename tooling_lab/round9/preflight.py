#!/usr/bin/env python3
"""Compare moment compiler with literal subset products on complete tiny shells."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
import importlib.util
import json
from pathlib import Path
import numpy as np
from two_insertions import from_matrix, rat

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round7/local_edits/review_results.py'
spec=importlib.util.spec_from_file_location('literal_review',OLD)
literal=importlib.util.module_from_spec(spec);spec.loader.exec_module(literal)


def check(S,c,a,d):
    result=from_matrix(S,c,a,d)
    outside=sorted(set(range(len(S)))-set(c));base=set(c)-set(a)
    values=[literal.target(S,sorted(base|{b,e}),d) for b,e in combinations(outside,2)]
    assert result['target_sum']==2*sum(values)
    assert result['target_square_sum']==2*sum(v*v for v in values)
    assert result['distinct_final_sets']==len(values)
    mean=F(sum(values),len(values));var=F(sum(v*v for v in values),len(values))-mean*mean
    assert result['mean']==rat(mean) and result['variance']==rat(var)
    return result,values


def main():
    cases=[]; evaluations=0
    for q,c in [(5,[0,1]),(5,[0,1,2]),(13,list(range(6))),(17,list(range(7))),
                (13,list(range(11))),(29,[0,1,3,4,7,11,18,20])]:
        S=literal.literal_signs(q)
        for d in range(min(6,len(c))+1):
            for a in list(combinations(c,2))[:3]:
                r,v=check(S,c,a,d);evaluations+=len(v)
                cases.append({'q':q,'selected':c,'deleted':a,'degree':d,'mean':r['mean'],'variance':r['variance']})
    twins=[[[0,1,2,3,4,6,10],[0,1,2,3,4,6,15]],
           [[0,1,2,4,7,12,14],[0,1,2,5,6,7,12]]]
    S=literal.literal_signs(17);pairs=[]
    for pair in twins:
        members=[]
        for c in pair:
            values=[];fibres=[]
            for a in combinations(c,2):
                r,v=check(S,c,a,6);values.extend(v);evaluations+=len(v)
                fibres.append({'deleted':a,'histogram':sorted(Counter(v).items()),'result':r})
            members.append({'selected':c,'histogram':sorted(Counter(values).items()),
                            'mean':rat(F(sum(values),len(values))),
                            'variance':rat(F(sum(v*v for v in values),len(values))-F(sum(values),len(values))**2),
                            'fibres':fibres})
        pairs.append({'members':members,'equal_two_swap_histograms':members[0]['histogram']==members[1]['histogram'],
                      'equal_two_swap_variances':members[0]['variance']==members[1]['variance'],
                      'equal_unlabelled_fixed_deletion_histograms':sorted(x['histogram'] for x in members[0]['fibres'])==sorted(x['histogram'] for x in members[1]['fibres'])})
    out={'status':'passed','boundary_cases':cases,'literal_distinct_insertion_pair_evaluations':evaluations,'pairs':pairs,
         'source_sha256':{f.name:sha256(f.read_bytes()).hexdigest() for f in (Path(__file__),HERE/'two_insertions.py')},
         'input_sha256':{'../round7/local_edits/review_results.py':sha256(OLD.read_bytes()).hexdigest()}}
    (HERE/'preflight_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases),'literal_values':evaluations,
                      'pair_outcomes':[{k:v for k,v in p.items() if k!='members'} for p in pairs]}),flush=True)


if __name__=='__main__':main()
