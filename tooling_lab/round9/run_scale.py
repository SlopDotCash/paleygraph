#!/usr/bin/env python3
"""Run fixed deletion pair queries on the previously frozen eight actual sets."""
from hashlib import sha256
import json
from pathlib import Path
from two_insertions import query

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round7/local_edits/results.json'


def main():
    cases=[]
    for prev in json.loads(OLD.read_text())['scale_cases']:
        r=query(prev['q'],prev['selected'],prev['selected'][:2]);r['family']=prev['family']
        cases.append(r)
        (HERE/'scale_partial.json').write_text(json.dumps(cases,indent=2)+'\n')
        print(json.dumps({'q':prev['q'],'n':prev['n'],'family':r['family'],'seconds':r['seconds'],
                          'distinct_final_sets':r['distinct_final_sets'],'Q':r['statistics']['Q'],
                          'variance':r['variance'],'Q_variance_correction':r['Q_variance_correction']}),flush=True)
    names=('two_insertions.py','two_insertions_backend.cpp','two_insertions_backend','run_scale.py')
    out={'status':'passed','scope':'Exact moments for one fixed deletion pair in each actual input. No full two-swap shell census at large q.',
         'cases':cases,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in names},
         'input_sha256':{'../round7/local_edits/results.json':sha256(OLD.read_bytes()).hexdigest()}}
    (HERE/'scale_results.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'scale_partial.json').unlink()


if __name__=='__main__':main()
