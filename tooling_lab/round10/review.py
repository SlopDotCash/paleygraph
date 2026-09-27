#!/usr/bin/env python3
"""Small dense checks, then complete independent inward-target scale readback."""
from hashlib import sha256
import importlib.util
import json
from math import comb
from pathlib import Path
import subprocess
import time
from query import query
from deletion_pair_means import compile_means
from decomposition import components

HERE=Path(__file__).resolve().parent
OLD=HERE.parent/'round7/local_edits/review_results.py'
spec=importlib.util.spec_from_file_location('literal',OLD);literal=importlib.util.module_from_spec(spec);spec.loader.exec_module(literal)


def inward(q,c,d):
    payload=' '.join(map(str,[q,len(c),d,*c]))+'\n'
    p=subprocess.run([str(HERE/'inward_targets')],input=payload,text=True,capture_output=True,check=True)
    return json.loads(p.stdout)


def small():
    cases=[];entries=0;targets=0
    for q,c in [(5,[0,1]),(5,[0,1,2]),(13,list(range(11))),(17,[0,1,2,3,4,6,10]),
                (29,[0,1,3,4,7,11,18,20]),(257,sorted({x*x%257 for x in range(1,257)})[:64])]:
        S=literal.literal_signs(q)
        for d in range(min(6,len(c))+1):
            r=query(q,c,d);dense=compile_means(S,c,d)
            for k,v in dense.items():assert r[k]==v,(q,d,k)
            t=inward(q,c,d);decomp=components(q,c,d,t['induced'],t)
            for a,b in zip(r['pairs'],decomp['pairs']):assert a['target_sum']==b['target_sum'],(q,d,a,b)
            if q<=29:
                for record in t['pairs']:
                    A=[x for i,x in enumerate(c) if i not in record['indices']]
                    for j,v in enumerate(record['values']):
                        expected=literal.target(S,A,d-2*j) if d-2*j>=0 else 0
                        assert expected==v;targets+=1
            entries+=len(r['pairs']);cases.append({'q':q,'n':len(c),'degree':d,'pair_means':len(r['pairs'])})
    out={'status':'passed','cases':cases,'all_pair_means_compared':entries,'literal_inward_targets_compared':targets,
         'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('review.py','query.py','pair_means_backend.cpp','pair_means_backend','deletion_pair_means.py','decomposition.py','inward_targets.cpp','inward_targets')},
         'input_sha256':{'../round7/local_edits/review_results.py':sha256(OLD.read_bytes()).hexdigest()}}
    (HERE/'boundary_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps({'status':'passed','cases':len(cases),'pair_means':entries,'literal_targets':targets}),flush=True)


def scale():
    rows=json.loads((HERE/'scale_results.json').read_text())['cases'];out=[]
    for row in rows:
        begin=time.monotonic();t=inward(row['q'],row['selected'],row['degree']);decomp=components(row['q'],row['selected'],row['degree'],t['induced'],t)
        for a,b in zip(row['pairs'],decomp['pairs']):assert a['deleted']==b['deleted'] and a['target_sum']==b['target_sum']
        r={'q':row['q'],'n':len(row['selected']),'family':row['family'],'all_pair_means_checked':len(row['pairs']),
           'inward_targets':t,'decomposition':decomp,'seconds':time.monotonic()-begin}
        out.append(r);(HERE/'review_partial.json').write_text(json.dumps(out,indent=2)+'\n')
        print(json.dumps({k:v for k,v in r.items() if k not in ('inward_targets','decomposition')}),flush=True)
    result={'status':'passed','scope':'Complete independent readback of every pair mean through direct inward targets and induced boundary; no sampled large entries.',
            'cases':out,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('review.py','decomposition.py','inward_targets.cpp','inward_targets')},
            'input_sha256':{'scale_results.json':sha256((HERE/'scale_results.json').read_bytes()).hexdigest()}}
    (HERE/'scale_verification.json').write_text(json.dumps(result,indent=2)+'\n');(HERE/'review_partial.json').unlink()


if __name__=='__main__':
    import sys
    if sys.argv[1:] == ['scale']:scale()
    else:small()
