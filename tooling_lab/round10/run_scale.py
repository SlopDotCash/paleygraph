#!/usr/bin/env python3
from hashlib import sha256
import json
from pathlib import Path
from query import query
from deletion_pair_means import compile_means

HERE=Path(__file__).resolve().parent


def main():
    oldpath=HERE.parent/'round9/scale_results.json';old=json.loads(oldpath.read_text());cases=[]
    for prev in old['cases']:
        s=prev['statistics'];r=query(s['q'],s['selected']);r['family']=prev['family']
        first=next(p for p in r['pairs'] if p['deleted']==s['deleted'])
        assert first['mean']==prev['mean'] and first['target_sum']==prev['target_sum']
        cases.append(r);(HERE/'scale_partial.json').write_text(json.dumps(cases,indent=2)+'\n')
        print(json.dumps({'q':r['q'],'n':len(r['selected']),'family':r['family'],'pair_means':len(r['pairs']),'seconds':r['seconds']}),flush=True)
    out={'status':'passed','cases':cases,'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('run_scale.py','query.py','pair_means_backend.cpp','pair_means_backend')},
         'input_sha256':{'../round9/scale_results.json':sha256(oldpath.read_bytes()).hexdigest()}}
    (HERE/'scale_results.json').write_text(json.dumps(out,indent=2)+'\n');(HERE/'scale_partial.json').unlink()


if __name__=='__main__':main()
