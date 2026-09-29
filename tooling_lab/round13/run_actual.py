#!/usr/bin/env python3
"""Find and verify covers, then run actual received words with prior complete oracles."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from query_cover import find_partition
from deterministic_pruning import prepare,run

HERE=Path(__file__).resolve().parent


def main():
    previous=HERE.parent/'round12';prior=json.loads((previous/'pruning_results.json').read_text())
    cases=[];certificate_bindings={}
    for name,source,sizes,runs in (
        ('small','hard_rank3.certificate.json',[4,4,3,3,2],prior['small_runs']),
        ('full_length','full_length.certificate.json',[5]*201+[6]*3+[1],prior['large_runs'])):
        data=json.loads((previous/source).read_text())['input']
        start=perf_counter();search=find_partition(data,sizes,seed=1306);construction_seconds=perf_counter()-start
        assert search['status']=='found';c=search.pop('certificate')
        certificate=HERE/f'{name}.certificate.json';certificate.write_text(json.dumps(c,separators=(',',':'))+'\n')
        certificate_bindings[certificate.name]=sha256(certificate.read_bytes()).hexdigest()
        start=perf_counter();prepared=prepare(c);preparation_seconds=perf_counter()-start
        records=[]
        for old in runs:
            start=perf_counter();result=run(prepared,old['run']['received']);seconds=perf_counter()-start
            if name=='small':assert old['oracle_parameters']==[o['parameters'] for o in result['outputs']]
            else:assert sorted(old['oracle_polynomials'])==sorted(o['coefficients'] for o in result['outputs'])
            records.append({'scalar':old['scalar'],'run_seconds':seconds,'run':result})
        cases.append({'name':name,'source':f'../round12/{source}','certificate':certificate.name,
                      'search':search,'block_sizes':sizes,'construction_seconds':construction_seconds,
                      'verified_preparation_seconds':preparation_seconds,'coverage_review':prepared['review'],'runs':records})
        print(json.dumps({'case':name,'queries':c['query_count'],'capacity':c['no_query_capacity'],
                          'attempts':search['attempts'],'construction_seconds':construction_seconds,
                          'verified_preparation_seconds':preparation_seconds,
                          'run_seconds':[r['run_seconds'] for r in records],
                          'output_counts':[len(r['run']['outputs']) for r in records]}),flush=True)
    files=['query_cover.py','cover_verifier.py','deterministic_pruning.py','run_actual.py']
    out={'status':'passed','scope':'Deterministic query-cover certificates and actual complete-list controls for two declared affine spaces. Construction is bounded search, not a universal algorithm or minimum-cover proof.',
         'cases':cases,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in files},
         'certificate_sha256':certificate_bindings,
         'input_sha256':{f'../round12/{f}':sha256((previous/f).read_bytes()).hexdigest() for f in ['pruning_results.json','hard_rank3.certificate.json','full_length.certificate.json']}}
    (HERE/'actual_results.json').write_text(json.dumps(out,separators=(',',':'))+'\n')


if __name__=='__main__':main()
