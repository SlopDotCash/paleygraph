#!/usr/bin/env python3
"""Run the transcript-aware pruner on actual words and lower thresholds."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from arc_certificate import certify_blocks,find_cover
from arc_pruning import prepare,run

HERE=Path(__file__).resolve().parent


def main():
    sourcepath=HERE.parent/'round6/pencil_tracks/results.json';large=next(x['result'] for x in json.loads(sourcepath.read_text())['large'] if x['name']=='skew-two')
    smallpath=HERE.parent/'round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json';small=json.loads(smallpath.read_text())
    priorpath=HERE.parent/'round12/pruning_results.json';prior=json.loads(priorpath.read_text())
    configs=[('small',None,11)]+[(f'd{d}_s{s}',d,s) for d,s in [(3,410),(4,410),(5,410),(3,192),(4,192),(3,128)]]
    cases=[];artifacts={}
    for name,dimension,threshold in configs:
        start=perf_counter()
        if name=='small':
            zero=json.loads((HERE/'zero_results.json').read_text());data=zero['input'];c=certify_blocks(data,zero['blocks']);search={'status':'existing_blocks'};scalars=range(17);raw=small
        else:
            old=json.loads((HERE/f'dimension{dimension}.preflight.json').read_text());data={**old['input'],'s':threshold};scalars=(0,1);raw=large
            if threshold==410:c=certify_blocks(data,old['blocks']);search={'status':'existing_blocks'}
            else:
                search=find_cover(data,seed=1407,max_attempts=128);assert search['status']=='found',search
                c=search.pop('certificate')
        construction_seconds=perf_counter()-start
        cfile=HERE/f'{name}.certificate.json';cfile.write_text(json.dumps(c,separators=(',',':'))+'\n');artifacts[cfile.name]=sha256(cfile.read_bytes()).hexdigest()
        start=perf_counter();prepared=prepare(c);preparation_seconds=perf_counter()-start
        records=[];p=data['p'];k=data['k']
        for scalar in scalars:
            received=[(a+scalar*b)%p for a,b in zip(raw['u0'],raw['u1'])];result=run(prepared,received)
            if name=='small':
                truth=prior['small_runs'][scalar]['oracle_parameters'];assert truth==[r['parameters'] for r in result['filter']['outputs']]
            else:
                pieces=[]
                for track in large['verified_tracks']:
                    a=track['intercept_coefficients'];b=track['slope_coefficients']
                    pieces.append([((a[j] if j<len(a) else 0)+scalar*(b[j] if j<len(b) else 0))%p for j in range(k)])
                assert sorted(pieces)==sorted(r['coefficients'] for r in result['filter']['outputs'])
                assert 2*(k-1)<threshold
            path=HERE/f'{name}_scalar{scalar}.run.json';path.write_text(json.dumps(result,separators=(',',':'))+'\n');artifacts[path.name]=sha256(path.read_bytes()).hexdigest()
            record={'scalar':scalar,'run_file':path.name,'candidate_count':len(result['transcript']['candidates']),
                    'queries':result['transcript']['queries_executed'],'query_seconds':result['query_seconds'],'filter_seconds':result['filter_seconds'],
                    **{key:result['filter'][key] for key in ('candidate_coordinate_evaluations','output_word_coordinate_evaluations','total_explicit_coordinate_evaluations','full_scan_coordinate_evaluations')},
                    'output_count':len(result['filter']['outputs'])}
            records.append(record)
            print(json.dumps({'case':name,**record}),flush=True)
        cases.append({'name':name,'certificate':cfile.name,'search':search,'construction_seconds':construction_seconds,'verified_preparation_seconds':preparation_seconds,
                      'review':prepared['review'],'ambient_johnson_comparison':{'s_squared':data['s']**2,'k_minus_one_times_n':(k-1)*data['n'],'below_classical_agreement_threshold':data['s']**2<(k-1)*data['n']},
                      'runs':records})
    files=['arc_certificate.py','arc_profile.py','dimension_preflight.py','arc_verifier.py','arc_pruning.py','run_arc.py']
    out={'status':'passed','scope':'Actual words in declared spaces, including deliberately reduced agreements128 and192 below the ambient Johnson threshold. Supplied-space completeness only; these special words retain a complete two-piece root-bound oracle. Query-count reductions and coordinate-evaluation reductions do not establish new general decoding or historical novelty.',
         'cases':cases,'artifact_sha256':artifacts,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in files},
         'input_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['zero_results.json','dimension3.preflight.json','dimension4.preflight.json','dimension5.preflight.json','../round12/pruning_results.json','../round6/pencil_tracks/results.json','../round6/overlap_preflight/coset_small_characteristic_sample3.certificate.json']}}
    (HERE/'actual_arc_results.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(cases),'received_words':sum(len(r['runs']) for r in cases)}),flush=True)


if __name__=='__main__':main()
