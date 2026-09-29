#!/usr/bin/env python3
"""Check identical full decision traces with the shared-order implementation."""
from hashlib import sha256
import json
from pathlib import Path
from time import perf_counter
from arc_pruning import prepare
from arc_filter_fast import filter_candidates

HERE=Path(__file__).resolve().parent


def main():
    sourcepath=HERE/'actual_arc_results.json';source=json.loads(sourcepath.read_text());rows=[];bindings={sourcepath.name:sha256(sourcepath.read_bytes()).hexdigest()}
    for case in source['cases']:
        cpath=HERE/case['certificate'];prepared=prepare(json.loads(cpath.read_text()));bindings[cpath.name]=sha256(cpath.read_bytes()).hexdigest()
        for record in case['runs']:
            path=HERE/record['run_file'];result=json.loads(path.read_text());bindings[path.name]=sha256(path.read_bytes()).hexdigest()
            start=perf_counter();fast=filter_candidates(prepared,result['transcript']);seconds=perf_counter()-start
            assert fast==result['filter']
            row={'case':case['name'],'scalar':record['scalar'],'decisions_matched':len(fast['decisions']),
                 'observed_shared_order_filter_seconds':seconds,'earlier_reference_filter_seconds':result['filter_seconds']}
            rows.append(row);print(json.dumps(row),flush=True)
    controls=json.loads((HERE/'arc_controls.json').read_text());boundary=controls['boundary_absence_and_single_occurrence_control'];prepared=prepare(boundary['certificate'])
    assert filter_candidates(prepared,boundary['run']['transcript'])==boundary['run']['filter']
    bindings['arc_controls.json']=sha256((HERE/'arc_controls.json').read_bytes()).hexdigest()
    stress=json.loads((HERE/'stress_results.json').read_text());prepared=prepare(json.loads((HERE/stress['certificate']).read_text()))
    for row in stress['cases']:
        path=HERE/row['run_file'];result=json.loads(path.read_text());assert filter_candidates(prepared,result['transcript'])==result['filter'];bindings[path.name]=sha256(path.read_bytes()).hexdigest()
    bindings['stress_results.json']=sha256((HERE/'stress_results.json').read_bytes()).hexdigest()
    out={'status':'passed','scope':'Every decision, coordinate-test trace and output matches the reference filter on all actual runs plus boundary and stress controls. Timings are observations from separate runs on a shared machine; they are not a controlled speed comparison or a bound on total decoding cost.',
         'cases':rows,'actual_decisions_matched':sum(r['decisions_matched'] for r in rows),'boundary_and_stress_runs_matched':4,
         'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['compare_filter.py','arc_filter_fast.py','arc_pruning.py','arc_verifier.py']},
         'input_sha256':bindings}
    (HERE/'filter_comparison.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
