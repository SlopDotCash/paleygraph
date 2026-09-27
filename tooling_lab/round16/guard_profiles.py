#!/usr/bin/env python3
"""Reserve singleton tests in each conditional branch, then certify all queries.

This is a profile constraint, not an unconditional guard rejection rule. All
decisions still use the exact transcript bounds. A true candidate may disagree
with a guard and must not be discarded merely for that disagreement.
"""
from copy import deepcopy
from hashlib import sha256
import json
from pathlib import Path
import random
from conditional_cover import inputs,projective_groups
from conditional_verifier import validate_conditional
from capacity_search import audit_queries,certify_blocks,validate
from projective_profile import optimize,pack
from conditional_pruning import run

HERE=Path(__file__).resolve().parent


def guarded_cover(data,min_unqueried=1,seed=160907,max_attempts=64):
    if type(min_unqueried) is not int or min_unqueried<1:raise ValueError('positive singleton count required')
    columns=inputs(data);groups,zeros=projective_groups(columns,data['p'])
    shape=optimize(list(map(len,groups)),len(zeros),data['s'],data['dimension'])
    rows=[r for r in shape['candidate_profiles'] if r['admissible'] and r['unqueried_coordinates']>=min_unqueried]
    if not rows:return {'status':'no_guarded_profile','profile':shape}
    shape['best']=min(rows,key=lambda r:(r['queries'],r['queried_blocks']))
    if shape['best']['queries']>250000:return {'status':'query_budget_exceeded','profile':shape}
    rng=random.Random(seed);attempts=[]
    for attempt in range(max_attempts):
        if attempt:
            rng.shuffle(groups)
            for group in groups:rng.shuffle(group)
        blocks=pack(groups,zeros,shape);audit=audit_queries(columns,blocks,data['p'],data['dimension'])
        attempts.append({'attempt':attempt,**audit})
        if audit['dependent_queries']:continue
        cert=certify_blocks(data,blocks);checked=validate(cert)
        return {'status':'found','min_unqueried':min_unqueried,'seed':seed,'profile':shape,'attempts':attempts,'certificate':cert,'review':checked}
    return {'status':'not_found','min_unqueried':min_unqueried,'seed':seed,'profile':shape,'attempts':attempts}


def replace_covers(c):
    out=deepcopy(c);searches=[]
    for index,part in enumerate(out['parts']):
        result=guarded_cover(part['certificate']['input'],seed=160907+index)
        if result['status']!='found':return {'status':'guard_search_failed','failed_branch':index,'search':result}
        part['certificate']=result['certificate'];searches.append(result)
    out['searches']=searches;out['query_count']=sum(p['certificate']['query_count'] for p in out['parts'])
    out['profile_policy']={'min_unqueried_per_branch':1,'scope':'The histogram cutoff is retained; branch profiles reserve at least one singleton. The original plan costs remain lower-bound estimates, not the costs of these modified profiles.'}
    validate_conditional(out)
    return out


def main():
    rows=[];bindings={};artifacts={}
    for name in ['boundary','late_near_miss','true_guard_mismatch']:
        source=HERE/f'{name}.conditional.json';c=replace_covers(json.loads(source.read_text()));assert c['status']=='found';bindings[source.name]=sha256(source.read_bytes()).hexdigest()
        certpath=HERE/f'{name}.guarded.json';certpath.write_text(json.dumps(c,separators=(',',':'))+'\n');artifacts[certpath.name]=sha256(certpath.read_bytes()).hexdigest()
        print(json.dumps({'phase':'guard_cover_certified','case':name,'queries':c['query_count'],'singletons':[s['profile']['best']['unqueried_coordinates'] for s in c['searches']]}),flush=True)
        result=run(c);outpath=HERE/f'{name}.guarded.decoded.json';outpath.write_text(json.dumps(result,separators=(',',':'))+'\n');artifacts[outpath.name]=sha256(outpath.read_bytes()).hexdigest()
        row={'case':name,'certificate':certpath.name,'decoded':outpath.name,'metrics':result['metrics'],'outputs':len(result['outputs'])};rows.append(row);print(json.dumps(row),flush=True)
    out={'status':'produced','scope':'Guard-reserving conditional profiles certified and decoded. Separate complete branch replay is required before accepting numerical improvements.',
         'cases':rows,'source_sha256':{f:sha256((HERE/f).read_bytes()).hexdigest() for f in ['guard_profiles.py','conditional_pruning.py','conditional_verifier.py','conditional_cover.py','../round15/capacity_search.py','../round15/projective_profile.py','../round14/arc_pruning.py','../round14/arc_filter_fast.py','../round14/arc_verifier.py','../round14/arc_certificate.py','../round14/dimension_preflight.py','../round14/arc_profile.py']},'input_sha256':bindings,'artifact_sha256':artifacts}
    (HERE/'guard_results.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
