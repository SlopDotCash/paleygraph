#!/usr/bin/env python3
"""Validate a single signed convolution against full relations and direct sums."""
from collections import Counter
import hashlib,json,math,platform,resource,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent
ROUND=HERE.parent
sys.dont_write_bytecode=True
sys.path.insert(0,str(ROUND/'marked_moments'))
from conditional_moments import conditional_moments
from contraction_reduction import reduce_from_cells
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def pair_values(u):
    if 0 in u:return 0,0
    a,b,c=u;return a*b+a*c+b*c,a*b*c
def relation_Q(record):
    return sum(pair_values(r['left_pattern'])[0]*pair_values(r['right_pattern'])[1]*r['sign']*r['count'] for r in record['relations'])
def character(p):
    assert p%4==1 and all(p%d for d in range(2,math.isqrt(p)+1))
    return [0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
def main():
    inputs=[]
    for path in sorted((ROUND/'marked_counts').glob('counts_*.json')):
        record=json.loads(path.read_text())
        if len(record['marks'])==3:inputs.append((path,record))
    path=HERE/'witness.json';w=json.loads(path.read_text())
    for side in w['pair']:inputs.append((path,side['record']))
    rows=[];cache={}
    for input_path,record in inputs:
        p=record.get('p',record.get('q'));M=record['marks'];start=time.monotonic()
        raw=subprocess.check_output([str(HERE/'single_contraction'),str(p),*map(str,M)],text=True);out=json.loads(raw)
        assert out['cells']==record['cells'];expected=relation_Q(record);assert out['Q']==expected
        if p not in cache:cache[p]=character(p)
        chi=cache[p];h2=[];h3=[]
        for x in range(p):
            a,b=pair_values([chi[(x-m)%p] for m in M]);h2.append(a);h3.append(b)
        samples={r['x']:r['value'] for r in out['sampled_convolution_rows']};totalQ=0;directrows=0
        for x in (range(p) if p<=1297 else samples):
            value=sum(chi[(x-y)%p]*h3[y] for y in range(p))
            if x in samples:assert value==samples[x]
            totalQ+=h2[x]*value;directrows+=1
        if p<=1297:assert totalQ==out['Q']
        reduced=[]
        if p>=17:
            for n in [6,7,8]+([31] if p>=31 else []):
                r=reduce_from_cells(p,M,{tuple(c['pattern']):c['size'] for c in out['cells']},out['Q'],n)
                full=conditional_moments(record,n)
                assert all(r[key]==full[key] for key in ['mean','second_moment','variance'])
                reduced.append(r)
        out['reduced_moments']=reduced
        out['validation']={'Q_from_full_exact_relations':expected,'relation_Q_matched':True,
            'direct_rows':directrows,'direct_integer_terms':directrows*p,'full_independent_direct_Q':p<=1297,
            'all_character_values_checked_by_Euler_criterion':p}
        out['provenance']={'input_file':str(input_path.relative_to(ROUND)),'input_sha256':sha(input_path),
            'cpp_sha256':sha(HERE/'single_contraction.cpp'),'header_sha256':sha(HERE/'ntt_exact.hpp'),
            'binary_sha256':sha(HERE/'single_contraction'),'verifier_sha256':sha(__file__),
            'reduction_source_sha256':sha(ROUND/'marked_moments/contraction_reduction.py'),
            'conditional_compiler_source_sha256':sha(ROUND/'marked_moments/conditional_moments.py'),
            'wall_including_independent_validation_seconds':time.monotonic()-start,
            'runner_lifetime_max_child_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*(1 if sys.platform=='darwin' else 1024)}
        target=HERE/('contraction_p'+str(p)+'_marks_'+'_'.join(map(str,M))+'.json')
        target.write_text(json.dumps(out,indent=2)+'\n');rows.append({'file':target.name,'sha256':sha(target),'p':p,'marks':M,'Q':out['Q'],**out['metadata'],**out['validation']})
        print(json.dumps({'p':p,'marks':M,'Q':out['Q'],'full_relation_match':True,'backend_seconds':out['metadata']['backend_seconds'],'direct_rows':directrows}),flush=True)
    failures=[]
    for args in [[],['49','0','1','2'],['29','0','0','1'],['29','0','1','29'],['5242889','0','1','2']]:
        p=subprocess.run([str(HERE/'single_contraction'),*args],capture_output=True,text=True);assert p.returncode!=0 and not p.stdout
        failures.append({'args':args,'error':p.stderr.strip()})
    (HERE/'contraction_validation.json').write_text(json.dumps({'status':'single-convolution Q matches all exact full relation counts and independent direct checks',
        'cases':rows,'guard_rejections':failures,'runtime':{'python':sys.executable,'version':platform.python_version()},'source_sha256':sha(__file__)},indent=2)+'\n')
if __name__=='__main__':main()
