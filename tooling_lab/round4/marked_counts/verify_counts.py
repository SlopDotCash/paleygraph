#!/usr/bin/env python3
"""Independent standard-library oracle: Euler symbols and direct ordered pairs.
No NTT, Fourier transform, NumPy, or imports from the producer.
"""
import collections,hashlib,json,math,platform,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def prime(p):
    return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))

def char_euler(p):
    assert prime(p) and p%4==1
    ret=[]
    for x in range(p):
        a=pow(x,(p-1)//2,p);assert a in [0,1,p-1];ret.append(-1 if a==p-1 else a)
    return ret

def relations(data):
    ret={}
    for r in data['relations']:
        key=(tuple(r['left_pattern']),tuple(r['right_pattern']),r['sign'])
        assert key not in ret and r['sign'] in [-1,0,1] and isinstance(r['count'],int) and r['count']>0
        ret[key]=r['count']
    return ret

def verify(path,ch):
    start=time.monotonic();data=json.loads(path.read_text());p=data['p'];M=data['marks']
    patterns=[tuple(ch[(x-m)%p] for m in M) for x in range(p)]
    sizes=dict(collections.Counter(patterns));observed_sizes={tuple(c['pattern']):c['size'] for c in data['cells']}
    assert sizes==observed_sizes and list(observed_sizes)==sorted(sizes)
    expected=relations(data);assert sum(expected.values())==p*p
    for U,size in sizes.items():
        for sign in [-1,0,1]:
            assert sum(v for (u,w,s),v in expected.items() if u==U and s==sign)==size*(1 if sign==0 else (p-1)//2)
    matrix=None;direct_pairs=0
    if p<=1297:
        # Materialize actual direct sign matrices at the acceptance primes.
        if p<=101:matrix=[[ch[(x-y)%p] for y in range(p)] for x in range(p)]
        actual=collections.Counter()
        for x in range(p):
            for y in range(p):actual[(patterns[x],patterns[y],matrix[x][y] if matrix else ch[(x-y)%p])]+=1
        assert dict(actual)==expected;direct_pairs=p*p
    order=list(observed_sizes);idx={U:i for i,U in enumerate(order)};cell=[idx[t] for t in patterns]
    sample_terms=0
    for row in data['sampled_rows']:
        x=row['x'];assert 0<=x<p;actual=[0]*len(order)
        for y in range(p):actual[cell[y]]+=ch[(x-y)%p]
        assert actual==row['cell_character_sums'];sample_terms+=p
    prov=data['provenance']
    for key,source in [('backend_source_sha256','marked_counts.cpp'),('backend_executable_sha256','marked_counts'),('runner_source_sha256','run_counts.py')]:assert prov[key]==sha(HERE/source)
    assert prov['input_sha256']==sha(HERE/prov['input_file'])
    return {'file':path.name,'sha256':sha(path),'p':p,'marks':M,'cells':len(sizes),'relations':len(expected),
            'all_character_values_verified_by_Euler_criterion':p,'independent_direct_ordered_pairs':direct_pairs,
            'independent_sampled_rows':len(data['sampled_rows']),'independent_sampled_row_terms':sample_terms,
            'full_relation_counts_independently_recomputed':bool(direct_pairs),'seconds':time.monotonic()-start}

def main():
    start=time.monotonic();cache={};rows=[]
    paths=sorted(HERE.glob('counts_*.json'),key=lambda f:(json.loads(f.read_text())['p'],f.name))
    for path in paths:
        data=json.loads(path.read_text());p=data['p']
        if p not in cache:cache[p]=char_euler(p)
        row=verify(path,cache[p]);rows.append(row);print(path.name,'verified',flush=True)
    bymark={tuple(json.loads(p.read_text())['marks']):json.loads(p.read_text()) for p in paths if json.loads(p.read_text())['p']==101}
    # Exact affine covariance controls on coloured counts, including global sign.
    def transformed(data,epsilon):
        sizes={tuple(epsilon*t for t in c['pattern']):c['size'] for c in data['cells']}
        rel={(tuple(epsilon*t for t in r['left_pattern']),tuple(epsilon*t for t in r['right_pattern']),epsilon*r['sign']):r['count'] for r in data['relations']}
        return sizes,rel
    covariance=[]
    for left,right,epsilon,name in [((0,),(1,),1,'one-mark translation'),((0,1,2),(1,2,3),1,'three-mark translation'),
          ((0,1,2),(0,4,8),1,'square dilation'),((0,1,2),(0,2,4),-1,'nonsquare dilation/sign reversal'),
          ((0,1),(0,2),-1,'two-mark nonsquare dilation/sign reversal')]:
        assert transformed(bymark[left],epsilon)==transformed(bymark[right],1)
        covariance.append({'control':name,'left_marks':left,'right_marks':right,'sign_multiplier':epsilon,'passed':True})
    guard_inputs=[[],['9'],['15'],['5242889'],['4294967297'],['13','0','0'],['13','-1'],['13','13'],['13','0','1','2','3','4','5','6']]
    guards=[]
    for args in guard_inputs:
        result=subprocess.run([str(HERE/'marked_counts'),*args],capture_output=True,text=True)
        assert result.returncode!=0 and not result.stdout
        guards.append({'arguments':args,'exit_code':result.returncode,'error':result.stderr.strip()})
    out={'status':'all exact checks passed; large counts have sampled independent row validation, not a second full convolution',
         'cases':rows,'affine_covariance_controls':covariance,'rejection_guards':guards,
         'source_sha256':sha(__file__),'elapsed_seconds':time.monotonic()-start,
         'runtime':{'python':sys.executable,'version':platform.python_version()}}
    (HERE/'validation.json').write_text(json.dumps(out,indent=2)+'\n')

if __name__=='__main__':main()
