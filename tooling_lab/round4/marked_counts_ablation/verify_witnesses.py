#!/usr/bin/env python3
"""Run independent all-completion pointwise enumeration on actual witnesses."""
import hashlib,json,math,subprocess,time
from fractions import Fraction
from pathlib import Path
HERE=Path(__file__).resolve().parent
def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def main():
    out=[]
    for prefix in ['', 'twins_','cross_twins_']:
        path=HERE/(prefix+'witness.json')
        if not path.exists():continue
        w=json.loads(path.read_text());mp=HERE/(prefix+'witness_matrices.json');matrices=json.loads(mp.read_text())['graphs']
        payload='2\n'
        for row in matrices:
            S=row['S'];payload+=str(len(S))+'\n'+'\n'.join(' '.join(map(str,r)) for r in S)+'\n'
        queries=[(g,n) for g in range(2) for n in [6,7,8]];payload+=str(len(queries))+'\n'
        for g,n in queries:
            M=w['pair'][g]['marks'];payload+=f'{g} {n} {len(M)} '+' '.join(map(str,M))+'\n'
        start=time.monotonic();raw=subprocess.check_output([str(HERE/'direct_completions')],input=payload,text=True)
        rows=[]
        for line in raw.splitlines():
            r=json.loads(line);g=r['graph_index'];n=r['n'];M=w['pair'][g]['marks'];ref=w['pair'][g]['moments'][str(n)]
            assert r['completions']==math.comb(w['q']-len(M),n-len(M))
            mean=Fraction(r['sum_T6'],r['completions']);second=Fraction(r['sum_T6_squared'],r['completions'])
            assert [mean.numerator,mean.denominator]==ref['mean']
            assert [second.numerator,second.denominator]==ref['second_moment']
            r.update({'marks':M,'mean':[mean.numerator,mean.denominator],'second_moment':[second.numerator,second.denominator]});rows.append(r)
        left,right=w['pair'];assert left['record']['cells']==right['record']['cells']
        comparisons=[]
        for n in [6,7,8]:
            a=Fraction(*left['moments'][str(n)]['second_moment']);b=Fraction(*right['moments'][str(n)]['second_moment']);gap=b-a
            comparisons.append({'n':n,'right_minus_left_second_moment':[gap.numerator,gap.denominator],'cell_profile_identical':True})
        out.append({'witness':path.name,'witness_sha256':sha(path),'matrices_sha256':sha(mp),'all_direct_completions':rows,
                    'comparisons':comparisons,'elapsed_seconds':time.monotonic()-start})
        print(prefix or 'prime','all completion checks passed',flush=True)
    result={'status':'actual witnesses independently confirmed by exhaustive pointwise completion enumeration',
        'cases':out,'source_sha256':{'python':sha(__file__),'cpp':sha(HERE/'direct_completions.cpp'),'executable':sha(HERE/'direct_completions')}}
    (HERE/'direct_validation.json').write_text(json.dumps(result,indent=2)+'\n')
if __name__=='__main__':main()
