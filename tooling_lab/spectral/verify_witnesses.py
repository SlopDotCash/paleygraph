#!/usr/bin/env python3
"""Independent standard-library exact replay of exported spectral witnesses.

Does not import discovery code, NumPy, SciPy, or eigenvalue computations.
Recomputes all weighted differences and quadratic forms using Python ints.
"""
import json,math,hashlib
from pathlib import Path
ROOT=Path(__file__).resolve().parent

def isprime(p):return p>=2 and all(p%d for d in range(2,math.isqrt(p)+1))

def replay(w):
    p,d,q=w['characteristic'],w['degree'],w['q'];assert isprime(p) and q==p**d and d in (1,2)
    def cp(a):return 0 if not a%p else 1 if pow(a%p,(p-1)//2,p)==1 else -1
    nu=next(a for a in range(2,p) if cp(a)==-1)
    def sub(x,y):
        return (x-y)%p if d==1 else (x%p-y%p)%p+p*((x//p-y//p)%p)
    def ch(x):return cp(x) if d==1 else cp((x%p)**2-nu*(x//p)**2)
    C=[x for x in range(q) if ch(x)==1 and ch(sub(x,1))==1];assert C==w['C']
    z=w['integer_vector'];assert len(C)==len(z) and all(isinstance(v,int) for v in z)
    n=sum(v*v for v in z);s=sum(z);a=0;corr=[0]*q
    for x,v in zip(C,z):
        for y,u in zip(C,z):
            t=sub(y,x);prod=u*v;corr[t]+=prod;a+=prod*ch(t)
    assert (n,s,a)==(w['norm_squared'],w['sum_entries'],w['signed_quadratic_form'])
    assert corr[w['transport_shift']]==w['translation_inner_product']
    T=[t for t,c in enumerate(corr) if 5*abs(c)>=3*n]
    margin=4*a*a-3*q*n*n if s==0 else None
    return {'q':q,'weighted_difference_pairs':len(C)**2,'norm_squared':n,
            'large_overlap_translations':T,
            'mean_zero_edge_squared_margin':margin,
            'certified_two_anchor_edge_outlier':bool(margin is not None and margin>0),
            'explanation':'When sum entries is zero, 4 A^2 - 3 q n^2 > 0 exactly certifies |v^T H v|/n > sqrt(3)/2; null margin means this test does not apply.'}

def main():
    main=json.loads((ROOT/'results.json').read_text());refine=json.loads((ROOT/'refinement.json').read_text());expected={r['q']:r['signature']['large_overlap_translations'] for r in refine['cases']}
    rows=[]
    for r in main['cases']:
        w=json.loads((ROOT/r['certificate_file']).read_text());out=replay(w)
        assert out['large_overlap_translations']==expected[w['q']]
        rows.append(out);print(w['q'],'verified',out['certified_two_anchor_edge_outlier'],flush=True)
    special=next(r for r in rows if r['q']==401)
    assert special['certified_two_anchor_edge_outlier'] and special['large_overlap_translations']==[0]
    payload={'status':'all exported witnesses and translation signatures verified using Python integer arithmetic',
             'count':len(rows),'total_weighted_difference_pairs':sum(r['weighted_difference_pairs'] for r in rows),
             'cases':rows,'source_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    (ROOT/'validation.json').write_text(json.dumps(payload,indent=2)+'\n')
if __name__=='__main__':main()
