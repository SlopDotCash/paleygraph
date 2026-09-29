#!/usr/bin/env python3
"""Recover the entire distance-two mean from the frozen row histogram alone."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import comb
from pathlib import Path

HERE=Path(__file__).resolve().parent


def allocations(populations):
    for z in range(3):
        for p in range(3-z):
            m=2-z-p;counts=(z,p,m)
            if all(a<=b for a,b in zip(counts,populations)):
                yield counts,comb(populations[0],z)*comb(populations[1],p)*comb(populations[2],m)


def elementary(p,m,d):
    return sum((-1)**j*comb(m,j)*comb(p,d-j) for j in range(max(0,d-p),min(d,m)+1))


def main():
    oldpath=HERE.parent/'round7/local_edits/results.json';old=json.loads(oldpath.read_text())['scale_cases']
    new=json.loads((HERE/'scale_results.json').read_text())['cases'];rows=[]
    for prev,r in zip(old,new):
        q,n,d=prev['q'],prev['n'],prev['degree'];assert r['q']==q and r['selected']==prev['selected']
        total=0
        for z,p,count in prev['signed_row_count_histogram']:
            negative=n-z-p;outside=(1-z,(q-1)//2-p,(q-1)//2-negative)
            for (az,ap,am),w1 in allocations((z,p,negative)):
                for (bz,bp,bm),w2 in allocations(outside):
                    total+=count*w1*w2*elementary(p-ap+bp,negative-am+bm,d)
        assert 2*total==sum(x['target_sum'] for x in r['pairs'])
        count=comb(n,2)*comb(q-n,2);mean=F(total,count)
        rows.append({'q':q,'n':n,'family':r['family'],'distance_two_final_sets':count,'mean':[mean.numerator,mean.denominator]})
    out={'status':'passed','scope':'Exact full-shell mean using only the previously verified row-count histogram, compared with all new pair-mean numerators.',
         'cases':rows,'source_sha256':{'shell_mean_review.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{'../round7/local_edits/results.json':sha256(oldpath.read_bytes()).hexdigest(),
                         'scale_results.json':sha256((HERE/'scale_results.json').read_bytes()).hexdigest()}}
    (HERE/'shell_mean_verification.json').write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'status':'passed','cases':len(rows),'largest_distance_two_shell':rows[-1]['distance_two_final_sets']}),flush=True)


if __name__=='__main__':main()
