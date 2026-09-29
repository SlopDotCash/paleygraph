#!/usr/bin/env python3
from fractions import Fraction as F
from hashlib import sha256
from itertools import product
from pathlib import Path
import json
from orbit_intersection import count_pairs

HERE=Path(__file__).resolve().parent


def main():
    small=[]
    for p,g in [(17,9),(41,27),(97,64)]:
        records=[];m=p//2
        for delta in range(1,p):
            Fdelta=[(delta*pow(g,i,p)+m)%p-m for i in range(4)]
            for mask in [[0],[0,2],list(range(4))]:
                for H in product(*([Fdelta[i],Fdelta[i]-p*(1 if Fdelta[i]>0 else -1)] for i in mask)):
                    fixed=dict(zip(mask,H));out=count_pairs(p,g,delta,fixed,limit=p)
                    records.append({'delta':delta,'fixed':list(map(list,fixed.items())),'output':out})
        name=f'orbit_small_p{p}';path=HERE/(name+'.json');path.write_text(json.dumps({'p':p,'g':g,'N':4,'records':records},separators=(',', ':'))+'\n')
        small.append({'artifact':path.name,'queries':len(records)});print(json.dumps(small[-1]),flush=True)
    source=HERE/'consecutive40.json';c=json.loads(source.read_text())['certificate'];p,N,k=c['p'],c['N'],c['k'];adj=c['adj'];U=c['erased']
    coverpath=HERE/'cover40.json';cover=json.loads(coverpath.read_text())['cover']
    cut_index=cover['continuous_witness_cuts'][0];target=cover['cuts'][cut_index]['target']
    physical=[];free=[]
    for j,row in enumerate(c['basis']):
        numerator=[sum((adj[i-u] if i>=u else -adj[N+i-u])*x for u,x in zip(U,row)) for i in range(N)]
        assert all(v%k==0 for v in numerator);image=[v//k for v in numerator];physical.append(image)
        if j in c['invisible_directions']:
            support=[i for i,v in enumerate(image) if v]
            assert len(support)==1 and abs(image[support[0]])==p;free.append(support[0])
    assert len(free)==len(set(free));H=[sum(t*physical[j][i] for j,t in zip(c['visible_directions'],target)) for i in range(N)]
    delta=H[0]%p;fixed={i:H[i] for i in range(N) if i not in free}
    assert all(H[i]%p==delta*pow(c['g'],i,p)%p for i in range(N))
    out=count_pairs(p,c['g'],delta,fixed)
    large={'p':p,'g':c['g'],'N':N,'source_cover_cut':cut_index,'target':target,'delta':delta,
           'free_coordinates':sorted(free),'fixed':list(map(list,fixed.items())),'output':out}
    (HERE/'orbit_large.json').write_text(json.dumps(large,indent=2)+'\n');print(json.dumps({'large':out}),flush=True)
    report={'status':'produced','small_cases':small,'large_artifact':'orbit_large.json',
            'input_sha256':{path.name:sha256(path.read_bytes()).hexdigest() for path in [source,coverpath]},
            'source_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ['orbit_experiments.py','orbit_intersection.py']}}
    (HERE/'orbit_summary.json').write_text(json.dumps(report,indent=2)+'\n')


if __name__=='__main__':main()
