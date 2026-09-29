#!/usr/bin/env python3
"""Literal matrix histogram checks and full prime parameter-orbit census."""
from array import array
from collections import Counter
from hashlib import sha256
from itertools import combinations,permutations,product
import json
from pathlib import Path
import sys
import time
sys.dont_write_bytecode=True
from folded_moment import compiler,theta_family

HERE=Path(__file__).resolve().parent


def chi(p):
    return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]


def prime(p):
    return p>=2 and all(p%d for d in range(2,__import__('math').isqrt(p)+1))


def prime_orbits(p):
    start=time.perf_counter();assert prime(p) and p%4==1
    path=HERE.parent/'trace_backend'/f'tau_p{p}.bin'
    data=array('i');data.frombytes(path.read_bytes())
    if sys.byteorder!='little':data.byteswap()
    assert data[0]==p;taus=data[1:];assert len(taus)==p
    sign=chi(p);inverse=array('i',[0])*p;inverse[1]=1
    for t in range(2,p):inverse[t]=-(p//t)*inverse[p%t]%p
    seen=bytearray(p);census=Counter();folded=Counter();special=[]
    for t in range(2,p):
        theta=sign[t]*sign[t-1]*taus[t]
        mono=sign[t]==sign[t-1]==1
        assert theta_family(p,theta)==mono
        r=(1-t)%p;i=inverse[t]
        assert taus[r]==taus[t] and taus[i]==sign[t]*taus[t]
        assert sign[r]*sign[(r-1)%p]*taus[r]==theta
        assert sign[i]*sign[(i-1)%p]*taus[i]==theta
        assert (sign[r]==sign[(r-1)%p]==1)==mono
        assert (sign[i]==sign[(i-1)%p]==1)==mono
        folded[mono,theta]+=1
        if seen[t]:continue
        orbit=sorted({t,r,i,inverse[r],(t-1)*i%p,t*inverse[(t-1)%p]%p})
        assert len(orbit) in (2,3,6)
        for x in orbit:
            assert not seen[x];seen[x]=1
            assert sign[x]*sign[(x-1)%p]*taus[x]==theta
            assert (sign[x]==sign[(x-1)%p]==1)==mono
        census[len(orbit)]+=1
        if len(orbit)<6:special.append({'parameters':orbit,'orbit_size':len(orbit),
                                      'stabilizer_size':6//len(orbit),'monochromatic':mono,'theta':theta,
                                      'type':'equianharmonic' if len(orbit)==2 else 'harmonic'})
    assert sum(k*v for k,v in census.items())==p-2 and sum(seen)==p-2
    assert census[3]==1 and census[2]==(1 if p%3==1 else 0)
    assert next(x['parameters'] for x in special if x['orbit_size']==3)==sorted({p-1,2,inverse[2]})
    for item in special:
        if item['orbit_size']==2:assert all((t*t-t+1)%p==0 for t in item['parameters'])
    saved=json.loads((HERE/f'inventory_Paley{p}.json').read_text())
    assert Counter({(row['monochromatic'],row['theta']):row['count'] for row in saved['records']})==folded
    return {'p':p,'all_nonsingular_parameters_checked':p-2,'generator_checks':2*(p-2),
            'orbit_census':dict(sorted(census.items())),'special_orbits':special,
            'tau_sha256':sha256(path.read_bytes()).hexdigest(),'elapsed_seconds':time.perf_counter()-start}


def matrix_cases():
    for p in (13,17,29):
        sign=chi(p);yield f'Paley{p}',[[sign[(x-y)%p] for y in range(p)] for x in range(p)]
    def mul(a,b):return ((a%7)*(b%7)-(a//7)*(b//7))%7+7*(((a%7)*(b//7)+(a//7)*(b%7))%7)
    def sub(a,b):return (a%7-b%7)%7+7*((a//7-b//7)%7)
    powers=[];x=1
    for _ in range(48):powers.append(x);x=mul(x,9)
    assert x==1 and len(set(powers))==48
    logarithm={x:k for k,x in enumerate(powers)}
    for name,classes in (('Paley49',(0,2)),('Peisert49',(0,1))):
        sign=[0]+[1 if logarithm[x]%4 in classes else -1 for x in range(1,49)]
        yield name,[[sign[sub(x,y)] for y in range(49)] for x in range(49)]


def literal_matrix_symmetry(name,S):
    q=len(S);checked=0;folded=Counter();all_transforms=0
    for rows in combinations(range(q),3):
        a,b,c=rows;edges=(S[a][b],S[a][c],S[b][c]);mono=len(set(edges))==1
        tau=sum(S[a][x]*S[b][x]*S[c][x] for x in range(q))
        theta=__import__('math').prod(edges)*tau
        assert theta_family(q,theta)==mono
        target=(1,1,1) if mono else (1,1,-1)
        target_tau=theta if mono else -theta
        expected=compiler.distinct_histogram(q,target,target_tau)
        # Check every transform taking this triangle to the canonical one.
        matching=0
        for perm in permutations(rows):
            u,v,w=perm
            for sign in (-1,1):
                if tuple(sign*z for z in (S[u][v],S[u][w],S[v][w]))!=target:continue
                actual=Counter(tuple(sign*S[row][x] for row in perm) for x in range(q))
                assert actual==expected
                matching+=1
        assert matching==(6 if mono else 2)
        checked+=1;all_transforms+=matching;folded[mono,theta]+=6
    normalized=json.loads((HERE/f'inventory_{name}.json').read_text())
    assert folded==Counter({(r['monochromatic'],r['theta']):q*(q-1)*r['count'] for r in normalized['records']})
    return {'graph':name,'q':q,'unordered_distinct_triples_checked':checked,
            'literal_canonical_histograms_checked':all_transforms,
            'all_ordered_triples_match_normalized_weights':True}


def main():
    out={'prime_orbits':[],'literal_matrix_checks':[],'historical_novelty_claim':False}
    for p in (13,17,29,101,1297,65537,1000033):
        row=prime_orbits(p);out['prime_orbits'].append(row);print('prime',p,row['orbit_census'],flush=True)
    for name,S in matrix_cases():
        row=literal_matrix_symmetry(name,S);out['literal_matrix_checks'].append(row);print('matrix',name,row,flush=True)
    out['source_sha256']={p.name:sha256(p.read_bytes()).hexdigest() for p in (Path(__file__),HERE/'folded_moment.py')}
    (HERE/'symmetry_validation.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
