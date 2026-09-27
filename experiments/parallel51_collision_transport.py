#!/usr/bin/env python3
"""Exact cross-level primitive-fiber transport and a quartic counterexample.

These checks do not supply a uniform collision or Paley estimate.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT = Path(__file__).resolve().parents[1]


def pairs(A, B, p):
    return Counter((a+b) % p for a in A for b in B)


def energy(r):
    return sum(v*v for v in r.values())


def inspect(p, N, g, expected=None, full=False):
    assert (p-1) % N == 0
    assert pow(g,N,p) == 1 and pow(g,N//2,p) == p-1
    levels, transported, intrinsic, Tsum, loops = [], {}, 0, Fraction(0), 0
    all_e = Counter()
    s=4
    while s<=N:
        h=pow(g,N//s,p)
        H=[pow(h,j,p) for j in range(s)]
        K,L=H[::2],H[1::2]
        roots=[pow(1+pow(h,j,p),s,p) for j in range(1,s//2,2)]
        c=Counter(roots)
        b=Counter()
        for x,v in c.items():
            b[pow(x,N//s,p)]+=v
        assert sum(c.values()) == sum(b.values()) == s//4
        # Each labeled root has the same endpoint label when computed
        # directly from the endpoint generator.
        direct=[pow(1+pow(g,(N//s)*j,p),N,p) for j in range(1,s//2,2)]
        assert Counter(direct)==b
        a_s=sum(v*(v-1) for v in c.values())
        Y=sum(v*(v-2) for v in c.values() if v>=3)
        alias=sum(v*(v-1) for v in b.values())-a_s
        assert alias>=0 and alias%2==0 and a_s>=Y
        mixed,child,parent=pairs(K,L,p),pairs(K,K,p),pairs(H,H,p)
        loops+=s*s+2*(s//2)**2
        B,Ec,E=energy(mixed),energy(child),energy(parent)
        T=sum(v*mixed.get(x,0) for x,v in child.items())
        assert Counter(mixed.values())==Counter({v:s*sum(w==v for w in c.values()) for v in set(c.values())})
        assert E==2*Ec+6*B+8*T
        assert B==(s//2)**2+s*a_s
        intrinsic+=a_s
        Tsum+=Fraction(T,s)
        D=Fraction(E-(3*s*s-3*s),s)
        assert D==6*intrinsic+8*Tsum
        row=dict(s=s,Y=Y,pair_collision_mass=a_s,alias_ordered_mass=alias,
                 multiplicity_histogram=dict(sorted(Counter(c.values()).items())),
                 transported_histogram=dict(sorted(Counter(b.values()).items())),
                 E=E,T=T,D=str(D))
        if full:
            row['root_labels_and_multiplicities']=sorted(c.items())
            row['endpoint_labels_and_multiplicities']=sorted(b.items())
            if Y:
                x=min(x for x,v in mixed.items() if v>=3)
                Lset=set(L)
                representations=[(a,(x-a)%p) for a in K if (x-a)%p in Lset]
                assert len(representations)==mixed[x]>=3
                row['witness']=dict(x=x,representations=representations)
        levels.append(row);transported[s]=b;all_e.update(b)
        s*=2
    assert sum(all_e.values())==N//2-1
    H=[pow(g,j,p) for j in range(N)]
    kernel=Counter(pow(1+h,N,p) for h in H if h!=p-1)
    target=Counter({x:2*v for x,v in all_e.items()})
    target[pow(2,N,p)]+=1
    assert kernel==target
    literal_E=energy(pairs(H,H,p));loops+=N*N
    assert literal_E==N*N+N*sum(v*v for v in kernel.values())
    alias=sum(x['alias_ordered_mass'] for x in levels)
    cross=sum(sum(v*transported[t].get(x,0) for x,v in transported[s].items())
              for s in transported for t in transported if s<t)
    anchor=all_e[pow(2,N,p)]
    root_energy=sum(v*(v-1) for v in all_e.values())+anchor
    assert root_energy==intrinsic+alias+2*cross+anchor
    assert literal_E==3*N*N-3*N+4*N*root_energy
    assert alias+2*cross+anchor==Fraction(intrinsic,2)+2*Tsum
    Ytotal=sum(x['Y'] for x in levels)
    assert literal_E>=3*N*N-3*N+6*N*Ytotal
    if expected:
        assert literal_E==expected['parent_energy']
        assert levels[-1]['Y']==expected['triple_mass_Y']
    return dict(p=p,N=N,generator=g,quartic_window=N**4<=4*p<=4*N**4,
                levels=levels,total_triple_mass=Ytotal,intrinsic_pair_mass=intrinsic,
                alias_ordered_mass=alias,cross_level_product=cross,anchor_mass=anchor,
                scaled_unbalanced_sum=str(Tsum),E=literal_E,
                compatibility_left=str(alias+2*cross+anchor),
                compatibility_right=str(Fraction(intrinsic,2)+2*Tsum),
                literal_pair_checks=loops,all_passed=True)


def main():
    source=ROOT/'results/parallel48_collision_eliminants_2026_09_06.json'
    old=json.loads(source.read_text())
    towers=[inspect(c['p'],c['n'],c['generator'],c)
            for order in old['orders'] for c in order['fields']]
    p,N=2144280833,256
    old_child=next(c for order in old['orders'] if order['n']==128
                   for c in order['fields'] if c['p']==p)
    limit=isqrt(p)
    assert p>2 and p%2 and all(p%d for d in range(3,limit+1,2))
    for a in range(2,p):
        g=pow(a,(p-1)//N,p)
        if pow(g,N//2,p)==p-1:break
    witness=inspect(p,N,g,full=True)
    assert witness['quartic_window']
    assert witness['levels'][-2]['Y']==3 and witness['levels'][-1]['Y']==0
    assert witness['levels'][-2]['D']==witness['levels'][-1]['D']=='48'
    assert witness['E']==208128
    assert witness['total_triple_mass']==3
    assert witness['intrinsic_pair_mass']==8
    assert witness['levels'][-2]['Y']==old_child['triple_mass_Y']
    assert witness['levels'][-2]['E']==old_child['parent_energy']
    result=dict(scope='Cross-level compatibility and a quartic counterexample to using only the top primitive triple mass; no uniform estimate.',
                towers=towers,witness=witness,
                witness_prime_certificate=dict(method='trial division by all odd integers through floor(sqrt(p))',p=p,limit=limit),
                witness_child_in_prior_certificate=True,
                tower_cases=len(towers)+1,
                level_checks=sum(len(t['levels']) for t in towers)+len(witness['levels']),
                literal_pair_checks=sum(t['literal_pair_checks'] for t in towers)+witness['literal_pair_checks'],
                prior_certificate_sha256=sha256(source.read_bytes()).hexdigest(),all_passed=True)
    out=ROOT/'results/parallel51_collision_transport_2026_09_06.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({k:result[k] for k in ['tower_cases','level_checks','literal_pair_checks','all_passed']}))
    print(json.dumps({k:witness[k] for k in ['p','N','E','total_triple_mass','alias_ordered_mass','cross_level_product','anchor_mass']}))


if __name__=='__main__':
    main()
