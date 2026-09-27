#!/usr/bin/env python3
"""Exact unit transport on the centered cyclotomic coset section.

Multiplication by epsilon=1+X+X^-1 preserves the integral field norm. Returning
to the centered coefficient section introduces a ternary carry vector and can
change the norm defect. The tool retains that carry, rather than treating
coefficient centering as a norm-minimizing operation.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
import json
from pathlib import Path

HERE=Path(__file__).resolve().parent


def multiply(a,b):
    n=len(a);assert len(b)==n;out=[0]*n
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            if i+j<n:out[i+j]+=x*y
            else:out[i+j-n]-=x*y
    return out


def norm(a):
    """Classical quadratic resultant descent in Z[X]/(X^N+1)."""
    n=len(a)
    if n==1:return a[0]
    e,o=a[::2],a[1::2];ee=multiply(e,e);oo=multiply(o,o);down=ee[:]
    for j,v in enumerate(oo):
        if j+1<len(oo):down[j+1]-=v
        else:down[0]+=v
    return norm(down)


def unit_step(v):
    n=len(v)
    return [v[i]+(v[i-1] if i else -v[-1])+(v[i+1] if i+1<n else -v[0]) for i in range(n)]


def centered(a,p,n,g):return [((a*pow(g,j,p)+p//2)%p)-p//2 for j in range(n//2)]


def transition(p,n,g,a,cache):
    N=n//2;m=(1+g+pow(g,-1,p))%p;next_a=a*m%p;F=centered(a,p,n,g);next_F=centered(next_a,p,n,g);raw=unit_step(F)
    assert all((x-y)%p==0 for x,y in zip(raw,next_F))
    C=[(x-y)//p for x,y in zip(raw,next_F)];assert set(C)<=set((-1,0,1))
    for label,poly in [(a,F),(next_a,next_F)]:
        if label not in cache:
            value=norm(poly);assert value>0 and value%p**(N-1)==0;cache[label]=value//p**(N-1)
    E=sum(x*x for x in F);new_E=sum(x*x for x in next_F);raw_E=sum(x*x for x in raw);dot=sum(x*y for x,y in zip(raw,C));carry_E=sum(x*x for x in C)
    assert new_E==raw_E-2*p*dot+p*p*carry_E
    second_a=next_a*m%p;second_F=centered(second_a,p,n,g);second_raw=unit_step(next_F)
    next_C=[(x-y)//p for x,y in zip(second_raw,second_F)]
    aggregate=[x+y for x,y in zip(unit_step(C),next_C)]
    assert second_F==[x-p*y for x,y in zip(unit_step(raw),aggregate)]
    return {'a':a,'next_a':next_a,'F':F,'raw_unit_image':raw,'next_F':next_F,'carry':C,
            'coefficient_energy':E,'raw_energy':raw_E,'next_energy':new_E,'carry_energy':carry_E,'raw_carry_pairing':dot,
            'norm_defect':cache[a],'next_norm_defect':cache[next_a],'radius_down_norm_up':new_E<E and cache[next_a]>cache[a],
            'second_a':second_a,'second_F':second_F,'second_carry':next_C,'two_step_carry':aggregate}


def main():
    source=HERE.parent.parent/'results/parallel41_shell_real_recovery_2026_09_06.json';old=json.loads(source.read_text())
    assert (old['p'],old['n'],old['g'],old['tau'])==(6700417,64,2,1217)
    specifications=[('small_complete',257,16,pow(3,16,257),None,None),
                    ('medium_walk',65537,32,pow(3,2048,65537),1,64),
                    ('known_short_walk',6700417,64,2,1,64),
                    ('known_centering_walk',6700417,64,2,2922708,64),
                    ('quartic_window_walk',2013265921,128,pow(31,2013265920//128,2013265921),1,32)]
    cases=[]
    for name,p,n,g,start,steps in specifications:
        N=n//2;assert pow(g,N,p)==p-1 and pow(g,n,p)==1
        epsilon=[1,1]+[0]*(N-3)+[-1];assert norm(epsilon)==1
        multiplier=(1+g+pow(g,-1,p))%p;assert multiplier
        if steps is None:scalars=list(range(1,p))
        else:
            scalars=[];a=start
            for _ in range(steps):scalars.append(a);a=a*multiplier%p
        cache={};records=[transition(p,n,g,a,cache) for a in scalars]
        counter=[r for r in records if r['radius_down_norm_up']]
        selected=max(counter,key=lambda r:Fraction(r['next_norm_defect'],r['norm_defect'])) if counter else None
        case={'name':name,'p':p,'n':n,'N':N,'g':g,'multiplier':multiplier,'unit_coefficients':epsilon,
              'complete_nonzero_scalar_census':steps is None,'start':start,'steps':steps,'p_ge_n_to_fourth':p>=n**4,
              'records':records,'distinct_norms_evaluated':len(cache),'radius_down_norm_up_steps':len(counter),
              'carry_support_histogram':sorted(Counter(r['carry_energy'] for r in records).items()),
              'largest_observed_norm_ratio_during_radius_decrease':None if selected is None else {'a':selected['a'],'next_a':selected['next_a'],'norm_defect':selected['norm_defect'],'next_norm_defect':selected['next_norm_defect'],'energy':selected['coefficient_energy'],'next_energy':selected['next_energy'],'carry_support':selected['carry_energy']}}
        cases.append(case);print(json.dumps({k:v for k,v in case.items() if k not in ('records','unit_coefficients')}),flush=True)
    out={'status':'produced','scope':'Exact centered-section unit transport and carry cocycles on actual cyclotomic cosets, including a complete small scalar census and bounded walks through a quartic-window prime. No uniform norm/shell estimate or minimum-cofactor result is asserted; independent resultant review is required.',
         'cases':cases,'source_sha256':{'norm_carry.py':sha256((HERE/'norm_carry.py').read_bytes()).hexdigest()},
         'input_sha256':{'../../results/parallel41_shell_real_recovery_2026_09_06.json':sha256(source.read_bytes()).hexdigest()}}
    (HERE/'norm_carry.json').write_text(json.dumps(out,separators=(',',':'))+'\n')


if __name__=='__main__':main()
