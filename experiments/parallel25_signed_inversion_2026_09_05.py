#!/usr/bin/env python3
"""Actual-character sign halves, completed averages, and Mobius cocycles."""
from collections import Counter
from itertools import combinations, combinations_with_replacement, product
from math import comb, isqrt
from pathlib import Path
from hashlib import sha256
import json

ROOT = Path(__file__).resolve().parents[1]
COUNT = Counter()

def chi_table(p):
    return [0] + [1 if pow(x, (p-1)//2, p)==1 else -1 for x in range(1,p)]

def fvec(C,ch,w=None):
    if w is None:w=[1]*len(C)
    return [sum(a*ch[(x-c)%len(ch)] for c,a in zip(C,w)) for x in range(len(ch))]

def moment(F,r=6):return sum(v**r for v in F)

def radical_le(lhs,A,p,B):
    # lhs <= A sqrt(p) + B; no floating-point acceptance.
    return lhs<=B or (lhs-B)**2<=A*A*p

def sidon(C,p):
    sums=[(a+b)%p for a,b in combinations_with_replacement(C,2)]
    return len(sums)==len(set(sums))

def greedy_sidon(p,n):
    C=[]
    for c in range(p):
        if sidon(C+[c],p):C.append(c)
        if len(C)==n:return C
    raise AssertionError('greedy construction failed')

def cocycles(C,p,ch):
    n=len(C)
    poles=[z for z in range(p) if z not in C]
    for z in poles[:2]+poles[-1:]:
        D=[pow((c-z)%p,-1,p) for c in C]
        w=[ch[(c-z)%p] for c in C]
        choices=[u for u in range(p) if u not in D]
        for u in choices[:3]:
            E=[pow((d-u)%p,-1,p) for d in D]
            ww=[a*ch[(d-u)%p] for a,d in zip(w,D)]
            if u==0:
                assert E==[(c-z)%p for c in C] and ww==[1]*n
            else:
                v=(z+pow(u,-1,p))%p
                assert v not in C
                assert E==[(-pow(u,-1,p)-pow(u*u%p,-1,p)*pow((c-v)%p,-1,p))%p for c in C]
                assert ww==[ch[-u%p]*ch[(c-v)%p] for c in C]
            COUNT['two_inversion_point_and_weight_cocycles']+=1
    w=[1 if j%2==0 else -1 for j in range(n)]
    F=fvec(C,ch,w)
    for a,b,c,d in [(1,0,0,1),(2,3,1,1),(1,1,1,2),(0,1,1,0)]:
        det=(a*d-b*c)%p
        if not det or any((c*t+d)%p==0 for t in C):continue
        E=[(a*t+b)*pow((c*t+d)%p,-1,p)%p for t in C]
        ww=[v*ch[(c*t+d)%p] for t,v in zip(C,w)]
        FE=fvec(E,ch,ww)
        for x in range(p):
            if (c*x+d)%p:
                y=(a*x+b)*pow((c*x+d)%p,-1,p)%p
                assert FE[y]==ch[det]*ch[(c*x+d)%p]*F[x]
                COUNT['mobius_weighted_coordinates']+=1
        assert moment(FE)+sum(ww)**6==moment(F)+sum(w)**6
        COUNT['mobius_augmented_moment_invariance']+=1

def inspect(C,p,tensor=False):
    C=list(C);n=len(C);ch=chi_table(p);eps=ch[-1]
    F=fvec(C,ch);M=moment(F);M4=moment(F,4);M2=moment(F,2)
    assert M2==p*n-n*n
    T=[0]*7;full_plus=full_minus=outside_plus=outside_minus=0
    records=[]
    for z in range(p):
        Q=[sum(ch[(x-c)%p]*ch[(z-c)%p] for c in C) for x in range(p)]
        for j in range(7):T[j]+=sum(f**(6-j)*q**j for f,q in zip(F,Q))
        jp=sum((q+eps*f)**6 for f,q in zip(F,Q))
        jm=sum((q-eps*f)**6 for f,q in zip(F,Q))
        full_plus+=jp;full_minus+=jm
        if z in C:continue
        D=[pow((c-z)%p,-1,p) for c in C]
        plus=[d for d in D if ch[d]==1];minus=[d for d in D if ch[d]==-1]
        mp=moment(fvec(plus,ch));mm=moment(fvec(minus,ch))
        assert 64*mp==jp and 64*mm==jm
        COUNT['actual_sign_half_moment_identities']+=2
        assert 32*(mp+mm)>=M
        assert len(plus)-len(minus)==eps*F[z]
        COUNT['joint_positive_lower_and_size_identities']+=1
        outside_plus+=mp;outside_minus+=mm
        records.append((z,mp+mm,moment(fvec(D,ch)),sidon(D,p)))
    A4=sum(F[c]**4 for c in C)
    assert T[0]==p*M and T[1]==0
    assert T[2]==p*(n*M4-A4)-M
    assert full_plus+full_minus==2*(p*M+15*T[2]+15*T[4]+T[6])
    assert full_plus-full_minus==2*eps*(20*T[3]+6*T[5])
    COUNT['complete_even_and_odd_average_identities']+=1
    assert radical_le(abs(T[3]),30*p*n**3,p,10*p*n**6)
    assert radical_le(abs(T[5]),60*p*n**3,p,20*p*n**6)
    assert radical_le(M4,3*n**4,p,3*p*n*n)
    assert radical_le(T[4],3*n**4*M2,p,3*p*n*n*M2)
    assert T[6]<=15*p*p*n**3+25*p*n**6
    COUNT['mixed_character_sum_bounds']+=5
    # Completed J values have denominator 64.
    assert radical_le(abs(full_plus-full_minus),1920*p*n**3,p,640*p*n**6)
    actual_diff=2*abs(outside_plus-outside_minus)-n*M
    assert radical_le(actual_diff,60*p*n**3,p,21*p*n**6)
    COUNT['completed_and_actual_sign_asymmetry_bounds']+=2
    if tensor:
        direct=[0]*7
        for word in product(C,repeat=6):
            K6=sum(prod_ch(word,x,ch) for x in range(p))
            for j in range(7):
                Kj=p if j==0 else sum(prod_ch(word[6-j:],z,ch) for z in range(p))
                direct[j]+=K6*Kj
        assert direct==T
        COUNT['independent_ordered_tuple_expansions']+=1
    if p>=n**4 and n>=2:
        B=(p-n-15)*M+195*p*p*n**3+25*p*n**6
        assert 32*(outside_plus+outside_minus)<=B
        U=15*p*p*n**3+25*p*n**6
        Nh=comb(n+1,2);Qh=n+2*comb(Nh,2);H=p-Qh
        if H>0:
            good=[v for v in records if v[3]]
            assert len(good)>=H
            both=[v for v in good if v[2]*H<=3*U and 32*v[1]*H<=3*B]
            assert 3*len(both)>=H
            COUNT['simultaneous_good_pole_fraction']+=1
        COUNT['critical_slice_joint_upper']+=1
    cocycles(C,p,ch)
    return {'p':p,'C':C,'sidon':sidon(C,p),'M6':M,
            'sum_actual_plus_M6':outside_plus,'sum_actual_minus_M6':outside_minus,
            'complete_T_j':T}

def prod_ch(word,x,ch):
    v=1
    for c in word:v*=ch[(x-c)%len(ch)]
    return v

def main():
    records=[]
    for p in [5,7,11,13]:
        for n in [2,3]:
            for C in combinations(range(p),n):
                records.append(inspect(C,p,tensor=(C==(0,1,2) and p in [7,11])))
    for p,n in [(17,2),(19,2),(83,3),(89,3),(257,4),(1297,6)]:
        records.append(inspect(greedy_sidon(p,n),p))
    inputs=['research/parallel24-inversion-orbit-upper-2026-09-05.md',
            'research/parallel24-inversion-independent-review-2026-09-05.md',
            'research/parallel24-classical-exceptions-2026-09-05.md',
            'research/parallel20-inversion-moments-2026-09-05.md',
            'experiments/parallel25_signed_inversion_2026_09_05.py']
    out={'status':'PASS; exact arithmetic identities and feedback bounds, no uniform target saving',
         'counts':dict(COUNT),'field_cases':records,
         'input_sha256':{f:sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs},
         'arithmetic':'Python integers; radical inequalities checked by signed squaring'}
    dest=ROOT/'results/parallel25_signed_inversion_2026_09_05.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'output':str(dest),'counts':dict(COUNT),'field_cases':len(records)},indent=2))

if __name__=='__main__':main()
