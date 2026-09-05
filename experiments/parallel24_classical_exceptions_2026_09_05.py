#!/usr/bin/env python3
"""Exact bounded audit of local drift, shells, and transported signs."""
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations
from math import comb,isqrt
from pathlib import Path
from hashlib import sha256
import json

ROOT=Path(__file__).resolve().parents[1]
COUNT=Counter()

def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))
def chars(p):return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
def fvec(C,ch,weights=None):
    if weights is None:weights=[1]*len(C)
    return [sum(w*ch[(x-c)%len(ch)] for c,w in zip(C,weights)) for x in range(len(ch))]
def m6(F):return sum(x**6 for x in F)
def variance(vals):return Q(sum(x*x for x in vals),len(vals))-Q(sum(vals),len(vals))**2

def small():
    reports=[]
    for p in [5,7,11,13]:
        ch=chars(p)
        for n in range(2,min(4,p-1)+1):
            sets=list(combinations(range(p),n))
            fs={C:fvec(C,ch) for C in sets}
            ms={C:m6(fs[C]) for C in sets}
            den=n*(p-n)
            for C in sets:
                Cset=set(C);F=fs[C]
                assert sum(x*x for x in F)==p*n-n*n
                actual=sum(ms[tuple(sorted((Cset-{a})|{b}))] for a in C for b in range(p) if b not in Cset)
                rhs=(den-6*(p-5))*ms[C]
                for x,v in enumerate(F):
                    N=n-int(x in Cset)
                    S=n*(p-1)+(p-2*n)*N
                    rhs+=(15*S-80*p+180)*v**4
                    rhs+=(15*S+90*N*(p-1-N)-96*p+122)*v*v
                    rhs+=S+30*N*(p-1-N)
                    COUNT['row_drift_coefficients']+=1
                assert actual==rhs
                COUNT['complete_swap_drift_identities']+=1
            for C in [sets[0],sets[-1]]:
                F=fs[C]
                for t in range(n+1):
                    shell=[D for D in sets if len(set(C)&set(D))==t]
                    if not shell:continue
                    lam=Q(p*t-n*n,n*(p-n))
                    for x in range(p):
                        assert Q(sum(fs[D][x] for D in shell),len(shell))==lam*F[x]
                    assert Q(sum(ms[D] for D in shell),len(shell))>=lam**6*ms[C]
                    COUNT['exact_shell_mean_and_moment_checks']+=1
                # Local insertion variance uses actual A moments, not a global average.
                A=C[:-1];FA=fvec(A,ch)
                H=[6*x**5+20*x**3+6*x for x in FA]
                Z=[x**4+x*x for x in FA]
                QH=[sum(ch[(x-b)%p]*H[x] for x in range(p)) for b in range(p)]
                X=sum(x*x for x in H);Y=sum(x*x for x in Z)
                assert sum(x*x for x in QH)==p*X-sum(H)**2
                vals=[ms[tuple(sorted(A+(b,)))] for b in range(p) if b not in A]
                lhs=variance(vals)*len(vals)-p*X-225*Y
                assert lhs<=0 or lhs*lhs<=900*p*X*Y
                COUNT['local_full_transform_variance_checks']+=1
                strong=sum(2*abs(F[z])>n for z in range(p) if z not in C)
                assert strong*n*n<=4*(p*n-n*n)
                COUNT['unbalanced_pole_second_moment_checks']+=1
                for z in range(p):
                    if z in C:continue
                    D=[pow((c-z)%p,-1,p) for c in C]
                    weights=[ch[d] for d in D]
                    assert weights==[ch[(c-z)%p] for c in C]
                    plus=sum(w==1 for w in weights)
                    eta=1 if plus>=n-plus else -1
                    E=[d for d,w in zip(D,weights) if w!=eta]
                    GD=fvec(D,ch);GE=fvec(E,ch);GW=fvec(D,ch,weights)
                    assert all(w==eta*(d-2*e) for w,d,e in zip(GW,GD,GE))
                    assert len(E)==(n-abs(F[z]))//2
                    assert m6(GW)==ms[C]-abs(F[z])**6+n**6
                    COUNT['actual_transported_sign_split_identities']+=1
            reports.append({'p':p,'n':n,'complete_slice_sets':len(sets)})
    return reports

def volume_certificates():
    reports=[]
    for n in [6,8,16,32,64,128,256]:
        p=n**4+1
        if p%2==0:p+=1
        while not prime(p):p+=2
        assert n**4<=p<(n+1)**4
        t=next(t for t in range(1,n+1) if (15+5*n)*t**6>15*n**6)
        assert t**6>n**5
        mass=Q(sum(comb(n,j)*comb(p-n,n-j) for j in range(t,n+1)),comb(p,n))
        p_dependent=Q(2*n*n,p)**t
        coarse=Q(2,n*n)**t
        assert mass<=p_dependent<=coarse<Q(n*n,p)
        assert Q(n,p-n)**6*(15+5*n)<15
        # Positive lower-order drift terms on the actual thin slice.
        den=n*(p-n)
        for N in [n-1,n]:
            S=n*(p-1)+(p-2*n)*N
            assert 15*S-80*p+180>=0
            assert 15*S+90*N*(p-1-N)-96*p+122>=0
            assert S+30*N*(p-1-N)>=0
        COUNT['thin_slice_volume_and_drift_certificates']+=1
        reports.append({'p':p,'n':n,'necessary_retained_points':t,
                        'p_dependent_global_tail_bound':{'rational_base':str(Q(2*n*n,p)),'exponent':t},
                        'coarse_global_tail_bound':{'rational_base':str(Q(2,n*n)),'exponent':t},
                        'exact_tail_fraction_numerator_bits':mass.numerator.bit_length(),
                        'exact_tail_fraction_denominator_bits':mass.denominator.bit_length(),
                        'tail_below_n2_over_p':True,
                        'scope':'integer combinatorial certificate; no large-field character computation'})
    return reports

def main():
    result={'status':'exact identities and quantitative losses verified; no exceptional-set upper bound proved',
            'small_actual_fields':small(),'thin_slice_certificates':volume_certificates()}
    inputs=['research/parallel23-classical-upper-2026-09-05.md',
            'research/parallel23-classical-independent-review-2026-09-05.md',
            'research/parallel20-inversion-moments-2026-09-05.md',
            'experiments/parallel24_classical_exceptions_2026_09_05.py']
    result['input_sha256']={f:sha256((ROOT/f).read_bytes()).hexdigest() for f in inputs}
    result['counts']=dict(COUNT)
    result['arithmetic']='Python integers, exact rational shell probabilities and radical comparisons; no floating acceptance'
    out=ROOT/'results/parallel24_classical_exceptions_2026_09_05.json'
    out.write_text(json.dumps(result,indent=2)+'\n')
    print(json.dumps({'output':str(out),'counts':result['counts']},indent=2))

if __name__=='__main__':main()
