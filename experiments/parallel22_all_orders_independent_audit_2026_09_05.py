from collections import Counter
from itertools import product, combinations_with_replacement
from math import prod
from fractions import Fraction
from pathlib import Path
from hashlib import sha256
import json
import random

ROOT = Path(__file__).resolve().parents[1]
counts=Counter()

def rows(k,p):
    L=k*k; W=p-k-L; den=1<<k
    out=[]
    for mask in range(1<<k):
        F=k-2*mask.bit_count()
        e2=(F*F-k)//2
        w=W-(L+1)*e2
        assert 2*w>=W and 2*w<=3*W
        if mask in (0,(1<<k)-1):
            w+=L*(1<<(k-1))
        out.append((mask,0,w))
    for zero in range(k):
        for mask in range(1<<k):
            if not (mask>>zero)&1:
                out.append((mask,1<<zero,2))
    assert sum(w for _,_,w in out)==p*den
    counts['explicit_weighted_models']+=1
    counts['explicit_weighted_rows']+=len(out)
    return out,den

def us(s,m):
    if m==0:
        return 1 if s==0 else 0
    # Independent symmetric binomial distribution computation.
    from math import comb
    return Fraction(sum(comb(m,j)*(m-2*j)**(2*s) for j in range(m+1)),1<<m)

def moment(R,den,b,s):
    total=0
    for mask,zero,w in R:
        g=sum(a*(-1 if mask&(1<<i) else 1) for i,a in enumerate(b) if not zero&(1<<i))
        total+=w*g**(2*s)
    return Fraction(total,den)

for k in range(2,101):
    for p in [k**4,k**4+1,k**5+17]:
        W=p-k-k*k
        assert W>0
        assert W-(k*k+1)*k*(k-1)>=k*k*(k-2)
        assert W-(k*k+1)*k>=k*(k-2)*(k*k+k+1)
        for j in range(k+1):
            F=k-2*j;e2=(F*F-k)//2
            rho=Fraction(W-(k*k+1)*e2,W)
            assert Fraction(1,2)<=rho<=Fraction(3,2)
            counts['density_sign_sum_checks']+=1
        counts['positivity_parameter_checks']+=1

rng=random.Random(2219)
fixtures=[]
for r,k in [(3,j) for j in range(2,9)]+[(4,10),(5,12)]:
    p=k**(r+1)
    R,den=rows(k,p)
    L=k*k
    for i in range(k):
        zsum=sum(w for _,zero,w in R if zero&(1<<i))
        assert zsum==den
        counts['zero_mass_checks']+=1
    if k<=8:
        for D in range(1,1<<k):
            total=sum(w*(-1 if (mask&D).bit_count()%2 else 1) for mask,zero,w in R if not zero&D)
            d=D.bit_count()
            expected=0 if d%2 else -1 if d==2 else L
            assert total==expected*den
            counts['all_subset_correlations']+=1
    coefficients=[[1]*k,[(-1)**j for j in range(k)],[int(j==0) for j in range(k)],
                  [rng.randrange(-3,4) for _ in range(k)],[0]*k]
    for s in range(1,r):
        Ds=prod(range(1,2*s,2))
        for b in coefficients:
            m=moment(R,den,b,s)
            assert 2*m<=(3*Ds+2)*p*sum(a*a for a in b)**s
            counts['real_coefficient_integer_fixtures']+=1
        for m in range(1,k+1):
            a=moment(R,den,[1]*m+[0]*(k-m),s)
            V=(us(s+1,m)-m*us(s,m))/2
            assert V>=0
            formula=L*m**(2*s)+(p-L-m)*us(s,m)-(L+1)*V+m*us(s,m-1)
            assert a==formula
            assert a<=(Ds+1)*p*m**s
            if s==1: assert a==p*m-m*m
            counts['exact_unsigned_formula_checks']+=1
    target=moment(R,den,[1]*k,r)
    assert target>=k**(2*r+2)
    if k>=2*(r+1):
        assert 2*target>=p*k**(r+1)
        counts['next_moment_growth_checks']+=1
    fixtures.append({'r':r,'k':k,'p_integer_parameter':p,'target_ratio':str(target/(p*k**r))})

for r in range(3,31):
    q=r+1
    for k in [2*q,2*q+1,10*q]:
        assert (k+1)**q<2*k**q
        counts['size_ratio_checks']+=1
        for h in range(2,q//2+1):
            f=(k-1)+sum((k-1)**(2*h-j) for j in range(1,h+1))
            assert f<(h+1)*k**(2*h-1)<k**(2*h)<=k**q
            counts['relation_label_threshold_checks']+=1
result={'scope':'independent finite consistency audit; general claims reviewed algebraically; large integer fixtures are not prime-field kernels',
        'reviewed_proof_sha256':sha256((ROOT/'research/parallel22-all-orders-obstruction-2026-09-05.md').read_bytes()).hexdigest(),
        'counts':dict(sorted(counts.items())),'fixtures':fixtures}
(ROOT/'results/parallel22_all_orders_independent_audit_2026_09_05.json').write_text(json.dumps(result,indent=2)+'\n')
print(json.dumps(result,indent=2))
