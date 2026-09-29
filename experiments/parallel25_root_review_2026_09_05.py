#!/usr/bin/env python3
"""Independent root checks; no imports from the three worker verifiers."""
from collections import Counter
from hashlib import sha256
from itertools import combinations, combinations_with_replacement, product
from math import comb, factorial
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()
def check(v,k):
    assert v,k
    COUNTS[k]+=1
def chi(p):
    return [0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
def F(C,ch,weights=None):
    weights=weights or [1]*len(C)
    return [sum(w*ch[(x-c)%len(ch)] for c,w in zip(C,weights)) for x in range(len(ch))]
def M(v):return sum(x**6 for x in v)

def signed():
    report=[]
    for p,C in ((31,[0,2,9]),(41,[1,4,13,22])):
        ch=chi(p);n=len(C);fv=F(C,ch);m=M(fv)
        T=[0]*7; actual_plus=actual_minus=0
        for z in range(p):
            Q=[sum(ch[(x-c)%p]*ch[(z-c)%p] for c in C) for x in range(p)]
            for j in range(7):T[j]+=sum(fv[x]**(6-j)*Q[x]**j for x in range(p))
            if z in C:continue
            halves=[[pow((c-z)%p,-1,p) for c in C if ch[(c-z)%p]==sgn] for sgn in (1,-1)]
            ms=[M(F(D,ch)) for D in halves]
            actual_plus+=ms[0];actual_minus+=ms[1]
            for sgn,md in zip((1,-1),ms):
                check(64*md==sum((Q[x]+sgn*ch[-1]*fv[x])**6 for x in range(p)),
                      'independent_actual_half_moments')
            check(32*sum(ms)==m+15*sum(fv[x]**4*Q[x]**2+fv[x]**2*Q[x]**4 for x in range(p))+sum(v**6 for v in Q),
                  'independent_joint_half_identity')
        A4=sum(fv[c]**4 for c in C)
        check(T[1]==0 and T[2]==p*(n*sum(v**4 for v in fv)-A4)-m,'independent_mixed_complete_identity')
        for g in ((2,3,1,5),(1,4,0,2),(0,1,1,7)):
            a,b,c,d=g;det=(a*d-b*c)%p
            if not det or any((c*u+d)%p==0 for u in C):continue
            w=[1,-1,1,-1][:n]
            D=[(a*u+b)*pow((c*u+d)%p,-1,p)%p for u in C]
            wp=[ww*ch[(c*u+d)%p] for u,ww in zip(C,w)]
            check(M(F(D,ch,wp))+sum(wp)**6==M(F(C,ch,w))+sum(w)**6,'independent_augmented_Mobius_moment')
        report.append({'p':p,'C':C,'M6':m,'sum_actual_plus':actual_plus,'sum_actual_minus':actual_minus})
    return report

def subgroup():
    report=[]
    for p,n in ((97,8),(33713,16),(37201,16)):
        H=[x for x in range(1,p) if pow(x,n,p)==1]
        assert len(H)==n
        r2=Counter((a+b)%p for a in H for b in H)
        e2=sum(x*x for x in r2.values())
        r4=sum(sum(t)%p==2 for t in product(H,repeat=4))
        totals=Counter()
        for t in combinations_with_replacement(H,6):
            if sum(t)%p:continue
            freq=Counter(t);weight=factorial(6)
            for v in freq.values():weight//=factorial(v)
            nu=sum(comb(v,2) for v in freq.values())
            intr=all(freq[x]==freq[-x%p] for x in freq)
            free=all(-x%p not in freq for x in freq)
            totals['E3']+=weight;totals['nu']+=weight*nu
            if intr:totals['intrinsic']+=weight;totals['intrinsic_nu']+=weight*nu
            elif not free:totals['J6']+=weight
            if free:
                balanced=False
                for idx in combinations(range(6),3):
                    a=b=1
                    for j,x in enumerate(t):
                        if j in idx:a=a*x%p
                        else:b=b*x%p
                    balanced |= (a+b)%p==0
                if len(freq)<6:totals['repeated']+=weight
                elif balanced:totals['distinct_balanced']+=weight
                else:totals['D6']+=weight
        check(totals['nu']==15*n*r4,'independent_equal_pair_incidence')
        check(totals['intrinsic_nu']==15*n*(6*n-8),'independent_intrinsic_incidence')
        check(totals['repeated']<=15*n*(r4-6*n+8),'independent_repeated_upper')
        check(r4<=e2,'independent_four_sum_Cauchy')
        T6=15*n**3-45*n*n+40*n
        check(totals['intrinsic']==T6 and totals['E3']==T6+totals['J6']+totals['repeated']+totals['distinct_balanced']+totals['D6'],
              'independent_disjoint_sixth_partition')
        report.append({'p':p,'n':n,'quartic_window':n**4//4<=p<=n**4,'r4_at_2':r4,'E2':e2,'counts':dict(totals)})
    # Check the largest stored witness through a different resultant route:
    # integer Bareiss determinant, not recursive quadratic norm descent.
    prev=json.loads((ROOT/'results/parallel25_subgroup_growing_orders_2026_09_05.json').read_text())
    witness=prev['cases'][-1]['cyclic_norm_witness'];f=witness['primitive_coefficients'];d=len(f)
    A=[[0]*d for _ in range(d)]
    for col in range(d):
        for j,v in enumerate(f):
            q,r=divmod(col+j,d);A[r][col]+=(-1 if q%2 else 1)*v
    sign=1;old=1
    for k in range(d-1):
        pivot=next(i for i in range(k,d) if A[i][k])
        if pivot!=k:A[k],A[pivot]=A[pivot],A[k];sign=-sign
        a=A[k][k]
        for i in range(k+1,d):
            for j in range(k+1,d):
                numerator=a*A[i][j]-A[i][k]*A[k][j]
                assert numerator%old==0
                A[i][j]=numerator//old
            A[i][k]=0
        old=a
    det=abs(sign*A[-1][-1])
    check(det==witness['norm'],'independent_128_dim_Bareiss_norm')
    p=prev['cases'][-1]['p']
    check(det%p==0 and det%(p*p)!=0,'independent_single_prime_norm_factor')
    return report

def spectral():
    report=[]
    for p in (13,17,61,269):
        ch=np.array(chi(p),dtype=np.int64);allx=np.arange(p)
        S=ch[(allx[:,None]-allx)%p]
        C=np.array([x for x in range(p) if ch[x]==ch[(x-1)%p]==1]);m=len(C)
        # Bound all integer matrix products below by p^4 before array use.
        check(p**4<np.iinfo(np.int64).max,'independent_spectral_integer_safety')
        SC=S[np.ix_(C,C)];one=np.ones((m,m),dtype=np.int64);eye=np.eye(m,dtype=np.int64)
        K=(S*ch[None,:])@S;K0=K[np.ix_(C,C)]
        idx={int(x):i for i,x in enumerate(C)}
        rp=np.array([idx[(1-int(x))%p] for x in C]);ip=np.array([idx[pow(int(x),-1,p)] for x in C])
        K1=K0[np.ix_(rp,rp)];A3=K0+K1+K1[np.ix_(ip,ip)]
        check(np.array_equal(4*(SC@SC),p*eye+A3-6*one),'independent_full_square_kernel')
        # Compute both irrational and rational coefficients of the actual
        # off-block leakage directly from the full involution.
        O=np.array([x for x in range(p) if x not in set(C)])
        off=S[np.ix_(O,C)]
        rational=p*(off.T@off)+(p-m)*one
        irrational=-(off.sum(axis=0)[:,None]+off.sum(axis=0)[None,:])
        check(np.array_equal(4*rational,3*p*p*eye-p*A3+(6*p-4*m)*one),'independent_direct_leakage_rational')
        check(np.array_equal(irrational,SC@one+one@SC),'independent_direct_leakage_surd')
        R=eye[:,rp];Inv=eye[:,ip];U=R@Inv;U2=U@U
        Ns=[(eye+R)@(eye+U+U2),(eye-R)@(eye+U+U2),
            (eye-R)@(2*eye-U-U2),(eye+R)@(2*eye-U-U2)]
        check(np.array_equal(sum(Ns),6*eye),'independent_full_sector_projection_sum')
        for N in Ns:
            check(np.array_equal(N@N,6*N) and np.array_equal(N@SC,SC@N),'independent_sector_invariance')
        dims=[int(np.trace(N))//6 for N in Ns]
        check(sum(dims)==m and dims[2]==dims[3],'independent_sector_dimensions')
        report.append({'p':p,'m':m,'dimensions':dims,'outside_worker_range':p>257})
    return report

def main():
    report={'status':'passed','signed':signed(),'subgroup':subgroup(),'spectral':spectral()}
    paths=['research/parallel25-signed-inversion-2026-09-05.md',
           'research/parallel25-subgroup-growing-orders-2026-09-05.md',
           'research/parallel25-spectral-full-operator-2026-09-05.md',
           'results/parallel25_signed_inversion_2026_09_05.json',
           'results/parallel25_subgroup_growing_orders_2026_09_05.json',
           'results/parallel25_spectral_full_operator_2026_09_05.json',
           'sources/parallel25-mrss-1712.00410v1.html',
           'experiments/parallel25_root_review_2026_09_05.py']
    report['input_sha256']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in paths}
    report['counts']=dict(COUNTS)
    report['scope']='Root review of worker deductions; no separate-agent review of the root-authored prize theorem, no human or formal verification.'
    out=ROOT/'results/parallel25_root_review_2026_09_05.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':report['status'],'counts':report['counts']},indent=2))

if __name__=='__main__':main()
