#!/usr/bin/env python3
"""Exact matrix-energy and boundary checks; no claim to prove the open depth bound."""
from fractions import Fraction as F
from hashlib import sha256
from itertools import combinations
from math import isqrt
from pathlib import Path
import importlib.util
import json

import numpy as np

ROOT=Path(__file__).resolve().parents[1]
SPEC=importlib.util.spec_from_file_location('transfer',ROOT/'experiments/parallel2_spectral_transfer_2026_09_04.py')
OLD=importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(OLD)


def add(x,y): return x[0]+y[0],x[1]+y[1]
def scale(x,v): return x[0]*v,x[1]*v
def minus(x,y): return add(x,scale(y,-1))


def signq(x,p):
    a,b=x
    if not a: return (b>0)-(b<0)
    if not b: return (a>0)-(a<0)
    if (a>0)==(b>0): return (a>0)-(a<0)
    diff=a*a-p*b*b
    if not diff: return 0
    return ((a>0)-(a<0)) if diff>0 else ((b>0)-(b<0))


def sqrt_upper(x, precision=10**12):
    x=F(x)
    assert x>=0
    n=isqrt(x.numerator*precision*precision//x.denominator)
    if F(n*n,precision*precision)<x: n+=1
    out=F(n,precision)
    assert out*out>=x
    return out


def upperq(x,p):
    upper=sqrt_upper(p)
    lower=F(p)/upper
    value=x[0]+x[1]*(upper if x[1]>=0 else lower)
    assert signq(minus((value,F(0)),x),p)>=0
    return value


def frobenius_pair(x,p,den=1,rows=None):
    u,v=x
    if rows is not None: u,v=u[rows,:],v[rows,:]
    return F(int(np.sum(u*u)+p*np.sum(v*v)),den),F(int(2*np.sum(u*v)),den)


def budget(p,a,m,j):
    b,n=2**a,2**a-1
    lam=2**(a-1)*(a+1)-1
    rank_sum=(a*lam**j-n**j)//(a-1) if a>=2 else j+1
    if a>=2: assert (a*lam**j-n**j)%(a-1)==0
    fj=F(rank_sum**2,n**j)
    epsilon=(a*fj+1)*sqrt_upper(F(1,p))
    g=sqrt_upper(fj/p)
    v=1+epsilon+2*g*sqrt_upper(1+epsilon)+g*g
    u=p*v+4*(a+1)*fj+(a+1)*n**j
    d0=F(b,2)*sqrt_upper(F(a*(p-1),p*n))+sqrt_upper(F(m*n*n+p-m,p*n))
    bound=(sqrt_upper(u)+j*d0*sqrt_upper(n**(j-1)))**2
    return dict(rank_sum=rank_sum,f=fj,epsilon=epsilon,g=g,v=v,u=u,d0=d0,bound=bound)


def full_case(p,anchors,chi,s,eye,jmat,depth=6):
    a=len(anchors); b=2**a; n=b-1
    cell=[x for x in range(p) if all(chi[(x-z)%p]==1 for z in anchors)]
    m=len(cell)
    assert not set(cell)&set(anchors) and m<=(p-1)//2
    d=np.array([b*int(x in cell)-1 for x in range(p)],dtype=object)
    masks={}
    for mask in range(1,1<<a):
        masks[mask]=np.array([int(np.prod([chi[(x-z)%p] for i,z in enumerate(anchors) if mask>>i&1]))
                             for x in range(p)],dtype=object)
    soft=sum(masks.values(),np.zeros(p,dtype=object))
    restored=soft.copy()
    for x in anchors: restored[x]-=b//2
    assert np.array_equal(restored,d)
    assert max(abs(int(x)) for x in soft)<=n
    w=(-jmat*d,s*d)
    powers=[(eye,np.zeros((p,p),dtype=object))]
    for k in range(1,2*depth+1): powers.append(OLD.mul(powers[-1],w,p))
    energies=[frobenius_pair(powers[j],p,p**(2*j)*n**j) for j in range(depth+1)]
    rowenergies=[frobenius_pair(powers[j],p,p**(2*j)*n**j,cell) for j in range(depth)]
    assert energies[0]==(F(p),F(0))
    assert energies[1]==(F(m*n*n+p-m,n),F(0))
    recurrence=curvature=trace_bounds=converse_bounds=0
    bias_signs=[]
    for j in range(depth):
        predicted=add(scale(energies[j],F(1,n)),scale(rowenergies[j],F(n*n-1,n)))
        assert predicted==energies[j+1]
        recurrence+=1
        bias_signs.append(signq(minus(scale(rowenergies[j],b),energies[j]),p))
    for j in range(1,depth+1):
        tr=scale(OLD.trace(powers[2*j]),F(1,p**(2*j)*n**j))
        assert signq(minus(energies[j],tr),p)>=0 and signq(add(energies[j],tr),p)>=0
        trace_bounds+=1
        absolute_tr=tr if signq(tr,p)>=0 else scale(tr,-1)
        converse=add((F(p),F(0)),scale(add(scale(absolute_tr,F(1,2)),(F(p),F(0))),F((n-1)**2*j*j,n)))
        assert signq(minus(converse,energies[j]),p)>=0
        converse_bounds+=1
        if j<depth:
            delta=add(minus(energies[j+1],scale(energies[j],2)),energies[j-1])
            assert delta==scale(tr,F((n-1)**2,n))
            curvature+=1
    if a==1: assert all(h==(F(p),F(0)) for h in energies)
    zmat=s*soft
    zpower=eye
    wordpaths={():s}
    specific=pointwise=missing=boundary_checks=0
    bounds=[]
    for j in range(1,depth+1):
        zpower=zpower@zmat
        zwiths=zpower@s
        budgetj=budget(p,a,m,j)
        bhs=F(int(np.sum(zpower*zpower)),(p*n)**j)
        kernelhs=F(int(np.sum(zwiths*zwiths)),p**(j+1)*n**j)
        colconstant=zpower@np.ones(p,dtype=object)
        consths=F(int(colconstant@colconstant),p**(j+1)*n**j)
        assert bhs==kernelhs+consths and consths<=n**j
        missing+=1
        # Decompose the original integer kernel into the exact open and omitted fibers.
        generic=[y for y in range(p) if y not in anchors]
        open_num=boundary_num=0
        for y in generic:
            for x in range(p):
                if x in anchors or x==y: boundary_num+=int(zwiths[x,y])**2
                else: open_num+=int(zwiths[x,y])**2
        exceptional_num=sum(int(zwiths[x,y])**2 for y in anchors for x in range(p))
        den=p**(j+1)*n**j
        assert F(open_num,den)<=p*budgetj['v']
        assert F(boundary_num,den)<=4*(a+1)*budgetj['f']
        assert F(exceptional_num,den)<=a*n**j
        assert F(open_num+boundary_num+exceptional_num,den)==kernelhs
        assert bhs<=budgetj['u']
        boundary_checks+=4
        # Exact Frobenius perturbation: sqrt(p*N)^j B^j = Z^j,
        # while p^j sqrt(N)^j T^j = W^j. Represent sqrt(p)^j exactly.
        if j%2: bp=(np.zeros_like(zpower),p**((j-1)//2)*zpower)
        else: bp=(p**(j//2)*zpower,np.zeros_like(zpower))
        diff=(powers[j][0]-bp[0],powers[j][1]-bp[1])
        diffnorm=frobenius_pair(diff,p,p**(2*j)*n**j)
        assert upperq(diffnorm,p)<=j*j*budgetj['d0']**2*n**(j-1)
        assert upperq(energies[j],p)<=budgetj['bound']
        bounds.append(dict(j=j,energy=[str(q) for q in energies[j]],upper_bound=str(budgetj['bound']),
                           below_trivial=budgetj['bound']<p*n**j))
        if p in (5,13) and j<=3:
            wordpaths={word+(mask,):s@(dm[:,None]*path) for word,path in wordpaths.items() for mask,dm in masks.items()}
            assert np.array_equal(sum(wordpaths.values(),np.zeros_like(s)),zwiths)
            specific+=1
            for word,path in wordpaths.items():
                cs=[0]*a
                for mask in word:
                    r=1+sum(cs)
                    cs=[r if mask>>v&1 else cs[v] for v in range(a)]
                r=1+sum(cs)
                for y in generic:
                    for x in range(p):
                        assert int(path[x,y])**2<=4*r*r*p**j
                        pointwise+=1
    return dict(p=p,anchors=anchors,cell=cell,depth=depth,recurrence_checks=recurrence,
                curvature_checks=curvature,trace_bounds=trace_bounds,converse_bounds=converse_bounds,missing_direction_checks=missing,
                boundary_checks=boundary_checks,specific_word_sums=specific,pointwise_checks=pointwise,
                bias_signs=bias_signs,full_bounds=bounds)


def jordan_edge():
    t=np.array([[F(8,5),F(-3,10)],[F(6,5),F(2,5)]],dtype=object)
    eye=np.eye(2,dtype=object)
    assert np.array_equal((t-eye)@(t-eye),np.zeros((2,2),dtype=object))
    power=eye
    for j in range(41):
        assert np.trace(power)==2
        assert np.sum(power*power)==2+F(9,4)*j*j
        power=power@t
    return dict(N=4,dimension=2,checks=41,meaning='An exact abstract projection block at the spectral edge, not a Paley example; H_j=2+9j^2/4.')


def large_first_step():
    p=10009
    assert all(p%d for d in range(2,isqrt(p)+1))
    chi=[0]+[1 if pow(x,(p-1)//2,p)==1 else -1 for x in range(1,p)]
    m=sum(chi[x]==1 and chi[(x-1)%p]==1 for x in range(p))
    value=F(9*m+p-m,3)
    bd=budget(p,2,m,1)['bound']
    assert value<=bd<2*p
    return dict(p=p,anchors=[0,1],cell_size=m,j=1,exact_energy=str(value),upper_bound=str(bd),
                upper_bound_below_twice_p=True,meaning='A nonvacuous finite instance of the full bound; not asymptotic evidence.')


def main():
    cases=[]
    for p in (5,13,17,29,41):
        chi,s,eye,jmat=OLD.field(p)
        for a in (1,2,3):
            count=0
            for anchors in combinations(range(p),a):
                if all(chi[(x-y)%p]==1 for x,y in combinations(anchors,2)):
                    cases.append(full_case(p,list(anchors),chi,s,eye,jmat))
                    count+=1
                    if count==2: break
        print(json.dumps(dict(p=p,cases_completed=sum(c['p']==p for c in cases))),flush=True)
    inputs=['research/parallel11-full-energy-2026-09-05.md','experiments/parallel11_full_energy_2026_09_05.py',
            'research/parallel10-word-aggregate-2026-09-05.md','research/parallel9-all-degrees-2026-09-05.md',
            'research/parallel2-spectral-transfer-2026-09-04.md','experiments/parallel2_spectral_transfer_2026_09_04.py',
            'sources/kunisky-2303.16475v1.html','sources/katz-rigid-local-systems.pdf',
            'sources/bbd-faisceaux-pervers.pdf','sources/katz-gauss-kloosterman-monodromy.pdf']
    report=dict(status='All finite exact checks passed; the uniform logarithmic range uses the written source-dependent proof. The longer-depth target remains open.',
                arithmetic='Integer object matrices, rational pairs u+v sqrt(p), exact sign comparisons, and rational upper square roots.',
                scope='Complete spectral trace including J, anchor corrections, exceptional fibers and the constant direction.',
                cases=cases,jordan_edge=jordan_edge(),large_first_step=large_first_step(),
                totals={k:sum(c[k] for c in cases) for k in ('recurrence_checks','curvature_checks','trace_bounds','converse_bounds','missing_direction_checks','boundary_checks','specific_word_sums','pointwise_checks')},
                input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    report['totals'].update(case_count=len(cases),full_upper_bounds=sum(len(c['full_bounds']) for c in cases),
                            positive_biases=sum(s>0 for c in cases for s in c['bias_signs']),negative_biases=sum(s<0 for c in cases for s in c['bias_signs']))
    out=ROOT/'results/parallel11_full_energy_2026_09_05.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps(dict(written=str(out),totals=report['totals'])),flush=True)


if __name__=='__main__': main()
