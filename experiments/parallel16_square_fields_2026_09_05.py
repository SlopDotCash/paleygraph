#!/usr/bin/env python3
"""Actual square-field Paley matrices: sharp bounds coexist with outliers.

These are NOT counterexamples over prime fields. All acceptance is exact.
"""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from math import isqrt,prod
from pathlib import Path
import json
import numpy as np

ROOT=Path(__file__).resolve().parents[1]


class QuadraticField:
    def __init__(self,ell):
        assert ell>2 and all(ell%d for d in range(2,isqrt(ell)+1))
        self.ell,self.q=ell,ell*ell
        cb=[0]+[1 if pow(x,(ell-1)//2,ell)==1 else -1 for x in range(1,ell)]
        self.nu=next(x for x in range(2,ell) if cb[x]==-1)
        self.coords=[(x%ell,x//ell) for x in range(self.q)]
        self.chi=[cb[(a*a-self.nu*b*b)%ell] for a,b in self.coords]
        for x in range(1,self.q):
            assert self.mul(x,self.inv(x))==1
            expected=1 if self.chi[x]==1 else ell-1
            assert self.power(x,(self.q-1)//2)==expected
        assert self.chi[0]==0 and all(self.chi[x]==1 for x in range(1,ell))
    def add(self,x,y):
        a,b=self.coords[x];c,d=self.coords[y];p=self.ell
        return (a+c)%p+p*((b+d)%p)
    def neg(self,x):
        a,b=self.coords[x];p=self.ell
        return -a%p+p*(-b%p)
    def sub(self,x,y):
        return self.add(x,self.neg(y))
    def mul(self,x,y):
        a,b=self.coords[x];c,d=self.coords[y];p=self.ell
        return (a*c+self.nu*b*d)%p+p*((a*d+b*c)%p)
    def inv(self,x):
        assert x
        a,b=self.coords[x];p=self.ell
        inv=pow((a*a-self.nu*b*b)%p,-1,p)
        return a*inv%p+p*(-b*inv%p)
    def div(self,x,y):
        return self.mul(x,self.inv(y))
    def power(self,x,n):
        out=1
        while n:
            if n&1:out=self.mul(out,x)
            x=self.mul(x,x);n//=2
        return out


class Curve:
    def __init__(self,F,r=2):
        self.F,self.r=F,r
        roots={}
        for y in range(F.q):roots.setdefault(F.mul(y,y),[]).append(y)
        self.points=[None]+[(x,y) for x in range(F.q)
                            for y in roots.get(F.mul(F.mul(x,F.sub(x,1)),F.sub(x,r)),[])]
        self.index={P:i for i,P in enumerate(self.points)}
    def add(self,P,Q):
        if P is None:return Q
        if Q is None:return P
        F,r=self.F,self.r;x,y=P;u,v=Q
        if x==u and F.add(y,v)==0:return None
        if P==Q:
            numerator=F.add(F.sub(F.mul(3,F.mul(x,x)),F.mul(F.mul(2,F.add(1,r)),x)),r)
            slope=F.div(numerator,F.mul(2,y))
        else:
            slope=F.div(F.sub(v,y),F.sub(u,x))
        z=F.sub(F.sub(F.add(F.mul(slope,slope),F.add(1,r)),x),u)
        return z,F.sub(F.mul(slope,F.sub(x,z)),y)
    def neg(self,P):return None if P is None else (P[0],self.F.neg(P[1]))
    def rho(self,P):
        Q=self.add(P,P)
        return self.F.q if Q is None else Q[0]
    def f(self,P):return 0 if P is None else self.F.chi[P[1]]


def psd_certificate(matrix):
    """Fraction-free symmetric elimination; zero pivot requires a zero row.

This also handles genuine equality in the square-field sharp bound.
"""
    b=[[int(x) for x in row] for row in matrix]
    assert b==[list(row) for row in zip(*b)]
    previous=1;rank=0;zero=0;bits=0
    for k in range(len(b)):
        pivot=b[k][k]
        assert pivot>=0,('negative pivot',k,pivot)
        if pivot==0:
            assert all(b[k][j]==0 for j in range(k+1,len(b))),('zero diagonal with nonzero row',k)
            zero+=1
            continue
        rank+=1;bits=max(bits,pivot.bit_length())
        for i in range(k+1,len(b)):
            for j in range(i,len(b)):
                value,rem=divmod(pivot*b[i][j]-b[i][k]*b[k][j],previous)
                assert rem==0
                b[i][j]=b[j][i]=value
        previous=pivot
    return dict(rank=rank,nullity=zero,largest_pivot_bits=bits)


def elliptic_case(F,S,C):
    E=Curve(F);G=E.points;q=F.q
    size=len(G)
    plus=[[E.index[E.add(P,Q)] for Q in G] for P in G]
    neg=[E.index[E.neg(P)] for P in G]
    fs=[E.f(P) for P in G];rho=[E.rho(P) for P in G]
    two=[i for i in range(size) if plus[i][i]==0]
    four=[i for i in range(size) if plus[i][i] in two]
    assert len(two)==4 and len(four)==16 and size==8*len(C)+16
    fibers=Counter(rho)
    assert set(fibers)==set(C)|{0,1,2,q}
    assert all(fibers[x]==(8 if x in C else 4) for x in fibers)
    assert all(fs[plus[i][t]]==fs[i] for i in range(size) for t in two)
    kernel_entries=0
    for i in range(size):
        for j in range(size):
            x,y=rho[i],rho[j]
            expected=0 if x==y else 1 if q in (x,y) else int(S[x,y])
            assert fs[plus[i][j]]*fs[plus[i][neg[j]]]==expected
            kernel_entries+=1
    ids,reps={},[]
    for i in range(size):
        if i not in ids:
            h=len(reps);reps.append(i)
            for t in two:ids[plus[i][t]]=h
    L=len(reps)
    table=[[ids[plus[i][j]] for j in reps] for i in reps]
    minus=[next(k for k in range(L) if table[h][k]==0) for h in range(L)]
    f=[fs[i] for i in reps]
    assert f[0]==0 and all(abs(v)==1 for v in f[1:])
    special=[h for h in range(L) if table[h][h]==0]
    M=np.array([[f[table[h][k]]*f[table[h][minus[k]]]
                 for k in range(L)] for h in range(L)],dtype=np.int64)
    assert all(M[h,k]==(0 if h==k else 1) for h in special for k in range(L))
    # Complete normalized quotient, checked without irrational square roots:
    # ordinary fibers have size 2, exceptional fibers size 1 on H.
    rhos=[rho[i] for i in reps];labels=C+[0,1,2,q]
    incidence=np.array([[int(x==y) for y in labels] for x in rhos],dtype=np.int64)
    sizes=np.diag([2]*len(C)+[1]*4)
    K=np.array([[0 if x==y else 1 if q in (x,y) else S[x,y]
                 for y in labels] for x in labels],dtype=np.int64)
    assert np.array_equal(incidence.T@incidence,sizes)
    assert np.array_equal(M,incidence@K@incidence.T)
    assert np.array_equal(M@incidence,incidence@K@sizes)
    good=[h for h in range(L) if h not in special]
    U=incidence[np.ix_(good,range(len(C)))]
    Mg=M[np.ix_(good,good)]
    assert np.array_equal(U.T@Mg@U,4*S[np.ix_(C,C)])
    assert np.array_equal(U.T@np.ones((len(good),len(good)),dtype=np.int64)@U,
                          4*np.ones((len(C),len(C)),dtype=np.int64))
    certs=[]
    supports=[(0,), (0,1), (0,1,2), (0,1,2,3)]
    for shifts in supports:
        g=[prod(f[table[h][t]] for t in shifts) for h in range(L)]
        A=np.array([[g[table[h][minus[k]]] for k in range(L)] for h in range(L)],dtype=object)
        cert=psd_certificate(len(shifts)**2*q*np.eye(L,dtype=object)-A@A.T)
        assert sum(g)**2<=len(shifts)**2*q
        certs.append(dict(shifts=shifts,sharp_coefficient=len(shifts),**cert))
    # All repeated four-shift products and M^2 boundary classifications.
    squared=M@M;strata=Counter()
    for u in range(L):
        for v in range(L):
            counts=Counter([u,minus[u],v,minus[v]])
            odd=[t for t,c in counts.items() if c%2]
            even=[t for t,c in counts.items() if c%2==0]
            g=[prod(f[table[h][t]] for t in odd) for h in range(L)]
            expected=sum(g)-sum(g[minus[t]] for t in even)
            assert int(squared[u,v])==expected
            if odd:assert sum(g)**2<=len(odd)**2*q
            else:assert expected==L-len(even)
            strata[f'odd_{len(odd)}_even_{len(even)}']+=1
    return dict(group_order=size,quotient_order=L,kernel_entries=kernel_entries,
                scalar_zero_fourier=sum(f),sharp_twisted_certificates=certs,
                square_entries=L*L,square_strata=dict(strata),
                full_border_and_centering_checked=True)


def field_case(ell,do_elliptic=False,do_powers=False):
    F=QuadraticField(ell);q=F.q
    S=np.array([[F.chi[F.sub(x,y)] for y in range(q)] for x in range(q)],dtype=np.int64)
    J=np.ones((q,q),dtype=np.int64);I=np.eye(q,dtype=np.int64)
    assert np.array_equal(S,S.T) and np.array_equal(S@S,q*I-J)
    assert np.array_equal(S@np.ones(q,dtype=np.int64),np.zeros(q,dtype=np.int64))
    Pnum=q*I-J+ell*S
    assert np.array_equal(Pnum@Pnum,2*q*Pnum)
    assert int(np.trace(Pnum))==q*(q-1)
    assert np.array_equal(S[:ell,:ell],np.ones((ell,ell),dtype=np.int64)-np.eye(ell,dtype=np.int64))
    subfield=np.array([int(x<ell) for x in range(q)],dtype=np.int64)
    assert np.array_equal(S@subfield,ell*subfield-np.ones(q,dtype=np.int64))
    centered_subfield=ell*subfield-np.ones(q,dtype=np.int64)
    assert np.array_equal(S@centered_subfield,ell*centered_subfield)
    assert np.array_equal(Pnum@centered_subfield,2*q*centered_subfield)
    cases=[];C3=None;power_checks=0
    for a in range(1,min(5,ell-1)+1):
        anchors=list(range(a))
        C=[x for x in range(q) if all(S[x,z]==1 for z in anchors)]
        B=list(range(a,ell));m=len(C);n=len(B)
        assert set(B)<=set(C)
        actual=Fraction(int(Pnum[np.ix_(B,B)].sum()),2*q*n)
        expected=1-Fraction(a+2,2*ell)+Fraction(a,2*q)
        assert actual==expected
        b=2**a;N=b-1;beta=1-Fraction(a+2,ell)+Fraction(a,q)
        zsq=Fraction(b*b,4*N)*beta*beta
        row=dict(anchors=a,cell_size=m,subfield_support=n,projection_rayleigh=str(actual),
                 normalized_rayleigh_squared=str(zsq),outside_unit_interval=zsq>1,
                 limiting_squared_rayleigh=str(Fraction(b*b,4*N)))
        if a==3:
            C3=C
            if ell>=29:
                margin=96*(q-5*ell+3)-79*q
                assert margin>0
                row['transfer_threshold_margin']=margin
                row['energy_lower_base']='63/16'
        if do_powers and a in (2,3):
            D=np.array([b-1 if x in C else -1 for x in range(q)],dtype=object)
            W=(ell*S.astype(object)-J.astype(object))*D[np.newaxis,:]
            Z=ell*S[np.ix_(C,C)].astype(object)-np.ones((m,m),dtype=object)
            Hprev=2*np.eye(m,dtype=object);H=b*Z
            power=np.eye(q,dtype=object);rank=(q-1)//2
            for k in range(1,4):
                power=power@W
                if k>1:Hprev,H=H,b*Z@H-q*q*N*Hprev
                residual=q**k*((q-rank-m)+(-1)**k*(rank-m))
                assert int(np.trace(power))==int(np.trace(H))+residual
                power_checks+=1
        cases.append(row)
    result=dict(characteristic=ell,field_order=q,nonresidue=F.nu,
                ambient_square_entries=q*q,ambient_projection_entries=q*q,
                subfield_clique_size=ell,subfield_eigenvector_entries=3*q,
                localizations=cases,power_trace_checks=power_checks)
    if do_elliptic:result['elliptic']=elliptic_case(F,S,C3)
    return result


def main():
    # Unit checks for singular and indefinite pivot cases.
    assert psd_certificate([[1,1],[1,1]])['nullity']==1
    assert psd_certificate([[2,0,1],[0,0,0],[1,0,1]])['nullity']==1
    for A in [[[0,1],[1,0]],[[1,2],[2,1]]]:
        try:psd_certificate(A)
        except AssertionError:pass
        else:raise AssertionError('Indefinite test matrix was accepted')
    cases=[]
    for ell in [7,11,17,19,29]:
        c=field_case(ell,ell in (7,17,19),ell in (7,11))
        cases.append(c)
        print(json.dumps(dict(characteristic=ell,field_order=ell*ell,
                              elliptic_quotient=c.get('elliptic',{}).get('quotient_order'),
                              three_anchor_outlier=c['localizations'][2]['outside_unit_interval'])),flush=True)
    totals=dict(fields=len(cases),ambient_square_entries=sum(c['ambient_square_entries'] for c in cases),
                ambient_projection_entries=sum(c['ambient_projection_entries'] for c in cases),
                localizations=sum(len(c['localizations']) for c in cases),
                subfield_eigenvector_entries=sum(c['subfield_eigenvector_entries'] for c in cases),
                exact_power_traces=sum(c['power_trace_checks'] for c in cases),
                elliptic_kernel_entries=sum(c.get('elliptic',{}).get('kernel_entries',0) for c in cases),
                squared_kernel_entries=sum(c.get('elliptic',{}).get('square_entries',0) for c in cases),
                sharp_twisted_certificates=sum(len(c.get('elliptic',{}).get('sharp_twisted_certificates',[])) for c in cases))
    inputs=['research/parallel16-square-field-barrier-2026-09-05.md',
            'experiments/parallel16_square_fields_2026_09_05.py',
            'research/parallel14-elliptic-model-2026-09-05.md',
            'research/parallel15-elliptic-correlations-2026-09-05.md',
            'research/parallel2-spectral-transfer-2026-09-04.md',
            'sources/katz-gauss-kloosterman-monodromy.pdf',
            'sources/kunisky-2303.16475v1.html']
    out=dict(status='Exact actual square-field Paley examples passed; prime-field conjecture remains unproved.',
             cases=cases,totals=totals,
             arithmetic='Exact quadratic finite fields, integer matrices, rational Rayleigh quotients, fraction-free positive-semidefinite certificates with zero-pivot checks. No floating-point acceptance.',
             input_sha256={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs})
    dest=ROOT/'results/parallel16_square_fields_2026_09_05.json'
    dest.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps(dict(written=str(dest),totals=totals)),flush=True)


if __name__=='__main__':main()
