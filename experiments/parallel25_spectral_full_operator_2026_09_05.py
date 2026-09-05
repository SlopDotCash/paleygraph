#!/usr/bin/env python3
"""Exact full-sector decomposition and elliptic identities in actual prime fields."""
from collections import Counter
from fractions import Fraction as F
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]
COUNTS=Counter()

def check(ok, name):
    assert ok, name
    COUNTS[name]+=1

def prime(p):
    return p>=2 and all(p%d for d in range(2,isqrt(p)+1))

def add(*args):
    return [[sum(A[i][j] for A in args) for j in range(len(args[0]))] for i in range(len(args[0]))]

def scale(c,A): return [[c*x for x in row] for row in A]
def mul(A,B):
    Bt=list(zip(*B))
    return [[sum(a*b for a,b in zip(row,col)) for col in Bt] for row in A]
def tr(A): return sum(A[i][i] for i in range(len(A)))
def transpose(A): return [list(x) for x in zip(*A)]
def conjugate(A,perm): return [[A[perm[i]][perm[j]] for j in range(len(A))] for i in range(len(A))]
def trace_product(A,B): return sum(A[i][j]*B[j][i] for i in range(len(A)) for j in range(len(A)))

def field(p):
    ch=[0]+[1 if pow(t,(p-1)//2,p)==1 else -1 for t in range(1,p)]
    C=[x for x in range(2,p) if ch[x]==ch[x-1]==1]; m=len(C); idx={x:i for i,x in enumerate(C)}
    E=[[int(i==j) for j in range(m)] for i in range(m)]; J=[[1]*m for _ in C]
    S=[[ch[(x-y)%p] for y in C] for x in C]
    rp=[idx[(1-x)%p] for x in C]; ip=[idx[pow(x,-1,p)] for x in C]
    R=[[int(i==rp[j]) for j in range(m)] for i in range(m)]
    Inv=[[int(i==ip[j]) for j in range(m)] for i in range(m)]
    U=mul(R,Inv); U2=mul(U,U)
    check(m==(p-5)//4 and mul(U2,U)==E,'actual_S3_action')
    check(conjugate(S,rp)==S and conjugate(S,ip)==S,'actual_sign_matrix_equivariance')
    Acyc=add(E,U,U2); Dst=add(scale(2,E),scale(-1,U),scale(-1,U2))
    Ns=[mul(add(E,R),Acyc),mul(add(E,scale(-1,R)),Acyc),
        mul(add(E,scale(-1,R)),Dst),mul(add(E,R),Dst)]
    check(add(*Ns)==scale(6,E),'complete_projection_sum')
    for i,N in enumerate(Ns):
        check(transpose(N)==N and mul(N,N)==scale(6,N),'rational_orthogonal_projection')
        for M in Ns[:i]:check(mul(N,M)==scale(0,E),'different_sector_orthogonality')
    V=add(U,scale(-1,U2))
    check(mul(transpose(V),V)==Dst and add(mul(V,R),mul(R,V))==scale(0,E),
          'standard_intertwiner')
    eps2=int(ch[2]==1);eps3=int(p%3==1)
    g=F(m-3*eps2-2*eps3,6)
    dims=[F(tr(N),6) for N in Ns]
    check(g.denominator==1 and dims==[g+eps2+eps3,g+eps3,2*g+eps2,2*g+eps2],
          'all_sector_dimensions')
    unseen=set(C); sizes=Counter()
    while unseen:
        x=min(unseen)
        orbit={x,(1-x)%p,pow(x,-1,p),(1-pow(x,-1,p))%p,
               pow((1-x)%p,-1,p),x*pow((x-1)%p,-1,p)%p}
        check(orbit<=set(C),'actual_orbit_membership')
        unseen-=orbit;sizes[len(orbit)]+=1
    check(sizes[2]==eps3 and sizes[3]==eps2 and sizes[6]==g,'exceptional_orbit_counts')
    L=[sum(ch[y*(y-1)*(y-t)%p] for y in range(p)) for t in range(p)]
    K0=[[L[y*pow(x,-1,p)%p] for y in C] for x in C]
    K1=conjugate(K0,rp)
    Kthird=conjugate(K1,ip)
    A3=add(K0,K1,Kthird)
    S2=mul(S,S)
    check(scale(4,S2)==add(scale(p,E),A3,scale(-6,J)),'exact_full_elliptic_square_identity')
    check(conjugate(A3,rp)==A3 and conjugate(A3,ip)==A3,'elliptic_conjugacy_average_equivariance')
    if p<=53:
        for i,x in enumerate(C):
            for j,y in enumerate(C):
                direct0=sum(ch[(x-t)%p]*ch[t]*ch[(t-y)%p] for t in range(p))
                direct01=sum(ch[(x-t)%p]*ch[t]*ch[(t-1)%p]*ch[(t-y)%p] for t in range(p))
                check(direct0==K0[i][j],'direct_elliptic_kernel')
                check(direct01==Kthird[i][j]-1,'inversion_boundary_minus_one')
    # Rational part of 4 p^2 [I-(S/sqrt(p)-J/p)^2].
    direct=add(scale(4*p*p,E),scale(-4*p,S2),scale(-4*m,J))
    derived=add(scale(3*p*p,E),scale(-p,A3),scale(6*p-4*m,J))
    check(direct==derived,'exact_full_leakage_identity')
    rows=[sum(row) for row in S];ell=F(sum(L[x] for x in C),m)
    a0=F(sum(rows),m);gv=[F(L[x])-ell for x in C]
    for i in range(m):
        for j in range(m):
            check(rows[i]+rows[j]==2*a0+(gv[i]+gv[j])/4,'retained_rank_two_border_entries')
    if p<=101:
        power=E;vec=[1]*m
        for k in range(1,7):
            power=mul(power,S)
            ts=[F(trace_product(N,power),6) for N in Ns]
            check(ts[2]==ts[3] and tr(power)==ts[0]+ts[1]+2*ts[2],
                  'complete_sector_power_traces')
            vec=[sum(S[i][j]*vec[j] for j in range(m)) for i in range(m)]
            check([sum(Ns[0][i][j]*vec[j] for j in range(m)) for i in range(m)]==[6*x for x in vec],
                  'uniform_seed_Krylov_in_trivial_sector')
    if p in (13,17):check(S==add(E,scale(-1,J)),'nontrivial_sector_extremum_small_example')
    return {'p':p,'m':m,'orbit_counts':dict(sizes),
            'block_dimensions':dict(zip(('trivial','sign','standard_one_copy'),map(int,dims[:3]))),
            'uniform_seed_omitted_dimension_at_least':m-int(dims[0])}

def main():
    cases=[field(p) for p in range(13,258,4) if prime(p)]
    inputs=['research/parallel25-spectral-full-operator-2026-09-05.md',
            'experiments/parallel25_spectral_full_operator_2026_09_05.py',
            'research/parallel24-spectral-operator-2026-09-05.md',
            'research/parallel6-seeded-kernels-2026-09-04.md',
            'research/parallel13-principal-budget-2026-09-05.md',
            'research/parallel2-spectral-transfer-2026-09-04.md']
    report={'status':'PASS','arithmetic':'Python integers and exact Fraction; no floating point.',
            'scope':'Full-sector algebra and exact elliptic reduction; no sector norm saving or full edge.',
            'counts':dict(COUNTS),'fields':cases,
            'input_sha256':{p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs},
            'limitations':['Finite identities do not supply the all-vector elliptic estimate A_sigma <= (2/3+o(1))p.',
                           'No growing-depth character correlation estimate is proved.']}
    out=ROOT/'results/parallel25_spectral_full_operator_2026_09_05.json'
    out.write_text(json.dumps(report,indent=2)+'\n')
    print(json.dumps({'status':'PASS','counts':dict(COUNTS),'fields':len(cases)},indent=2))

if __name__=='__main__':main()
