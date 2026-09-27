#!/usr/bin/env python3
"""Exact ambient elliptic convolution, Jacobi identities, and odd eigenvector.

The asymptotic obstruction uses the cited discrepancy theorem separately.
This checker does not prove the desired compressed spectral bound.
"""
from fractions import Fraction
from hashlib import sha256
from math import isqrt
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def prime(p):return p>=2 and all(p%d for d in range(2,isqrt(p)+1))


def primitive_root(p):
    m=p-1;factors=[];d=2
    while d*d<=m:
        if m%d==0:
            factors.append(d)
            while m%d==0:m//=d
        d+=1
    if m>1:factors.append(m)
    return next(g for g in range(2,p) if all(pow(g,(p-1)//d,p)!=1 for d in factors))


def verify(p, matrix=False):
    assert prime(p) and p%4==1
    g=primitive_root(p);N=p-1;q=N//2
    chi=[0]+[1 if pow(a,q,p)==1 else -1 for a in range(1,p)]
    G=[pow(g,j,p) for j in range(N)];Q=G[::2]
    f=[chi[(1-t)%p] for t in G]
    assert all(f[-r%N]==((-1)**r)*f[r] for r in range(N))
    L=[sum(chi[y*(y-1)*(y-t)%p] for y in range(p)) for t in G]
    convolution=[sum(f[j]*f[(r-j)%N] for j in range(N)) for r in range(N)]
    assert L==convolution
    autocorrelation=[sum(f[j]*f[(j+r)%N] for j in range(N)) for r in range(N)]
    assert autocorrelation==[(p if r==0 else 0)-1-(-1)**r for r in range(N)]
    k=L[::2]
    assert all(k[r]==k[-r%q] for r in range(q)) and sum(k)==1
    result=dict(p=p,generator=g,Q_size=q,
                convolution_coefficients_checked=N,
                autocorrelation_coefficients_checked=N,
                constant_mode_eigenvalue=1,all_passed=True)
    if p==89:
        v=[0,1,0,-1]*(q//4)
        assert len(v)==q and all(v[-r%q]==-v[r] for r in range(q))
        Kv=[sum(k[(j-i)%q]*v[j] for j in range(q)) for i in range(q)]
        assert Kv==[73*x for x in v]
        norm=sum(x*x for x in v);form=sum(x*y for x,y in zip(v,Kv))
        assert 3*form>2*p*norm
        outside=[Q[r] for r in range(q) if v[r] and chi[(Q[r]-1)%p]!=1]
        assert outside
        result['witness']=dict(Q_in_generator_order=Q,vector=v,eigenvalue=73,
                               norm_squared=norm,quadratic_form=form,
                               ratio_to_p=str(Fraction(73,p)),inversion_odd=True,
                               nonzero_coordinates_outside_common_neighbors=outside,
                               supported_on_common_neighbors=False)
    if matrix:
        C=[x for x in Q if chi[(x-1)%p]==1]
        m=len(C);index={x:j for j,x in enumerate(C)}
        assert m==(p-5)//4
        R=[index[(1-x)%p] for x in C]
        I=[index[pow(x,-1,p)] for x in C]
        Lmap=dict(zip(G,L))
        K=[[Lmap[y*pow(x,-1,p)%p] for y in C] for x in C]
        S=[[chi[(x-y)%p] for y in C] for x in C]
        for i,x in enumerate(C):
            for j,y in enumerate(C):
                ambient=sum(chi[(x-t)%p]*chi[t]*chi[(t-y)%p] for t in range(p))
                assert K[i][j]==ambient
                average3=K[i][j]+K[R[i]][R[j]]+K[R[I[i]]][R[I[j]]]
                square=sum(S[i][a]*S[a][j] for a in range(m))
                assert 4*square==p*(i==j)+average3-6
        result['common_neighbor_dimension']=m
        result['restricted_identity_entries']=m*m
    return result


def main():
    cases=[verify(p,matrix=p<=113) for p in [13,17,29,37,41,53,61,73,89,97,109,113,193,241,257,641]]
    # Parameter audit for the sourced discrepancy consequence. At
    # p>=2^20, sqrt(p)>=1024 and pi<22/7 make the gap less than1/3.
    gap_bound=32*Fraction(22,7)**2/1024
    assert gap_bound<Fraction(1,3)
    p0=2**20
    assert 2*(p0-3)>2*32  # 2(p-3)p^(-1/4)>2, so remove two quartic modes.
    source=ROOT/'sources/lu-zheng-zheng-1305.3405v3.html'
    out=dict(scope='Ambient quadratic-residue elliptic operator: exact diagonalization ingredients and an inversion-odd eigenvalue73 at p89. No bound on the required compressed and averaged operator.',
             cases=cases,prime_cases=len(cases),
             convolution_coefficients=sum(x['convolution_coefficients_checked'] for x in cases),
             autocorrelation_coefficients=sum(x['autocorrelation_coefficients_checked'] for x in cases),
             restricted_identity_entries=sum(x.get('restricted_identity_entries',0) for x in cases),
             parameter_audit=dict(p_threshold=p0,relative_gap_upper_at_threshold=str(gap_bound),
                                  required_upper_gap='1/3',uses_sourced_discrepancy=True,
                                  asymptotic_claim_is_not_a_finite_test_conclusion=True),
             discrepancy_source_sha256=sha256(source.read_bytes()).hexdigest(),all_passed=True)
    path=ROOT/'results/parallel53_ambient_elliptic_2026_09_06.json'
    path.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({k:out[k] for k in ['prime_cases','convolution_coefficients','autocorrelation_coefficients','restricted_identity_entries','all_passed']}))
    print(json.dumps(next(x['witness'] for x in cases if 'witness' in x)))


if __name__=='__main__':main()
