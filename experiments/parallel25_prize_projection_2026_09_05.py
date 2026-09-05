#!/usr/bin/env python3
"""Exact projection, interleaving and extension-field checks for pass25."""
from collections import Counter
from fractions import Fraction
from hashlib import sha256
from itertools import combinations, product
from math import comb
from pathlib import Path
from random import Random
import json

ROOT = Path(__file__).resolve().parents[1]
COUNTS = Counter()
RNG = Random(2026090525)

def check(ok, name):
    assert ok, name
    COUNTS[name] += 1

def scalar_code(p, n, k):
    return [tuple(sum(c[j]*pow(x,j,p) for j in range(k)) % p for x in range(n))
            for c in product(range(p), repeat=k)]

def exact_model(p, n, k, r, s):
    code = scalar_code(p, n, k)
    subsets = [T for size in range(s, n+1) for T in combinations(range(n), size)]
    restrictions = [{tuple(c[x] for x in T) for c in code} for T in subsets]
    words = list(product(range(p), repeat=n*r))
    flags = {}
    for W in words:
        flags[W] = sum(1 << i for i, T in enumerate(subsets)
                       if all(tuple(W[row*n+x] for x in T) in restrictions[i] for row in range(r)))
    return words, flags

def bad_fast(W, Z, p, flags):
    return {g for g in range(p)
            if flags[tuple((w+g*z) % p for w,z in zip(W,Z))] & ~flags[Z]}

def exhaustive_mca():
    records = []
    for p,n,k,s in ((2,3,1,2), (3,2,1,2)):
        maxima = []
        for r in (1,2):
            words, flags = exact_model(p,n,k,r,s)
            counts = [len(bad_fast(W,Z,p,flags)) for W in words for Z in words]
            maxima.append(max(counts))
            COUNTS['exhaustive_MCA_word_pairs'] += len(counts)
        check(maxima[0] == maxima[1], 'exact_MCA_interleaving_maxima')
        records.append({'p':p,'n':n,'k':k,'agreement':s,'scalar_max':maxima[0],'two_row_max':maxima[1],
                        'code_family':'Linear repetition code; domain labels need not be distinct field elements.'})
    # Endpoint B_1=q-1, with an actual short RS code.
    p,n,k,s = 5,4,2,3
    words, flags = exact_model(p,n,k,1,s)
    K = max(len(bad_fast(W,Z,p,flags)) for W in words for Z in words)
    COUNTS['exhaustive_MCA_word_pairs'] += len(words)**2
    check(K == p-1, 'nontrivial_q_minus_one_MCA_endpoint')
    # Restriction masks factor rowwise, so no large two-row table is needed.
    def rflags(W):
        a,b = W[:n],W[n:]
        return flags[a] & flags[b]
    projections = [a for a in product(range(p), repeat=2) if any(a)]
    project = lambda a,W: tuple((a[0]*W[x]+a[1]*W[n+x]) % p for x in range(n))
    max_seen = 0
    for _ in range(80):
        W,Z = [tuple(RNG.randrange(p) for _ in range(2*n)) for _ in range(2)]
        bad = {g for g in range(p)
               if rflags(tuple((w+g*z)%p for w,z in zip(W,Z))) & ~rflags(Z)}
        survived = Counter()
        for a in projections:
            bp = bad_fast(project(a,W),project(a,Z),p,flags)
            check(len(bp) <= K, 'actual_projected_MCA_scalar_bound')
            survived.update(bad & bp)
        for g in bad:
            check(survived[g] >= p*(p-1), 'nonzero_projection_witness_survival')
        check(len(bad) <= K, 'sampled_MCA_interleaving_bound')
        max_seen = max(max_seen,len(bad))
    records.append({'p':p,'n':n,'k':k,'agreement':s,'scalar_max':K,
                    'two_row_scope':'80 sampled pairs, every nonzero projection; not an exhaustive two-row maximum.',
                    'sample_two_row_max':max_seen})
    return records

def interp_line(row, a, b, p):
    slope = (row[b]-row[a])*pow((b-a)%p,-1,p) % p
    intercept = (row[a]-slope*a) % p
    return tuple((intercept+slope*x)%p for x in range(len(row)))

def list_lines(rows,p,s):
    n=len(rows[0])
    candidates = {tuple(interp_line(row,a,b,p) for row in rows) for a,b in combinations(range(n),2)}
    return [F for F in candidates if sum(all(F[j][x]==rows[j][x] for j in range(len(rows)))
                                       for x in range(n)) >= s]

def list_projection_checks():
    records=[]
    for p,r in ((5,2),(29,2),(29,3)):
        n=4;s=2
        rows=[tuple(pow(x,j+2,p) for x in range(n)) for j in range(r)]
        ls=list_lines(rows,p,s)
        check(len(ls)==6, 'six_distinct_interleaved_explanations')
        # The scalar maximum is exactly six: at most one distinct degree<2
        # interpolant per two-subset, with six attained by the quadratic word.
        scalar=list_lines([tuple(x*x%p for x in range(n))],p,s)
        check(len(scalar)==comb(n,2)==6,'scalar_list_maximum_certificate')
        nonzero=[a for a in product(range(p),repeat=r) if any(a)]
        collisions=0; injective=0; largest=0
        for a in nonzero:
            projected=[tuple(sum(a[j]*F[j][x] for j in range(r))%p for x in range(n)) for F in ls]
            Wp=tuple(sum(a[j]*rows[j][x] for j in range(r))%p for x in range(n))
            available={F[0] for F in list_lines([Wp],p,s)}
            check(set(projected)<=available, 'projected_list_agreement_preservation')
            col=sum(v*(v-1)//2 for v in Counter(projected).values())
            collisions+=col;injective+=col==0;largest=max(largest,len(set(projected)))
        check(collisions <= comb(len(ls),2)*(p**(r-1)-1),'projected_pair_collision_count')
        if comb(6+1,2)<p:
            check(injective>0,'large_alphabet_injective_projection')
        if 6<p:
            check(len(ls)*(p-6)<=6*(p-1),'dimension_free_list_inequality')
        records.append({'p':p,'rows':r,'list_size':len(ls),'nonzero_projections':len(nonzero),
                        'injective_projections':injective,'colliding_pairs_total':collisions,
                        'largest_distinct_projected_list':largest})
    return records

class QuadraticField:
    def __init__(self,p,nu):
        self.p=p;self.nu=nu;self.q=p*p
        assert pow(nu,(p-1)//2,p)==p-1
    def add(self,a,b):
        return (a%self.p+b%self.p)%self.p+self.p*((a//self.p+b//self.p)%self.p)
    def neg(self,a):
        return (-a%self.p)%self.p+self.p*((- (a//self.p))%self.p)
    def sub(self,a,b): return self.add(a,self.neg(b))
    def mul(self,a,b):
        p=self.p
        x,y=a%p,a//p;u,v=b%p,b//p
        return (x*u+self.nu*y*v)%p+p*((x*v+y*u)%p)
    def inv(self,a):
        p=self.p;x,y=a%p,a//p
        d=pow((x*x-self.nu*y*y)%p,-1,p)
        return x*d%p+p*((-y*d)%p)
    def div(self,a,b): return self.mul(a,self.inv(b))

def residual(row,T,F):
    a,b=T[:2]
    slope=F.div(F.sub(row[b],row[a]),(b-a)%F.p)
    return tuple(F.sub(row[x],F.add(row[a],F.mul(slope,(x-a)%F.p))) for x in T[2:])

def bad_extension(W,Z,F,s):
    n=len(W)
    subsets=[T for k in range(s,n+1) for T in combinations(range(n),k)]
    out=set()
    for T in subsets:
        rw,rz=residual(W,T,F),residual(Z,T,F)
        if not any(rz):continue
        j=next(j for j,z in enumerate(rz) if z)
        g=F.div(F.neg(rw[j]),rz[j])
        if all(F.add(w,F.mul(g,z))==0 for w,z in zip(rw,rz)):out.add(g)
    # Independent direct folded-word check, without residual ratios.
    direct=set()
    for g in range(F.q):
        folded=tuple(F.add(w,F.mul(g,z)) for w,z in zip(W,Z))
        if any(any(residual(Z,T,F)) and not any(residual(folded,T,F)) for T in subsets):direct.add(g)
    check(out==direct,'extension_remainder_ratio_vs_direct_fold')
    return out

def extension_checks():
    F=QuadraticField(5,2);p=F.p
    for a in range(1,F.q):check(F.mul(a,F.inv(a))==1,'quadratic_field_inverse')
    lines={tuple(sorted({F.add(alpha,F.mul(beta,t)) for t in range(p)}))
           for alpha in range(F.q) for beta in range(1,F.q)}
    check(len(lines)==p*(p+1),'all_affine_base_field_lines')
    # The base-field scalar MCA maximum at n4,k2,s3 was exhaustively four.
    for _ in range(50):
        W,Z=[tuple(RNG.randrange(F.q) for _ in range(4)) for _ in range(2)]
        bad=bad_extension(W,Z,F,3)
        for line in lines:
            check(len(set(line)&bad)<=4,'extension_bad_set_affine_line_bound')
        # Expand one arbitrary line and compare exact bad parameters, using
        # both base-field coefficient rows. Field multiplication need not
        # commute with a scalar coordinate projection.
        alpha,beta=RNG.randrange(F.q),RNG.randrange(1,F.q)
        U=tuple(F.add(w,F.mul(alpha,z)) for w,z in zip(W,Z))
        V=tuple(F.mul(beta,z) for z in Z)
        descended=set()
        for t in range(p):
            folded=tuple(F.add(u,F.mul(t,v)) for u,v in zip(U,V))
            for size in (3,4):
                for T in combinations(range(4),size):
                    rv=residual(V,T,F);rf=residual(folded,T,F)
                    vrows=[tuple((v//(p**j))%p for v in rv) for j in (0,1)]
                    frows=[tuple((v//(p**j))%p for v in rf) for j in (0,1)]
                    if any(any(row) for row in vrows) and not any(any(row) for row in frows):descended.add(t)
        check({F.add(alpha,F.mul(beta,t)) for t in descended}==bad & {F.add(alpha,F.mul(beta,t)) for t in range(p)},
              'exact_affine_line_descent')
    return {'field':'F_5[u]/(u^2-2)','received_pairs':50,'base_field_lines':len(lines),
            'scope':'Actual RS n4,k2,s3; all challenges and witnesses for every sampled pair.'}

def official_check():
    p=2130706433;n=262144;k=131072;q=p**6;R=q//2**128
    check(R==274980728111395087,'official_combined_integer_budget')
    check(comb(R+1,2)<q,'official_large_alphabet_list_hypothesis')
    A=p-18
    bound=A*(p-1)//(p-A)
    check(bound==252217214619120071<R,'conditional_base_list_projection_example')
    # Pure integer endpoint in the interleaving theorem, no asymptotic floor.
    for q0 in (2,3,5,9,25):
        for r in (1,2,3,4):
            alpha=Fraction(q0**(r-1)-1,q0**r-1)
            for K in range(q0):
                check(Fraction(K,1)/(1-alpha)<K+1,'MCA_exact_integer_rounding')
    return {'p':p,'q':q,'n':n,'k':k,'budget':R,
            'large_alphabet_slack':2*q-R*(R+1),
            'conditional_base_list_A':A,'conditional_interleaved_upper':bound,
            'scope':'Pinned contract arithmetic only; neither scalar maximum nor spot-check obligation proved.'}

def main():
    d={'status':'passed','mca':exhaustive_mca(),'lists':list_projection_checks(),
       'extension':extension_checks(),'official':official_check()}
    inputs=['research/parallel25-prize-projection-2026-09-05.md',
            'experiments/parallel25_prize_projection_2026_09_05.py',
            'research/parallel23-prize-bridge-2026-09-05.md',
            'results/parallel23_prize_bridge_2026_09_05.json',
            'research/parallel-prize-subfield-2026-09-04.md',
            'sources/official-prize-2026-09-04/ProximityPrize/Benchmark/IRSProfile.lean',
            'sources/official-prize-2026-09-04/ProximityPrize/Benchmark/TargetLower.lean',
            'sources/official-prize-2026-09-04/dependencies/ArkLib/ProofSystem/ToyProblem/Impl/IRS.lean',
            'sources/gopalan-guruswami-raghavendra-0811.4395-abstract.html']
    d['input_sha256']={p:sha256((ROOT/p).read_bytes()).hexdigest() for p in inputs}
    d['counts']=dict(COUNTS)
    d['limitations']='Projection identities and exact conditional reduction; no production maxima, full Paley/prize proof, or formal certification.'
    out=ROOT/'results/parallel25_prize_projection_2026_09_05.json'
    out.write_text(json.dumps(d,indent=2)+'\n')
    print(json.dumps({'status':d['status'],'counts':d['counts'],'output':str(out)},indent=2))

if __name__=='__main__':main()
