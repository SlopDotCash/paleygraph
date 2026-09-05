#!/usr/bin/env python3
"""Exact list-to-winning-set checks and a finite IRS unsafe-radius bound.

This script does not prove Paley, does not build a Lean proof, and does not
submit a benchmark entry. All acceptance conditions are integer comparisons.
"""
from collections import Counter,defaultdict
from fractions import Fraction
from hashlib import sha256
from itertools import combinations,product
from math import comb
from pathlib import Path
import json

ROOT=Path(__file__).resolve().parents[1]


def poly_roots(roots,p):
    out=[1]
    for x in roots:
        following=[0]*(len(out)+1)
        for i,c in enumerate(out):
            following[i]=(following[i]-x*c)%p
            following[i+1]=(following[i+1]+c)%p
        out=following
    return out


def evaluate(coeffs,x,p):
    out=0
    for c in reversed(coeffs):
        out=(out*x+c)%p
    return out


def dot(a,b,p):
    return sum(x*y for x,y in zip(a,b))%p


def finite_case(p,domain,k,r):
    n=len(domain)
    assert k<=r<n<p and len(set(domain))==n
    families=defaultdict(list)
    for roots in combinations(domain,r):
        poly=poly_roots(roots,p)
        key=tuple(poly[k:r])
        families[key].append((roots,tuple((-x)%p for x in poly[:k])))
    key=max(families,key=lambda x:(len(families[x]),x))
    family=families[key]
    minimum=(comb(n,r)+p**(r-k)-1)//p**(r-k)
    assert len(family)>=minimum
    received_poly=[0]*k+list(key)+[1]
    received=tuple(evaluate(received_poly,x,p) for x in domain)
    listed=[]
    minimum_weight=n
    for message in product(range(p),repeat=k):
        encoded=tuple(evaluate(message,x,p) for x in domain)
        if any(message):
            minimum_weight=min(minimum_weight,sum(x!=0 for x in encoded))
        if sum(a!=b for a,b in zip(encoded,received))<=n-r:
            listed.append(message)
    assert minimum_weight==n-k+1 and n-r<minimum_weight
    assert set(listed)=={a for roots,a in family}
    L=len(listed)
    best_v=None
    best_image=set()
    total_collisions=0
    for v in product(range(p),repeat=k):
        counts=Counter(dot(a,v,p) for a in listed)
        total_collisions+=sum(c*c for c in counts.values())
        if len(counts)>len(best_image):
            best_v=v
            best_image=set(counts)
    assert total_collisions==p**k*L+p**(k-1)*L*(L-1)
    assert len(best_image)*(p+L-1)>=p*L
    if comb(L,2)<p:
        assert len(best_image)==L
    # Independently enumerate all output messages and winning challenges.
    winning=set()
    least_invalid_second_weight=n+1
    for message in product(range(p),repeat=k):
        encoded=tuple(evaluate(message,x,p) for x in domain)
        value=dot(message,best_v,p)
        if value==1:
            least_invalid_second_weight=min(least_invalid_second_weight,sum(x!=0 for x in encoded))
        if sum(a!=b for a,b in zip(encoded,received))<=n-r:
            winning.add(value)
    assert winning==best_image and least_invalid_second_weight>=minimum_weight
    return {'p':p,'domain':domain,'n':n,'k':k,'r':r,
            'coefficient_labels':list(key),'received_polynomial_ascending':received_poly,
            'pigeonhole_list_lower_bound':minimum,'exact_list_size':L,
            'message_vectors':listed,'v':list(best_v),'mu1':0,'mu2':1,
            'winning_challenges':sorted(winning),'winning_count':len(winning),
            'hash_bound_numerator':L,'hash_bound_denominator':p+L-1,
            'minimum_codeword_weight':minimum_weight,
            'least_second_word_weight_for_claim_one':least_invalid_second_weight,
            'whole_admissible_suffix_by_minimum_distance':True,
            'projection_vectors_exhausted':p**k,
            'collision_total':total_collisions,
            'large_list_regime':L>=p}


def int_digest(x):
    return {'bits':x.bit_length(),
            'sha256_big_endian':sha256(x.to_bytes((x.bit_length()+7)//8 or 1,'big')).hexdigest()}


def profile():
    p=2130706433
    q=p**6
    n=2**18
    k=2**17
    required=q//2**128+1
    assert required*2**128>q and (required-1)*2**128<=q
    assert comb(required,2)<q
    # Walk exact binomial ratios from the central coefficient. There is
    # one initial math.comb call, then checked exact divisions.
    r=k
    numerator=comb(n,r)
    denominator=1
    last=None
    while numerator>(required-1)*denominator:
        last=(r,numerator,denominator)
        next_numerator=numerator*(n-r)
        assert next_numerator%(r+1)==0
        numerator=next_numerator//(r+1)
        denominator*=p
        r+=1
    best_r,best_numerator,best_denominator=last
    lower=(best_numerator+best_denominator-1)//best_denominator
    next_lower=(numerator+denominator-1)//denominator
    assert lower>=required and next_lower<required
    assert comb(lower,2)<q
    assert best_numerator==comb(n,best_r)
    assert numerator==comb(n,best_r+1)
    assert best_denominator==p**(best_r-k)
    unsafe=n-best_r
    assert 0<unsafe<n-k+1
    # Exact centibit score: raise the rational inequality to the 100th
    # power and use integers, with no transcendental acceptance condition.
    score_numerator=best_r**12800
    score_denominator=n**12800
    centibits=max(0,score_denominator.bit_length()-score_numerator.bit_length())
    while (score_numerator<<centibits)<score_denominator:
        centibits+=1
    while centibits and (score_numerator<<(centibits-1))>=score_denominator:
        centibits-=1
    assert (score_numerator<<centibits)>=score_denominator
    assert centibits==0 or (score_numerator<<(centibits-1))<score_denominator
    # Existing root-lift lower bound upgrades to a near-one winning bound
    # under the general scalar-list lemma, without any L<q requirement.
    old=json.loads((ROOT/'results/subset_sum_list_certificate.json').read_text())
    previous=old['official_profile']
    assert (previous['p'],previous['q'],previous['n'],previous['scalar_dimension'])==(p,q,n,k)
    selected=next(x for x in previous['records'] if x['a']==32)
    bound=selected['subset_bound']
    assert (bound['n'],bound['s'],bound['t'],bound['M_upper'])==(8192,4097,4095,2691)
    assert bound['lower_log2_floor']==8154
    old_L=(comb(8192,4097)-(p-1)*comb(6785,4095)+p-1)//p
    assert old_L>=2**8154
    assert q<2**186
    near_one={'radius':'4095/8192','list_at_least_power_of_two':8154,
              'field_less_than_power_of_two':186,
              'winning_density_strictly_above':'1 - 2^(-7968)'}
    # A simple, explicit common-sum antipodal family from July ABF C.7.
    antipodal=comb(63,32)
    assert antipodal>=required
    return {'p':p,'q':q,'domain_size':n,'scalar_dimension':k,'interleaving':8,
            'total_dimension':2**20,'minimum_list_size_to_beat_epsilon':required,
            'selected_list_size':lower,
            'injectivity_union_bound_pairs':comb(lower,2),
            'injectivity_union_bound_less_than_q':True,
            'best_agreement_count_in_pigeonhole_family':best_r,
            'unsafe_index':unsafe,'unsafe_radius':str(Fraction(unsafe,n)),
            'minimum_distance_numerator':n-k+1,'minimum_distance_denominator':n,
            'coefficient_constraints':best_r-k,
            'guaranteed_list_size':lower,'next_agreement_guaranteed_list_size':next_lower,
            'best_binomial':int_digest(best_numerator),'best_field_power':int_digest(best_denominator),
            'centibits':centibits,'score_verified_by_power_12800':True,
            'centibits_minimal_for_this_radius':True,
            'scope':'Optimal only within this scalar coefficient-pigeonhole bound at fixed profile',
            'antipodal_family':{'image_order':128,'subset_size':65,'fixed_sum':1,
                                'list_size':antipodal,'radius':'63/128'},
            'earlier_large_list_consequence':near_one}


def main():
    finite=[]
    for args in [(5,[1,2,3,4],2,2),(7,list(range(1,7)),3,3),
                 (13,list(range(1,7)),3,4),(17,sorted({pow(2,j,17) for j in range(8)}),4,5)]:
        finite.append(finite_case(*args))
        print(f'finite p={args[0]}, n={len(args[1])}: exact lists and winning sets passed',flush=True)
    result=profile()
    paths=['experiments/list_to_winning_set.py','results/subset_sum_list_certificate.json',
           'sources/abf-2026-680-july.pdf',
           'sources/prize-attack-2026-09-04/SoundnessBounds.lean',
           'sources/prize-attack-2026-09-04/SimplifiedIOR.lean',
           'sources/official-prize-2026-09-04/ProximityPrize/Benchmark/TargetUpper.lean',
           'sources/official-prize-2026-09-04/ProximityPrize/Benchmark/IRSProfile.lean']
    output={'status':'Ordinary theorem and exact finite checks; not Lean formalized or submitted',
            'finite_cases':finite,'profile':result,
            'source_sha256':{name:sha256((ROOT/name).read_bytes()).hexdigest() for name in paths}}
    (ROOT/'results/list_to_winning_set.json').write_text(json.dumps(output,indent=2)+'\n')
    print(json.dumps(result,indent=2),flush=True)


if __name__=='__main__':
    main()
