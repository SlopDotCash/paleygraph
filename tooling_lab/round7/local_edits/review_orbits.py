#!/usr/bin/env python3
"""Independent Burnside and representative-feature check of the saturation audit."""
from collections import defaultdict
from hashlib import sha256
import json
from pathlib import Path
from review_results import literal_signs, elementary, target

HERE=Path(__file__).resolve().parent


def canonical(c,q,multipliers):
    return min(tuple(sorted((a*x+b)%q for x in c)) for a in multipliers for b in range(q))


def burnside(q,n,multipliers):
    fixed_sum=0
    for a in multipliers:
        for b in range(q):
            seen=set();cycles=[]
            for start in range(q):
                if start in seen:continue
                x=start;length=0
                while x not in seen:seen.add(x);length+=1;x=(a*x+b)%q
                cycles.append(length)
            coefficients=[1]+[0]*n
            for length in cycles:
                for j in range(n,length-1,-1):coefficients[j]+=coefficients[j-length]
            fixed_sum+=coefficients[n]
    assert fixed_sum%(q*len(multipliers))==0
    return fixed_sum//(q*len(multipliers))


def feature(S,c):
    q=len(S);n=len(c);hist=defaultdict(int);summaries=[]
    for row in S:
        signs=[row[a] for a in c];hist[(signs.count(0),signs.count(1))]+=1
    for a in c:
        r=[elementary([row[b] for b in c if b!=a],5) for row in S]
        internal=[sum(S[x][b]*r[x] for x in range(q)) for b in c]
        summaries.append((sum(r),sum(x*x for x in r),internal[c.index(a)],sum(internal)))
    return (tuple((*key,value) for key,value in sorted(hist.items())),tuple(sorted(summaries)))


def main():
    data=json.loads((HERE/'orbit_results.json').read_text());rows=[]
    for name,value in data['source_sha256'].items():assert sha256((HERE/name).read_bytes()).hexdigest()==value
    for case in data['cases']:
        q,n=case['q'],case['n'];squares={x*x%q for x in range(1,q)};S=literal_signs(q)
        reps=[tuple(x) for x in case['canonical_orbit_representatives']]
        assert len(reps)==len(set(reps))==burnside(q,n,squares)==case['square_affine_orbits']
        groups=defaultdict(set);full={}
        for c in reps:
            assert c==canonical(c,q,squares)
            groups[feature(S,c)].add(c);full[c]=canonical(c,q,range(1,q))
        assert len(set(full.values()))==burnside(q,n,range(1,q))==case['full_affine_orbits']
        counts=defaultdict(int);full_counts=defaultdict(int);non_affine_pairs=[]
        for group in groups.values():counts[len(group)]+=1;full_counts[len({full[c] for c in group})]+=1
        for group in groups.values():
            if len({full[c] for c in group})<=1:continue
            profiles=[]
            for c in sorted(group):
                histogram=defaultdict(int)
                for a in c:
                    for b in set(range(q))-set(c):histogram[target(S,sorted((set(c)-{a})|{b}),6)]+=1
                profiles.append({'selected':list(c),'neighbour_histogram':[[v,k] for v,k in sorted(histogram.items())],
                                 'raw_sums_through_four':[sum(k*v**j for v,k in histogram.items()) for j in range(5)]})
            non_affine_pairs.append({'profiles':profiles,'all_neighbour_histograms_equal':all(x['neighbour_histogram']==profiles[0]['neighbour_histogram'] for x in profiles)})
        assert len(groups)==case['deletion_summary_fibres']
        assert dict(counts)=={int(k):v for k,v in case['fibres_by_number_of_distinct_orbits'].items()}
        assert dict(full_counts)=={int(k):v for k,v in case['fibres_by_number_of_full_affine_orbits'].items()}
        rows.append({'q':q,'n':n,'square_affine_orbits':len(reps),'full_affine_orbits':len(set(full.values())),
                     'summary_fibres':len(groups),'summary_fibres_by_full_affine_orbits':dict(full_counts),
                     'non_affine_pairs':non_affine_pairs})
    out={'status':'passed','scope':'Independent all-translation canonicalization, Burnside cycle counting and literal derivative profiles on every square-affine representative; no production import.',
         'cases':rows,'source_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('review_orbits.py','review_results.py')},
         'input_sha256':{'orbit_results.json':sha256((HERE/'orbit_results.json').read_bytes()).hexdigest()}}
    (HERE/'orbit_verification.json').write_text(json.dumps(out,indent=2)+'\n');print(json.dumps(rows),flush=True)


if __name__=='__main__':main()
