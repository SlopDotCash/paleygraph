#!/usr/bin/env python3
"""Independent Euler-symbol and literal row-product oracle; no NTT imports."""
from array import array
from collections import Counter,defaultdict
from fractions import Fraction
from itertools import combinations,product
import hashlib,json,math,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def character(p):
    assert p%4==1 and all(p%d for d in range(2,math.isqrt(p)+1))
    return [0 if x==0 else (1 if pow(x,(p-1)//2,p)==1 else -1) for x in range(p)]
def e6(values):
    c=[1,0,0,0,0,0,0]
    for value in values:
        for k in range(6,0,-1):c[k]+=value*c[k-1]
    return c[6]

def expected_types(p,edges,tau):
    e01,e02,e12=edges;zero=[(0,e01,e02),(e01,0,e12),(e02,e12,0)];ret=Counter(zero)
    for a,b,c in product((-1,1),repeat=3):
        terms=(p-3-a*(e01+e02)-b*(e01+e12)-c*(e02+e12)
               -a*b*(1+e02*e12)-a*c*(1+e01*e12)-b*c*(1+e01*e02)+a*b*c*tau)
        assert terms>=0 and terms%8==0
        if terms:ret[a,b,c]+=terms//8
    return ret

def direct_row_moment(p,chi):
    S=[[chi[(x-y)%p] for y in range(p)] for x in range(p)]
    equal=sum(e6(S[x]) for x in range(p))
    double=3*sum(e6(S[a][x]**2*S[b][x] for x in range(p)) for a in range(p) for b in range(p) if a!=b)
    distinct=6*sum(e6(S[a][x]*S[b][x]*S[c][x] for x in range(p)) for a,b,c in combinations(range(p),3))
    return equal,double,distinct

def main():
    rows=[];start=time.monotonic()
    preflight_path=HERE.parent/'preflight/third_moment_preflight.json';preflight=json.loads(preflight_path.read_text())
    for path in sorted(HERE.glob('inventory_p*.json'),key=lambda p:json.loads(p.read_text())['p']):
        began=time.monotonic();data=json.loads(path.read_text());p=data['p'];chi=character(p)
        raw=(HERE/data['provenance']['tau_file']).read_bytes();assert len(raw)==4*(p+1) and int.from_bytes(raw[:4],'little')==p
        tau=array('i');assert tau.itemsize==4;tau.frombytes(raw[4:])
        if sys.byteorder!='little':tau.byteswap()
        assert len(tau)==p and tau[0]==tau[1]==-1 and sum(tau)==0
        histogram=Counter((1,chi[t],chi[t-1],int(tau[t])) for t in range(2,p))
        actual={(tuple(r['edges']),r['tau']):r['count'] for r in data['records']}
        assert actual=={(key[:3],key[3]):count for key,count in histogram.items()}
        power=defaultdict(lambda:[0]*7)
        for (edges,t),count in actual.items():
            for k in range(7):power[edges][k]+=count*t**k
        assert data['edge_power_sums']==[{'edges':list(edges),'powers':values} for edges,values in sorted(power.items())]
        assert data['total_power_sums']==[sum(v[k] for v in power.values()) for k in range(7)]
        assert sum(r['count'] for r in data['records'])==p-2
        assert data['ordered_distinct_row_weight']==p*(p-1)
        selected=range(p) if p<=1297 else [r['t'] for r in data['sampled_tau_rows']]
        terms=0
        for t in selected:
            value=sum(chi[x]*chi[(x-1)%p]*chi[(x-t)%p] for x in range(p));assert value==tau[t];terms+=p
        types=0
        for t in (range(2,p) if p<=101 else [2,p//2,p-1]):
            actual_types=Counter((chi[(-x)%p],chi[(1-x)%p],chi[(t-x)%p]) for x in range(p))
            edges=(chi[p-1],chi[p-t],chi[(1-t)%p]);assert edges==(1,chi[t],chi[t-1])
            assert actual_types==expected_types(p,edges,int(tau[t]));types+=1
        # Complete parameter involution controls, with signed inversion.
        inverse=array('I',[0])*p;inverse[1]=1
        for t in range(2,p):inverse[t]=(p-(p//t)*inverse[p%t]%p)%p
        for t in range(2,p):
            assert tau[(1-t)%p]==tau[t]
            assert tau[inverse[t]]==chi[t]*tau[t]
        # Exact affine/sign-complement checks, including inverses.
        affine=0
        for d in (range(1,p) if p<=101 else [1,2,3,p//2,p-1]):
            inv=pow(d,p-2,p);assert inv*d%p==1
            for x in (range(p) if p<=101 else [0,1,2,3,p//2,p-2,p-1]):
                assert chi[(x*inv)%p]==chi[d]*chi[x];affine+=1
        direct_moment=None
        if p<=29:
            equal,double,distinct=direct_row_moment(p,chi);ref=data['n6_global_third_moment']
            assert [equal,double,distinct]==[ref['all_equal_numerator'],ref['exactly_two_equal_numerator'],ref['all_distinct_numerator']]
            val=Fraction(equal+double+distinct,math.comb(p,6));assert [val.numerator,val.denominator]==ref['third_moment']
            direct_moment={'all_unordered_distinct_row_triples':math.comb(p,3),'all_ordered_two_equal_pairs':p*(p-1),'third_moment':ref['third_moment']}
        if p in [13,17]:
            old=next(r for r in preflight['cases'] if r['graph']==f'Paley{p}')
            assert old['joint_edges_tau_histogram']==data['records']
            assert old['third_moment']==data['n6_global_third_moment']['third_moment']
        prov=data['provenance']
        for key,src in [('input_sha256',prov['input_file']),('tau_sha256',prov['tau_file']),('cpp_sha256','trace_inventory.cpp'),('ntt_header_sha256','ntt_exact.hpp'),('executable_sha256','trace_inventory'),('runner_sha256','run_inventory.py')]:assert prov[key]==sha(HERE/src)
        assert prov['copied_ntt_source_sha256']==prov['ntt_header_sha256']
        row={'p':p,'inventory_file':path.name,'inventory_sha256':sha(path),'all_character_values_verified':p,
            'direct_tau_rows':len(selected),'direct_integer_terms':terms,'full_tau_array_independently_recomputed':p<=1297,
            'literal_triple_type_checks':types,'reflection_and_signed_inversion_parameters_checked':p-2,
            'affine_sign_map_point_checks':affine,'direct_row_third_moment':direct_moment,'seconds':time.monotonic()-began}
        rows.append(row);print(json.dumps(row),flush=True)
    rejected=[]
    # Failing commands must not create the requested output file.
    target=HERE/'rejected_should_not_exist.bin'
    for args in [[],['9',str(target)],['49',str(target)],['15',str(target)],['5242889',str(target)]]:
        r=subprocess.run([str(HERE/'trace_inventory'),*args],capture_output=True,text=True);assert r.returncode!=0 and not r.stdout and not target.exists()
        rejected.append({'args':args,'error':r.stderr.strip()})
    (HERE/'validation.json').write_text(json.dumps({'status':'exact inventory and normalization validation passed; large tau arrays independently sampled',
        'cases':rows,'guard_rejections':rejected,'preflight_sha256':sha(preflight_path),'elapsed_seconds':time.monotonic()-start,'source_sha256':sha(__file__)},indent=2)+'\n')
if __name__=='__main__':main()
