#!/usr/bin/env python3
"""Attach big-integer power sums, the n6 third moment, and exact provenance."""
from collections import defaultdict
from fractions import Fraction
import argparse,hashlib,json,math,platform,resource,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(p):return hashlib.sha256(Path(p).read_bytes()).hexdigest()
def choose(n,k):return math.comb(n,k) if 0<=k<=n else 0

def k6(m,tau):
    assert abs(tau)<=m and (m+tau)%2==0
    positive,negative=(m+tau)//2,(m-tau)//2
    direct=sum((-1)**j*choose(negative,j)*choose(positive,6-j) for j in range(7))
    numerator=tau**6-(15*m-40)*tau**4+(45*m*m-210*m+184)*tau*tau-15*m*(m-2)*(m-4)
    assert numerator==720*direct
    return direct

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--p',type=int,nargs='+',default=[13,17,29,101,1297,65537,1000033]);args=parser.parse_args()
    rows=[]
    for p in args.p:
        source=HERE/f'input_p{p}.json';source.write_text(json.dumps({'p':p,'family':'prime Paley','degree':6,'parameter_normalization':'ordered rows (0,1,t), t!=0,1','trace_sign':'tau=-FrobeniusTrace'},indent=2)+'\n')
        binary=HERE/f'tau_p{p}.bin';start=time.monotonic()
        out=json.loads(subprocess.check_output([str(HERE/'trace_inventory'),str(p),str(binary)],text=True))
        sums=defaultdict(lambda:[0]*7)
        for record in out['records']:
            edges=tuple(record['edges']);tau=record['tau'];count=record['count']
            for k in range(7):sums[edges][k]+=count*tau**k
        assert sum(v[0] for v in sums.values())==p-2 and sum(v[1] for v in sums.values())==2
        out['edge_power_sums']=[{'edges':list(edges),'powers':values} for edges,values in sorted(sums.items())]
        powers=[sum(v[k] for v in sums.values()) for k in range(7)];out['total_power_sums']=powers
        all_equal=-p*choose((p-1)//2,3);two_equal=-3*p*(p-1)*choose((p-3)//2,3)
        distinct=p*(p-1)*sum(r['count']*k6(p-3,r['tau']) for r in out['records'])
        numerator=all_equal+two_equal+distinct;third=Fraction(numerator,choose(p,6))
        # Independent aggregation of the same Krawtchouk formula via power sums.
        m=p-3;num=powers[6]-(15*m-40)*powers[4]+(45*m*m-210*m+184)*powers[2]-15*m*(m-2)*(m-4)*powers[0]
        assert num%720==0 and p*(p-1)*(num//720)==distinct
        out['n6_global_third_moment']={'n':6,'degree':6,'all_equal_numerator':all_equal,'exactly_two_equal_numerator':two_equal,
            'all_distinct_numerator':distinct,'sum_T6_cubed':numerator,'number_of_six_sets':choose(p,6),'third_moment':[third.numerator,third.denominator],
            'third_moment_display':float(third),'six_column_subsets_enumerated':0}
        out['normalization']={'family':'prime Paley only','row_order':['0','1','t'],
            'edge_order':['S[0,1]','S[0,t]','S[1,t]'],'edge_formula':['1','chi(t)','chi(t-1)'],
            'trace_definition':'tau(t)=sum_x chi(x(x-1)(x-t))',
            'frobenius_trace_sign':'a_t(p)=-tau(t)','excluded_parameters':[0,1],
            'affine_map':'x -> (x-a)/(b-a); t=(c-a)/(b-a)',
            'sign_complement_multiplier':'chi(b-a)',
            'valid_target_symmetry':'Degree6 row coefficients and their third products are invariant under simultaneous sign complement.',
            'distinct_ordered_triples_per_parameter':p*(p-1),
            'weights_applied_to_records':False}
        out['provenance']={'input_file':source.name,'input_sha256':sha(source),'tau_file':binary.name,'tau_sha256':sha(binary),
            'tau_format':'little-endian uint32 p, followed by p little-endian signed int32 tau values including t0,1',
            'cpp_sha256':sha(HERE/'trace_inventory.cpp'),'ntt_header_sha256':sha(HERE/'ntt_exact.hpp'),
            'copied_ntt_source':'round4/marked_counts_ablation/ntt_exact.hpp',
            'copied_ntt_source_sha256':sha(HERE.parents[1]/'round4/marked_counts_ablation/ntt_exact.hpp'),
            'executable_sha256':sha(HERE/'trace_inventory'),'runner_sha256':sha(__file__),
            'wall_seconds':time.monotonic()-start,
            'maximum_child_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
            'rss_scope':'maximum over child processes in this runner lifetime, not isolated per input',
            'runtime':{'python':sys.executable,'version':platform.python_version()}}
        target=HERE/f'inventory_p{p}.json';target.write_text(json.dumps(out,indent=2)+'\n')
        row={'p':p,'file':target.name,'bins':len(out['records']),'third_moment':out['n6_global_third_moment']['third_moment'],
             'backend_seconds':out['metadata']['backend_seconds'],'wall_seconds':out['provenance']['wall_seconds']}
        rows.append(row);print(json.dumps(row),flush=True)
    (HERE/'run_index.json').write_text(json.dumps({'cases':rows,'runner_sha256':sha(__file__)},indent=2)+'\n')
if __name__=='__main__':main()
