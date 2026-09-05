#!/usr/bin/env python3
"""Run the exact backend; retain input, source, executable and schema provenance."""
import argparse,hashlib,json,platform,resource,subprocess,sys,time
from pathlib import Path
HERE=Path(__file__).resolve().parent

def sha(path):return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def run(p,marks):
    stem=f'p{p}_marks_'+('_'.join(map(str,marks)) if marks else 'none')
    inp={'p':p,'marks':marks,'matrix':'S[x,y]=LegendreSymbol(x-y,p)','ordered_row_pairs':True}
    path=HERE/f'input_{stem}.json';path.write_text(json.dumps(inp,sort_keys=True,indent=2)+'\n')
    start=time.monotonic();proc=subprocess.run([str(HERE/'marked_counts'),str(p),*map(str,marks)],capture_output=True,text=True,check=True)
    out=json.loads(proc.stdout)
    out['provenance']={'input_sha256':sha(path),'input_file':path.name,'backend_source_sha256':sha(HERE/'marked_counts.cpp'),
        'backend_executable_sha256':sha(HERE/'marked_counts'),'runner_source_sha256':sha(__file__),
        'wall_seconds':time.monotonic()-start,
        'children_peak_rss_bytes':resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss*(1 if sys.platform=='darwin' else 1024),
        'children_peak_rss_scope':'maximum over this runner process lifetime, not isolated per case',
        'python_executable':sys.executable,'python_version':platform.python_version()}
    target=HERE/f'counts_{stem}.json';target.write_text(json.dumps(out,indent=2)+'\n')
    print(json.dumps({'file':target.name,'p':p,'marks':marks,'cells':len(out['cells']),'relations':len(out['relations']),
        'wall_seconds':out['provenance']['wall_seconds'],'ntt_transforms':out['metadata']['ntt_transforms']}),flush=True)
    return target

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--p',type=int);parser.add_argument('--marks',type=int,nargs='*',default=[0,1,2]);parser.add_argument('--million',action='store_true');args=parser.parse_args()
    if args.p:run(args.p,args.marks);return
    for p in [13,17,101,1297,65537]:run(p,[0,1,2])
    # Extra finite inputs let the conditional compiler compare arithmetic marks.
    for marks in [[],[0],[1],[0,1],[0,2],[0,1,3],[0,1,4],[1,2,3],[0,4,8],[0,2,4]]:run(101,marks)
    if args.million:run(1000033,[0,1,2])

if __name__=='__main__':main()
