#!/usr/bin/env python3
"""Plot exact saved round6 results and their explicit limitations."""
from fractions import Fraction as F
from hashlib import sha256
import json
from math import prod,sqrt
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    names=('critical_third/scaling_results.json','critical_third/inventory_free_results.json',
           'trace_envelopes/results.json','overlap_preflight/results.json','pencil_tracks/results.json')
    data={name:json.loads((HERE/name).read_text()) for name in names}
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':11,'axes.spines.top':False,
                         'axes.spines.right':False,'axes.titleweight':'bold','text.color':'#263343',
                         'axes.labelcolor':'#263343','xtick.color':'#263343','ytick.color':'#263343'})
    fig,axs=plt.subplots(2,2,figsize=(14,8.9))
    fig.subplots_adjust(left=.08,right=.97,top=.87,bottom=.13,wspace=.30,hspace=.51)
    fig.suptitle('Round 6 · quantify a limitation, then change the certificate',x=.08,y=.974,ha='left',fontsize=20,weight='bold')
    fig.text(.08,.927,'Exact global-moment bounds and complete finite codeword lists; no exceptional-set or prize theorem.',fontsize=11)
    blue,orange,gray='#37879b','#b66939','#74859c'
    ax=axs[0,0];ns=list(range(9,250))
    proxy=[40*sqrt(5)*prod(range(n-8,n+1))/prod(range(n-5,n+1))**1.5 for n in ns]
    ax.plot(ns,proxy,color=gray,label='Finite-n leading proxy',lw=1.6)
    ax.axhline(40*sqrt(5),color=gray,ls='--',lw=1,label='Limit 40√5')
    old=[r for r in data['critical_third/scaling_results.json']['finite_checks'] if r['q']>10000]
    ax.scatter([r['n'] for r in old],[sqrt(r['q'])*r['actual_standardized_float'] for r in old],s=55,color=blue,zorder=3,label='Exact measured moment')
    new=data['critical_third/inventory_free_results.json']['cases']
    mids=[sqrt(r['q'])*sum(r['standardized_interval_float'])/2 for r in new]
    ax.scatter([r['n'] for r in new],mids,s=60,marker='s',facecolors='white',edgecolors=orange,zorder=3,label='Certified enclosure midpoint')
    ax.set_xscale('log');ax.set_xlim(13,270);ax.set_ylim(0,101)
    ax.set_xlabel('Subset size n (log scale)');ax.set_ylabel('√q × standardized third moment')
    ax.set_title('A  A finite-size trend becomes an exact law',loc='left',pad=12)
    ax.legend(frameon=False,fontsize=9,loc='lower right')

    ax=axs[0,1];records=data['trace_envelopes/results.json']['cases']
    for model,color,label in (('conference_only',blue,'Conference + Hasse support'),
                               ('ordinary_powers_through_six',orange,'Also ordinary moments 0–6')):
        rows=[r for r in records if r['model']==model and r['q']>=1297]
        ax.plot([r['q'] for r in rows],[float(F(*r['width_relative_to_absolute_raw_third_moment'])) for r in rows],marker='o',color=color,label=label)
    ax.set_xscale('log');ax.set_yscale('log');ax.set_ylim(1e-20,2)
    ax.set_xlabel('Prime order q (log scale)');ax.set_ylabel('Certified raw-moment relative interval width')
    ax.set_title('B  Lost arithmetic information has a small effect',loc='left',pad=12)
    ax.legend(frameon=False,fontsize=9,loc='upper right')

    ax=axs[1,0];pencils=data['pencil_tracks/results.json']['large']
    # Read actual symbolic certificate totals instead of hardcoding them.
    certs=[r['result']['certificate'] for r in pencils]
    def total(cert):
        for key in ('total_node_count','complete_node_count','complete_symbolic_node_count','node_count'):
            if key in cert:return cert[key]
        raise KeyError(list(cert))
    counts=[total(cert) for cert in certs]
    bars=ax.bar(range(len(counts)),counts,color=[blue,orange],width=.52)
    ax.bar_label(bars,padding=4,labels=[f'{c:,}' for c in counts]);ax.set_yscale('log');ax.set_ylim(1,1e7)
    ax.set_xticks(range(2),['Two skew tracks','Three colliding tracks'])
    ax.set_ylabel('Complete scalar/codeword nodes (log scale)')
    ax.set_title('C  Exact lists without a codeword anchor',loc='left',pad=12)
    ax.text(.03,.95,'F65,537 · n = 1024 · k = 64 · s = 410',transform=ax.transAxes,fontsize=10)

    ax=axs[1,1];rows=data['overlap_preflight/results.json']['results'];x=np.arange(2);w=.31
    b=ax.bar(x-w/2,[r['old_anchored_bases'] for r in rows],w,color=gray,label='Previous anchored bases')
    c=ax.bar(x+w/2,[r['after_redundancy_deletion'] for r in rows],w,color=orange,label='Selected overlapping tracks')
    ax.bar_label(b,padding=3);ax.bar_label(c,padding=3)
    ax.set_yscale('log');ax.set_ylim(20,10000);ax.set_xticks(x,['F1009 · sample 0','F17 · sample 3'])
    ax.set_ylabel('Verified cover size (log scale)')
    ax.set_title('D  Overlap bypasses the partition-cap barrier',loc='left',pad=12)
    ax.legend(frameon=False,fontsize=9,loc='upper right')
    for ax in axs.flat:ax.grid(axis='y',alpha=.13);ax.set_axisbelow(True)
    fig.text(.08,.045,'A: enclosure widths are too small to display; new q = 6,700,417 and 2,013,265,921 require no character array.\nB: nonnegative histogram relaxations, not constructed graph pairs. D: both toy covers are complete; pure k-support cover size is exponential at fixed rate.',fontsize=10,color='#526173')
    for ext in ('png','pdf'):fig.savefig(HERE/f'overview.{ext}',dpi=180,facecolor='white')
    (HERE/'plot_metadata.json').write_text(json.dumps({'source_sha256':sha256(Path(__file__).read_bytes()).hexdigest(),
         'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in names},
         'scope':'Saved exact results; plot uses rounded real display values and labeled enclosure midpoints.'},indent=2)+'\n')


if __name__=='__main__':main()
