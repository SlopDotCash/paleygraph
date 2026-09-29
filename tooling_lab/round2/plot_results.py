#!/usr/bin/env python3
"""Standalone research figure generated only from saved exact/numeric results."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def read(path):return json.loads((HERE/path).read_text())


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.titlesize':12,
                         'axes.spines.top':False,'axes.spines.right':False,'figure.facecolor':'#fafaf7',
                         'axes.facecolor':'#fafaf7','savefig.facecolor':'#fafaf7'})
    fig,axs=plt.subplots(2,2,figsize=(12.5,9.0))
    fig.suptitle('Mathematical tooling lab · second iteration',fontsize=19,x=.06,ha='left',y=.97)
    fig.text(.06,.927,'Exact diagnostics, complete finite censuses, and a scaled counterexample to an overstrong explanation.',fontsize=11,color='#4c555c')
    blue='#285c8e';orange='#be612f';green='#36826b'
    specs=read('observables/review_exchange_results.json')['spectra']
    x=np.arange(len(specs));low=[r['lower_energy_ratio'][0]/r['lower_energy_ratio'][1] for r in specs]
    ax=axs[0,0];ax.bar(x,low,color=[blue]*3+[green]*3,width=.62)
    ax.set_yscale('log');ax.set_ylim(1e-6,1);ax.set_xticks(x,[f"{r['p']:,}\nn={r['n']}" for r in specs])
    ax.set_ylabel('Share of centered L² energy in levels 1–5');ax.set_xlabel('Field size p · exact all-input computation')
    ax.set_title('A  Sixth-order variation is mostly at level six',loc='left',pad=12)
    ax.grid(axis='y',alpha=.15);ax.set_axisbelow(True)
    rows=read('observables/cohort_review.json')['cases'];ax=axs[0,1];x=np.arange(len(rows))
    for j,(key,label,color) in enumerate([('baseline','Baseline',blue),('actual','True incidence',orange),('shuffled','Shuffled incidence','#89969e')]):
        ys=[r['model_comparison'][key]['random_holdout']['r2_against_cohort_mean'] for r in rows]
        ax.bar(x+(j-1)*.23,ys,width=.21,color=color,label=label)
    ax.axhline(0,color='#56616a',linewidth=.8);ax.set_xticks(x,[f"{r['p']:,}\nn={r['n']}" for r in rows])
    ax.set_ylabel('Affine-orbit-disjoint holdout R²');ax.set_xlabel('Toy field at left; three critical-size cases at right')
    ax.set_title('B  Toy prediction gain fails to transfer',loc='left',pad=12);ax.legend(fontsize=8,frameon=False)
    results=read('proximity/results.json')['scale_experiments'];ax=axs[1,0];x=np.arange(len(results))
    ax.bar(x-.18,[float(r['full_subset_count_not_enumerated']) for r in results],width=.35,color='#becad3',label='Full agreement subsets')
    ax.bar(x+.18,[r['result']['cover_certificate']['basis_count'] for r in results],width=.35,color=blue,label='Certified interpolation bases')
    ax.set_yscale('log');ax.set_ylim(1,1e38);ax.set_xticks(x,[f"n={r['result']['n']}\nk={r['result']['k']}, s={r['result']['s']}" for r in results]);ax.set_ylabel('Combinatorial candidates · logarithmic scale')
    ax.set_title('C  Complete per-stack census becomes practical',loc='left',pad=12);ax.legend(fontsize=8,frameon=False)
    profiles=read('novelty/arithmetic_lift_results.json')['profiles'];ax=axs[1,1]
    ax.plot(range(1,5),profiles[0]['hensel_survival'],'o-',color=orange,label='F41 example · lawful lift',linewidth=2)
    ax.plot(range(1,5),profiles[1]['hensel_survival'],'s-',color=green,label='Coset controls · lawful lift',linewidth=2)
    ax.plot([1,2],[1,0],'x--',color='#88949b',label='All controls · naive residues',markersize=9)
    ax.set_xticks(range(1,5),['p','p²','p³','p⁴']);ax.set_yticks([0,1]);ax.set_ylim(-.12,1.25)
    ax.set_xlabel('Arithmetic precision');ax.set_ylabel('Surviving initial kernel dimension')
    ax.set_title('D  Preserving the domain changes the diagnosis',loc='left',pad=12);ax.legend(fontsize=8,frameon=False,loc='center right',bbox_to_anchor=(1,.57))
    fig.text(.06,.035,'Separate spectral result: p=1,000,033, exact Rayleigh >0.867; every nonzero translation overlap <0.011210.\nThese are finite diagnostics. They establish neither conjecture, historical novelty, nor a worst-case bound.',fontsize=10,color='#4c555c')
    fig.subplots_adjust(left=.075,right=.975,bottom=.135,top=.855,hspace=.55,wspace=.3)
    fig.savefig(HERE/'overview.png',dpi=190)
    fig.savefig(HERE/'overview.pdf')
    print('Saved round2 overview.png and overview.pdf')


if __name__=='__main__':main()
