#!/usr/bin/env python3
"""Standalone research figure from exact round4 result artifacts."""
from fractions import Fraction
import json
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent


def main():
    data = json.loads((HERE/'marked_moments/contraction_envelopes.json').read_text())['cases']
    fibers = json.loads((HERE/'scalar_fibers/results.json').read_text())['inherited_fixtures']
    aligned = next(x for x in fibers if x['name'] == 'aligned-n48-k6-s28')
    interleaved = next(x for x in fibers if x['name'] == 'interleaved-n48-k6-s28')
    operations = [aligned['previous_track_interpolations'], interleaved['previous_track_interpolations'],
                  interleaved['certificate']['static_direction_certificate']['distinct_piece_polynomials']]
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False,
                         'axes.titlesize':12,'axes.labelsize':10})
    fig, axes = plt.subplots(1,3,figsize=(13.0,4.8), gridspec_kw={'width_ratios':[1,1.15,1]})
    fig.suptitle('Retain the missing information, then measure what it can change',
                 fontsize=16,x=.065,ha='left',y=.99)
    ax=axes[0]
    corrections=[float(Fraction(*x['Q_coefficient'])*x['Q']) for x in data[:2]]
    ax.bar([0,1],corrections,color=['#16697a','#dc6c3e'],width=.56)
    ax.axhline(0,color='#444444',linewidth=.8)
    ax.set_xticks([0,1],['(1, 4, 0)\nQ = −42','(1, 9, 0)\nQ = 86'])
    ax.set_ylabel('Second moment minus the common cell-only term')
    ax.set_title('A. Real inputs with equal cell data',loc='left',pad=16)
    ax.set_ylim(-.031,.024)
    ax.text(.5,.91,'Paley29 · n = 7, d = 6\nExact gap = 256 / 7475',
            transform=ax.transAxes,ha='center',va='top',fontsize=10)
    ax.set_xlabel('Ordered marked triple')

    ax=axes[1]
    rows=[x for x in data if (x['q'],x['n']) in [(1297,8),(65537,16),(1000033,31)]]
    radii=[x['relative_radius_upper_float'] for x in rows]
    ax.plot(range(3),radii,'o-',color='#16697a',linewidth=2,markersize=7)
    ax.set_yscale('log');ax.set_ylim(1e-16,1e-6)
    ax.set_xticks(range(3),['q = 1,297\nn = 8*','q = 65,537\nn = 16','q = 1,000,033\nn = 31'])
    ax.set_ylabel('Certified correction radius / |cell-only term|')
    ax.set_title('B. The extra correction is small',loc='left',pad=16)
    ax.grid(axis='y',alpha=.2)
    for i,v in enumerate(radii):
        ax.annotate(f'{v:.2e}',(i,v),xytext=((5,-16) if i == 0 else (0,12)),
                    textcoords='offset points',ha=('left' if i == 0 else 'center'),fontsize=9)
    ax.set_xlabel('Marks (0, 1, 2); *off critical size')

    ax=axes[2]
    labels=['Aligned\ntrack batching','Interleaved\ntrack batching','Interleaved\npiece certificate']
    ax.bar(range(3),operations,color=['#16697a','#dc6c3e','#427a50'],width=.58)
    ax.set_yscale('log');ax.set_ylim(1,12000)
    ax.set_xticks(range(3),labels,fontsize=9)
    ax.set_ylabel('Tracks interpolated / supplied pieces checked')
    ax.set_title('C. A failed batching rule was refined',loc='left',pad=16)
    for i,v in enumerate(operations):ax.text(i,v*1.18,str(v),ha='center',va='bottom')
    ax.text(.5,.95,'n = 48, k = 6, s = 28\nSame single decoded node',
            transform=ax.transAxes,ha='center',va='top',fontsize=10)
    ax.set_xlabel('Different operations; not a timing comparison')
    fig.subplots_adjust(left=.07,right=.99,top=.82,bottom=.28,wspace=.48)
    fig.text(.07,.07,'Exact finite diagnostics. Panel B bounds only the relation correction to a conditional average.\n'
             'Panel C verifies supplied piece witnesses; it does not discover them or bound arbitrary inputs.',fontsize=10,color='#444444')
    fig.savefig(HERE/'overview.png',dpi=180,bbox_inches='tight')
    fig.savefig(HERE/'overview.pdf',bbox_inches='tight')
    plt.close(fig)


if __name__=='__main__':main()
