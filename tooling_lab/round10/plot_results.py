#!/usr/bin/env python3
from fractions import Fraction as F
from hashlib import sha256
import json
from math import sqrt
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent


def main():
    scales=json.loads((HERE/'projection_scale.json').read_text())['cases']
    toy=json.loads((HERE/'projection_toy.json').read_text())['cases']
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,ax=plt.subplots(1,2,figsize=(12,4.5),layout='constrained')
    for family,color in [('arithmetic_progression','#286c9d'),('seeded_uniform','#b44936')]:
        rows=[r for r in scales if r['family']==family and r['remainder_relative_to_high_norm_squared'] is not None]
        ax[0].loglog([r['q'] for r in rows],[sqrt(float(F(*r['remainder_relative_to_high_norm_squared']))) for r in rows],'o-',color=color,label=family.replace('_',' '))
    ax[0].set_xlabel('Prime q');ax[0].set_ylabel('Relative Euclidean norm of the remainder')
    ax[0].set_title('Forward means mostly reproduce inward T6 incidence')
    ax[0].grid(alpha=.2);ax[0].legend(frameon=False)
    x=np.arange(3)
    for i,(key,label,color) in enumerate([('full_affine_orbits','Full affine orbits','#777777'),('mean_deck_fibres','Distinct mean collections','#286c9d'),('projected_norm_fibres','Distinct projected norms','#b44936')]):
        counts=[r[key] for r in toy];bars=ax[1].bar(x+(i-1)*.25,counts,width=.24,label=label,color=color)
        ax[1].bar_label(bars,padding=2,fontsize=9)
    ax[1].set_xticks(x,[f'n = {r["n"]}' for r in toy]);ax[1].set_ylabel('Count at q = 17');ax[1].set_ylim(0,120)
    ax[1].set_title('At n = 8, apparent separation identifies every orbit')
    ax[1].legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('Round 10 — a faster compiler and an exact information audit',fontsize=14)
    for ext in ('png','pdf'):fig.savefig(HERE/f'overview.{ext}',dpi=180)
    plt.close(fig)
    out={'source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{name:sha256((HERE/name).read_bytes()).hexdigest() for name in ('projection_scale.json','projection_toy.json')},
         'scope':'Exact fractions converted to floating point only for plotting; chosen-input diagnostics, not uniform limits.'}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
