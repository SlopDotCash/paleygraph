#!/usr/bin/env python3
"""Publication-style static overview from frozen exact arithmetic results."""
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
    twins=json.loads((HERE/'twin_analysis.json').read_text())['pairs'][0]['members']
    scales=json.loads((HERE/'scale_results.json').read_text())['cases']
    coding=json.loads((HERE/'coding_block_results.json').read_text())
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12,8),layout='constrained')
    colors=['#286c9d','#b44936']
    values=[x for x,n in twins[0]['histogram']];x=np.arange(len(values))
    for i,m in enumerate(twins):axs[0,0].bar(x+(i-.5)*.34,[n for v,n in m['histogram']],width=.34,color=colors[i],label=f'Set {i+1}')
    axs[0,0].set_xticks(x,values);axs[0,0].set_xlabel('T6 after two swaps');axs[0,0].set_ylabel('Final sets out of 945')
    axs[0,0].set_title('Actual second neighbourhoods differ');axs[0,0].legend(frameon=False)
    for i,m in enumerate(twins):axs[0,1].plot(range(1,22),[float(F(*v)) for v in m['fixed_deletion_mean_deck']],'.-',color=colors[i],label=f'Set {i+1}')
    axs[0,1].set_xlabel('Rank among 21 deletion pairs');axs[0,1].set_ylabel('Mean T6 over 45 insertion pairs')
    axs[0,1].set_title('Conditioning reveals a difference in the mean')
    for family,color in zip(('arithmetic_progression','seeded_uniform'),colors):
        rows=[r for r in scales if r['family']==family and r['statistics']['q']>1297]
        qs=[r['statistics']['q'] for r in rows]
        radius=[sqrt(float(F(*r['Q_free_variance_radius_squared'])/F(*r['variance'])**2)) for r in rows]
        actual=[abs(float(F(*r['Q_variance_correction'])/F(*r['variance']))) for r in rows]
        axs[1,0].loglog(qs,radius,'o-',color=color,label=family.replace('_',' ')+' bound')
        axs[1,0].loglog(qs,actual,'x--',color=color,alpha=.8)
    axs[1,0].set_xlabel('Prime q');axs[1,0].set_ylabel('Size relative to exact variance')
    axs[1,0].set_title('The global Q correction is small on these inputs')
    axs[1,0].legend(fontsize=8,frameon=False);axs[1,0].grid(alpha=.2)
    axs[1,0].text(.02,.03,'Solid: Cauchy radius; dashed: actual |4Q/N|',transform=axs[1,0].transAxes,fontsize=8)
    maxima=[coding['max_all_basepoints'],coding['max_disjoint_basepoints']]
    axs[1,1].bar(['Overlapping blocks','Disjoint blocks'],maxima,color=colors,width=.6)
    axs[1,1].axhline(1.5,color='#555555',ls='--',label='Displayed all-basepoints bound: 1.5')
    axs[1,1].set_ylim(0,2.5);axs[1,1].set_ylabel('Maximum vanishing blocks')
    axs[1,1].set_title('Reusing a root changes the coding count')
    axs[1,1].legend(fontsize=8,frameon=False)
    axs[1,1].text(.04,.05,'All 5,220 monic degree < 4 polynomials over F17\nBlock length 2, multiplier 3',transform=axs[1,1].transAxes,fontsize=9,bbox={'facecolor':'white','alpha':.9,'edgecolor':'none'})
    fig.suptitle('Round 9 — actual two-insertion outcomes and block-incidence checks',fontsize=15)
    for suffix in ('png','pdf'):fig.savefig(HERE/f'overview.{suffix}',dpi=180)
    plt.close(fig)
    out={'source_sha256':{'plot_results.py':sha256(Path(__file__).read_bytes()).hexdigest()},
         'input_sha256':{n:sha256((HERE/n).read_bytes()).hexdigest() for n in ('twin_analysis.json','scale_results.json','coding_block_results.json')},
         'scope':'Plot only: exact rational results converted to floating point for display.'}
    (HERE/'plot_metadata.json').write_text(json.dumps(out,indent=2)+'\n')


if __name__=='__main__':main()
