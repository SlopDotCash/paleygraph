#!/usr/bin/env python3
"""Export a static, source-grounded research figure, not an asymptotic claim."""
import json
from pathlib import Path
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

ROOT=Path(__file__).resolve().parent


def read(p):return json.loads((ROOT/p).read_text())


def main():
    plt.rcParams.update({'font.family':'DejaVu Sans','font.size':10,
                         'axes.spines.top':False,'axes.spines.right':False})
    fig,axs=plt.subplots(2,2,figsize=(12.8,9))
    fig.subplots_adjust(top=.86,bottom=.11,hspace=.51,wspace=.30)
    fig.suptitle('Tools that expose missing information',x=.065,y=.97,ha='left',fontsize=22,fontweight='bold')
    fig.text(.065,.915,'Exact finite experiments and numerical discovery · September 5, 2026',fontsize=12,color='#555555')

    ax=axs[0,0]
    counts=[63,1,0]
    ax.bar(range(3),counts,color=['#6c7e95','#dc8e32','#288d78'],width=.57)
    for i,y in enumerate(counts):ax.text(i,y+1.8,str(y),ha='center',fontweight='bold',fontsize=13)
    ax.set_ylim(0,73);ax.set_xticks(range(3),['Moments + boundary\n+ energy + quartic deck','Retain rooted\nquartic incidence','Also retain\ntriangle contraction'])
    ax.set_ylabel('Fibers with different sixth moments')
    ax.set_title('A   Refine the measurements',loc='left',fontweight='bold',pad=15)
    ax.text(0,-.31,'All 455,126 normalized six-sets in F₆₁.\nZero collisions here is not a universal determination theorem.',transform=ax.transAxes,fontsize=9,color='#555555')

    ax=axs[0,1]
    for degree,label,color,marker in [(1,'Prime fields','#356cb3','o'),(2,'Square fields','#288d78','s')]:
        rows=sorted([r for r in read('spectral/results.json')['cases'] if r['degree']==degree],key=lambda r:r['q'])
        ax.plot([r['q'] for r in rows],[r['tail_subspace']['max_abs_overlap'] for r in rows],marker=marker,color=color,label=label)
        ax.plot([r['q'] for r in rows],[r['same_spectrum_relabelled_tail']['max_abs_overlap'] for r in rows],linestyle='--',color=color,alpha=.5)
    ax.set_xscale('log');ax.set_ylim(0,1);ax.set_xlabel('Field size q (log scale)');ax.set_ylabel('Maximum near-edge translation overlap')
    ax.legend(frameon=False,fontsize=9,loc='center right')
    ax.set_title('B   Keep arithmetic coordinates',loc='left',fontweight='bold',pad=15)
    ax.text(0,-.31,'Solid: original operator. Dashed: same-spectrum relabelling.\nFloating-point discovery; integer witnesses replayed separately.',transform=ax.transAxes,fontsize=9,color='#555555')

    ax=axs[1,0]
    profile=read('proximity/deformation_results.json')['profiles'][0]
    for key,label,color,style in [('translation','Translation','#288d78','-'),('tangent_basis_0','Tangent path','#dc8e32','--'),('single_root','Single root','#6c7e95',':')]:
        ax.plot(range(4),profile['jet_survival'][key],marker='o',color=color,linestyle=style,label=label)
    ax.set_xticks(range(4));ax.set_yticks([0,1]);ax.set_ylim(-.12,1.2);ax.set_xlabel('Order of truncated root motion');ax.set_ylabel('Surviving initial syzygy dimension')
    ax.legend(frameon=False,fontsize=9,loc='center right')
    ax.set_title('C   First order can miss failure',loc='left',fontweight='bold',pad=15)
    ax.text(0,-.31,'Known F₄₁ witness. All 24 tested tangent paths die at order two.\nAmbient straight paths; subgroup admissibility is tested separately.',transform=ax.transAxes,fontsize=9,color='#555555')

    ax=axs[1,1]
    records=read('proximity/stack_results.json')['records']
    for i,p in enumerate([17,1009,2147483713]):
        rows=[r for r in records if r['p']==p]
        values=[r['finite_bad_scalar_count'] for r in rows]
        ax.scatter(i+np.array([-.12,-.04,.04,.12]),values,color=['#356cb3','#288d78','#6c7e95'][i],s=45)
        for x,y in zip(i+np.array([-.12,-.04,.04,.12]),values):
            if p==17:ax.text(x,y+.25,str(y),ha='center',fontsize=9)
    ax.set_xlim(-.4,2.4);ax.set_ylim(0,15);ax.set_xticks(range(3),['17','1,009','2,147,483,713']);ax.set_xlabel('Field characteristic p');ax.set_ylabel('Exact close-scalar count per sampled pencil')
    ax.set_title('D   Realize the coding configuration',loc='left',fontweight='bold',pad=15)
    ax.text(0,-.31,'Same coset construction; n=16, k=8, agreement threshold=11.\nFour seeded pencils per field. No worst-case or prize count claimed.',transform=ax.transAxes,fontsize=9,color='#555555')
    fig.savefig(ROOT/'overview.png',dpi=180,facecolor='white')
    fig.savefig(ROOT/'overview.pdf',facecolor='white')
    print(ROOT/'overview.png')


if __name__=='__main__':main()
