"""Six proof-linked panels, drawn at journal width, including grayscale exports."""
from pathlib import Path
from fractions import Fraction as F
import csv,json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager
ROOT=Path(__file__).resolve().parents[1];DATA=ROOT/'results'/'revision_analysis';OUT=ROOT/'figures'
OUT.mkdir(parents=True,exist_ok=True)
rows=list(csv.DictReader((DATA/'all_candidates.csv').open()))
FONT='Arial' if any(f.name=='Arial' for f in font_manager.fontManager.ttflist) else 'DejaVu Sans'
plt.rcParams.update({'font.family':FONT,'font.size':9,'font.weight':'normal',
 'axes.labelweight':'normal','axes.titleweight':'normal','axes.linewidth':1.2,
 'axes.grid':False,'xtick.major.width':1,'ytick.major.width':1,
 'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','svg.hashsalt':'soft-east-figures-v14','mathtext.default':'regular'})
def select(k,method,N=None,s=None):
    return sorted([x for x in rows if F(x['kappa'])==F(k) and x['method']==method
          and (s is None or F(x['s'])==F(s)) and (N is None or int(x['N'] or 0)==N)],key=lambda x:int(x['R']))
def style(ax,letter=None):
    ax.grid(False,which='both');ax.tick_params(direction='out',length=3,pad=2,labelsize=8)
    for sp in ax.spines.values():sp.set_visible(True);sp.set_linewidth(1.2)
    if letter:ax.text(-.20,1.10,letter,transform=ax.transAxes,fontsize=12,fontweight='bold',va='bottom')
def export(fig,stem,gray):
    for ext in ['pdf','svg','png']:
        metadata={'CreationDate':None,'ModDate':None} if ext=='pdf' else {'Date':None} if ext=='svg' else {}
        fig.savefig(OUT/(stem+('_grayscale' if gray else '')+'.'+ext),dpi=400,facecolor='white',metadata=metadata)
    plt.close(fig)
for gray in [False,True]:
    per,var=('#111111','#707070') if gray else ('#006DA4','#9B5428')
    fig=plt.figure(figsize=(5.125,5.2));gs=fig.add_gridspec(2,2,left=.13,right=.985,bottom=.10,top=.87,wspace=.36,hspace=.52)
    handles=[]
    for i,(k,s) in enumerate([('1','.02'),('.1','.01'),('.01','.001')]):
        ax=fig.add_subplot(gs[i//2,i%2]);style(ax,'ABC'[i]);ax.set_yscale('log')
        for method,color,marker,ls,face in [('perron',per,'o','-','white'),('variational',var,'v','--',var)]:
            rr=select(k,method,18 if method=='perron' else None,s)
            h,=ax.plot([int(x['R']) for x in rr],[float(F(x['exact_width'])) for x in rr],
                       color=color,marker=marker,mfc=face,mew=1.1,ms=4,lw=1.35,ls=ls,
                       label='Perron, N=18' if method=='perron' else 'Variational')
            if i==0:handles.append(h)
        const=json.loads((ROOT/'results'/'reference_sensitivity'/f'constant_k{float(k):g}_s{float(s):g}_verified.json').read_text())
        h=ax.axhline(float(F(const['exact_width'])),color='#222222',lw=1.35,ls=':',label='Constant')
        if i==0:handles.append(h)
        ax.set_xticks(range(3,10));ax.set_xlabel('Range R',labelpad=3);ax.set_ylabel('Certified width',labelpad=3)
        ax.set_title(f'κ={float(k):g}, s={float(s):g}',fontsize=9,pad=6)
    fig.legend(handles=handles,loc='upper center',bbox_to_anchor=(.55,.985),ncol=3,frameon=False,fontsize=9,
               columnspacing=.9,handlelength=2.1,handletextpad=.45)
    sub=gs[1,1].subgridspec(2,1,height_ratios=[1.0,1.0],hspace=.63)
    ax=fig.add_subplot(sub[0]);style(ax,'D')
    rr=list(csv.DictReader((DATA/'reference_zoom.csv').open()))
    ax.plot([int(x['N']) for x in rr],[float(F(x['exact_width']))*1e3 for x in rr],
            color=per,marker='o',mfc='white',ms=4,mew=1.1,lw=1.35,label='Perron, R=7')
    ax.set_title('κ=0.01, s=0.001',fontsize=9,pad=6)
    ax.text(.97,.95,'Perron, R=7',transform=ax.transAxes,ha='right',va='top',fontsize=8.5)
    ax.set_xticks([14,16,18,20]);ax.set_ylabel('W × 10³',labelpad=3)
    ax=fig.add_subplot(sub[1]);style(ax)
    rr=[x for x in rr if int(x['N'])>=17]
    ax.plot([int(x['N']) for x in rr],[float(F(x['delta_width_from_N18_exact']))*1e6 for x in rr],
            color=per,marker='s',mfc='white',ms=4,mew=1.1,lw=1.35)
    ax.axhline(0,color='#777777',lw=1.05,ls=':');ax.set_xticks([17,18,19,20]);ax.set_xlabel('Reference size N',labelpad=3)
    ax.set_ylabel('ΔW × 10⁶',labelpad=3)
    export(fig,'width_reference_sensitivity',gray)
    fig,axes=plt.subplots(1,2,figsize=(5.125,2.65));fig.subplots_adjust(left=.13,right=.985,bottom=.21,top=.87,wspace=.37)
    ax=axes[0];style(ax,'A');ax.set_yscale('log')
    tails=json.loads((DATA/'spatial_tail_bounds.json').read_text())
    ax.plot([x['K'] for x in tails],[float(x['flip_tail_upper']['upper']) for x in tails],color=per,lw=1.4,label='Certified tail bound')
    ax.set_xlabel('Tail diameter K');ax.set_ylabel('Flip-tail upper bound');ax.set_xticks([10,30,50,70])
    ax.set_title('Fixed R=7 candidate',fontsize=9,pad=6)
    ax.legend(loc='upper right',frameon=False,fontsize=8.5,handlelength=1.4,handletextpad=.3)
    ax=axes[1];style(ax,'B');ax.set_yscale('log')
    bounds=json.loads((DATA/'reference_bounds.json').read_text())
    for key,label,color,marker,ls,face in [('prior','A priori',var,'s','--','white'),('posterior','A posteriori',per,'^','-',per)]:
        ax.plot([x['R'] for x in bounds],[float(x[key]['upper']) for x in bounds],color=color,
                marker=marker,ls=ls,mfc=face,mew=1.1,ms=4,lw=1.35,label=label)
    ax.set_xlabel('Range R');ax.set_ylabel('Reference error upper bound');ax.set_xticks(range(3,10))
    ax.set_title('Fixed reference N=18',fontsize=9,pad=6)
    ax.legend(loc='lower left',bbox_to_anchor=(.00,.025),frameon=False,fontsize=8.5,handlelength=1.4,handletextpad=.3)
    export(fig,'practical_overlap',gray)
print('Two vector figures and direct grayscale variants generated from exact records.')
