"""Plot actual frozen widths and measured verification times; no fitted curve."""
from pathlib import Path
from fractions import Fraction as F
import argparse,json,hashlib,statistics
import numpy as np
import matplotlib
matplotlib.use('Agg')
from matplotlib import pyplot as plt,font_manager,ticker
ROOT=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--data',type=Path,required=True);a.add_argument('--output',type=Path,required=True);args=a.parse_args()
rows=json.loads(args.data.read_text(encoding='utf-8'));family='Arial' if any(f.name=='Arial' for f in font_manager.fontManager.ttflist) else 'DejaVu Sans';font=font_manager.findfont(family,fallback_to_default=False)
plt.rcParams.update({'font.family':family,'font.size':9,'axes.labelsize':9.5,'axes.labelweight':'bold','axes.linewidth':1.0,'xtick.labelsize':8.5,'ytick.labelsize':8.5,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none','mathtext.fontset':'custom','mathtext.rm':family,'mathtext.it':family+':italic','mathtext.bf':family+':bold'})
fig,axs=plt.subplots(3,2,figsize=(5.25,6.15));fig.subplots_adjust(left=.135,right=.975,bottom=.085,top=.895,wspace=.52,hspace=.70)
styles=[dict(color='#1261D9',marker='o',linestyle='-',label='Graph'),dict(color='#E54B18',marker='s',linestyle='--',label='Quadratic limit')]
letters='ABCDEF';table=[]
for i,k in enumerate(['1','0.1','0.01']):
 r=sorted([x for x in rows if x['kappa']==k and x['N']==18],key=lambda x:x['R']);xs=[x['R'] for x in r]
 widths=[[float(F(x['graph']['upper_density'])-F(x['graph']['lower_density']))for x in r],[float(F(x['infinite']['spectral_width_upper']))for x in r]]
 for j,st in enumerate(styles):
  axs[i,0].plot(xs,widths[j],**st,lw=1.45,ms=4,mew=.6)
  key=['graph_seconds','quadratic_seconds'][j];times=[x['timing'][key] for x in r];med=np.array([statistics.median(t)for t in times]);lo=med-np.array([min(t)for t in times]);hi=np.array([max(t)for t in times])-med
  axs[i,1].errorbar(xs,med,yerr=np.array([lo,hi]),**st,lw=1.45,ms=4,mew=.6,capsize=2,elinewidth=1)
 for j in (0,1):
  ax=axs[i,j];ax.set_yscale('log');ax.set_xlim(2.75,9.25);ax.set_xticks([3,5,7,8,9]);ax.set_xlabel('Range R');ax.set_ylabel('Interval width' if j==0 else 'Verification time (s)')
  ax.grid(False,which='both');ax.yaxis.set_major_locator(ticker.LogLocator(base=10,numticks=4));ax.yaxis.set_minor_locator(ticker.NullLocator());ax.tick_params(width=1,length=3)
  for sp in ax.spines.values():sp.set_visible(True)
  ax.text(-.27,1.12,letters[2*i+j],transform=ax.transAxes,fontweight='bold',fontsize=13,va='bottom')
  ax.set_title('κ = '+k+'; s = '+r[0]['s'],fontsize=9.5,fontweight='bold',pad=5)
 for x in r:table.append(dict(candidate=x['candidate'],R=x['R'],kappa=k,graph_width=x['graph']['exact_width'],quadratic_width=x['infinite']['spectral_width_upper'],timing=x['timing']))
handles,labels=axs[0,0].get_legend_handles_labels();fig.legend(handles,labels,loc='upper center',bbox_to_anchor=(.54,.991),ncol=2,frameon=False,fontsize=10)
args.output.parent.mkdir(parents=True,exist_ok=True)
for ext in ['pdf','svg','png']:fig.savefig(args.output.with_suffix('.'+ext),dpi=300)
meta=dict(data_sha256=hashlib.sha256(args.data.read_bytes()).hexdigest(),font_path=font,font_family=family,figure_inches=[5.25,6.15],panels_inches=[list(np.array(ax.get_position().size)*fig.get_size_inches())for ax in axs.flat],matplotlib=matplotlib.__version__,series_styles=styles,records=table)
args.output.with_suffix('.json').write_text(json.dumps(meta,indent=2)+'\n');print(json.dumps({k:v for k,v in meta.items() if k not in ['records','series_styles']}))
