#!/usr/bin/env python3
"""Supplement S1: cell identity and CNV evidence, FigureYa adapted.

Run 13_supplement_S1_prepare.R once, then this script with system Python.
No normalization, cell classification, UMAP, or CNV inference is rerun.
"""
from pathlib import Path
import os, json, hashlib, shutil, colorsys, io, tempfile
ROOT = Path(__file__).resolve().parents[1]
SOURCE_ROOT = Path(os.environ.get("EMT_SOURCE_ROOT", "/Users/sodajay/Desktop/EMT_reproduce_final_20260910"))
OUT = Path(os.environ.get("EMT_S1_OUTPUT", str(ROOT / "figures/Supplement_S1")))
INP = OUT / 'inputs'
OUT.mkdir(parents=True, exist_ok=True)
os.environ.setdefault('MPLCONFIGDIR', str(OUT / '.mplconfig'))
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import font_manager as fm
from matplotlib.colors import LinearSegmentedColormap, TwoSlopeNorm
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import linkage, leaves_list

FONTDIR = Path('/System/Library/Fonts/Supplemental')
for name in ['Arial.ttf','Arial Bold.ttf','Arial Italic.ttf','Arial Bold Italic.ttf']:
    fm.fontManager.addfont(str(FONTDIR/name))
assert Path(fm.findfont(fm.FontProperties(family='Arial'))).name == 'Arial.ttf'
plt.rcParams.update({'font.family':'Arial','font.size':10,'axes.titlesize':14,
    'axes.labelsize':10,'xtick.labelsize':10,'ytick.labelsize':10,
    'legend.fontsize':10,'pdf.fonttype':42,'ps.fonttype':42,'svg.fonttype':'none',
    'axes.linewidth':.65,'axes.edgecolor':'#666666','text.color':'#262626',
    'axes.labelcolor':'#262626','xtick.color':'#444444','ytick.color':'#444444',
    'savefig.facecolor':'white','figure.facecolor':'white'})
PAL = {'Hepatocyte':'#E89B86','Epithelial (cholangiocyte)':'#7895B2',
       'Fibroblast/CAF':'#9CCDBD','Endothelial':'#F2B39A','Myeloid':'#72B9A8',
       'T cell':'#83C4D8','B cell':'#A89078','Plasma cell':'#A3ADCA'}
ROLECOL = {'qHSC':'#83C4D8','Activated HSC':'#72B9A8',
           'Portal fibroblast-like':'#9CCDBD','Fibrotic MFB':'#A3ADCA'}
CNVCOL={'malignant_epi':'#F2B39A','MFB':'#A3ADCA','other_CAF':'#9CCDBD','immune_ref':'#ADADAD'}
CNVLAB={'malignant_epi':'Malignant epithelium','MFB':'MFB','other_CAF':'Other CAF','immune_ref':'Immune reference'}
CMAP=LinearSegmentedColormap.from_list('blue_white_coral',['#7895B2','#FFFFFF','#E89B86'])
W,H=9.8,20.15
fig=plt.figure(figsize=(W,H))
M=.20; GAP=.20; CW=(W-2*M-GAP)/2
X=[M,M+CW+GAP]
rects={}

def axbox(x,top,w,h): return fig.add_axes([x/W,1-(top+h)/H,w/W,h/H])
def txt(x,top,s,**kw):
    return fig.text(x/W,1-top/H,s,va='top',fontsize=kw.pop('fontsize',10),**kw)
def header(letter,title,x,top,w):
    txt(x,top,letter,fontsize=16,fontweight='bold')
    txt(x+.27,top+.015,title,fontsize=14,fontweight='bold')
def panel(letter,title,x,top,w,h,subtitle=None):
    rects[letter]={'x':x,'top':top,'w':w,'h':h}
    header(letter,title,x,top,w)
    if subtitle: txt(x+.27,top+.30,subtitle,color='#555555')
def mapax(x,top,w,h,frame):
    ax=axbox(x,top,w,h)
    xx=frame.UMAP_1; yy=frame.UMAP_2
    dx=xx.max()-xx.min();dy=yy.max()-yy.min()
    ax.set_xlim(xx.min()-.045*dx,xx.max()+.045*dx)
    ax.set_ylim(yy.min()-.045*dy,yy.max()+.045*dy)
    ax.set_aspect('equal');ax.axis('off')
    return ax
def categorical(ax,d,col,palette,size=1.3,alpha=.9,rasterized=True):
    z=d.sample(frac=1,random_state=20260911)
    ax.scatter(z.UMAP_1,z.UMAP_2,s=size,c=[palette[v] for v in z[col]],
        linewidths=0,alpha=alpha,rasterized=rasterized)
def number_labels(ax,d,col):
    med=d.groupby(col)[['UMAP_1','UMAP_2']].median()
    # Spread only labels, never data coordinates; deterministic local repulsion.
    xy=med.to_numpy().copy();orig=xy.copy()
    xr=np.ptp(ax.get_xlim());yr=np.ptp(ax.get_ylim())
    for _ in range(110):
        for i in range(len(xy)):
            for j in range(i):
                delta=(xy[i]-xy[j])/np.array([xr,yr])
                if abs(delta[0])<.048 and abs(delta[1])<.06:
                    v=np.array([delta[0],delta[1]])
                    if np.linalg.norm(v)<1e-5:v=np.array([.01,-.01])
                    v=v/np.linalg.norm(v)*.0028*np.array([xr,yr])
                    xy[i]+=v;xy[j]-=v
        xy+=(orig-xy)*.006
    for key,pos,original in zip(med.index,xy,orig):
        ax.annotate(str(key),xy=original,xytext=pos,fontsize=10,ha='center',va='center',
             bbox=dict(facecolor='white',edgecolor='none',pad=.12,alpha=.86),
             arrowprops=dict(arrowstyle='-',lw=.35,color='#999999') if np.linalg.norm(pos-original)>.4 else None)
def handles(colors,labels=None,ms=6):
    return [Line2D([],[],ls='',marker='o',markersize=ms,markerfacecolor=c,
        markeredgewidth=0,label=labels.get(k,k) if labels else k) for k,c in colors.items()]

global_source=INP/'global_umap_plot_data.csv.gz'
if not global_source.exists():
    shutil.copy2(SOURCE_ROOT/'EMT_Resutls/data/global_umap_plot_data.csv.gz',global_source)
g=pd.read_csv(global_source)
assert len(g)==95022 and g.barcode.nunique()==95022
g['cluster']=g['cluster'].astype(int)
st=pd.read_csv(INP/'stromal_umap_frozen.csv')
st['cluster']=st['cluster'].astype(int)
assert len(st)==3018 and st.barcode.nunique()==3018

# A-C: shared pre-existing coordinates, cluster and annotation coloring only.
top=.20; RH=3.90
panel('A','Before cell-type annotation',X[0],top,CW,RH,'95,022 cells | 40 cell clusters')
a=mapax(X[0],top+.51,CW,RH-.96,g)
clcol={i:colorsys.hsv_to_rgb((i*.61803398875)%1,.40,.78) for i in sorted(g.cluster.unique())}
categorical(a,g,'cluster',clcol,size=.45,alpha=.77);number_labels(a,g,'cluster')
panel('B','Major cell-type annotation',X[1],top,CW,RH,'95,022 cells | eight major cell types')
b=mapax(X[1],top+.51,CW,RH-.96,g)
categorical(b,g,'celltype',PAL,size=.45,alpha=.82)
labs={k:('Cholangiocyte' if k=='Epithelial (cholangiocyte)' else k) for k in PAL}
legax=axbox(X[1],top+RH-.43,CW,.39);legax.axis('off')
legax.legend(handles=handles(PAL,labs),loc='center',ncol=4,frameon=False,
    columnspacing=.68,handletextpad=.27,handlelength=.6,labelspacing=.28)
# D/E: role overlay and the stored qHSC module score, without new calls.
top=4.30; RH=3.65
mapping={'qHSC_root':'qHSC','activated_HSC_intermediate':'Activated HSC',
 'portal_fibroblast_root':'Portal fibroblast-like','portal_fibroblast_like_root':'Portal fibroblast-like',
 'fibrotic_MFB_terminal':'Fibrotic MFB'}
st['role']=st.frozen_role.map(mapping); selected=st[st.role.notna()]
assert len(selected)==757
panel('C','Stromal cell types',X[0],top,CW,RH,'757 classified cells; other stromal cells in gray')
d=mapax(X[0],top+.51,CW,RH-1.0,st)
d.scatter(st.UMAP_1,st.UMAP_2,s=2.3,c='#E7E7E7',lw=0,rasterized=True)
categorical(d,selected,'role',ROLECOL,size=7,alpha=.92)
rlab={k:k.replace('Portal fibroblast-like','Portal fibroblast(-like)')+f' (n={sum(selected.role==k)})' for k in ROLECOL}
legax=axbox(X[0],top+RH-.45,CW,.42);legax.axis('off')
legax.legend(handles=handles(ROLECOL,rlab),loc='center',ncol=2,frameon=False,
    columnspacing=.6,handletextpad=.3,handlelength=.6,labelspacing=.45)
panel('D','qHSC identity program',X[1],top,CW,RH,'qHSC program score in 3,018 stromal cells')
e=mapax(X[1],top+.51,CW,RH-1.0,st)
qlo,qhi=np.quantile(st.qHSC_score,[.01,.99]);norm=TwoSlopeNorm(vmin=qlo,vcenter=0,vmax=qhi)
so=st.sort_values('qHSC_score')
im=e.scatter(so.UMAP_1,so.UMAP_2,c=so.qHSC_score,s=3.8,lw=0,cmap=CMAP,norm=norm,rasterized=True)
cbax=axbox(X[1]+1.04,top+RH-.34,CW-2.08,.10)
cb=fig.colorbar(im,cax=cbax,orientation='horizontal',extend='both');cb.set_ticks([qlo,0,qhi]);cb.set_ticklabels([f'{qlo:.2f}','0',f'{qhi:.2f}'])
cb.ax.tick_params(length=2,pad=2);cb.outline.set_linewidth(.4)

fields={'malignant_epi':'malig_LOO_corr','MFB':'MFB_corr','other_CAF':'CAF_corr','immune_ref':'immune_corr'}

# G/H: arm-level input summaries, preserving genomic order.
arm=pd.read_csv(INP/'arm_cnv_by_group.csv')
arms=arm.arm.drop_duplicates().tolist();assert len(arms)==39
cohorts=['China-HCC','USA-HCC'];groups=['malignant_epi','MFB']
matrix=[];rowlab=[];devout=[]
for cohort in cohorts:
    sub=arm[arm.dataset==cohort].pivot(index='arm',columns='celltype',values='cnv').loc[arms]
    for role in groups:
        vals=sub[role]-sub['other_CAF'];matrix.append(vals.to_numpy())
        rowlab.append(f'{cohort} | '+('malignant' if role=='malignant_epi' else 'MFB'))
        for aa,v in vals.items():devout.append({'cohort':cohort,'role':role,'arm':aa,'deviation_vs_other_CAF':v})
matrix=np.array(matrix);pd.DataFrame(devout).to_csv(INP/'arm_deviation_plot_values.csv',index=False)
top=8.15;RH=1.31
panel('E','Chromosome-arm CNV relative to other CAF',M,top,W-2*M,RH)
gx=M+2.01;gw=W-gx-M-.84
ga=axbox(gx,top+.40,gw,.80)
vmax=max(.2,min(.3,float(abs(matrix).max())*1.05))
hm=ga.imshow(matrix,aspect='auto',cmap=CMAP,vmin=-vmax,vmax=vmax,interpolation='nearest')
ga.set_yticks(np.arange(4),rowlab);ga.tick_params(axis='both',length=0);ga.set_xticks([])
ga.set_yticks(np.arange(-.5,4,1),minor=True);ga.grid(which='minor',axis='y',color='white',lw=2)
ga.tick_params(which='minor',length=0)
for ii,role in enumerate(groups*2):
    ga.add_patch(Rectangle((-.73,ii-.5),.18,1,facecolor=CNVCOL[role],edgecolor='none',clip_on=False))
for spine in ga.spines.values():spine.set_visible(False)
cbax=axbox(gx+gw+.17,top+.40,.12,.8)
cb=fig.colorbar(hm,cax=cbax);cb.set_ticks([-vmax,0,vmax]);cb.set_ticklabels([f'{-vmax:.2f}','0',f'{vmax:.2f}']);cb.ax.tick_params(length=2,pad=3);cb.outline.set_linewidth(.4)

top=9.66;RH=2.08
panel('F','Chromosome-arm CNV profiles',M,top,W-2*M,RH)
linesax=[]
for i,cohort in enumerate(cohorts):
    hh=axbox(gx,top+.43+i*.61,gw,.54)
    sub=arm[(arm.dataset==cohort)&arm.celltype.isin(groups)]
    for role in groups:
        yy=sub[sub.celltype==role].set_index('arm').loc[arms,'cnv'].to_numpy()
        hh.plot(np.arange(39),yy,color=CNVCOL[role],lw=1.25,marker='o',ms=2.25,label=CNVLAB[role])
    hh.axhline(0,lw=.6,color='#BDBDBD');hh.set_xlim(-.5,38.5);hh.set_ylim(-.15,.29)
    hh.set_yticks([0,.2]);hh.tick_params(length=2,pad=3)
    hh.set_ylabel(cohort,rotation=0,ha='right',va='center',labelpad=14)
    hh.spines[['top','right']].set_visible(False)
    hh.set_xticks(range(39))
    if i==1:hh.set_xticklabels([a.replace('chr','') for a in arms],rotation=90,va='top');hh.set_xlabel('Chromosome arm',labelpad=4)
    else:hh.set_xticklabels([]);hh.tick_params(axis='x',length=0)
    linesax.append(hh)
linesax[0].legend(loc='upper right',bbox_to_anchor=(1,1.69),frameon=False,ncol=2,
    handlelength=1.5,columnspacing=1.4,handletextpad=.5)

# G/H: original display cells; only within-group column order is clustered.
top=11.94;RH=7.67;HW=(W-2*M-GAP)/2
cnvnorms=[];display=[]
for sid in ['HCC06T','S042']:
    mat=pd.read_csv(INP/f'cellarm_{sid}.csv',index_col=0)
    audit=pd.read_csv(INP/f'display_cells_{sid}.csv')
    assert mat.columns.tolist()==audit.cell.tolist() and mat.shape[0]==39
    lim=max(.25,min(.40,float(np.quantile(np.abs(mat.to_numpy()),.995))))
    cnvnorms.append(lim);display.append((sid,mat,audit))
clim=max(cnvnorms)
counts=pd.read_csv(INP/'cnv_display_counts.csv')
allselected=[]
for index,(sid,mat,audit) in enumerate(display):
    x=M+index*(HW+GAP);letter=['G','H'][index];cohort=cohorts[index]
    panel(letter,f'{cohort}: {sid}',x,top,HW,RH,'Cell-level chromosome-arm CNV')
    # 39 chromosome-arm labels at fixed 10 pt require a sufficiently tall matrix.
    heatx=x+.44;heatw=HW-.59;heattop=top+1.09;heath=6.30
    ha=axbox(heatx,heattop,heatw,heath)
    order=[];boundaries=[];bounds=[];cursor=0
    for role in fields:
        cellids=audit.loc[audit.celltype==role,'cell'].tolist()
        mm=mat[cellids].to_numpy()
        if len(cellids)>1:
            leaf=leaves_list(linkage(mm.T,method='complete',metric='euclidean'))
            cellids=[cellids[j] for j in leaf]
        order+=cellids
        bounds.append((cursor,cursor+len(cellids),role));cursor+=len(cellids);boundaries.append(cursor)
    plot=mat[order].to_numpy()
    im=ha.imshow(plot,aspect='auto',cmap=CMAP,vmin=-clim,vmax=clim,interpolation='nearest',rasterized=True)
    ha.set_yticks(range(39),[a.replace('chr','') for a in mat.index]);ha.tick_params(length=0,pad=3);ha.set_xticks([])
    for bd in boundaries[:-1]:ha.axvline(bd-.5,color='white',lw=3.2)
    for spine in ha.spines.values():spine.set_visible(False)
    for start,end,role in bounds:
        xx=heatx+heatw*start/len(order);ww=heatw*(end-start)/len(order)
        # Equal widths encode only balanced display selection, not abundance.
        fig.patches.append(Rectangle((xx/W,1-(heattop-.03)/H),ww/W,.085/H,
            transform=fig.transFigure,facecolor=CNVCOL[role],edgecolor='none'))
        lab={'malignant_epi':'Malignant\nepithelium','MFB':'MFB','other_CAF':'Other CAF','immune_ref':'Immune\nreference'}[role]
        txt(xx+ww/2,heattop-.62,lab,ha='center',color='#333333',linespacing=1.3)
        txt(xx+ww/2,heattop-.28,f'n = {end-start}',ha='center',color='#555555')
    allselected.append({'sample':sid,'cells_displayed':len(order),'cells_per_group':[b-a for a,b,_ in bounds],
         'arm_count':39,'individual_original_color_limit':cnvnorms[index],
         'shared_new_color_limit':clim,'within_group_clustering':'Euclidean complete linkage; display only'})
    pd.DataFrame({'display_order':np.arange(len(order))+1,'cell':order}).to_csv(INP/f'cell_order_{sid}.csv',index=False)
    txt(x+.45,top+RH-.20,'Columns: cells; group widths do not indicate abundance.',color='#555555')

# One shared heatmap legend beneath both representative samples.
cbax=axbox(W/2-1.3,19.80,2.6,.09)
cb=fig.colorbar(im,cax=cbax,orientation='horizontal');cb.set_ticks([-clim,0,clim]);cb.set_ticklabels([f'{-clim:.2f}','0',f'{clim:.2f}']);cb.ax.tick_params(length=2,pad=2);cb.outline.set_linewidth(.4)
txt(W/2-1.56,19.80,'Relative CNV',ha='right')

stem=OUT/'Supplement_S1_cell_identity_CNV_Arial'
for ext in ['pdf','png','svg']:
    # Render to memory, then replace atomically: an interrupted write cannot
    # leave a truncated deliverable on a synced Desktop filesystem.
    buffer=io.BytesIO();fig.savefig(buffer,format=ext,dpi=400)
    with tempfile.NamedTemporaryFile(dir=OUT,prefix='.export_',suffix='.'+ext,delete=False) as tmp:
        tmp.write(buffer.getvalue());tmp_path=Path(tmp.name)
    tmp_path.replace(Path(str(stem)+'.'+ext))
plt.close(fig)

sources={
 'global_umap':'EMT_Resutls/data/global_umap_plot_data.csv.gz',
 'stromal_object':'results_v4/S3_stroma_reclustered.rds',
 'china_frozen_roles':'results_v4/tables/S6_HBV_stromal_frozen_roles_by_cell.csv',
 'usa_frozen_roles':'results_v4/tables/S6_USAHCC_stromal_frozen_roles_by_cell.csv',
 'cca_frozen_roles':'results_v4/tables/S6_CCAcombined_stromal_frozen_roles_by_cell.csv',
 'cnv_matching':'results_v4/tables/S4_cohort_cnv_matching.csv',
 'arm_cnv':'results_v4/tables/S4_armcnv_by_group_celltype.csv',
 'cellarm':'results_v4/cnv_matching/cohort/cellarm_list.rds',
 'HCC06T_selection':'results_v4/tables/S4_review_cellCNV_ChinaHCC_HCC06T_v4_audit.csv',
 'S042_selection':'results_v4/tables/S4_review_cellCNV_USAHCC_S042_v4_audit.csv'}
source_rows=[]
for role,path in sources.items():
    p=SOURCE_ROOT/path
    source_rows.append({'role':role,'path':path,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()})
pd.DataFrame(source_rows).to_csv(OUT/'source_manifest.csv',index=False)
summary={'panels':8,'previous_to_current_panel_mapping':{'A':'A','B':'B','D':'C','E':'D','G':'E','H':'F','I':'G','J':'H'},'removed_previous_panels':['C','F'],'page_inches':[W,H],'PNG_dpi':400,'font':str(FONTDIR/'Arial.ttf'),
 'title_pt':14,'label_pt':10,'panel_letter_pt':16,'global_n':len(g),'stromal_n':len(st),
 'stromal_celltype_n':len(selected),'stromal_celltype_counts':selected.role.value_counts().to_dict(),
 'qHSC_color_clip_percentiles':[1,99],'qHSC_color_limits':[qlo,qhi],
 'cnv_scope':'historical matching-specific MFB groups, distinct from stromal cell-type assignments',
 'representative_panels':allselected,'panel_rectangles_inches':rects}
(OUT/'plot_audit.json').write_text(json.dumps(summary,indent=2))
print('Wrote',stem,'PDF + 400 dpi PNG + editable SVG')
print(json.dumps(summary,indent=2))
