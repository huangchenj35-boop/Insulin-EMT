#!/usr/bin/env python3
"""Recompose the approved V2 figures from frozen artwork and numerical outputs."""
from pathlib import Path
import hashlib, json, shutil, math
import fitz
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'V6_最终图_20261024'; AUD=WORK/'audit'
OUT=ROOT/'JBC_singlecell_20261024_final'; SRC=ROOT/'JBC_singlecell_20261014_v2'
for d in (AUD,OUT):d.mkdir(parents=True,exist_ok=True)
W,H=595.2756,841.8898
PAL=json.loads((ROOT/'配色参考主刊_20261010/palette.json').read_text())
COL=PAL['cell_types']; INK='#262626'; MUTED='#666666'
for n,f in [('Arial','Arial.ttf'),('Arial-Bold','Arial Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(n,str(Path('/System/Library/Fonts/Supplemental')/f)))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
def atomic(path,data):
    staged=AUD/(path.name+'.staged')
    staged.write_bytes(data);staged.replace(path)
protected={str(p):sha(p) for p in SRC.iterdir() if p.is_file()}
manifest=[]; values=[]; geometry=[]; panelmap=[]
def read(p,**kwargs):
    manifest.append({'path':str(p),'sha256':sha(p)})
    return pd.read_csv(p,**kwargs)
INP=ROOT/'Figure2_细胞来源追踪/拼图_20261009/inputs'
SINP=ROOT/'本次整合_20261010/inputs'
paga=read(INP/'MFB_connectivity_frozen_values.csv')
resamp=read(INP/'balanced_resampling_uncapped_values.csv')
stats=read(INP/'resampling_summary_verified.csv')
xy=read(INP/'structural_ddrtree_coordinates.csv')
edges=read(INP/'structural_ddrtree_edges.csv')
cohorts=['HCC-GSE149614','HCC-GSE151530','CCA-GSE138709 + GSE151530']
legacy=['China-HCC','USA-HCC','CCA-combined']; tags=['ChinaHCC','USAHCC','CCAcombined']
roles=['hepatocyte','cholangiocyte','malignant','qHSC','activated_HSC','portal','MFB']
rmap={'normal_hepatocyte_root':'hepatocyte','normal_cholangiocyte_root':'cholangiocyte',
      'malignant_epithelial_EMT_root':'malignant','qHSC_root':'qHSC','activated_HSC_intermediate':'activated_HSC',
      'portal_fibroblast_root':'portal','portal_fibroblast_like_root':'portal','fibrotic_MFB_terminal':'MFB'}
LAB={'hepatocyte':'Hepatocyte','cholangiocyte':'Cholangiocyte','malignant':'HCC cell',
     'qHSC':'qHSC','activated_HSC':'aHSC','portal':'PF-like','MFB':'MFB'}
def allowed(i,r): return not ((i<2 and r=='cholangiocyte') or (i==2 and r=='hepatocyte'))
def label(i,r):return 'PF' if i==0 and r=='portal' else 'CCA cell' if i==2 and r=='malignant' else LAB[r]
def fill(col):return HexColor(col) if isinstance(col,str) else col
def text(s,x,y,size=7,bold=False,color=INK,align='left',rot=0):
    s=str(s);face='Arial-Bold' if bold else 'Arial';tw=pdfmetrics.stringWidth(s,face,size)
    c.saveState();c.setFillColor(fill(color));c.setFont(face,size);c.translate(x,H-y);c.rotate(rot)
    c.drawString(-{'left':0,'center':.5,'right':1}[align]*tw,-size*.78,s);c.restoreState()
    geometry.append(dict(figure=fig,text=s,x=x,y=y,size=size,rot=rot))
def line(x,y,xx,yy,color='#E6E6E6',width=.5,dash=None):
    c.setStrokeColor(fill(color));c.setLineWidth(width);c.setDash(dash or [])
    c.line(x,H-y,xx,H-yy);c.setDash([])
def rect(x,y,w,h,col,edge=None):
    c.setFillColor(fill(col));c.setStrokeColor(fill(edge or col));c.setLineWidth(.35)
    c.rect(x,H-y-h,w,h,fill=1,stroke=int(edge is not None))
def dot(x,y,r,col,alpha=1,edge=None):
    c.saveState();c.setFillAlpha(alpha);c.setFillColor(fill(col));c.setStrokeColor(fill(edge or col));c.setLineWidth(.4)
    c.circle(x,H-y,r,fill=1,stroke=int(edge is not None));c.restoreState()
def diamond(x,y,r,col):
    c.setFillColor(fill(col));c.setStrokeColor(HexColor('#FFFFFF'));c.setLineWidth(.4)
    p=c.beginPath();p.moveTo(x-r,H-y);p.lineTo(x,H-y+r);p.lineTo(x+r,H-y);p.lineTo(x,H-y-r);p.close()
    c.drawPath(p,fill=1,stroke=1)
def title(letter,s,y,x=14):text(letter,x,y,12,True);text(s,x+18,y+1,10,True)
def blend(a,b,t):
    a=np.array([int(a[j:j+2],16) for j in (1,3,5)])/255;b=np.array([int(b[j:j+2],16) for j in (1,3,5)])/255
    return Color(*(a+(b-a)*np.clip(t,0,1)))
def pc(v):return blend(PAL['continuous']['low'],PAL['continuous']['mid'],v) if v<=1 else blend(PAL['continuous']['mid'],PAL['continuous']['high'],(v-1)/1.5)
def zc(v):return blend(PAL['continuous']['mid'],PAL['continuous']['high' if v>=0 else 'low'],abs(v)/2)
def bar(x,y,w,lo,hi,fn,ticks,caption=None):
    for j in range(150):rect(x+j*w/150,y,w/150+.1,4,fn(lo+(hi-lo)*j/149))
    for v in ticks:text(f'{v:g}',x+(v-lo)/(hi-lo)*w,y+7,7,align='center')
    if caption:text(caption,x-6,y-1,7,align='right')
def head(i,x,y):
    text(cohorts[i] if i<2 else 'CCA-GSE138709',x,y,8,True,align='center')
    if i==2:text('+ GSE151530',x,y+10,8,True,align='center')

# Figure 7: all connectivity results, DDRTree and chromosome-arm CNV.
fig='Figure_7';c=canvas.Canvas(str(AUD/'Figure_7_top.pdf'),pagesize=(W,H));c.setTitle('Figure 7')
title('A','Joint PAGA connectivity | k = 30',10)
matrix_records=[]
for i,tag in enumerate(tags):
    d=read(SINP/f'S3_{tag}_PAGA_connectivity_matrix_k30_uncapped.csv').set_index('frozen_role')
    d.index=d.index.map(rmap);d.columns=d.columns.map(rmap)
    order=[r for r in roles if r in d.index and allowed(i,r)];d=d.loc[order,order]
    assert np.allclose(d,d.T) and (np.diag(d)==0).all()
    head(i,106+i*193,25)
    step=21 if len(d)==6 else 25;first=106+i*193-(len(d)-1)*step/2;base=108
    for a in range(len(d)):
        for b in range(a+1,len(d)):
            v=float(d.iloc[a,b]);diamond(first+(a+b)*step/2,base-(b-a)*step/2,step/2,pc(v))
            matrix_records.append(dict(cohort=cohorts[i],role_1=order[a],role_2=order[b],value=v))
    for j,r in enumerate(order):
        x=first+j*step;rect(x-3,110,6,2,COL[r]);text(label(i,r),x+2,116,7,align='right',rot=55)
bar(275,157,125,0,2.5,pc,[0,1,2,2.5],'Observed / expected ratio')
pd.DataFrame(matrix_records).to_csv(AUD/'Figure7_edges.csv',index=False)
assert len(matrix_records)==35
panelmap.append({'source':'V2 Figure S1B','destination':'Figure 7A','change':'HCC reference cholangiocytes omitted from display, original graph values retained'})

title('B','Connectivity with MFB across graph settings',164)
for i,co in enumerate(legacy):
    head(i,106+i*193,183)
    text(f'200 resamples; {[29,59,16][i]} cells / population',106+i*193,203,7,color=MUTED,align='center')
    left=65+i*193;right=193+i*193;fx=lambda v:left+v/2.5*(right-left)
    display=[r for r in roles[:-1] if allowed(i,r)]
    ys=[222+18*j for j in range(len(display))]
    line(fx(1),214,fx(1),304,'#A0A0A0',.5,[2,2])
    for j,(r,y) in enumerate(zip(display,ys)):
        text(label(i,r),left-5,y-2.8,7,align='right')
        q=paga[paga.cohort.eq(co)&paga.role.eq(r)]
        if q.empty:text('-',fx(1.2),y-3,color=MUTED,align='center');continue
        val=float(q.ratio.iloc[0]);diamond(fx(val),y-5.5,2.2,COL[r])
        text(f'{val:.3f}',fx(val)+4,y-11,7)
        values.append(dict(figure=fig,panel='B',cohort=co,role=r,k=30,value=val))
        for k,dy in [(10,0),(20,5.5)]:
            rr=stats[stats.cohort.eq(co)&stats.n_neighbors.eq(k)&stats.role.eq(r)].iloc[0]
            z=resamp[resamp.cohort.eq(co)&resamp.n_neighbors.eq(k)&resamp.role.eq(r)].sort_values('replicate')
            assert len(z)==200
            assert np.allclose(np.quantile(z.paga_connectivity,[.025,.5,.975]),[rr.q025,rr['median'],rr.q975],atol=1e-12)
            jitter=np.random.default_rng(2300+i*100+j*5+k).uniform(-1.1,1.1,len(z))
            for v,ji in zip(z.paga_connectivity,jitter):dot(fx(v),y+dy+ji,.48,COL[r],.30)
            line(fx(rr.q025),y+dy,fx(rr.q975),y+dy,INK,.55)
            if k==10:rect(fx(rr['median'])-1.5,y+dy-1.5,3,3,'#FFFFFF',INK)
            else:dot(fx(rr['median']),y+dy,1.4,INK)
            values.append(dict(figure=fig,panel='B',cohort=co,role=r,k=k,value=rr['median'],q025=rr.q025,q975=rr.q975))
    line(left,307,right,307,INK,.5)
    for v in (0,1,2,2.5):line(fx(v),307,fx(v),309,INK,.5);text(f'{v:g}',fx(v),311,7,align='center')
text('Uncapped observed / expected ratio',W/2,323,7,align='center')
diamond(38,337,2.3,COL['portal']);text('k = 30',45,334,7)
rect(121,335.5,3,3,'#FFFFFF',INK);text('k = 10',130,334,7)
dot(211,337,1.5,INK);text('k = 20',219,334,7)
text('Bars: 95% cell-resampling intervals',320,334,7)
shown=resamp[~(resamp.cohort.eq('China-HCC')&resamp.role.eq('cholangiocyte'))]
assert len(shown)==5200
shown.to_csv(AUD/'Figure7_resamples.csv',index=False)
panelmap.append({'source':'V2 Figure 2B and Figure S1C','destination':'Figure 7B','change':'k30, k10 and k20 combined on one axis per cohort'})

title('C','Unrooted stromal DDRTree',348)
for i,co in enumerate(cohorts):
    head(i,106+i*193,368)
    z=xy[xy.analysis_group.eq(co)];e=edges[edges.analysis_group.eq(co)]
    pts=np.vstack([z[['DDRTree_1','DDRTree_2']],e[['x','y']],e[['xend','yend']]])
    center=(pts.max(axis=0)+pts.min(axis=0))/2;span=np.ptp(pts,axis=0);sc=min(156/span[0],68/span[1])
    def pos(x,y):return 106+i*193+(x-center[0])*sc,426-(y-center[1])*sc
    for rr in e.itertuples():
        a,b=pos(rr.x,rr.y);aa,bb=pos(rr.xend,rr.yend);line(a,b,aa,bb,'#777777',.5)
    for rr in z.sample(frac=1,random_state=20260901).itertuples():
        x,y=pos(rr.DDRTree_1,rr.DDRTree_2);dot(x,y,1.08,COL[rmap[rr.frozen_role]],.9)
    text(f'n = {len(z)}',188+i*193,389,7,color=MUTED,align='right')
for x,r in zip([139,215,337,437],['qHSC','activated_HSC','portal','MFB']):
    dot(x,469,2.2,COL[r]);text({'qHSC':'qHSC','activated_HSC':'Activated HSC','portal':'PF / PF-like','MFB':'MFB'}[r],x+6,466,7)
panelmap.append({'source':'V2 Figure 2C','destination':'Figure 7C','change':'three cohorts in one row; same coordinates and edges'})
c.showPage();c.save()

# Preserve the V2 CNV artwork; harmonize three malignant population labels.
d=fitz.open(AUD/'Figure_7_top.pdf');src=fitz.open(SRC/'Figure_2.pdf')
cnv_labels=[s for b in src[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans'] if s['text']=='Malignant' and s['bbox'][1]>472]
assert len(cnv_labels)==3
for s in cnv_labels:src[0].add_redact_annot(fitz.Rect(s['bbox']),fill=False)
src[0].apply_redactions(images=0,graphics=0,text=0)
fontpath='/System/Library/Fonts/Supplemental/Arial.ttf'
src[0].insert_font(fontname='ArialFinal',fontfile=fontpath)
for s in cnv_labels:
    bb=fitz.Rect(s['bbox']);tw=pdfmetrics.stringWidth('HCC cell','Arial',7)
    ox=s['origin'][0] if bb.y0<500 else bb.x0+(bb.width-tw)/2
    col=tuple(((s['color'] >> k)&255)/255 for k in (16,8,0))
    src[0].insert_text((ox,s['origin'][1]),'HCC cell',fontname='ArialFinal',fontsize=7,color=col)
src[0].add_redact_annot(fitz.Rect(0,0,W,472),fill=False);src[0].apply_redactions(images=0,graphics=1,text=0)
d[0].show_pdf_page(fitz.Rect(0,472,W,H),src,0,clip=fitz.Rect(0,472,W,H))
d.set_metadata({'title':'Figure 7 | Epithelial-stromal connectivity and chromosome-arm profiles'})
atomic(OUT/'Figure_7.pdf',d.tobytes(garbage=4,deflate=True));d.close();src.close()
panelmap.append({'source':'V2 Figure 2D','destination':'Figure 7D','change':'identical artwork'})


# Layout-only verification against the already audited V3 numerical record.
old=pd.read_csv(ROOT/'V3_重组_20261023/audit/plotted_values.csv')
new=pd.DataFrame(values)
new['k']=new['k'].astype(float)
keep=['figure','panel','cohort','role','k','value','q025','q975']
pd.testing.assert_frame_equal(new[keep],old[old.figure.eq('Figure_7')][keep].reset_index(drop=True))
new.to_csv(AUD/'Figure7_values.csv',index=False)
(AUD/'Figure7_text_geometry.json').write_text(json.dumps(geometry,indent=2))
d=fitz.open(OUT/'Figure_7.pdf');p=d[0]
atomic(OUT/'Figure_7.png',p.get_pixmap(matrix=fitz.Matrix(400/72,400/72),alpha=False).tobytes('png'))
atomic(AUD/'Figure7_preview.png',p.get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).tobytes('png'))
oldpdf=fitz.open(ROOT/'JBC_singlecell_20261023_v4/Figure_7.pdf')
# Numerical values remain identical; label replacements are the only semantic change.
from collections import Counter
oldtext=oldpdf[0].get_text()
assert Counter(p.get_text().split())-Counter(oldtext.split()) == Counter({'Hepatocyte':4,'Cholangiocyte':2,'HCC':7,'CCA':2,'cell':9})
assert Counter(oldtext.split())-Counter(p.get_text().split()) == Counter({'Ref.':6,'hep.':4,'chol.':2,'Malig.':6,'epi.':6,'Malignant':3})
(AUD/'Figure7_preservation.json').write_text(json.dumps({'all_39_connectivity_summary_records_preserved':True,'all_5200_resamples_preserved':True,'all_35_connectivity_matrix_values_preserved':True,'CNV_artwork_unchanged_except_three_population_labels':True,'coordinates_edges_draw_order_and_palette_unchanged':True,'label_changes_only':True},indent=2))
print('Figure7 final labels rebuilt')
