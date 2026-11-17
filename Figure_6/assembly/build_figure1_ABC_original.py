#!/usr/bin/env python3
"""Rebuild Figure 1 from saved coordinates and audited numerical tables.

Run with the Codex bundled Python (numpy, pandas, Pillow, reportlab, pypdf).
The first run snapshots the small numerical inputs; subsequent runs use them.
No clustering, cell selection, expression scoring, or GSVA is rerun.
"""
from pathlib import Path
import colorsys
import hashlib
import json
import math
import shutil
import subprocess

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color, white
from reportlab.lib.utils import ImageReader
from pypdf import PdfReader

WENG = Path(__file__).resolve().parents[2]
PALETTE_ROOT = WENG / '配色参考主刊_20261010'
PALETTE = json.loads((PALETTE_ROOT / 'palette.json').read_text())
REV = WENG / 'Figure1_静态表达/拼图_20261009'
WORKSPACE = WENG.parent
PROJECT = WORKSPACE.parent
UPSTREAM = PROJECT / 'EMT_reproduce_final_20260910'
TABLES = UPSTREAM / 'emt_invivo_bodyfig_standalone_20260901/tables'
INPUTS = REV / 'inputs'
OUTPUT = WENG / 'JBC_singlecell_20261013_v2'
AUDIT = WENG / 'Figure1_修订_20261013_v2/audit'
for directory in [INPUTS, AUDIT, OUTPUT]:
    directory.mkdir(parents=True, exist_ok=True)

KEYS = ['ChinaHCC', 'USAHCC', 'CCAcombined']
COHORTS = ['HCC-GSE149614', 'HCC-GSE151530', 'CCA-GSE138709 + GSE151530']
ROLEMAP = {'Normal hepatocyte': 'hepatocyte', 'Normal cholangiocyte': 'cholangiocyte',
           'Malignant epithelial': 'malignant', 'qHSC': 'qHSC',
           'Activated HSC': 'activated_HSC', 'Portal fibroblast': 'portal',
           'Portal fibroblast-like': 'portal', 'Fibrotic MFB': 'MFB'}
RAWMAP = {'normal_hepatocyte_root': 'hepatocyte', 'normal_cholangiocyte_root': 'cholangiocyte',
          'malignant_epithelial_EMT_root': 'malignant', 'qHSC_root': 'qHSC',
          'activated_HSC_intermediate': 'activated_HSC', 'portal_fibroblast_root': 'portal',
          'portal_fibroblast_like_root': 'portal', 'fibrotic_MFB_terminal': 'MFB'}
ORDER = ['hepatocyte', 'cholangiocyte', 'malignant', 'qHSC', 'activated_HSC', 'portal', 'MFB']
COLORS = {role: PALETTE['cell_types'][role] for role in ORDER}
LABELS = {'hepatocyte': 'Ref. hepatocyte', 'cholangiocyte': 'Ref. cholangiocyte',
          'malignant': 'Malignant epithelial', 'qHSC': 'qHSC',
          'activated_HSC': 'Activated HSC', 'portal': 'PF / PF-like', 'MFB': 'MFB'}
MAJOR = {'Hepatocyte': COLORS['hepatocyte'], 'Epithelial (cholangiocyte)': COLORS['cholangiocyte'],
         'Fibroblast/CAF': COLORS['portal'], 'Endothelial': PALETTE['major_extra']['Endothelial'],
         'Myeloid': COLORS['activated_HSC'], 'T cell': COLORS['qHSC'],
         'B cell': PALETTE['major_extra']['B cell'], 'Plasma cell': COLORS['MFB']}
PHENO = ['EPCAM','KRT8','KRT18','KRT19','KRT7','CDH1','OCLN','TJP1','DSP','JUP','CLDN3','CLDN4','CLDN7','TACSTD2',
         'VIM','FN1','CDH2','ITGA5','ITGB1','LGALS1','SPARC','COL1A1','COL1A2','COL3A1','DCN','LUM','ACTA2','TAGLN','FAP','POSTN','CTHRC1','COL11A1']
TFS = ['GRHL2','OVOL1','OVOL2','ELF3','HNF4A','HNF1A','HNF1B','SOX9','KLF5','GATA6',
       'SNAI1','SNAI2','ZEB1','ZEB2','TWIST1','TWIST2','PRRX1','PRRX2','FOXC2','TCF4','FOSL1','RUNX2']
PATHWAYS = ['TGF', 'INSULIN-AKT', 'HGF-AKT', 'EGF-RAS']
SIGNAL_GENES = {
    'TGF': ['TGFB1','TGFB2','TGFBR1','TGFBR2','SMAD2','SMAD3','SMAD4','SMAD7','SERPINE1','CTGF','COL1A1'],
    'INSULIN-AKT': ['INSR','IRS1','IRS2','PIK3CA','PIK3CB','AKT1','AKT2','MTOR','FOXO1','GSK3B','RPS6KB1'],
    'HGF-AKT': ['HGF','MET','GAB1','PIK3CA','AKT1','MAPK1','MAPK3','SOS1','GRB2'],
    'EGF-RAS': ['EGF','EGFR','GRB2','SOS1','KRAS','NRAS','HRAS','RAF1','BRAF','MAP2K1','MAPK1','ELK1']}
QHSC_GENES = ['LRAT','RELN','NGFR','GFAP','VIPR1','ANGPTL6','RBP1','CYGB','DES','HGF','COLEC11','PTH1R']
FONTDIR = Path('/System/Library/Fonts/Supplemental')
for face, filename in [('Arial','Arial.ttf'), ('Arial-Bold','Arial Bold.ttf'), ('Arial-Italic','Arial Italic.ttf')]:
    pdfmetrics.registerFont(TTFont(face, str(FONTDIR/filename)))


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


sources = {'global_umap_plot_data.csv.gz': PROJECT/'figures/Supplement_S1/inputs/global_umap_plot_data.csv.gz',
           'stromal_umap_frozen.csv': PROJECT/'figures/Supplement_S1/inputs/stromal_umap_frozen.csv'}
for key in KEYS:
    for suffix in ['phenotype_optimized_zscore', 'TF_expression_optimized_zscore',
                   'current_gene_mean_by_role', 'signaling_GSVA_activity_by_role', 'signaling_GSVA_activity_zscore']:
        name = f'fig2_{key}_{suffix}.csv'
        sources[name] = TABLES/name
    name = f'S2_{key}_joint_with_normal_roles_input.csv'
    sources[name] = UPSTREAM/'emt_invivo_bodyfig_20260901/tables'/name
manifest = []
for name, source in sources.items():
    target = INPUTS/name
    if not target.exists():
        shutil.copy2(source, target)
    assert digest(source) == digest(target), f'Snapshot differs from source: {name}'
    manifest.append({'input': name, 'source': str(source), 'sha256': digest(target), 'bytes': target.stat().st_size})
(AUDIT/'input_manifest.json').write_text(json.dumps(manifest, indent=2)+'\n')

global_cells = pd.read_csv(INPUTS/'global_umap_plot_data.csv.gz')
stroma = pd.read_csv(INPUTS/'stromal_umap_frozen.csv')
global_cells['cluster'] = global_cells.cluster.astype(int)
assert len(global_cells) == global_cells.barcode.nunique() == 95022
assert global_cells.cluster.nunique() == 40 and global_cells.celltype.nunique() == 8
assert len(stroma) == stroma.barcode.nunique() == 3018
assert set(stroma.barcode).issubset(set(global_cells.barcode))
stroma['role'] = stroma.frozen_role.map(RAWMAP)
selected = stroma[stroma.role.isin(['qHSC','activated_HSC','portal','MFB'])]
assert selected.role.value_counts().to_dict() == {'portal':289, 'qHSC':247, 'MFB':159, 'activated_HSC':62}
qlo, qhi = np.quantile(stroma.qHSC_score, [.01,.99])
assert np.isfinite(stroma.qHSC_score).all()

values = []
counts = []
matrices = {}
roles = {}
errors = {}
for i, key in enumerate(KEYS):
    means = pd.read_csv(INPUTS/f'fig2_{key}_current_gene_mean_by_role.csv')
    means['role'] = means.group.map(RAWMAP)
    cells = pd.read_csv(INPUTS/f'S2_{key}_joint_with_normal_roles_input.csv')
    cells['role'] = cells.frozen_role.map(RAWMAP)
    for kind, suffix, genes in [('phenotype','phenotype_optimized_zscore', PHENO),
                                 ('TF','TF_expression_optimized_zscore', TFS),
                                 ('GSVA','signaling_GSVA_activity_zscore',PATHWAYS)]:
        frame = pd.read_csv(INPUTS/f'fig2_{key}_{suffix}.csv', index_col=0)
        assert frame.index.tolist() == genes
        frame = frame.rename(columns=ROLEMAP)
        keep = [r for r in ORDER if r in frame.columns]
        roles[key] = keep
        frame = frame[keep]
        raw = (pd.read_csv(INPUTS/f'fig2_{key}_signaling_GSVA_activity_by_role.csv',index_col=0).rename(columns=ROLEMAP)[keep]
               if kind == 'GSVA' else means.pivot(index='gene',columns='role',values='mean_expression').loc[genes,keep])
        expected = raw.sub(raw.mean(axis=1),axis=0).div(raw.std(axis=1,ddof=1).replace(0,np.nan),axis=0).fillna(0)
        error = float(np.max(np.abs(expected.to_numpy()-frame.to_numpy())))
        assert error < 1e-12, (kind,key,error)
        assert np.isfinite(frame.to_numpy()).all() and abs(frame.to_numpy()).max() < 2.5
        errors[f'{kind}_{key}'] = error
        matrices[(kind,key)] = (frame, raw)
        for gene in genes:
            for role in keep:
                values.append({'panel':{'phenotype':'B','TF':'C','GSVA':'D'}[kind], 'cohort':COHORTS[i],
                               'feature':gene,'role':role,'source_value':float(frame.loc[gene,role]),
                               'raw_mean':float(raw.loc[gene,role]),'all_zero_row':bool((raw.loc[gene]==0).all()) if kind!='GSVA' else False})
    cc = cells[cells.role.isin(roles[key])].role.value_counts()
    assert cc.sum() == [705,454,369][i]
    counts.extend({'cohort':COHORTS[i],'role':role,'cells':int(cc[role])} for role in roles[key])
# Retain the exact original four signaling matrices. Display labels only change.
rename_signal = {'TGF':'TGF','INSULIN-AKT':'PI3K-AKT','HGF-AKT':'HGF','EGF-RAS':'EGF'}
for key in KEYS:
    frame,raw=matrices[('GSVA',key)]
    matrices[('GSVA',key)]=(frame.rename(index=rename_signal),raw.rename(index=rename_signal))
PATHWAYS=['TGF','PI3K-AKT','HGF','EGF']
SIGNAL_GENES={rename_signal[k]:v for k,v in SIGNAL_GENES.items()}
for v in values:
    v['feature']=rename_signal.get(v['feature'],v['feature'])
pd.DataFrame(values).to_csv(AUDIT/'source_values.csv',index=False)
pd.DataFrame(counts).to_csv(AUDIT/'cell_counts.csv',index=False)
gene_sets = [('A','qHSC identity',QHSC_GENES),('B','Epithelial phenotype',PHENO[:14]),
             ('B','Mesenchymal phenotype',PHENO[14:]),('C','Epithelial TF',TFS[:10]),
             ('C','Mesenchymal TF',TFS[10:])]
gene_sets += [('D',name,genes) for name,genes in SIGNAL_GENES.items()]
pd.DataFrame([{'panel':p,'gene_set':name,'gene':gene,'order':i+1}
              for p,name,genes in gene_sets for i,gene in enumerate(genes)]).to_csv(AUDIT/'Figure_1_gene_sets.csv',index=False)

W,H = 612,840
pdf_path = OUTPUT/'Figure_1.pdf'
A4W,A4H=595.2756,841.8898
INK='#262626'
MUTED=PALETTE['neutral']['muted']
TYPE_SIZES={'panel_letter':12, 'panel_title':10, 'cohort_title':8, 'body':7}
c = canvas.Canvas(str(pdf_path),pagesize=(A4W,A4H),pageCompression=1)
PAGE_SCALE=A4W/W
c.translate(0,(A4H-H*PAGE_SCALE)/2)
c.scale(PAGE_SCALE,PAGE_SCALE)
c.setTitle('Figure 1 | Cell identities and transcriptional programs')
c.setAuthor('')
text_records=[]
tile_records=[]
scatter_records=[]


def text(x, top, value, size=7, face='Arial', color=INK, align='left', rotation=0):
    value=str(value)
    if size < 7.9: size=7
    elif size < 9: size=8
    size=size/PAGE_SCALE  # Keep requested type size in physical A4 points.
    width=pdfmetrics.stringWidth(value,face,size)
    anchor={'left':0,'center':.5,'right':1}[align]
    c.saveState();c.setFillColor(HexColor(color));c.setFont(face,size)
    if rotation:
        c.translate(x,H-top);c.rotate(rotation)
        c.drawString(-anchor*width,-size*.78,value)
    else:
        c.drawString(x-anchor*width,H-top-size*.78,value)
    c.restoreState()
    if rotation == 0:
        bounds=[x-anchor*width,top,x+(1-anchor)*width,top+size]
    elif rotation == 90:
        bounds=[x-size*.22,top-(1-anchor)*width,x+size*.78,top+anchor*width]
    else:
        raise ValueError(rotation)
    text_records.append({'text':value,'size':size,'bbox':bounds,'rotation':rotation})


def rect(x, top, w, h, color, border=None):
    c.setFillColor(HexColor(color) if isinstance(color,str) else color)
    if border:
        c.setStrokeColor(HexColor(border));c.setLineWidth(.25)
    c.rect(x,H-top-h,w,h,stroke=int(border is not None),fill=1)


def circle(x,top,r,color):
    c.setFillColor(HexColor(color));c.circle(x,H-top,r,stroke=0,fill=1)


def title(letter, value, x, y):
    text(x,y,letter,12,'Arial-Bold')
    text(x+20,y+1,value,10,'Arial-Bold')


def legend_item(x,y,label,color,size=7.4):
    circle(x+2,y+3.6,2.35,color)
    text(x+8,y,label,size)


def palette(value, low=-2.5, high=2.5):
    # Blue-white-coral endpoints; exact values and the existing scale are preserved.
    stops=np.array([[int(PALETTE['continuous'][name][j:j+2],16) for j in (1,3,5)]
                    for name in ('low','mid','high')],float)
    t=np.clip((value-low)/(high-low),0,1)
    rgb=(stops[0]+t*2*(stops[1]-stops[0])) if t<=.5 else (stops[1]+(t-.5)*2*(stops[2]-stops[1]))
    return tuple(np.clip(np.rint(rgb),0,255).astype(int))


def scatter(frame,x,top,w,h,colors,radius=.36,alpha=.8,labels=False):
    # Original UMAP coordinates, rigid scale, no coordinate refit or point subsampling.
    extent=frame[['UMAP_1','UMAP_2']].agg(['min','max'])
    lo=extent.loc['min'].to_numpy();hi=extent.loc['max'].to_numpy()
    delta=hi-lo;lo=lo-.045*delta;hi=hi+.045*delta
    scale=min(w/(hi[0]-lo[0]),h/(hi[1]-lo[1]))
    pw=(hi[0]-lo[0])*scale;ph=(hi[1]-lo[1])*scale
    x=x+(w-pw)/2;top=top+(h-ph)/2
    scatter_records.append({'cells':len(frame),'x':x,'top':top,'width':pw,'height':ph,
                            'original_coordinate_scale':float(scale)})
    uv=frame[['UMAP_1','UMAP_2']].to_numpy()
    xy=np.column_stack(((uv[:,0]-lo[0])*scale,(hi[1]-uv[:,1])*scale))
    pixels=7
    im=Image.new('RGB',(math.ceil(pw*pixels),math.ceil(ph*pixels)),(255,255,255))
    draw=ImageDraw.Draw(im,'RGBA')
    rr=radius*pixels
    for (xx,yy),color in zip(xy,colors):
        rgb=tuple(int(color[j:j+2],16) for j in (1,3,5)) if isinstance(color,str) else tuple(color)
        draw.ellipse((xx*pixels-rr,yy*pixels-rr,xx*pixels+rr,yy*pixels+rr),fill=(*rgb,int(alpha*255)))
    c.drawImage(ImageReader(im),x,H-top-ph,pw,ph,mask='auto')
    if labels:
        med=frame.groupby('cluster')[['UMAP_1','UMAP_2']].median()
        anchors=np.column_stack(((med.UMAP_1-lo[0])*scale,(hi[1]-med.UMAP_2)*scale))
        loc=anchors.copy()
        for _ in range(260):
            for ii in range(len(loc)):
                for jj in range(ii):
                    d=loc[ii]-loc[jj]
                    if abs(d[0])<12 and abs(d[1])<8:
                        if np.linalg.norm(d)<1e-5:d=np.array([1.,-1.])
                        shift=d/np.linalg.norm(d)*.28
                        loc[ii]+=shift;loc[jj]-=shift
            loc+=(anchors-loc)*.005
        for key,p,a in zip(med.index,loc,anchors):
            s=str(key);tw=pdfmetrics.stringWidth(s,'Arial',7/PAGE_SCALE)
            if np.linalg.norm(p-a)>3:
                c.setStrokeColor(HexColor(PALETTE['neutral']['edge']));c.setLineWidth(.25)
                c.line(x+p[0],H-top-p[1],x+a[0],H-top-a[1])
            rect(x+p[0]-tw/2-.4,top+p[1]-3.4,tw+.8,6.8,'#FFFFFF')
            text(x+p[0],top+p[1]-3.1,s,7,align='center')


# A-B share one compact band. Both B maps use the available width side by side.
title('A','Cell clustering and stromal subtypes',16,14)
text(87.5,35,'40 clusters',8,'Arial-Bold',align='center')
text(235.5,35,'Eight major cell types',8,'Arial-Bold',align='center')
g=global_cells.sample(frac=1,random_state=20260911)
# Golden-angle hue spacing balances cool and warm colors while keeping moderate
# chroma. Varying lightness also separates nearby hues across all 40 clusters.
cluster_colors={i:tuple(round(v*255) for v in colorsys.hls_to_rgb(
    (i*.61803398875)%1, [.61,.69,.56,.64][i//10], .44))
    for i in sorted(g.cluster.unique())}
scatter(g,18,49,139,139,[cluster_colors[i] for i in g.cluster],radius=.30,alpha=.85,labels=True)
scatter(g,166,49,139,139,[MAJOR[i] for i in g.celltype],radius=.30,alpha=.9)
# Explicit epithelial/tumor labels, with a compact two-line legend.
for x,label,role in [(23,'Hepatocytes /','Hepatocyte'),
                     (166,'Cholangiocytes /','Epithelial (cholangiocyte)')]:
    legend_item(x,190,label,MAJOR[role],7)
    text(x+8,200,'tumor epithelial cells',7)
for x,y,label,role in [(23,215,'Fibroblasts/CAFs','Fibroblast/CAF'),
                        (145,215,'Endothelial','Endothelial'),
                        (235,215,'Myeloid','Myeloid'),
                        (23,228,'T cells','T cell'),
                        (145,228,'B cells','B cell'),
                        (235,228,'Plasma cells','Plasma cell')]:
    legend_item(x,y,label,MAJOR[role],7)
text(161,245,'Integrated HCC / CCA atlas; 95,022 cells',7,color=MUTED,align='center')
text(387.5,35,'Stromal subtypes',8,'Arial-Bold',align='center')
text(528.5,35,'qHSC identity program',8,'Arial-Bold',align='center')
ordered=pd.concat([stroma[~stroma.index.isin(selected.index)],selected.sample(frac=1,random_state=20260911)])
st_colors=[COLORS[r] if r in ['qHSC','activated_HSC','portal','MFB'] else PALETTE['neutral']['grid'] for r in ordered.role]
scatter(ordered,319,49,137,139,st_colors,radius=.65,alpha=.95)
for index,(role,label) in enumerate([('qHSC','qHSC (247)'),('activated_HSC','Activated HSC (62)'),('portal','PF / PF-like (289)'),('MFB','MFB (159)')]):
    legend_item(344,190+index*13,label,COLORS[role],7)
q=stroma.sort_values('qHSC_score')
normalized=np.where(q.qHSC_score<0,q.qHSC_score/abs(qlo),q.qHSC_score/qhi)
scatter(q,460,49,137,139,[palette(v,-1,1) for v in normalized],radius=.6,alpha=.95)
for j in range(120):
    rgb=palette(-1+2*j/119,-1,1)
    rect(468.5+j,199,1.1,5,Color(*[v/255 for v in rgb]))
for xx,label in [(468.5,f'{qlo:.2f}'),(528.5,'0'),(588.5,f'{qhi:.2f}')]:text(xx,210,label,7,align='center')
text(528.5,230,'All 3,018 stromal cells',7,color=MUTED,align='center')

# B-D. Independent column labels for each panel, with original matrix values.
LEFT=16;RIGHT=318;CELL_W=10.55;GAP=10.5;LABEL_W=61
MATRIX_TOP=319


def displayed_roles(kind,key):
    return [r for r in roles[key]
            if not (key == 'ChinaHCC' and kind in ('phenotype','TF','GSVA') and r == 'cholangiocyte')]


def cohort_header(x, top, kind):
    offset=x+LABEL_W
    for i,key in enumerate(KEYS):
        width=len(displayed_roles(kind,key))*CELL_W
        labels=['HCC\nGSE149614','HCC\nGSE151530','CCA\nGSE138709\n+ GSE151530'][i].split('\n')
        for line,s in enumerate(labels):
            text(offset+width/2,top+line*9,s,8,'Arial-Bold',align='center')
        offset+=width+GAP


def draw_matrix(kind,genes,x,top,split=None,show_roles=False,role_top=718,row_h=12.5):
    offset=x+LABEL_W
    ycenters=[]
    for row,gene in enumerate(genes):
        y=top+row*row_h+(5.5 if split is not None and row>=split else 0)
        ycenters.append(y+row_h/2)
        text(x+LABEL_W-7,y+(row_h-7/PAGE_SCALE)/2,gene,7,'Arial-Italic' if kind!='GSVA' else 'Arial',align='right')
        if kind!='GSVA':
            rect(x+LABEL_W-4.5,y,2,row_h,PALETTE['emt_states']['epithelial'] if row<split else PALETTE['emt_states']['mesenchymal'])
    for i,key in enumerate(KEYS):
        frame,raw=matrices[(kind,key)]
        for column,role in enumerate(displayed_roles(kind,key)):
            xx=offset+column*CELL_W
            rect(xx,top-6,CELL_W-.5,3,COLORS[role])
            for row,gene in enumerate(genes):
                y=ycenters[row]-row_h/2
                val=float(frame.loc[gene,role]);rgb=palette(val)
                rect(xx,y,CELL_W,row_h,Color(*[v/255 for v in rgb]),'#FFFFFF')
                if kind!='GSVA' and bool((raw.loc[gene]==0).all()):
                    circle(xx+CELL_W/2,y+row_h/2,.55,PALETTE['neutral']['edge'])
                tile_records.append({'panel':{'phenotype':'B','TF':'C','GSVA':'D'}[kind],
                                     'cohort':COHORTS[i],'feature':gene,'role':role,'value':val,
                                     'x':xx,'top':y,'width':CELL_W,'height':row_h,'rgb':rgb})
            if show_roles:
                text(xx+CELL_W/2-2,role_top,LABELS[role],7,align='right',rotation=90)
        offset+=len(displayed_roles(kind,key))*CELL_W+GAP


title('B','EMT phenotype marker RNA',LEFT,264)
title('C','EMT transcription factor RNA',RIGHT,264)
cohort_header(LEFT,285,'phenotype')
cohort_header(RIGHT,285,'TF')
draw_matrix('phenotype',PHENO,LEFT,MATRIX_TOP,14,True,row_h=12)
draw_matrix('TF',TFS,RIGHT,MATRIX_TOP,10,True,role_top=553,row_h=10)
title('D','Signaling programs (GSVA)',RIGHT,625)
cohort_header(RIGHT,645,'GSVA')
draw_matrix('GSVA',PATHWAYS,RIGHT,678,None,True,role_top=727,row_h=10.55)

# One compact legend for all three heatmap groups.
scale_y=807
for j in range(100):
    rgb=palette(-2.5+5*j/99)
    rect(48+j,scale_y,1.1,5,Color(*[v/255 for v in rgb]))
text(98,795,'Within-cohort row z-score',7,align='center')
for xx,label in [(48,'-2.5'),(98,'0'),(148,'2.5')]:text(xx,816,label,7,align='center')
rect(209,805,5,5,PALETTE['emt_states']['epithelial']);text(219,804,'Epithelial-associated genes',7.2)
rect(407,805,5,5,PALETTE['emt_states']['mesenchymal']);text(417,804,'Mesenchymal-associated genes',7.2)
text(401,819,'Ref., reference group; PF, portal fibroblast; dots, all-zero RNA rows.',7,align='center',color=MUTED)
c.showPage();c.save()

pd.DataFrame(tile_records).to_csv(AUDIT/'heatmap_tile_geometry.csv',index=False)
pd.DataFrame(tile_records)[['panel','cohort','feature','role','value']].to_csv(AUDIT/'plotted_values.csv',index=False)
(AUDIT/'text_geometry.json').write_text(json.dumps(text_records,indent=2)+'\n')
(AUDIT/'figure1_scatter_geometry.json').write_text(json.dumps(scatter_records,indent=2)+'\n')
outside=[r for r in text_records if r['bbox'][0]<0 or r['bbox'][1]<0 or r['bbox'][2]>W or r['bbox'][3]>H]
assert not outside, outside
assert len(tile_records)==928
for panel,n in [('B',512),('C',352),('D',64)]:
    assert sum(r['panel']==panel for r in tile_records)==n
reader=PdfReader(pdf_path)
extracted=reader.pages[0].extract_text()
assert len(reader.pages)==1
for forbidden in ['China','USA','metabolism','Anabolic','Catabolic','EMT-high']:
    assert forbidden not in extracted, forbidden
for wanted in ['GSE149614','GSE151530','GSE138709','PI3K-AKT','qHSC','GSVA']:
    assert wanted in extracted,wanted
assert all(extracted.splitlines().count(letter)==1 for letter in 'ABCD')
(AUDIT/'figure_text.txt').write_text(extracted)
# PyMuPDF is installed in system Python; authoring uses the bundled runtime.
check = json.loads(subprocess.check_output(['/usr/bin/python3',str(Path(__file__).with_name('validate_pdf.py')),str(pdf_path)],text=True))
assert set(check['type_sizes_pt']) == set(TYPE_SIZES.values())
assert not check['physical_out_of_bounds_text']
assert not check['text_bbox_collisions'], check['text_bbox_collisions']
report={'figure':str(pdf_path),'sha256':digest(pdf_path),'global_cells':95022,'global_clusters':40,
        'stromal_cells':3018,'selected_stromal_cells':757,'matrix_values':928,
        'analysis_rerun_in_this_revision':False,'display_filter':'HCC-GSE149614 reference cholangiocytes omitted in B, C and D; original z-scores retained','requested_PI3K_AKT_source':'Original custom 11-gene INSULIN-AKT, display label only; original 705/454/369 cells','max_zscore_reconstruction_error':max(errors.values()),
        'verification_by_input':errors,'qHSC_color_limits':[qlo,qhi],
        'page_points':[A4W,A4H],'page_mm':[210,297], 'out_of_bounds_text':outside,'panels':list('ABCD'),
        'fonts':check['fonts'], 'type_sizes_pt':TYPE_SIZES,
        'shared_cell_type_colors':COLORS, 'diverging_palette':[PALETTE['continuous'][n] for n in ('low','mid','high')],
        'journal_reference_palette':PALETTE['name'],
        'major_cell_type_colors':MAJOR,'cluster_colors':{str(k):v for k,v in cluster_colors.items()},
        'text_colors':{'main':INK, 'secondary':MUTED},
        'physical_out_of_bounds_text':check['physical_out_of_bounds_text'], 'text_bbox_collisions':check['text_bbox_collisions'],
        'signaling_gene_sets':'Original custom TGF, 11-gene INSULIN-AKT, HGF-AKT and EGF-RAS, with display labels TGF, PI3K-AKT, HGF and EGF',
        'source_files_unchanged':all(digest(Path(m['source']))==m['sha256'] for m in manifest)}
(AUDIT/'verification.json').write_text(json.dumps(report,indent=2)+'\n')
subprocess.run(['pdftoppm','-r','400','-png','-singlefile',str(pdf_path),str(OUTPUT/pdf_path.stem)],check=True)
subprocess.run(['pdftoppm','-scale-to','1800','-png','-singlefile',str(pdf_path),str(AUDIT/'figure_preview')],check=True)
print(json.dumps(report,indent=2))
