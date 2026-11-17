#!/usr/bin/env python3
"""Compact A4 landscape layout; replay frozen scientific artwork and exact values."""
from pathlib import Path
import hashlib,json,math,shutil
from collections import Counter
import fitz
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor,Color

ROOT=Path(__file__).resolve().parents[2]
WORK=ROOT/'V10_按正文修图_20261024';AUD=WORK/'audit';OUT=WORK/'figures'
AUD.mkdir(parents=True,exist_ok=True)
BACK=WORK/'before';BACK.mkdir(exist_ok=True)
for name in []:
    if not (BACK/name).exists():shutil.copy2(OUT/name,BACK/name)
W,H=841.8898,595.2756
PAL=json.loads((ROOT/'配色参考主刊_20261010/palette.json').read_text())
INK='#262626';MUTED='#555555'
for n,f in [('Arial','Arial.ttf'),('Arial-Bold','Arial Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(n,str(Path('/System/Library/Fonts/Supplemental')/f)))
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
manifest=[];geometry=[];heat_records=[]
def read(p):manifest.append({'path':str(p),'sha256':sha(p)});return pd.read_csv(p)
def fill(col):return HexColor(col) if isinstance(col,str) else col
c=canvas.Canvas(str(AUD/'S1_text.pdf'),pagesize=(W,H))
def text(s,x,y,size=8,bold=False,color=INK,align='left',rot=0):
    s=str(s);face='Arial-Bold' if bold else 'Arial';tw=pdfmetrics.stringWidth(s,face,size)
    c.saveState();c.setFillColor(fill(color));c.setFont(face,size);c.translate(x,H-y);c.rotate(rot)
    c.drawString(-{'left':0,'center':.5,'right':1}[align]*tw,-size*.78,s);c.restoreState()
    geometry.append(dict(text=s,x=x,y=y,size=size,rot=rot))
def rect(x,y,w,h,col,edge=None):
    c.setFillColor(fill(col));c.setStrokeColor(fill(edge or col));c.setLineWidth(.35)
    c.rect(x,H-y-h,w,h,fill=1,stroke=int(edge is not None))
def title(letter,s,y):text(letter,18,y,12,True);text(s,37,y+1,10,True)
def blend(a,b,t):
    a=np.array([int(a[j:j+2],16) for j in (1,3,5)])/255
    b=np.array([int(b[j:j+2],16) for j in (1,3,5)])/255
    return Color(*(a+(b-a)*np.clip(t,0,1)))

title('B','Malignant epithelial expression states',276)
source_path=ROOT/'JBC_singlecell_20261014_v2/Figure_2.pdf'
source=fitz.open(source_path);manifest.append({'path':str(source_path),'sha256':sha(source_path)})
p=source[0]
spans=[dict(s,direction=l['dir']) for b in p.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
cohorts=['HCC-GSE149614','HCC-GSE151530','CCA-GSE138709']
art=[];a_records=[];scale=1.50;atop=343
for i,left in enumerate([9,286,563]):
    clip=fitz.Rect(10+193*i,52,203+193*i,185)
    dest=fitz.Rect(left,atop,left+193*scale,atop+133*scale);art.append((clip,dest))
    center=left+90.75*scale
    text(cohorts[i],center,303,8,True,align='center')
    if i==2:text('+ GSE151530',center,314,8,True,align='center')
    text(f'{["13,430","10,914","5,304"][i]} {"HCC cells" if i<2 else "CCA cells"}',center,328,8,align='center')
    for s in spans:
        bb=fitz.Rect(s['bbox'])
        if not fitz.Rect(10+193*i,25,203+193*i,185).contains(bb):continue
        a_records.append({'cohort':cohorts[i],'text':s['text']})
        if 'malignant epithelial cells' in s['text'] or s['size']>=7.9 or s['text']=='GSE151530' or s['text'].startswith('+ GSE'):continue
        if s['text']=='Hallmark EMT score':
            text(s['text'],left+1,atop+102,8,align='center',rot=90);continue
        if s['text']=='Epithelial-associated score':
            text(s['text'],center,529,8,align='center');continue
        x=left+((bb.x0+bb.x1)/2-clip.x0)*scale
        y=atop+((bb.y0+bb.y1)/2-clip.y0)*scale
        text(s['text'],x,y-4,8,'Bold' in s['font'],align='center')

title('A','Epithelial states in paired HCC tissues',16)
text('HCC-GSE149614 | 8 paired patients',W/2,39,8,True,align='center')
paired=read(ROOT/'S1A_HCC配对热图_FigureYa311_20261010/inputs/patient_percentages.csv')
stats=read(ROOT/'S1A_HCC配对热图_FigureYa311_20261010/inputs/paired_statistics.csv')
pats=list(paired.patient.drop_duplicates());cell=42;top=110
for metric,key,left in [('E/M high','EM_high',80),('M-high total','M_high',488)]:
    center=left+4*cell
    text(metric,center,62,8,True,align='center')
    text('Both scores high' if key=='EM_high' else 'Mesenchymal high + E/M high',center,76,8,align='center')
    text('Patient',left-8,97,8,align='right')
    for j,patient in enumerate(pats):text(str(patient).replace('HCC',''),left+(j+.5)*cell,97,8,align='center')
    for row,pop in enumerate(['Non-tumor','Malignant']):
        y=top+row*cell;mid=y+cell/2
        if row==0:
            text('Non-tumor',left-8,mid-13,8,align='right')
            text('hepatocyte',left-8,mid-3,8,align='right')
            text('n = 1,575',left-8,mid+8,7,align='right')
        else:
            text('HCC cell',left-8,mid-7,8,align='right')
            text('n = 11,153',left-8,mid+5,7,align='right')
        rect(left-4,y,2.5,cell,PAL['cell_types']['hepatocyte' if row==0 else 'malignant'])
        for j,patient in enumerate(pats):
            x=left+j*cell
            rr=paired[(paired.metric==key)&(paired.patient==patient)&(paired.population==pop)].iloc[0]
            v=float(rr.percentage);assert abs(v-100*rr.positive_cells/rr.cells)<1e-10
            rect(x,y,cell,cell,blend('#FFFFFF',PAL['continuous']['high'],v/100),'#FFFFFF')
            text(f'{v:.1f}%',x+cell/2,y+17,8,align='center')
            heat_records.append(dict(metric=key,patient=patient,population=pop,percentage=v,x=x,y=y,width=cell,height=cell))
    pval=float(stats.loc[stats.metric.eq(key),'p_holm'].iloc[0])
    text(f'Paired P (Holm) = {pval:.3f}',center,210,8,align='center')

text('Cells (%)',342,239,8,align='right')
for j in range(150):rect(352+j,240,1.1,5,blend('#FFFFFF',PAL['continuous']['high'],j/149))
for v in [0,50,100]:text(str(v),352+1.5*v,251,8,align='center')
text('A: Percentages are calculated within each patient and cell population.',18,554,8)
text('B: Dashed lines indicate pooled malignant-cell medians; percentages use all malignant cells in each cohort.',18,568,8)
c.showPage();c.save()

# Replay the original, complete vector artwork, removing only source text.
source[0].add_redact_annot(source[0].rect,fill=False)
source[0].apply_redactions(images=0,graphics=0,text=0)
assert not source[0].get_text().strip()
out=fitz.open(AUD/'S1_text.pdf')
for clip,dest in art:out[0].show_pdf_page(dest,source,0,clip=clip,overlay=False)
out.set_metadata({'title':'Figure S1 | Epithelial expression states'})
(OUT/'Figure_S1.pdf').write_bytes(out.tobytes(garbage=4,deflate=True));out.close();source.close()
pd.DataFrame(heat_records).to_csv(AUD/'S1A_percentages.csv',index=False)
pd.DataFrame(a_records).to_csv(AUD/'S1B_labels.csv',index=False)
(AUD/'source_manifest.json').write_text(json.dumps(manifest,indent=2,ensure_ascii=False))
(AUD/'text_geometry.json').write_text(json.dumps(geometry,indent=2))
(AUD/'art_geometry.json').write_text(json.dumps([{'source_clip':list(a),'destination':list(b),'scale':scale} for a,b in art],indent=2))
pdf=fitz.open(OUT/'Figure_S1.pdf');page=pdf[0]
assert len(pdf)==1 and abs(page.rect.width-W)<.01 and abs(page.rect.height-H)<.01
old=pd.read_csv(ROOT/'V6_最终图_20261024/audit/S1B_percentages.csv')
key=['metric','patient','population'];cols=key+['percentage']
pd.testing.assert_frame_equal(pd.DataFrame(heat_records)[cols].sort_values(key).reset_index(drop=True),old[cols].sort_values(key).reset_index(drop=True))
assert len(heat_records)==32 and all(r['width']==r['height']==42 for r in heat_records)
assert all(sha(Path(r['path']))==r['sha256'] for r in manifest)
sp=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
collisions=[]
for i,a in enumerate(sp):
    for b in sp[i+1:]:
        inter=fitz.Rect(a['bbox'])&fitz.Rect(b['bbox'])
        if not inter.is_empty and inter.width>.4 and inter.height>.4:collisions.append({'a':a['text'],'b':b['text'],'rect':list(inter)})
report={'A4_landscape':True,'Arial_only':all('Arial' in s['font'] for s in sp),
        'font_sizes':dict(Counter(round(s['size'],2) for s in sp)),
        'outside_page':[s['text'] for s in sp if not page.rect.contains(fitz.Rect(s['bbox']))],
        'text_collisions':collisions,'unchanged_percentage_records':32,'square_cells_pt':42,
        'new_pdf_sha256':sha(OUT/'Figure_S1.pdf')}
(AUD/'verification.json').write_text(json.dumps(report,indent=2))
# Render before assertions to facilitate visual review of any genuine collision.
page.get_pixmap(dpi=120,alpha=False).save(AUD/'S1_preview.png')
page.get_pixmap(dpi=400,alpha=False).save(OUT/'Figure_S1.png')
assert not report['outside_page'] and not collisions,report
assert report['Arial_only']
print(json.dumps(report,indent=2))
