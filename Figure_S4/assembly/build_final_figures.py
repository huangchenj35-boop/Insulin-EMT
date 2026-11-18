#!/usr/bin/env python3
"""Relayout S1 only; preserve the existing percentage data and the complete B panel."""
from pathlib import Path
import csv, hashlib, json, shutil
from collections import Counter
import fitz
import numpy as np
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.lib.colors import HexColor, Color

ROOT = Path(__file__).resolve().parents[2]
WORK = Path(__file__).resolve().parents[1]
OUT, AUD = WORK/'figures', WORK/'audit'
OUT.mkdir(parents=True, exist_ok=True)
AUD.mkdir(parents=True, exist_ok=True)
SOURCE = ROOT/'JBC_singlecell_20261024_final'
W, H = 841.8898, 528
PAL = json.loads((ROOT/'配色参考主刊_20261010/palette.json').read_text())
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
inputs = [SOURCE/f'Figure_{n}.{ext}' for n in ['6','7','S1'] for ext in ['pdf','png']]
data_path = ROOT/'S1A_HCC配对热图_FigureYa311_20261010/inputs/patient_percentages.csv'
stats_path = data_path.with_name('paired_statistics.csv')
# Non-figure manuscript hash omitted in public release (not an analysis input).
inputs += [data_path, stats_path]
before = {str(p):sha(p) for p in inputs}
for n in ['6','7']:
    for ext in ['pdf','png']:
        shutil.copy2(SOURCE/f'Figure_{n}.{ext}', OUT/f'Figure_{n}.{ext}')

for name, filename in [('Arial','Arial.ttf'),('Arial-Bold','Arial Bold.ttf')]:
    pdfmetrics.registerFont(TTFont(name,str(Path('/System/Library/Fonts/Supplemental')/filename)))
c = canvas.Canvas(str(AUD/'S1_labels.pdf'), pagesize=(W,H))
def text(s,x,y,size=9,bold=False,align='left'):
    face='Arial-Bold' if bold else 'Arial'
    c.setFont(face,size);c.setFillColor(HexColor('#262626'))
    width=pdfmetrics.stringWidth(str(s),face,size)
    c.drawString(x-{'left':0,'center':.5,'right':1}[align]*width,H-y-size*.78,str(s))
def rect(x,y,w,h,col,edge=None):
    c.setFillColor(col if isinstance(col,Color) else HexColor(col))
    c.setLineWidth(.35)
    if edge:c.setStrokeColor(HexColor(edge))
    c.rect(x,H-y-h,w,h,fill=1,stroke=bool(edge))
def blend(t):
    low=np.ones(3)
    high=np.array([int(PAL['continuous']['high'][i:i+2],16) for i in (1,3,5)])/255
    return Color(*(low+(high-low)*t))
paired=list(csv.DictReader(data_path.open()))
stats={r['metric']:float(r['p_holm']) for r in csv.DictReader(stats_path.open())}
patients=list(dict.fromkeys(r['patient'] for r in paired))
assert len(patients)==8 and len(paired)==32
text('A',18,12,12,True)
text('Epithelial states in paired HCC tissues',37,13,11,True)
text('HCC-GSE149614 | 8 paired patients',W/2,34,9,True,'center')
records=[]
for metric,key,left,subtitle in [
    ('E/M high','EM_high',80,'Both scores high'),
    ('M-high total','M_high',488,'Mesenchymal high + E/M high')]:
    center=left+168
    headw=pdfmetrics.stringWidth(metric,'Arial-Bold',9.5)
    subw=pdfmetrics.stringWidth(subtitle,'Arial',8.5)
    start=center-(headw+10+subw)/2
    text(metric,start,54,9.5,True)
    text(subtitle,start+headw+10,55,8.5)
    text('Patient',left-8,72,8.5,align='right')
    for j,p in enumerate(patients):text(p.replace('HCC',''),left+(j+.5)*42,72,9,align='center')
    for i,pop in enumerate(['Non-tumor','Malignant']):
        y=84+i*42
        if i==0:
            text('Non-tumor',left-8,y+7,9,align='right')
            text('hepatocyte',left-8,y+18,9,align='right')
            text('n = 1,575',left-8,y+30,8,align='right')
        else:
            text('HCC cell',left-8,y+12,9,align='right')
            text('n = 11,153',left-8,y+25,8,align='right')
        rect(left-4,y,2.5,42,PAL['cell_types']['hepatocyte' if i==0 else 'malignant'])
        for j,p in enumerate(patients):
            rr=[r for r in paired if r['metric']==key and r['population']==pop and r['patient']==p]
            assert len(rr)==1
            r=rr[0];v=float(r['percentage'])
            assert abs(v-100*int(r['positive_cells'])/int(r['cells']))<1e-10
            rect(left+j*42,y,42,42,blend(v/100),'#FFFFFF')
            text(f'{v:.1f}%',left+(j+.5)*42,y+17,10,align='center')
            records.append({**r,'label':f'{v:.1f}%'})
    text(f'Paired P (Holm) = {stats[key]:.3f}',center,179,9,align='center')

# Compact shared scale retains its original tick order and endpoints.
text('Cells (%)',342,191,8.5,align='right')
for j in range(150):rect(352+j,191,1.1,5,blend(j/149))
for value in [0,50,100]:text(str(value),352+1.5*value,201,8,align='center')
text('A: Percentages are calculated within each patient and cell population.',18,495,8)
text('B: Dashed lines indicate pooled malignant-cell medians; percentages use all malignant cells in each cohort.',18,509,8)
c.showPage();c.save()

source=fitz.open(SOURCE/'Figure_S1.pdf')
out=fitz.open(AUD/'S1_labels.pdf')
clip=fitz.Rect(0,269,W,542)
dest=fitz.Rect(0,209,W,482)
out[0].show_pdf_page(dest,source,0,clip=clip,keep_proportion=True,overlay=False)
out.set_metadata({'title':'Figure S4 | Compact layout | 2026-10-26',
                  'subject':'Existing rectangular heatmaps; unchanged scientific data and B panel'})
# Serialize before writing and replace atomically, avoiding an in-place PDF write.
pending=OUT/'Figure_S1.pending.pdf'
pending.write_bytes(out.tobytes(garbage=4,deflate=True))
pending.replace(OUT/'Figure_S1.pdf')
out.close()
new=fitz.open(OUT/'Figure_S1.pdf')
new[0].get_pixmap(dpi=400,alpha=False).save(OUT/'Figure_S1.png')
new[0].get_pixmap(dpi=144,alpha=False).save(AUD/'S1_preview.png')
source[0].get_pixmap(dpi=100,alpha=False).save(AUD/'S1_before.png')
for n in ['6','7']:
    d=fitz.open(OUT/f'Figure_{n}.pdf');d[0].get_pixmap(dpi=110,alpha=False).save(AUD/f'Figure_{n}_preview.png');d.close()

# Compare B at the same raster scale after its integer-point translation.
orig_b=source[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=clip,alpha=False)
new_b=new[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=dest,alpha=False)
assert (orig_b.width,orig_b.height)==(new_b.width,new_b.height)
diff=np.abs(np.frombuffer(orig_b.samples,dtype=np.uint8).astype(int)-np.frombuffer(new_b.samples,dtype=np.uint8).astype(int))
old_text=source[0].get_text('text',clip=clip).split()
new_text=new[0].get_text('text',clip=dest).split()
assert Counter(old_text)==Counter(new_text)
spans=[s for b in new[0].get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
new_labels=Counter(s['text'] for s in spans if '%' in s['text'] and s['bbox'][1]<170)
assert new_labels==Counter(r['label'] for r in records)
collisions=[]
for i,a in enumerate(spans):
    for b in spans[i+1:]:
        inter=fitz.Rect(a['bbox'])&fitz.Rect(b['bbox'])
        if not inter.is_empty and inter.width>.4 and inter.height>.4:
            collisions.append([a['text'],b['text']])
outside=[s['text'] for s in spans if not new[0].rect.contains(fitz.Rect(s['bbox']))]
assert not outside and not collisions,(outside,collisions)
assert diff.max()==0, int(diff.max())
assert all(sha(Path(p))==h for p,h in before.items())
report={'source_sha256':before,'original_size_pt':list(source[0].rect)[2:],
        'new_size_pt':[W,H],'canvas_height_reduction_percent':100*(1-H/source[0].rect.height),
        'cell_label_size_pt_before':8,'cell_label_size_pt_after':10,
        'square_cell_size_pt':42,'unchanged_percentage_cells':32,
        'unchanged_p_holm':stats,'B_translation_y_pt':-60,'B_pixel_diff_max_144dpi':int(diff.max()),
        'B_words_unchanged':True,'text_collisions':collisions,'outside_page':outside,
        'figures_6_7_byte_identical':all(sha(OUT/f'Figure_{n}.{ext}')==before[str(SOURCE/f'Figure_{n}.{ext}')] for n in ['6','7'] for ext in ['pdf','png']),
        'source_files_unchanged':True,'S1_pdf_sha256':sha(OUT/'Figure_S1.pdf')}
(AUD/'figure_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2))
(AUD/'plotted_percentages.json').write_text(json.dumps(records,indent=2))
print(json.dumps({k:v for k,v in report.items() if k!='source_sha256'},ensure_ascii=False,indent=2))
