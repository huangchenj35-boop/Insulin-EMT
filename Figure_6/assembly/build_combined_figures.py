#!/usr/bin/env python3
"""Render complete curated pathway scores; retain the original Figure 6 A-C."""
from pathlib import Path
from io import BytesIO
import hashlib, json, subprocess, sys
import fitz
import numpy as np
import pandas as pd
from reportlab.pdfgen import canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

WORK = Path(__file__).resolve().parents[1]
ROOT = WORK.parent
OUTDIR = WORK / 'final_combination'
AUD = OUTDIR / 'audit'
SOURCE = ROOT / 'V16_Insulin与PI3K分开_20261026/Figure_6_Insulin_third.pdf'
OUT = OUTDIR / 'Figure_6_combined.pdf'
COMPARE = OUTDIR / 'Figure_6D_combined_comparison.pdf'
EXPECTED_SOURCE = 'a935806bb123376939ba28aef6ab8507f5fa30a72ceb1c06e0844ce2ea730d1b'
ROWS = [r['label'] for r in json.loads((AUD/'final_combination_selection.json').read_text())['selected']]
OLDROWS = ['TGF-β', 'PI3K-AKT', 'Insulin', 'HGF', 'EGF', 'IGF']
COHORTS = ['ChinaHCC', 'USAHCC', 'CCAcombined']
TITLES = ['HCC | GSE149614', 'HCC | GSE151530', 'CCA | GSE138709 + GSE151530']
ROLES = [
    ['Normal hepatocyte', 'Malignant epithelial', 'qHSC', 'Activated HSC', 'Portal fibroblast', 'Fibrotic MFB'],
    ['Normal hepatocyte', 'Malignant epithelial', 'qHSC', 'Portal fibroblast-like', 'Fibrotic MFB'],
    ['Normal cholangiocyte', 'Malignant epithelial', 'qHSC', 'Portal fibroblast-like', 'Fibrotic MFB']]
SHORT_ROLES = [['Hepatocyte','HCC cell','qHSC','Activated HSC','PF','MFB'],
               ['Hepatocyte','HCC cell','qHSC','PF-like','MFB'],
               ['Cholangiocyte','CCA cell','qHSC','PF-like','MFB']]
XSTARTS = [368.64306640625, 448.04306640625, 515.54306640625]
CW, Y0, RH = 11.9, 671.89389, 9.4
sha = lambda p: hashlib.sha256(p.read_bytes()).hexdigest()
pal = json.loads((ROOT / '配色参考主刊_20261010/palette.json').read_text())
STOPS = np.array([[int(pal['continuous'][n][j:j+2],16) for j in (1,3,5)] for n in ('low','mid','high')],float)

def color(value):
    t = np.clip((value + 2.5)/5, 0, 1)
    rgb = STOPS[0] + t*2*(STOPS[1]-STOPS[0]) if t <= .5 else STOPS[1]+(t-.5)*2*(STOPS[2]-STOPS[1])
    return tuple(int(v) for v in np.clip(np.rint(rgb),0,255))

def load_frames(folder, rows):
    frames = {k: pd.read_csv(folder / f'fig6_{k}_signaling_GSVA_activity_6pathways_zscore_displayed.csv',index_col=0) for k in COHORTS}
    for k, rr in zip(COHORTS,ROLES):
        assert list(frames[k].index) == rows and list(frames[k].columns) == rr
        assert np.isfinite(frames[k].values).all()
    return frames

def register_fonts():
    pdfmetrics.registerFont(TTFont('Arial','/System/Library/Fonts/Supplemental/Arial.ttf'))
    pdfmetrics.registerFont(TTFont('Arial-Bold','/System/Library/Fonts/Supplemental/Arial Bold.ttf'))

def build_main(frames):
    source = fitz.open(SOURCE)
    w,h = source[0].rect.width,source[0].rect.height
    out = fitz.open(); page = out.new_page(width=w,height=h)
    # Source vector art retained outside the replacement data rows and labels.
    for coords in [(0,0,w,246),(0,246,305,h),(305,246,w,610),(305,790,w,h),
                   (305,633,w,Y0),(367,Y0+6*RH,w,790)]:
        r = fitz.Rect(coords); page.show_pdf_page(r,source,0,clip=r)
    buf=BytesIO(); c=canvas.Canvas(buf,pagesize=(w,h),pageCompression=1)
    c.setFillColorRGB(.149,.149,.149)
    c.setFont('Arial-Bold',12);c.drawString(309.31,h-629.7,'D')
    c.setFont('Arial-Bold',10);c.drawString(328.76,h-629.12,'Signaling programs (RNA-based GSVA)')
    c.setFont('Arial',7)
    for i,name in enumerate(ROWS):
        c.setFillColorRGB(.149,.149,.149)
        c.drawRightString(361.84,h-(Y0+(i+.5)*RH+1.96),name)
        for key,x,rr in zip(COHORTS,XSTARTS,ROLES):
            for j,role in enumerate(rr):
                c.setFillColorRGB(*(v/255 for v in color(frames[key].loc[name,role])))
                c.setStrokeColorRGB(1,1,1); c.setLineWidth(.243168)
                c.rect(x+j*CW,h-(Y0+(i+1)*RH),CW,RH,stroke=1,fill=1)
    c.showPage(); c.save()
    overlay=fitz.open(stream=buf.getvalue(),filetype='pdf'); page.show_pdf_page(page.rect,overlay,0)
    out.set_metadata({'title':'Figure 6 | Original-based exploratory signaling combination',
                      'subject':'Four original V16 rows retained; two replacements: complete Reactome R-HSA-1257604 and KEGG Medicus HGF-MET-RAS-ERK. Post-score exploratory selection; unchanged RNA GSVA pipeline; insulin third; original panels A-C retained'})
    OUT.write_bytes(out.tobytes(garbage=4,deflate=True))
    out.close(); overlay.close(); source.close()


def draw_candidate_page(c):
    w,h=841.8898,595.2756
    def txt(x,y,t,size=9,bold=False):
        c.setFillColorRGB(.15,.15,.15);c.setFont('Arial-Bold' if bold else 'Arial',size);c.drawString(x,h-y,t)
    txt(35,32,'Figure 6D | 20 new candidates + prior Reactome reference',18,True)
    txt(35,51,'Complete curated sets; unchanged scoring. Asterisks mark the two replacement rows; four original rows stay unchanged.',10)
    topics=['TGF-beta receptor signaling','generic PI3K-AKT downstream signaling','insulin receptor signaling','HGF-MET','EGF/EGFR','IGF1/IGF1R']
    defs=json.loads((WORK/'inputs/candidate_specifications.json').read_text())
    defs.append({'standard_name':'R-HSA-1257604','topic':'generic PI3K-AKT downstream signaling','source':'Reactome','expected_members':284})
    ordered=sorted(enumerate(defs),key=lambda z:(topics.index(z[1]['topic']),z[0]));defs=[d for _,d in ordered]
    selected={x['standard_name'] for x in json.loads((AUD/'final_combination_selection.json').read_text())['selected']}
    names={
      'R-HSA-1257604':'PI3K-AKT | Reactome PIP3-AKT',
      'PID_TGFBR_PATHWAY':'TGF-beta | PID receptor network',
      'WP_TGFBETA_RECEPTOR_SIGNALING':'TGF-beta | WikiPathways receptor',
      'PID_PI3KCI_AKT_PATHWAY':'PI3K-AKT | PID AKT downstream',
      'WP_PI3KAKT_SIGNALING':'PI3K-AKT | WikiPathways broad',
      'PID_INSULIN_PATHWAY':'Insulin | PID receptor network',
      'WP_INSULIN_SIGNALING':'Insulin | WikiPathways broad',
      'PID_MET_PATHWAY':'HGF | PID MET network',
      'WP_HEPATOCYTE_GROWTH_FACTOR_RECEPTOR_SIGNALING':'HGF | WikiPathways receptor',
      'PID_ERBB1_RECEPTOR_PROXIMAL_PATHWAY':'EGF | PID receptor-proximal',
      'PID_ERBB1_DOWNSTREAM_PATHWAY':'EGF | PID broad downstream',
      'WP_EGFEGFR_SIGNALING':'EGF | WikiPathways broad',
      'PID_IGF1_PATHWAY':'IGF1 | PID receptor network'}
    med={'TGF_BETA_SIGNALING_PATHWAY':'TGF-beta | Medicus SMAD core','INS_AKT_SIGNALING_PATHWAY':'Insulin | Medicus IRS-PI3K','HGF_MET_PI3K_SIGNALING_PATHWAY':'HGF | Medicus PI3K-AKT-BAD','HGF_MET_RAS_ERK_SIGNALING_PATHWAY':'HGF | Medicus RAS-ERK','EGF_EGFR_PI3K_SIGNALING_PATHWAY':'EGF | Medicus PI3K-AKT-BAD','EGF_EGFR_RAS_ERK_SIGNALING_PATHWAY':'EGF | Medicus RAS-ERK','IGF_IGF1R_PI3K_SIGNALING_PATHWAY':'IGF1 | Medicus PI3K-AKT-mTOR','IGF_IGF1R_RAS_ERK_SIGNALING_PATHWAY':'IGF | Medicus RAS-ERK'}
    names.update({'KEGG_MEDICUS_REFERENCE_'+k:v for k,v in med.items()})
    frames={k:pd.read_csv(WORK/'tables'/f'{k}_all_candidate_role_zscores.csv',index_col=0) for k in COHORTS}
    for k in COHORTS:
        r=pd.read_csv(OUTDIR/'tables'/f'fig6_{k}_signaling_GSVA_activity_6pathways_zscore.csv',index_col=0).loc['PI3K-AKT']
        frames[k].loc['R-HSA-1257604']=r
    top,rh,cw=106,15.2,27;starts=[276,456,609];drawn=[]
    for key,title,x,rr in zip(COHORTS,TITLES,starts,ROLES):
        c.setFont('Arial-Bold',9);c.setFillColorRGB(.15,.15,.15);c.drawCentredString(x+len(rr)*cw/2,h-(top-12),title)
    previous=None
    for i,d in enumerate(defs):
        name=d['standard_name'];y=top+i*rh
        if previous and previous!=d['topic']:
            c.setStrokeColorRGB(.45,.45,.45);c.setLineWidth(.6);c.line(35,h-y,744,h-y)
        label=(' * ' if name in selected else '    ')+names[name]+f" ({d['expected_members']})"
        txt(35,y+rh*.5+2.8,label,8.4,name in selected)
        for key,x,rr in zip(COHORTS,starts,ROLES):
            for j,role in enumerate(rr):
                value=float(frames[key].loc[name,role]);rgb=color(value)
                c.setFillColorRGB(*(v/255 for v in rgb));c.setStrokeColorRGB(1,1,1);c.setLineWidth(.5)
                c.rect(x+j*cw,h-y-rh,cw,rh,fill=1,stroke=1)
                c.setFont('Arial',7.2);c.setFillColorRGB(.16,.16,.16)
                c.drawCentredString(x+(j+.5)*cw,h-y-rh*.5-2.4,f'{value:.2f}')
                drawn.append({'standard_name':name,'cohort':key,'role':role,'zscore':value,'sample_x':x+(j+.12)*cw,'sample_y':y+rh*.16,'expected_rgb':rgb})
        previous=d['topic']
    for x,short in zip(starts,SHORT_ROLES):
        for j,role in enumerate(short):
            c.saveState();c.translate(x+(j+.57)*cw,h-(top+len(defs)*rh+10));c.rotate(55)
            c.setFillColorRGB(.2,.2,.2);c.setFont('Arial',8.5);c.drawRightString(0,0,role);c.restoreState()
    txt(35,521,'Numbers in parentheses are complete gene-set sizes; detected member counts are provided in the source table.',9)
    txt(35,537,'All rows use the same -2.5 to +2.5 color scale. Compare cell groups within each row; row ranges represent different definitions.',9)
    txt(35,553,'Replacement rows were chosen after viewing results. This panel retains all alternatives, including discordant results.',9)
    txt(35,569,'China HCC standardization retains the hidden normal-cholangiocyte reference, matching the original pipeline.',9)
    (AUD/'candidate_plotted_cells.json').write_text(json.dumps(drawn,ensure_ascii=False,indent=2)+'\n')


def build_comparison(old,new):
    w,h = 841.8898,595.2756
    c=canvas.Canvas(str(COMPARE),pagesize=(w,h),pageCompression=1)
    c.setTitle('Figure 6D | Original and combined module definitions')
    def txt(x,y,s,size=10,bold=False):
        c.setFillColorRGB(.15,.15,.15); c.setFont('Arial-Bold' if bold else 'Arial',size); c.drawString(x,h-y,s)
    txt(35,32,'Figure 6D | Combination on the original figure',18,True)
    txt(35,51,'Same cells, normalization, GSVA settings and within-cohort row standardization',10)
    # Shared, fixed color scale, identical to the main figure.
    lx,ly,lw=655,30,145
    for j in range(145):
        c.setFillColorRGB(*(v/255 for v in color(-2.5+5*j/144)))
        c.rect(lx+j,h-ly-8,1,8,fill=1,stroke=0)
    txt(lx,ly+21,'-2.5',8); txt(lx+lw/2-2,ly+21,'0',8); txt(lx+lw-12,ly+21,'2.5',8)
    txt(lx,ly-5,'Row z-score',9)
    starts=[139,395,609]; cw=34; rh=20
    drawn=[]
    for version,frames,rows,top,label in [
        ('V16',old,OLDROWS,109,'A   Original definitions | custom sets + KEGG + Reactome'),
        ('Combined',new,ROWS,344,'B   Combined version | four original rows retained; PI3K-AKT and HGF replaced')]:
        txt(35,top-31,label,12,True)
        for key,title,x,rr,short in zip(COHORTS,TITLES,starts,ROLES,SHORT_ROLES):
            c.setFont('Arial-Bold',9.5); c.setFillColorRGB(.15,.15,.15)
            c.drawCentredString(x+len(rr)*cw/2,h-(top-9),title)
            for i,name in enumerate(rows):
                if x==starts[0]:
                    c.setFont('Arial',9.5); c.setFillColorRGB(.15,.15,.15)
                    c.drawRightString(x-10,h-(top+(i+.5)*rh+3.3),name)
                for j,role in enumerate(rr):
                    value=float(frames[key].loc[name,role]);rgb=color(value)
                    c.setFillColorRGB(*(v/255 for v in rgb));c.setStrokeColorRGB(1,1,1);c.setLineWidth(.65)
                    c.rect(x+j*cw,h-(top+(i+1)*rh),cw,rh,fill=1,stroke=1)
                    c.setFont('Arial',8);c.setFillColorRGB(.16,.16,.16)
                    c.drawCentredString(x+(j+.5)*cw,h-(top+(i+.5)*rh+2.65),f'{value:.2f}')
                    drawn.append({'version':version,'cohort':key,'pathway':name,'role':role,'zscore':value,
                                  'sample_x':x+(j+.12)*cw,'sample_y':top+(i+.16)*rh,'expected_rgb':rgb})
            for j,s in enumerate(short):
                c.saveState();c.translate(x+(j+.57)*cw,h-(top+6*rh+10));c.rotate(55)
                c.setFillColorRGB(.2,.2,.2);c.setFont('Arial',9);c.drawRightString(0,0,s);c.restoreState()
    txt(35,550,'Four original rows are unchanged: TGF-beta, Insulin, EGF and IGF. Insulin remains third.',9)
    txt(35,565,'Replacements: Reactome R-HSA-1257604 for PI3K-AKT; complete KEGG Medicus HGF-MET-RAS-ERK branch for HGF.',9)
    txt(35,580,'Exploratory selection after candidate comparison. Higher than the stromal mean does not mean higher than every stromal subtype.',8.5)
    c.showPage();draw_candidate_page(c);c.showPage();c.save()
    (AUD/'comparison_plotted_cells.json').write_text(json.dumps(drawn,ensure_ascii=False,indent=2)+'\n')

def verify():
    assert sha(SOURCE)==EXPECTED_SOURCE
    frames=load_frames(OUTDIR/'tables',ROWS)
    old=fitz.open(SOURCE);new=fitz.open(OUT);p=new[0];w,h=p.rect.width,p.rect.height
    assert len(new)==1
    identical=[]
    for name,coords in {'A':(0,0,w,246),'B':(0,246,305,h),'C':(305,246,w,610),'footer':(305,790,w,h)}.items():
        r=fitz.Rect(coords)
        a=old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=r,alpha=False)
        b=p.get_pixmap(matrix=fitz.Matrix(2,2),clip=r,alpha=False)
        assert (a.width,a.height)==(b.width,b.height) and a.samples==b.samples,name
        identical.append(name)
    spans=[s for b in p.get_text('dict')['blocks'] for line in b.get('lines',[]) for s in line['spans']
           if s['text'] in ROWS and 305<=s['bbox'][0]<362 and Y0<=s['bbox'][1]<Y0+6*RH]
    assert [s['text'] for s in sorted(spans,key=lambda s:s['bbox'][1])]==ROWS,spans
    assert min(s['size'] for s in spans)>=6.99
    pix=p.get_pixmap(matrix=fitz.Matrix(4,4),alpha=False);records=[]
    for key,x,rr in zip(COHORTS,XSTARTS,ROLES):
        for i,name in enumerate(ROWS):
            for j,role in enumerate(rr):
                value=float(frames[key].loc[name,role]);expected=color(value)
                actual=pix.pixel(int((x+(j+.5)*CW)*4),int((Y0+(i+.5)*RH)*4))
                assert max(abs(a-b) for a,b in zip(expected,actual))<=1,(key,name,role,expected,actual)
                records.append({'cohort':key,'pathway':name,'role':role,'zscore':value,'expected_rgb':expected,'rendered_rgb':actual})
    assert len(records)==96
    oldpix=old[0].get_pixmap(matrix=fitz.Matrix(4,4),alpha=False);unchanged_cells=0
    for i in [0,2,4,5]:
        for x,rr in zip(XSTARTS,ROLES):
            for j,_ in enumerate(rr):
                xx,yy=int((x+(j+.5)*CW)*4),int((Y0+(i+.5)*RH)*4)
                assert max(abs(a-b) for a,b in zip(oldpix.pixel(xx,yy),pix.pixel(xx,yy)))<=1
                unchanged_cells+=1
    assert unchanged_cells==64
    p.get_pixmap(dpi=400,alpha=False).save(OUT.with_suffix('.png'))
    p.get_pixmap(dpi=144,alpha=False).save(AUD/'Figure_6_preview.png')
    p.get_pixmap(matrix=fitz.Matrix(4,4),clip=fitz.Rect(305,610,w,790),alpha=False).save(OUTDIR/'Figure_6D_combined.png')
    comp=fitz.open(COMPARE);assert len(comp)==2
    cp=comp[0];cpix=cp.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    compare_cells=json.loads((AUD/'comparison_plotted_cells.json').read_text())
    assert len(compare_cells)==192
    for r in compare_cells:
        actual=cpix.pixel(int(r['sample_x']*2),int(r['sample_y']*2))
        assert max(abs(a-b) for a,b in zip(actual,r['expected_rgb']))<=1,r
    cpix.save(OUTDIR/'Figure_6D_combined_comparison.png')
    candpix=comp[1].get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    candcells=json.loads((AUD/'candidate_plotted_cells.json').read_text());assert len(candcells)==336
    for r in candcells:
        actual=candpix.pixel(int(r['sample_x']*2),int(r['sample_y']*2))
        assert max(abs(a-b) for a,b in zip(actual,r['expected_rgb']))<=1,r
    candpix.save(OUTDIR/'Figure_6D_all_candidates.png')
    report={'status':'verified','unchanged_original_row_cells_checked':unchanged_cells,'candidate_cells_checked':len(candcells),'original_source_sha256':sha(SOURCE),'updated_figure_sha256':sha(OUT),
            'comparison_pdf_sha256':sha(COMPARE),'unchanged_regions_pixel_identical_144dpi':identical,
            'row_order':ROWS,'main_heatmap_cells_checked':len(records),'comparison_cells_checked':len(compare_cells),
            'rgb_tolerance':1,'all_cell_colors_match_numeric_scores':True,'source_unchanged':True}
    (AUD/'figure_verification.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
    (AUD/'main_plotted_cells.json').write_text(json.dumps(records,ensure_ascii=False,indent=2)+'\n')
    print(json.dumps(report,ensure_ascii=False,indent=2))

if __name__=='__main__':
    if '--verify' in sys.argv: verify()
    else:
        assert sha(SOURCE)==EXPECTED_SOURCE
        OUTDIR.mkdir(exist_ok=True);AUD.mkdir(exist_ok=True)
        calc=json.loads((AUD/'calculation_notes.json').read_text()); assert calc['status']=='completed'
        new=load_frames(OUTDIR/'tables',ROWS)
        old=load_frames(ROOT/'V16_Insulin与PI3K分开_20261026/tables',OLDROWS)
        register_fonts();build_main(new);build_comparison(old,new)
        subprocess.run([sys.executable,__file__,'--verify'],check=True)
