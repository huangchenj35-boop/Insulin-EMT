#!/usr/bin/env python3
"""Replace only S1A in a review copy; preserve B-D at original coordinates."""
from pathlib import Path
import hashlib, json, shutil, subprocess
import fitz
import numpy as np
import pandas as pd

WORK=Path(__file__).resolve().parents[1]
ROOT=WORK.parent
for name in ('inputs','audit','reference'):
    (WORK/name).mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
score_path=ROOT/'Figure2_细胞来源追踪/拼图_20261009/inputs/S5_Cook2020_two_axis_EMT_scores_by_cell.csv'
normal_path=ROOT.parent/'revision_07_normal_EMT/output/normal_epithelium_scores_GSE149614.csv'
source_pdf=ROOT/'当前主图汇总_20261010/Figure_S1.pdf'
palette=ROOT/'配色参考主刊_20261010/palette.json'
protected={str(p):sha(p) for p in (ROOT/'当前主图汇总_20261010').glob('*') if p.is_file()}
s=pd.read_csv(score_path); h=pd.read_csv(normal_path)
e=float(s.epithelial_score.median()); m=float(s.hallmark_EMT_score.median())
h=h[h.subtype.eq('Hepatocyte-like')].copy();h['population']='Non-tumor'
a=s[s.cohort.eq('China-HCC')].copy()
a['patient']=a.patient.str.replace('T$','',regex=True)
a['population']='Malignant';a=a[a.patient.isin(h.patient)]
assert set(a.patient)==set(h.patient)=={f'HCC{i:02d}' for i in range(3,11)}
assert h.barcode.is_unique and a.barcode.is_unique and not set(a.barcode)&set(h.barcode)
assert len(h)==1575 and len(a)==11153
d=pd.concat([h,a],ignore_index=True)
d['EM_high']=(d.epithelial_score>e)&(d.hallmark_EMT_score>m)
d['M_high']=d.hallmark_EMT_score>m
rows=[]
for (patient,pop),g in d.groupby(['patient','population']):
    for metric in ('EM_high','M_high'):
        rows.append(dict(patient=patient,population=pop,metric=metric,cells=len(g),
                         positive_cells=int(g[metric].sum()),percentage=float(100*g[metric].mean())))
p=pd.DataFrame(rows).sort_values(['metric','patient','population'])
assert len(p)==32 and not p[['patient','population','metric']].duplicated().any()
p.to_csv(WORK/'inputs/patient_percentages.csv',index=False)
d[['barcode','patient','population','epithelial_score','hallmark_EMT_score','EM_high','M_high']].to_csv(WORK/'audit/included_cells.csv',index=False)
definitions=dict(E_cut=e,M_cut=m,reference_cells=len(s),patients=sorted(d.patient.unique().tolist()),
                 normal_cells=len(h),malignant_cells=len(a),
                 analysis='Saved Figure 2A scoring; no rescoring, decontX or depth downsampling',
                 EM_high='E > E_cut and M > M_cut',M_high='M > M_cut; Mesenchymal high plus E/M high',
                 unit='patient',paired_test='two-sided exact Wilcoxon signed-rank; Holm adjustment across two metrics')
(WORK/'audit/definitions.json').write_text(json.dumps(definitions,indent=2))
template=WORK/'reference/FigureYa195PanPaire.Rmd'
if not template.exists():shutil.copyfile('/private/tmp/emt_FigureYa195PanPaire.Rmd',template)
manifest=[dict(path=str(f.resolve()),sha256=sha(f)) for f in (score_path,normal_path,palette,source_pdf,template)]
(WORK/'audit/source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
subprocess.run(['Rscript',str(WORK/'code/draw_paired_panel.R')],check=True)
point=pd.read_csv(WORK/'audit/plotted_points.csv')
q=point.merge(p,on=['metric','patient','population'],validate='one_to_one',suffixes=('_plot','_source'))
assert len(q)==32 and np.allclose(q.percentage_plot,q.percentage_source,atol=1e-12,rtol=0)
assert np.allclose(q.y,174-1.04*q.percentage_source,atol=1e-12,rtol=0)
old=fitz.open(source_pdf); layer=fitz.open(WORK/'audit/S1A_paired_layer.pdf')
W,H=old[0].rect.width,old[0].rect.height
out=fitz.open();page=out.new_page(width=W,height=H)
page.show_pdf_page(fitz.Rect(0,0,W,215),layer,0)
lower=fitz.Rect(0,221,W,H)
page.show_pdf_page(lower,old,0,clip=lower)
out.set_metadata({'title':'Figure S1 | Paired HCC epithelial states and stromal connectivity'})
pdf=WORK/'Figure_S1_HCC_paired_preview.pdf'
out.save(pdf,garbage=4,deflate=True);out.close()
doc=fitz.open(pdf);page=doc[0]
assert len(doc)==1 and abs(W-595.2756)<.01 and abs(H-841.8898)<.01
text=page.get_text()
spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
assert all('Arial' in s['font'] for s in spans)
assert all(fitz.Rect(s['bbox']).x0>=4 and fitz.Rect(s['bbox']).x1<=W-4 and fitz.Rect(s['bbox']).y0>=4 and fitz.Rect(s['bbox']).y1<=H-4 for s in spans)
assert not any(x in text for x in ['decontX','KRT8','CDH1','2,000 UMI'])
assert '8 paired patients' in text and '1,575 cells' in text and '11,153 cells' in text
# Exact rendered equality protects every B-D mark, label, and color.
before=old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=lower,alpha=False)
after=page.get_pixmap(matrix=fitz.Matrix(2,2),clip=lower,alpha=False)
assert before.width==after.width and before.height==after.height
assert before.samples==after.samples
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(WORK/'Figure_S1_HCC_paired_preview.png')
page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).save(WORK/'audit/Figure_S1_review.png')
page.get_pixmap(matrix=fitz.Matrix(2.4,2.4),clip=fitz.Rect(0,0,W,215),alpha=False).save(WORK/'audit/S1A_review.png')
assert all(sha(Path(p))==v for p,v in protected.items())
checks=dict(A4=True,Arial_only=True,points_verified=32,patients=8,paired_lines=16,
            B_to_D_pixels_identical=True,current_final_files_unchanged=True,
            figure_sha256=sha(pdf))
(WORK/'audit/verification.json').write_text(json.dumps(checks,indent=2))
oldlegend=(ROOT/'当前主图汇总_20261010/Figure_S1_图注.txt').read_text()
tail=oldlegend[oldlegend.index('(B)'):]
legend='''Figure S1 Paired HCC epithelial states and stromal connectivity

(A) Patient-level proportions of E/M-high cells (left) and cells with a high mesenchymal score (right) in paired tissues from eight patients with HCC (GSE149614; HCC03-HCC10). The analysis includes 1,575 non-tumor hepatocyte-like cells and 11,153 malignant epithelial cells from the same eight patients. Each point represents the proportion within one patient's indicated population, and each line connects that patient's non-tumor and malignant populations. All eligible malignant cells from these patients are included before EMT-state selection. Scores use the same Cook-derived 17-gene epithelial signature, 200-gene Hallmark EMT signature, saved control genes, log-normalization and fixed pooled-malignant thresholds as Figure 2A (E=0.2013776236661785; M=0.0039257136053211). E/M high denotes scores above both thresholds. M-high total includes Mesenchymal high and E/M high. Displayed P values are from two-sided exact paired Wilcoxon signed-rank tests, adjusted across the two displayed metrics using Holm's method (n=8 patient pairs). Non-tumor hepatocyte-like cells are identified from non-tumor tissue; this label is distinct from the source reference populations used in the PAGA analysis. '''+tail
(WORK/'Figure_S1_legend_draft.txt').write_text(legend)
stat=pd.read_csv(WORK/'audit/paired_statistics.csv')
readme='''FigureYa195PanPaire adaptation: HCC paired S1A preview

Template: https://github.com/ying-ge/FigureYa/blob/main/FigureYa195PanPaire/FigureYa195PanPaire.Rmd
Adaptation: the paired segments/points arrangement; each unit is one patient.
Source code is archived under reference/. The template's TCGA and ssGSEA analyses were not run.
FigureYa license: CC BY-NC-SA 4.0. Cite Lu et al. (2025), FigureYa, iMetaMed 1:e70005, doi:10.1002/imm3.70005.

Reproduce: /usr/bin/python3 code/build_preview.py
The current finalized files were not overwritten. B-D retain identical rendered pixels.
HCC01 and HCC02 have no non-tumor hepatocyte-like counterpart and are excluded from this paired comparison.
The published malignant score values are unchanged. Non-tumor cells use the previously saved identical-score analysis.
No single-cell-level inferential tests are performed. Two related outcomes are retained and Holm-adjusted together.

'''+stat.to_string(index=False)+'\n'
(WORK/'README.txt').write_text(readme)
print(json.dumps({'pdf':str(pdf),'checks':checks,'statistics':stat.to_dict(orient='records')},ensure_ascii=False,indent=2))
