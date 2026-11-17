#!/usr/bin/env python3
"""Reformat the approved HCC patient comparison as a percentage heatmap."""
from pathlib import Path
import hashlib,json,shutil,subprocess
import fitz
import numpy as np
import pandas as pd

WORK=Path(__file__).resolve().parents[1];ROOT=WORK.parent
PREVIOUS=ROOT/'S1A_HCC配对_FigureYa195_20261010'
for name in ('inputs','audit','reference'):(WORK/name).mkdir(exist_ok=True)
sha=lambda p:hashlib.sha256(p.read_bytes()).hexdigest()
sources={
    'patient_percentages.csv':PREVIOUS/'inputs/patient_percentages.csv',
    'paired_statistics.csv':PREVIOUS/'audit/paired_statistics.csv',
    'definitions.json':PREVIOUS/'audit/definitions.json',
}
manifest=[]
for name,p in sources.items():
    target=WORK/'inputs'/name
    shutil.copyfile(p,target)
    assert sha(p)==sha(target)
    manifest.append(dict(source=str(p),copy=str(target),sha256=sha(p)))
template=ROOT/'本次整合_20261010/template_sources/FigureYa311PAM50Heatmap.Rmd'
shutil.copyfile(template,WORK/'reference'/template.name)
manifest.append(dict(source=str(template),sha256=sha(template)))
protected={str(p):sha(p) for p in (ROOT/'当前主图汇总_20261010').glob('*') if p.is_file()}
src_pdf=ROOT/'当前主图汇总_20261010/Figure_S1.pdf'
manifest.append(dict(source=str(src_pdf),sha256=sha(src_pdf)))
subprocess.run(['Rscript',str(WORK/'code/draw_heatmap.R')],check=True)
d=pd.read_csv(WORK/'inputs/patient_percentages.csv')
ann=pd.read_csv(WORK/'audit/heatmap_annotations.csv')
joined=ann.merge(d,on=['metric','patient','population'],validate='one_to_one',suffixes=('_plot','_source'))
assert len(joined)==32
assert np.allclose(joined.percentage_plot,joined.percentage_source,atol=1e-12,rtol=0)
assert np.allclose(joined.percentage_plot,100*joined.positive_cells/joined.cells,atol=1e-12,rtol=0)
assert ann.width.eq(36).all() and ann.height.eq(36).all()
old=fitz.open(src_pdf);layer=fitz.open(WORK/'audit/S1A_heatmap_layer.pdf')
W,H=old[0].rect.width,old[0].rect.height
# Remove old A content in an in-memory copy, not just visually clip it.
# A clipped Form XObject can still return hidden text in get_textbox().
lower_source=fitz.open(src_pdf)
lower_source[0].add_redact_annot(fitz.Rect(0,0,W,220),fill=False)
lower_source[0].apply_redactions(images=0,graphics=1,text=0)
out=fitz.open();page=out.new_page(width=W,height=H)
page.show_pdf_page(fitz.Rect(0,0,W,215),layer,0)
lower=fitz.Rect(0,221,W,H)
page.show_pdf_page(lower,lower_source,0,clip=lower)
out.set_metadata({'title':'Figure S1 | Paired HCC state proportions and stromal connectivity'})
pdf=WORK/'Figure_S1_HCC_heatmap_preview.pdf'
out.save(pdf,garbage=4,deflate=True);out.close()
doc=fitz.open(pdf);page=doc[0]
assert len(doc)==1 and abs(W-595.2756)<.01 and abs(H-841.8898)<.01
spans=[s for b in page.get_text('dict')['blocks'] for l in b.get('lines',[]) for s in l['spans']]
assert all('Arial' in s['font'] for s in spans)
assert all(fitz.Rect(s['bbox']).x0>=4 and fitz.Rect(s['bbox']).x1<=W-4 and fitz.Rect(s['bbox']).y0>=4 and fitz.Rect(s['bbox']).y1<=H-4 for s in spans)
for row in ann.itertuples():
    actual=page.get_textbox(fitz.Rect(row.x-17.5,row.y-8,row.x+17.5,row.y+8)).strip()
    assert actual==row.label,(row.patient,row.metric,row.population,actual,row.label)
drawings=page.get_drawings()
tiles=[]
for drawing in drawings:
    r=drawing['rect']
    if r.y1<=211 and abs(r.width-36)<.02 and abs(r.height-36)<.02:
        tiles.append(list(r))
assert len(tiles)==32,len(tiles)
before=old[0].get_pixmap(matrix=fitz.Matrix(2,2),clip=lower,alpha=False)
after=page.get_pixmap(matrix=fitz.Matrix(2,2),clip=lower,alpha=False)
assert before.samples==after.samples
assert all(sha(Path(p))==v for p,v in protected.items())
page.get_pixmap(matrix=fitz.Matrix(300/72,300/72),alpha=False).save(WORK/'Figure_S1_HCC_heatmap_preview.png')
page.get_pixmap(matrix=fitz.Matrix(1.6,1.6),alpha=False).save(WORK/'audit/Figure_S1_review.png')
page.get_pixmap(matrix=fitz.Matrix(2.4,2.4),clip=fitz.Rect(0,0,W,215),alpha=False).save(WORK/'audit/S1A_review.png')
checks=dict(A4=True,Arial_only=True,patient_pairs=8,percentage_tiles=32,tile_size_points=[36,36],
            all_numbers_verified_at_pdf_positions=True,common_scale=[0,100],
            patient_order_preserved=True,row_standardization=False,clustering=False,
            old_A_hidden_text_removed=True,B_to_D_pixels_identical=True,
            current_final_files_unchanged=True,figure_sha256=sha(pdf))
(WORK/'audit/verification.json').write_text(json.dumps(checks,indent=2))
(WORK/'audit/source_manifest.json').write_text(json.dumps(manifest,ensure_ascii=False,indent=2))
oldlegend=(PREVIOUS/'Figure_S1_legend_draft.txt').read_text()
legend=oldlegend.replace('E/M-high cells (left) and cells with a high mesenchymal score (right)',
                        'E/M-high cells (upper block) and cells with a high mesenchymal score (lower block)')
legend=legend.replace("Each point represents the proportion within one patient's indicated population, and each line connects that patient's non-tumor and malignant populations.",
    'Columns identify patients in the same order in both blocks. Within each block, the two rows show non-tumor hepatocyte-like cells and malignant epithelial cells. Row-label cell counts are pooled across the eight patients; each tile percentage uses that patient\'s population-specific cell count as its denominator. Heatmap color represents the cell proportion on a shared 0-100% scale; numbers give percentages to one decimal place. Values are neither row-standardized nor clustered. Narrow blue and coral strips identify the two populations; the numerical color scale is shared by both.')
(WORK/'Figure_S1_legend_draft.txt').write_text(legend)
(WORK/'README.txt').write_text('''FigureYa311 paired HCC percentage heatmap preview

Template: https://github.com/ying-ge/FigureYa/blob/main/FigureYa311PAM50Heatmap/FigureYa311PAM50Heatmap.Rmd
Adapted element: ComplexHeatmap with annotated sample groups and numeric cell annotations.
The PAM50 and GSVA computations in the original template are not used.
FigureYa license: CC BY-NC-SA 4.0. Cite Lu et al. (2025), FigureYa, iMetaMed 1:e70005, doi:10.1002/imm3.70005.

Reproduce: /usr/bin/python3 code/build_preview.py
Saved patient proportions and paired statistics are copied unchanged from the previous paired-plot preview.
Eight patient columns, two state blocks, two tissue rows per block: 32 square 36-point cells.
Shared 0-100% scale. No row standardization, clustering, rescoring, or new statistical analysis.
Patient-level paired Wilcoxon P values retain the previous two-metric Holm adjustment.
All existing final files remain unchanged; B-D retain identical rendered pixels.
''')
print(json.dumps({'pdf':str(pdf),'checks':checks},ensure_ascii=False,indent=2))
