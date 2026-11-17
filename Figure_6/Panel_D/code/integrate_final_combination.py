#!/usr/bin/env python3
"""Preserve four V16 rows and integrate two cached replacements; no GSVA rerun."""
from pathlib import Path
import csv
import hashlib
import io
import json
import zipfile
import numpy as np
import pandas as pd

WORK = Path(__file__).resolve().parents[1]
ROOT = WORK.parent
OUT = WORK / 'final_combination'
TABLES, AUDIT = OUT / 'tables', OUT / 'audit'
TABLES.mkdir(parents=True, exist_ok=True)
AUDIT.mkdir(parents=True, exist_ok=True)
V16 = ROOT / 'V16_Insulin与PI3K分开_20261026'
REACTOME = ROOT / 'V17_PI3K基因集复核_20261027/reactome_six'
V9 = ROOT / 'V9_严格逐条复核_20261024/audit/fig6_recomputed'
V6 = ROOT / 'V6_最终图_20261024/audit'
LABELS = ['TGF-β', 'PI3K-AKT', 'Insulin', 'HGF-RAS/ERK', 'EGF', 'IGF']
PRESERVED = ['TGF-β', 'Insulin', 'EGF', 'IGF']
HGF = 'KEGG_MEDICUS_REFERENCE_HGF_MET_RAS_ERK_SIGNALING_PATHWAY'
INITIAL_SHA = '175eff84d3f9b524ff62c34f1acfefd26138f48b409d7466cc5d4b3245ec2b9e'
INITIAL_SELECTION = WORK / 'audit/selection_before_scores.json'
FINAL_SELECTION = WORK / 'audit/final_combination_selection.json'
EXPECTED = {'ChinaHCC': ('06910d7e500ac6da7526787396966e96', 705),
            'USAHCC': ('d621adf1a8b5bc081a0ad938cdbb9329', 454),
            'CCAcombined': ('e19b2ca07f12f2985b37cdc256b48d50', 369)}
sha = lambda p: hashlib.sha256(Path(p).read_bytes()).hexdigest()
md5 = lambda p: hashlib.md5(Path(p).read_bytes()).hexdigest()
read = lambda p, **kwargs: pd.read_csv(p, float_precision='round_trip', **kwargs)
selection = json.loads(FINAL_SELECTION.read_text())
assert sha(INITIAL_SELECTION) == INITIAL_SHA
assert selection['status'] == 'post_score_exploratory_selection'
assert selection['selection_includes_observed_direction'] and not selection['independent_validation']
assert [x['label'] for x in selection['selected']] == LABELS
assert selection['preserved_V16_rows'] == PRESERVED
assert selection['replacement_rows'] == ['PI3K-AKT', 'HGF-RAS/ERK']
assert [x['standard_name'] for x in selection['selected']] == [
    'V16_ORIGINAL_TGFB', 'R-HSA-1257604', 'V16_KEGG_hsa04910', HGF,
    'V16_ORIGINAL_EGF', 'V16_REACTOME_R-HSA-2404192']
notes20 = json.loads((WORK / 'audit/calculation_notes.json').read_text())
notes17 = json.loads((REACTOME / 'audit/calculation_notes.json').read_text())
for note in [notes20, notes17]:
    assert note['status'] == 'completed' and all(note['checks'].values())
    for key, value in {'kcdf':'Gaussian', 'minSize':3, 'maxSize':'Inf', 'seed':1,
                       'BPPARAM':'BiocParallel::SerialParam', 'package_version':'2.6.6'}.items():
        assert note['gsva'][key] == value

legacy_coverage_path = V9 / 'signature_gene_coverage.csv'
legacy_coverage = read(legacy_coverage_path)
archive = V6 / 'inputs/ReactomePathways_current_20261024.gmt.zip'
assert sha(archive) == '8c1dbc8578431da5d2d5118262718c60b553a9be3398e93658daa069e4a9afd4'
with zipfile.ZipFile(archive) as z:
    archived_rows = [x.split('\t') for x in z.read('ReactomePathways.gmt').decode().splitlines()]
reactome_members = {}
for identifier in ['R-HSA-1257604', 'R-HSA-2404192']:
    found = [x for x in archived_rows if x[1] == identifier]
    assert len(found) == 1
    reactome_members[identifier] = found[0][2:]
insulin_gmt = V16 / 'inputs/KEGG_hsa04910_Insulin.gmt'
hgf_gmt = WORK / f'inputs/{HGF}.gmt'
hgf_json = WORK / f'inputs/{HGF}.json'
hgf_members = hgf_gmt.read_text().strip().split('\t')[2:]
assert hgf_members == json.loads(hgf_json.read_text())[HGF]['geneSymbols']
members = {
    'TGF-β': legacy_coverage[(legacy_coverage.cohort == 'ChinaHCC') & (legacy_coverage.pathway == 'TGF')].gene.tolist(),
    'PI3K-AKT': reactome_members['R-HSA-1257604'],
    'Insulin': insulin_gmt.read_text().strip().split('\t')[2:],
    'HGF-RAS/ERK': hgf_members,
    'EGF': legacy_coverage[(legacy_coverage.cohort == 'ChinaHCC') & (legacy_coverage.pathway == 'EGF-RAS')].gene.tolist(),
    'IGF': reactome_members['R-HSA-2404192']}
assert [len(members[x]) for x in LABELS] == [11, 284, 140, 15, 12, 54]
assert all(len(v) == len(set(v)) for v in members.values())
definitions = pd.DataFrame(selection['selected'])
definitions['complete_member_count'] = [len(members[x]) for x in LABELS]
definitions['treatment'] = ['preserve_original_V16' if x in PRESERVED else 'replace_from_cached_scores' for x in LABELS]
definitions['source_type'] = ['original_custom_signature', 'complete_archived_Reactome', 'complete_official_KEGG',
                              'complete_official_KEGG_MEDICUS', 'original_custom_signature', 'complete_archived_Reactome']
definitions['per_cell_availability'] = ['not_archived_in_this_integration', 'cached_and_verified', 'cached_and_verified',
                                       'cached_and_verified', 'not_archived_in_this_integration', 'cached_and_verified']
definitions['selection_status'] = 'POST-SCORE EXPLORATORY; outcome-informed; not independent validation'
definitions.to_csv(TABLES / 'pathway_definitions.csv', index=False)
pd.DataFrame([dict(pathway=label, standard_name=selection['selected'][i]['standard_name'],
                   scoring_identifier=g, member_order=j+1, source_type=definitions.iloc[i].source_type)
              for i, label in enumerate(LABELS) for j, g in enumerate(members[label])]).to_csv(
    TABLES / 'complete_gene_set_membership.csv', index=False)

protected = {INITIAL_SELECTION, FINAL_SELECTION, archive, legacy_coverage_path, insulin_gmt, hgf_gmt, hgf_json,
             WORK/'audit/calculation_notes.json', REACTOME/'audit/calculation_notes.json',
             WORK/'inputs/official_source_manifest.json'}
protected.update((WORK/'tables').glob('*.csv'))
protected.update((V16/'tables').glob('*.csv'))
protected.update((REACTOME/'tables').glob('*.csv'))
protected.update(V6/'IGF'/f'{k}_IGF_per_cell.csv' for k in EXPECTED)
before = {str(p):sha(p) for p in sorted(protected)}

def raw_csv(path):
    lines = Path(path).read_bytes().splitlines(keepends=True)
    header = next(csv.reader([lines[0].decode().rstrip('\r\n')]))
    rows = {next(csv.reader([line.decode().rstrip('\r\n')]))[0]: line for line in lines[1:]}
    return lines[0], header, rows

def reindex_raw_row(line, source_columns, wanted_columns, label):
    cells = next(csv.reader([line.decode().rstrip('\r\n')]))
    values = dict(zip(source_columns[1:], cells[1:]))
    result = [label] + [values[col] for col in wanted_columns]
    buf = io.StringIO(newline=''); csv.writer(buf, lineterminator='\n').writerow(result)
    return buf.getvalue().encode()

checks, preservation, long_scores, equal, pairwise, comparisons, coverage_records = [], [], [], [], [], [], []
cov20 = read(WORK/'tables/gene_coverage.csv')
cov17 = read(REACTOME/'tables/gene_coverage.csv')
cov16 = read(V16/'tables/gene_coverage.csv')
for key, (expected_hash, expected_n) in EXPECTED.items():
    original = Path('/Users/sodajay/Desktop/Deutsch_Project/03_EMT/archive/EMT_reproduce_final_20260910/emt_invivo_bodyfig_20260901') / f'S2_{key}_joint_with_normal_roles_frozen.rds'
    assert md5(original) == expected_hash
    p20 = read(WORK/f'tables/{key}_all_candidate_per_cell_scores.csv').set_index('barcode')
    p17 = read(REACTOME/f'tables/{key}_six_pathway_per_cell_scores.csv').set_index('barcode')
    p16 = read(V16/f'tables/{key}_two_pathway_per_cell_scores.csv').set_index('barcode')
    pigf = read(V6/f'IGF/{key}_IGF_per_cell.csv').set_index('barcode')
    assert len(p20) == expected_n and not p20.index.has_duplicates
    for source in [p17,p16,pigf]:
        assert len(source) == expected_n and not source.index.has_duplicates and set(source.index) == set(p20.index)
        assert source.loc[p20.index].role_id.equals(p20.role_id)
    for source in [p20,p17]:
        assert np.max(np.abs(source.loc[p20.index,'KEGG_PI3K_AKT_control']-p16.loc[p20.index,'PI3K-AKT'])) <= 1e-10
    available = p20[['cohort_key','cohort','role_id','role_label','displayed_in_figure6']].copy()
    available['PI3K-AKT'] = p17.loc[p20.index,'PI3K-AKT']
    available['HGF-RAS/ERK'] = p20[HGF]
    available['Insulin'] = p16.loc[p20.index,'Insulin']
    available['IGF'] = pigf.loc[p20.index,'IGF']
    available.to_csv(TABLES/f'{key}_available_per_cell_scores.csv', float_format='%.17g')
    available[['cohort_key','cohort','role_id','role_label','displayed_in_figure6','PI3K-AKT','HGF-RAS/ERK']].to_csv(
        TABLES/f'{key}_replacement_per_cell_scores.csv', float_format='%.17g')
    source_paths = {
      'by_role': (REACTOME/f'tables/fig6_{key}_signaling_GSVA_activity_6pathways_by_role.csv', WORK/f'tables/{key}_all_candidate_role_means.csv'),
      'zscore': (REACTOME/f'tables/fig6_{key}_signaling_GSVA_activity_6pathways_zscore.csv', WORK/f'tables/{key}_all_candidate_role_zscores.csv')}
    final = {}; previous = {}
    for kind in ['by_role','zscore']:
        oldpath = V16/f'tables/fig6_{key}_signaling_GSVA_activity_6pathways_{kind}.csv'
        previous[kind] = read(oldpath,index_col=0)
        _, a_columns, a_rows = raw_csv(source_paths[kind][0])
        _, b_columns, b_rows = raw_csv(source_paths[kind][1])
        for suffix in ['', '_displayed']:
            oldvariant = V16/f'tables/fig6_{key}_signaling_GSVA_activity_6pathways_{kind}{suffix}.csv'
            rawheader, columns, oldrows = raw_csv(oldvariant)
            newrows = {label:oldrows[label] for label in PRESERVED}
            newrows['PI3K-AKT'] = reindex_raw_row(a_rows['PI3K-AKT'],a_columns,columns[1:],'PI3K-AKT')
            newrows['HGF-RAS/ERK'] = reindex_raw_row(b_rows[HGF],b_columns,columns[1:],'HGF-RAS/ERK')
            dest = TABLES/f'fig6_{key}_signaling_GSVA_activity_6pathways_{kind}{suffix}.csv'
            dest.write_bytes(rawheader+b''.join(newrows[label] for label in LABELS))
            _, _, rereadrows = raw_csv(dest)
            oldvalues, newvalues = read(oldvariant,index_col=0),read(dest,index_col=0)
            for label in PRESERVED:
                assert rereadrows[label] == oldrows[label]
                assert np.array_equal(oldvalues.loc[label].to_numpy(),newvalues.loc[label].to_numpy())
                preservation.append(dict(cohort_key=key,kind=kind,displayed=bool(suffix),pathway=label,
                    source_file=str(oldvariant),output_file=str(dest),source_row_bytes_identical=True,numeric_maximum_difference=0.0))
            if not suffix: final[kind] = newvalues
    means,zscores = final['by_role'],final['zscore']
    displayed_columns = read(TABLES/f'fig6_{key}_signaling_GSVA_activity_6pathways_zscore_displayed.csv',index_col=0).columns.tolist()
    mean_check = available.groupby('role_label')[['PI3K-AKT','HGF-RAS/ERK','Insulin','IGF']].mean().T.reindex(columns=means.columns)
    mean_error = np.abs(mean_check-means.loc[mean_check.index]).to_numpy().max()
    zcheck = means.sub(means.mean(axis=1),axis=0).div(means.std(axis=1,ddof=1),axis=0).fillna(0)
    z_error = np.abs(zcheck-zscores).to_numpy().max()
    assert mean_error <= 1e-10 and z_error <= 1e-10
    role_ids = p20[['role_id','role_label']].drop_duplicates().set_index('role_label').role_id.to_dict()
    role_counts = p20.role_label.value_counts()
    epithelial = [r for r in means.columns if role_ids[r].startswith('normal_') or role_ids[r]=='malignant_epithelial_EMT_root']
    stromal = [r for r in means.columns if r not in epithelial]
    assert len(stromal) == (4 if key=='ChinaHCC' else 3)
    for label in LABELS:
        oldlabel = selection['selected'][LABELS.index(label)]['V16_label']
        for role in means.columns:
            long_scores.append(dict(cohort_key=key,pathway=label,role_id=role_ids[role],role_label=role,n_cells=int(role_counts[role]),
                mean_gsva=means.loc[label,role],within_cohort_role_zscore=zscores.loc[label,role],displayed_in_figure6=role in displayed_columns,
                preserved_original_row=label in PRESERVED,selection_status='post_score_exploratory'))
            comparisons.append(dict(cohort_key=key,pathway=label,previous_pathway=oldlabel,role_label=role,
                previous_mean_GSVA=previous['by_role'].loc[oldlabel,role],new_mean_GSVA=means.loc[label,role],
                previous_Z=previous['zscore'].loc[oldlabel,role],new_Z=zscores.loc[label,role],
                preserved_original_row=label in PRESERVED,displayed_in_figure6=role in displayed_columns))
        for ep in epithelial:
            shared=dict(cohort_key=key,pathway=label,epithelial_role_id=role_ids[ep],epithelial_role_label=ep,
                epithelial_type='malignant' if role_ids[ep]=='malignant_epithelial_EMT_root' else 'normal',
                displayed_in_figure6=ep in displayed_columns,non_TGF_target=label!='TGF-β')
            diff=float(means.loc[label,ep]-means.loc[label,stromal].mean())
            differences=[float(means.loc[label,ep]-means.loc[label,s]) for s in stromal]
            equal.append(dict(**shared,epithelial_mean_GSVA=means.loc[label,ep],equal_weight_stromal_role_mean=means.loc[label,stromal].mean(),
                epithelial_minus_equal_weight_stroma=diff,positive_difference=diff>1e-12,
                epithelial_Z=zscores.loc[label,ep],equal_weight_stromal_Z=zscores.loc[label,stromal].mean(),
                Z_difference=zscores.loc[label,ep]-zscores.loc[label,stromal].mean(),stromal_role_count=len(stromal),
                stromal_roles=';'.join(stromal),greater_than_every_stromal_role=all(d>1e-12 for d in differences),
                weighting='equal weight per stromal role, not cell-weighted',p_values_computed=False))
            for st,difference in zip(stromal,differences):
                pairwise.append(dict(**shared,stromal_role_id=role_ids[st],stromal_role_label=st,
                    epithelial_mean_GSVA=means.loc[label,ep],stromal_mean_GSVA=means.loc[label,st],
                    epithelial_minus_stromal_role=difference,positive_difference=difference>1e-12,
                    epithelial_Z=zscores.loc[label,ep],stromal_Z=zscores.loc[label,st],
                    Z_difference=zscores.loc[label,ep]-zscores.loc[label,st],p_values_computed=False))
    for label in LABELS:
        if label in ['TGF-β','EGF','IGF']:
            oldname={'TGF-β':'TGF','EGF':'EGF-RAS','IGF':'IGF'}[label]
            c=legacy_coverage[(legacy_coverage.cohort==key)&(legacy_coverage.pathway==oldname)].copy()
            c=c.rename(columns={'gene':'scoring_identifier','in_counts':'present_in_counts','nonzero':'detected_in_any_cell'})
            c['supplied_to_GSVA']=c.detected_in_any_cell;c['has_gene_symbol']=True;sourcefile=legacy_coverage_path
        elif label=='Insulin':
            c=cov16[(cov16.cohort_key==key)&(cov16.pathway=='Insulin')].copy();sourcefile=V16/'tables/gene_coverage.csv'
        elif label=='PI3K-AKT':
            c=cov17[(cov17.cohort_key==key)&(cov17.pathway=='PI3K-AKT')].copy();sourcefile=REACTOME/'tables/gene_coverage.csv'
        else:
            c=cov20[(cov20.cohort_key==key)&(cov20.standard_name==HGF)].copy();sourcefile=WORK/'tables/gene_coverage.csv'
        assert len(c)==len(members[label]) and set(c.scoring_identifier)==set(members[label])
        for x in c.to_dict('records'):
            coverage_records.append(dict(cohort_key=key,pathway=label,scoring_identifier=x['scoring_identifier'],
                has_gene_symbol=bool(x['has_gene_symbol']),present_in_counts=bool(x['present_in_counts']),
                detected_in_any_cell=bool(x['detected_in_any_cell']),variable=bool(x['variable']),
                supplied_to_GSVA=bool(x['supplied_to_GSVA']),source_coverage_file=str(sourcefile)))
    assert md5(original)==expected_hash
    checks.append(dict(cohort_key=key,cells=expected_n,frozen_RDS=str(original),frozen_MD5=expected_hash,
        cached_barcode_and_role_match=True,per_cell_to_mean_maximum_error=float(mean_error),row_Z_maximum_error=float(z_error),
        frozen_KEGG_control_pass=True,available_per_cell_rows=['PI3K-AKT','HGF-RAS/ERK','Insulin','IGF'],
        unavailable_per_cell_rows=['TGF-β','EGF']))

eq,pair,cov=pd.DataFrame(equal),pd.DataFrame(pairwise),pd.DataFrame(coverage_records)
assert len(eq)==42 and len(pair)==144 and len(cov)==3*516
pd.DataFrame(long_scores).to_csv(TABLES/'six_pathways_all_role_scores.csv',index=False,float_format='%.17g')
pd.DataFrame(comparisons).to_csv(TABLES/'V16_vs_final_role_comparison.csv',index=False,float_format='%.17g')
eq.to_csv(TABLES/'epithelial_vs_equal_weight_stromal_role_means.csv',index=False,float_format='%.17g')
pair.to_csv(TABLES/'epithelial_vs_each_stromal_role.csv',index=False,float_format='%.17g')
eq[eq.displayed_in_figure6].to_csv(TABLES/'displayed_epithelial_vs_equal_weight_stromal_role_means.csv',index=False,float_format='%.17g')
pair[pair.displayed_in_figure6].to_csv(TABLES/'displayed_epithelial_vs_each_stromal_role.csv',index=False,float_format='%.17g')
cov.to_csv(TABLES/'gene_coverage.csv',index=False)
coverage_summary=cov.groupby(['cohort_key','pathway'],sort=False).agg(complete_member_count=('scoring_identifier','size'),
    present_in_counts=('present_in_counts','sum'),detected_in_any_cell=('detected_in_any_cell','sum'),
    variable=('variable','sum'),supplied_to_GSVA=('supplied_to_GSVA','sum')).reset_index()
coverage_summary.to_csv(TABLES/'gene_coverage_summary.csv',index=False)
cov[~cov.supplied_to_GSVA].to_csv(TABLES/'unavailable_members_by_cohort.csv',index=False)
pd.DataFrame(preservation).to_csv(AUDIT/'four_original_rows_exact_preservation.csv',index=False)
assert len(preservation)==48 and all(x['source_row_bytes_identical'] for x in preservation)
after={p:sha(p) for p in before};assert before==after and sha(INITIAL_SELECTION)==INITIAL_SHA
pd.DataFrame([dict(path=p,SHA256_before=h,SHA256_after=after[p],unchanged=True) for p,h in before.items()]).to_csv(
    AUDIT/'provenance_and_preservation_hashes.csv',index=False)
visible=eq[eq.displayed_in_figure6]; target=visible[visible.non_TGF_target]
visible_pair=pair[pair.displayed_in_figure6];target_pair=visible_pair[visible_pair.non_TGF_target]
assert len(visible)==36 and len(target)==30 and len(visible_pair)==120 and len(target_pair)==100
summary={'all_six_displayed_vs_equal_stroma':{'positive':int(visible.positive_difference.sum()),'total':len(visible)},
         'non_TGF_displayed_vs_equal_stroma':{'positive':int(target.positive_difference.sum()),'total':len(target)},
         'all_six_displayed_greater_than_every_stromal_role':{'positive':int(visible.greater_than_every_stromal_role.sum()),'total':len(visible)},
         'non_TGF_displayed_greater_than_every_stromal_role':{'positive':int(target.greater_than_every_stromal_role.sum()),'total':len(target)},
         'all_six_displayed_pairwise_stromal_comparisons':{'positive':int(visible_pair.positive_difference.sum()),'total':len(visible_pair)},
         'non_TGF_displayed_pairwise_stromal_comparisons':{'positive':int(target_pair.positive_difference.sum()),'total':len(target_pair)},
         'positive_difference_threshold':1e-12}
(AUDIT/'contrast_direction_summary.json').write_text(json.dumps(summary,indent=2)+'\n')
report={'status':'completed','analysis_status':'POST-SCORE EXPLORATORY / outcome-informed selection / not independent validation',
    'selection_file':str(FINAL_SELECTION),'selection_SHA256':before[str(FINAL_SELECTION)],'initial_twenty_selection_SHA256':INITIAL_SHA,
    'row_order':LABELS,'preserved_original_rows':PRESERVED,'replacement_rows':['PI3K-AKT','HGF-RAS/ERK'],
    'definitions':definitions.to_dict('records'),'GSVA_rerun':False,'source_pipeline':'Same frozen RDS and cells; lnCP1e4 full expression background; GSVA2.6.6 Gaussian, minSize3, maxSizeInf, seed1, SerialParam',
    'assembly':'Four original V16 raw-mean and Z rows copied byte-for-byte in both full and displayed CSVs; two cached replacement rows independently reconstructed from source per-cell scores, then source summary values retained',
    'missing_per_cell_policy':'TGF-beta and EGF original role summaries are preserved exactly; their unavailable per-cell scores are neither invented nor regenerated',
    'source_membership_policy':'Original custom TGF-beta/EGF and original complete Insulin/IGF definitions preserved; complete Reactome284 and Medicus HGF15 used without union, subtraction, or gene editing',
    'cohorts':checks,'contrast_direction_summary':summary,'coverage_summary':coverage_summary.to_dict('records'),
    'limitations':['The two replacements were chosen after seeing results; this is exploratory selection, not validation.',
                   'Above equal-weight stromal-role mean does not imply above every stromal role.',
                   'RNA gene-set scores do not establish biochemical phosphorylation or independent ligand-specific activity.',
                   'The six rows are not all newly official or same-source; four original definitions remain.'],
    'checks':{'four_original_rows_exact_in_48_checks':True,'replacement_source_scores_reconstructed':True,
              'same_frozen_RDS_barcodes_roles_and_pipeline':True,'KEGG_source_controls_pass':True,
              'all_initial_twenty_and_previous_inputs_unchanged':True,'selection_manifest_unchanged':True,
              'complete_memberships_preserved':True,'no_GSVA_rerun':True,'no_P_values':True,'no_figures_written':True}}
(AUDIT/'calculation_notes.json').write_text(json.dumps(report,ensure_ascii=False,indent=2)+'\n')
print(json.dumps(summary,ensure_ascii=False,indent=2))
print(coverage_summary.to_string(index=False))
