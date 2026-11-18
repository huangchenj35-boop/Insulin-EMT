# Insulin depletion results in TGF-β-induced epithelial-to-mesenchymal transition in liver epithelial cells

This repository contains code and the figure-level data for **Figure 6, Figure 7 and Figure S4**, the three bioinformatics figures contributed by C.H. to the manuscript. Other authors' experimental figures, clinical data, manuscripts and historical working folders are not part of this repository.

## Repository contents

| Figure | Panel | Analysis | Directory |
| --- | --- | --- | --- |
| Figure 6 | A | 95,022-cell global atlas UMAP, stromal subtypes and qHSC score | Figure_6/Panel_A |
| Figure 6 | B | Epithelial and mesenchymal phenotype gene-expression heatmaps | Figure_6/Panel_B |
| Figure 6 | C | Transcription-factor expression heatmaps | Figure_6/Panel_C |
| Figure 6 | D | RNA-based GSVA of six signaling programs | Figure_6/Panel_D |
| Figure 7 | A | Joint PAGA observed/expected connectivity | Figure_7/Panel_A |
| Figure 7 | B | PAGA balanced-resampling results | Figure_7/Panel_B |
| Figure 7 | C | Unrooted stromal DDRTree embedding | Figure_7/Panel_C |
| Figure 7 | D | Chromosome-arm RNA-derived CNV | Figure_7/Panel_D |
| Figure S4 | A | Paired HCC epithelial-state percentages | Figure_S4/Panel_A |
| Figure S4 | B | Cook et al. two-axis E/M scores and proportions | Figure_S4/Panel_B |

Within each Figure directory:

- **Panel_A/, Panel_B/, etc.** contain the related historical upstream analysis scripts in **code/** and exported derived input data in **data/**. In several cases, the same original script supported multiple panels; its exact historical purpose is stated in the manifest instead of inventing a new de novo script.
- **assembly/** contains the original final multi-panel composition script and the necessary vector source components or frozen numeric tables.
- **final/** contains the three **approved, unchanged** final PDFs and corresponding PNG previews.

The files **PANEL_INDEX.tsv** and **SOURCE_PROVENANCE.tsv** document the provenance and SHA-256 of the original source files. The separate **runtime_map.json** provides paths for the historical assembly code in an isolated temporary workspace.

## Reproducing the final figure layouts

Install Python 3 with the packages in the provided requirements file:

    python3 -m pip install -r environment/requirements_plotting.txt

Run from the repository root:

    python3 rebuild_figures.py

Or rebuild one figure:

    python3 rebuild_figures.py --figure Figure_7

The entry point assembles the approved PDFs using the original Python scripts and frozen figure data/vector artwork, and checks each resulting PDF against the reference **pixel-by-pixel at 144 dpi**. It writes results to a newly created output folder. The approved files remain unchanged.

The scripts require a locally installed **Arial.ttf** and **Arial Bold.ttf**. Font files are not redistributed. By default they are read from the standard macOS Supplemental fonts folder; on other operating systems provide **--font-dir /path/to/fonts**.

### Figure 6

**A–C:** The prior multi-panel rendering script is stored at Figure_6/assembly/build_figure1_ABC_original.py; corresponding RNA expression and frozen UMAP inputs are separated by panel.

**D:** Historical candidate GSVA scoring and integration scripts are under Figure_6/Panel_D/code. The selected Figure 6D tables were used by Figure_6/assembly/build_combined_figures.py, which retains the previously approved vector artwork for A–C. Its vector source is included under Figure_6/assembly/assets.

### Figure 7

**A–B:** Original PAGA and balanced-resampling scripts and numerical results are in their respective panel directories.

**C:** The DDRTree script and frozen coordinates/edges are in Panel_C.

**D:** Chromosome-arm CNV preparation/export scripts and derived plotted data are in Panel_D. The approved composite embeds an audited pre-existing CNV vector panel supplied under assembly/assets. The published layout is not represented as a complete de novo CNV analysis.

The approved final assembly script is Figure_7/assembly/build_figure7.py.

### Figure S4

**A:** Original patient-level analysis, paired statistics, heatmap sources and the data displayed in the final paired panel are stored in Panel_A.

**B:** Original two-axis EMT-score analysis and the derived per-cell data are stored in Panel_B. Its approved image is retained as vector artwork, not silently regenerated with new values.

The final assembly script is Figure_S4/assembly/build_final_figures.py. The only alteration to this script for the public code package was the removal of a hash-only reference to an unrelated private manuscript DOCX. This document was not used to calculate or draw Figure S4; no analysis formula or plotted values were changed.

## Data sources and reproducibility scope

Published GEO datasets used in this analysis:

- [GSE149614](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614): HCC
- [GSE151530](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE151530): HCC / CCA
- [GSE138709](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE138709): CCA

The repository contains **derived, figure-level data**, including cell and sample pseudonymous identifiers, coordinates, gene expression summaries, module scores, GSVA matrices, PAGA connectivity and resampling values, patient-level proportions and CNV display inputs. These are the files found in the final archived figure workstation. They are not raw reads or a substitute for complete count matrices.

**Important reproducibility boundary:** The archive's final figure reconstruction from frozen numerical tables and artwork had previously been checked for pixel identity at 144 dpi (see environment/previous_rebuild_verification.json). The historical scripts for the upstream analyses are retained here, but the complete analysis from original counts was **not rerun** for this publication package; some intermediate RDS objects from the original analysis are missing. Therefore no claim of full count-to-figure reproducibility is made.

The preserved original script versions and derived inputs allow readers to inspect the actual analytical procedures and reproduce the approved final compositions. Any complete reanalysis should start from the original GEO accessions and use the methodological details in the accompanying manuscript.

## Version and release

- Figure 6: V20, 27 October 2026 — Figure_6_combined.pdf.
- Figure 7: V6, 24 October 2026 — Figure_7.pdf.
- Figure S4: V14, 26 October 2026 — Figure_S4.pdf.

This repository intentionally does not include unrelated figures, wet-lab data, raw FASTQ, historical redundant versions, private documents, bundled compiled libraries, or commercial font files. No reuse license is implied unless one is added by the authors.

The source workstation filed the corresponding supplement as Figure S1. The repository uses Figure S4, following the final manuscript. The archived source filenames within the original plotting scripts remain unchanged.
