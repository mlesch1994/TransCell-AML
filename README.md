# AML TransCell / CIBERSORTx analysis code

Code used to make the figures for the AML manuscript. Single-cell references from Guo et al. are used to deconvolve Beat AML bulk RNA-seq with CIBERSORTx (high-resolution mode). Drug sensitivity is then predicted with TransCell, and the predictions are compared with ex vivo drug response and clinical outcome.

## Folder layout

Run every script from this folder. All paths are relative to it.

- `data/`: input files (see below)
- `results/`: files written by one notebook and read by another
- `figures/`: all figures (created when the scripts run)

## Scripts

| Script | What it does |
| --- | --- |
| `01_reference_difference_matrix.ipynb` | Compares the AML progenitor fraction estimated with each single-cell reference to the bone-marrow blast percentage. |
| `02_single_cell_landscape_and_deconvolution_failure.ipynb` | UMAP of the CIBERSORTx high-resolution profiles across references, cell types and patients. Identifies failed deconvolution runs and relates failure to cell fractions and reference properties. Runs on the HPCC because of the size of the deconvolution output. |
| `03_cell_composition_reference_selection_and_virtual_screen.ipynb` | Cell-type composition across the single-cell references, selection of references by distance to the blast percentage (writes `results/dist_df.csv`), and the TransCell virtual drug screen: top predictions, MOA and target enrichment, and comparison with PBMC predictions. |
| `04_prediction_performance_bulk_vs_deconvolved.ipynb` | Correlation between TransCell scores and ex vivo AUC for bulk and deconvolved samples. Compares performance between low- and high-blast patients, between bulk and each reference, and between bulk and the best reference. |
| `05_cell_line_selection_kasumi1_thp1.ipynb` | Top predicted differentially sensitive compounds for the Kasumi-1 and THP-1 cell lines. |
| `06_tcga_cell_type_survival_analysis.R` | Scores TCGA AML samples for each single-cell cell-type signature (top 50 marker genes, `sigatlas`) and tests the association with overall survival (Kaplan-Meier, log-rank, Cox). |
| `07_supplementary_figures_S4_S8.ipynb` | Supplementary Figures S4 to S8 for the reviewer response. |
| `08_transcell_score_vs_clinical_outcome.ipynb` | Beat AML de novo patients treated with 7+3-type induction: TransCell scores for cytarabine, daunorubicin, idarubicin, venetoclax and azacitidine versus vital status, induction response and overall survival (Mann-Whitney U, Cox adjusted for age and ELN2017, Kaplan-Meier). |
| `09_exvivo_auc_vs_survival.ipynb` | Ex vivo drug AUC versus overall survival for every drug in Beat AML (Cox hazard ratio per SD adjusted for age and ELN2017, Kaplan-Meier by AUC thirds). |
| `10_venetoclax_score_vs_exvivo_auc_by_blast.ipynb` | Venetoclax TransCell score versus ex vivo AUC and IC50, in all patients and split by progenitor (blast) fraction. |

Notebooks 02 and 04 read `results/dist_df.csv`, which is included and is recreated when notebook 03 runs.

## Input files (`data/`)

| File | Used by |
| --- | --- |
| `cluster_info.csv`, `frac_df.csv` | 01, 03 |
| `AML_high_res_deconv/` (CIBERSORTx high-resolution output for all references, including `BULK_AML_CIBERSORTX.txt`), `num_of_cells.csv` | 02 (not included; HPCC only) |
| `total_df.csv`, `pbmc_transcell_results.csv`, `mix_transcell_results.csv`, `transcell_screen_fda_approved.csv` | 03 |
| `Repurposing_Hub_export.txt` (Drug Repurposing Hub annotations) | 03 |
| `Metadata.txt` (single-cell metadata) | 03 |
| `drug_exvivo.csv` | 03 |
| `beataml2_patients_drug.csv`, `beataml2_data.csv` (Beat AML ex vivo drug response and clinical data) | 04, 07, 08, 09, 10 |
| `DECONVALID_results.csv` | 04 |
| `ALLDECONVALID_results.csv` | 07, 08, 10 |
| `cell_line_predicitons.csv` | 05 |
| `phenoDF_withage.csv`, `tcga_pan-cancer_Survival.txt`, `markers/` | 06 |
| `tcga_panc_symbol.RData` (not included, see below) | 06 |


### TCGA expression data

`tcga_panc_symbol.RData` can be downloaded from TCGA pan-cancer gene expression data, for example from [UCSC Xena](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)). Save it as an R object named `exp_matrix_symbol`: a matrix with gene symbols as rows and TCGA sample IDs as columns, with `-` in the IDs replaced by `.`.

## Requirements

- Python: pandas, numpy, scipy, matplotlib, seaborn, lifelines, scikit-learn, umap-learn and tqdm.
- R: sigatlas, ggplot2, survival, survminer.
