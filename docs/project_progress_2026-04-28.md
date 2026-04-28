# MitoMet-POCD Project Progress Log

Date: 2026-04-28

## Project Positioning

This repository is positioned as a mechanism-driven, cross-dataset validation, reusable AI+Bio research framework for POCD-related mitochondrial and metabolic dysregulation.

The goal is not to perform a single ordinary GEO dataset analysis, but to build a structured proof-of-readiness portfolio for PhD applications in AI+Biology, computational biology, biomedical informatics, neuroscience-AI, and related fields.

## Current Core Goal

To generate strong supplemental evidence for PhD applications by demonstrating:

1. Ability to identify a biomedical mechanism-driven research question;
2. Ability to structure public biomedical datasets into a reusable research framework;
3. Ability to perform transcriptomic preprocessing and quality control;
4. Ability to move from raw GEO data toward gene-level and pathway/module-level interpretation;
5. Ability to design a cross-dataset validation plan rather than relying on a single dataset.

## Completed Repository Components

1. README professionalized.
2. `project_positioning.md` created.
3. `cross_dataset_validation_plan.csv` created.
4. `literature_matrix.csv` created.
5. `pocd_modules_v1.csv` created.
6. `GSE95426_series_matrix.txt.gz` uploaded.
7. `GSE95426_metadata.csv` formalized.
8. `dataset_registry.csv` updated to `metadata_verified`.
9. `GSE95426_metadata_notes.md` created.
10. `src/01_extract_GSE95426_matrix.py` created, executed, and pushed.
11. `data/processed/GSE95426_expression_matrix_raw.csv` generated and pushed.
12. `data/processed/GSE95426_metadata_verified.csv` generated and pushed.
13. `src/02_prepare_GSE95426_expression_matrix.py` created, executed, and pushed.
14. `data/processed/GSE95426_expression_matrix_clean_probe_level.csv` generated and pushed.
15. `results/tables/GSE95426_qc_summary.csv` generated and pushed.
16. `results/figures/GSE95426_sample_boxplot_probe_level.png` generated and pushed.

## Current Data Status

Dataset: GSE95426  
Species: Mouse  
Tissue: Hippocampus  
Condition: POCD model versus control  
Sample size: 12 total samples  
Current matrix level: Clean probe-level expression matrix  
Current QC status: Probe-level sample QC completed  

## Next Planned Stage

The next stage has not yet been started.

Planned next steps:

1. Download and inspect GPL22782 platform annotation file.
2. Identify reliable probe-to-gene-symbol mapping columns.
3. Create a cleaned GPL22782 annotation table.
4. Convert probe-level expression matrix to gene-level expression matrix.
5. Summarize probe-to-gene mapping quality.
6. Collapse multiple probes mapping to the same gene.
7. Generate gene-level QC summary.
8. Begin mechanism-module scoring for POCD-related biological processes.

## Planned Mechanism Modules

Initial mechanism-driven modules include:

1. Mitochondrial dysfunction;
2. Oxidative stress;
3. Neuroinflammation;
4. Synaptic dysfunction;
5. Blood-brain barrier / vascular dysfunction;
6. ER stress and proteostasis;
7. Apoptosis and cell death;
8. Metabolic remodeling.

## Important Project Principle

This project should continue to be framed as:

Mechanism-driven + cross-dataset validation + reusable framework

It should not be reduced to a simple GSE95426 preprocessing exercise.

## Next Script To Build

Expected next script:

`src/03_inspect_GPL22782_annotation.py`

Purpose:

To inspect the GPL22782 platform annotation file and determine the correct columns for probe ID, gene symbol, accession, transcript ID, or other annotation fields.

