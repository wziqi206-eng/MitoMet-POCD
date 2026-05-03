# GSE95426 Metadata Audit

## Purpose

This note freezes the sample-level metadata used for the GSE95426 anchor analysis in the MitoMet-POCD / NeuroMitoMap project.

## Local expression matrix audit

The local cleaned expression matrix was inspected:

- File: `data/processed/GSE95426_expression_matrix_clean_probe_level.csv`
- Shape: 97,904 probes × 13 columns
- Columns: `ID_REF` plus 12 GSM sample columns
- GSM range: GSM2510565–GSM2510576

Therefore, the current repository-level primary analysis uses 12 hippocampal samples.

## Frozen sample grouping

- Control: GSM2510565–GSM2510570, labeled Hippocampal C1–C6
- POCD/Surgery: GSM2510571–GSM2510576, labeled Hippocampal T1–T6

## Frozen sample count

- Control: n = 6
- POCD/Surgery: n = 6

## Dataset role

GSE95426 is retained as the Tier-1 hippocampal POCD anchor dataset.

## Biological metadata

- Tissue: hippocampus
- Species: Mus musculus
- Sex: male
- Age: 12–14 months
- Model: tibial fracture POCD model
- Platform: Arraystar mouse lncRNA + mRNA microarray

## Reason for this audit

Earlier secondary reports contain inconsistent sample-count descriptions for GSE95426, including n=3 vs 3 and n=6 vs 6. The local expression matrix used in this repository contains 12 hippocampal samples labeled C1–C6 and T1–T6. Therefore, this project freezes the actual repository-level sample manifest before cross-dataset concordance with GSE115440.

## Analysis rule

GSE95426 and GSE115440 must be processed independently. Raw expression matrices must not be merged across datasets. Cross-dataset comparison will be performed only at the module-summary / directionality level.

## Interpretation rule

GSE95426 serves as the Tier-1 anchor. GSE115440 serves as the Tier-2 validation dataset. Inflammatory modules are the primary reproducible signal. Mitochondrial/OXPHOS modules are treated as secondary, weak/mixed, and hypothesis-generating.
