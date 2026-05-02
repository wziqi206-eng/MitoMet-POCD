# GSE115440 Cross-Dataset Validation (Tier 2)

## Purpose

This folder contains the independent processing and validation analysis of GSE115440, used as a Tier 2 cross-dataset validator for the MitoMet-POCD / NeuroMitoMap project.

The Tier 1 anchor dataset is GSE95426. GSE115440 is not merged with GSE95426 at the raw expression level. It is only compared at the module-score / directionality level.

## Dataset Summary

- GEO ID: GSE115440
- Platform: GPL11533 Affymetrix Mouse Gene 1.1 ST Array
- Species: Mus musculus
- Strain: C57BL/6
- Tissue: brain hippocampus
- Age: 12-14 weeks

## Sample Groups

| Group | Samples | Use in this analysis |
|---|---|---|
| Control | GSM3178275, GSM3178276, GSM3178277 | Primary validation contrast |
| Surgery | GSM3178278, GSM3178279, GSM3178280 | Primary validation contrast |
| Surgery+Maresin1 | GSM3178281, GSM3178282, GSM3178283 | Auxiliary / rescue only; not in primary validation |

Primary validation contrast:

Surgery vs Control only.

The Surgery+Maresin1 group is treatment/rescue auxiliary and is excluded from the primary cross-dataset POCD validation.

## Scientific Rules

1. Do not merge raw expression data across GSE95426 and GSE115440.
2. Process GSE115440 independently.
3. Compare with GSE95426 only at the module-score / directionality level.
4. DEG analysis is secondary because n=3 vs n=3.
5. Module-level mitochondrial / inflammatory directionality is the main validation layer.
6. Be conservative and avoid causal overclaiming.

## Folder Structure

analysis/02_validation_GSE115440/
- README.md
- scripts/
- results/
  - figures/
  - tables/

## Planned Workflow

The R scripts will be generated incrementally, one at a time, in the scripts/ folder.

Planned stages:

1. Download and load raw CEL files from GEO.
2. Perform QC and RMA normalization independently for GSE115440.
3. Map probes to genes for GPL11533.
4. Run the primary Surgery vs Control limma analysis.
5. Treat DEG results as secondary evidence.
6. Compute mitochondrial and inflammatory module scores.
7. Run leave-one-out robustness checks for module scores.
8. Compare GSE115440 with GSE95426 only at module-score directionality level.
9. Use Surgery+Maresin1 only as an exploratory treatment/rescue auxiliary check.
10. Save figures and tables under results/.

## Relationship to Other Project Files

- data/metadata/GSE115440_metadata_verified.csv: verified GSE115440 sample metadata.
- data/metadata/dataset_registry.csv: central dataset registry.
- data/metadata/cross_dataset_validation_plan.csv: cross-dataset validation plan.
- analysis/01_*_GSE95426/: Tier 1 anchor analysis; do not overwrite.

## Expected Outputs

Expected tables:

- normalized expression matrix
- QC summary table
- limma DEG table
- module score matrix
- module group comparison table
- leave-one-out robustness table
- cross-dataset directionality comparison table

Expected figures:

- sample boxplot
- density plot
- PCA plot
- sample correlation heatmap
- module score boxplots
- cross-dataset directionality heatmap

## Status

- [x] Folder scaffold created
- [x] README.md created
- [ ] Script 01: download and load GSE115440 raw data
- [ ] Script 02: QC and RMA normalization
- [ ] Script 03: limma Surgery vs Control DEG analysis
- [ ] Script 04: GSVA / ssGSEA module scoring
- [ ] Script 05: leave-one-out module robustness
- [ ] Script 06: comparison with GSE95426 directionality
- [ ] Final figures and tables compiled

## Notes

This analysis is designed as part of a GitHub technical report / PhD application research portfolio.

Main claim level:

cross-dataset module-level directional consistency

Not main claim level:

single-gene causal mechanism
