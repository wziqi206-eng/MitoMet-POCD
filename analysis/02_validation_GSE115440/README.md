# GSE115440 Validation Analysis

## Role in project

GSE115440 is used as a Tier 2 validation dataset for the MitoMet-POCD / NeuroMitoMap project.

The main goal is not to chase single-gene DEGs, because the dataset is small. Instead, this analysis will test whether mitochondrial dysfunction and inflammatory activation modules show directionally consistent changes compared with the Tier 1 anchor dataset GSE95426.

## Dataset

- Accession: GSE115440
- Species: Mus musculus
- Tissue: hippocampus
- Disease/model context: perioperative neurocognitive disorder / postoperative neuroinflammation model
- Platform: Affymetrix Mouse Gene 1.0 ST Array / GPL11533
- Main use: independent validation of module-level mitochondrial and inflammatory signatures

## Analysis strategy

1. Download and inspect GEO metadata.
2. Normalize the microarray data independently.
3. Build a sample annotation table.
4. Run limma DEG analysis as secondary evidence.
5. Compute module-level scores for mitochondrial and inflammatory modules.
6. Compare module directionality with GSE95426.
7. Report concordance using direction tables and heatmaps.

## Important caveat

This dataset should not be merged directly with GSE95426 at the raw expression level because the platforms and study designs differ. Cross-dataset comparison should happen at the module-score level.

## Expected outputs

- processed expression matrix
- verified metadata table
- limma DEG table
- module score matrix
- group-level module summary
- GSE95426 vs GSE115440 concordance table
- validation figures
