# GSE115440 Metadata Check Plan

## Purpose

Before running differential expression or module scoring, GSE115440 sample metadata must be manually checked.

The key goal is to identify the cleanest comparison for validating the GSE95426 module-level mitochondrial and inflammatory signals.

## Questions to answer

1. Which samples are control / sham?
2. Which samples are surgery / POCD / PND model?
3. Are there Maresin 1 or other treatment groups?
4. Which samples should be used for the primary validation contrast?
5. Should treatment groups be excluded from the main validation?
6. Is hippocampus clearly confirmed for all selected samples?
7. What time point after surgery is represented?

## Expected primary contrast

Use the cleanest untreated disease-model contrast if available:

control_or_sham vs surgery_or_POCD_model

Treatment groups such as Maresin 1 should not be mixed into the primary validation contrast.

## Decision rule

- If GSE115440 has clean control vs POCD/surgery samples, use those as Tier 2 validation.
- If the dataset is mainly treatment-focused, use only the untreated model vs control samples.
- If metadata is ambiguous, pause DEG/module scoring and document the ambiguity.

## Outputs after metadata check

- data/processed/GSE115440_metadata_verified.csv
- analysis/02_validation_GSE115440/notes/GSE115440_metadata_notes.md
