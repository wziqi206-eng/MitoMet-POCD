# Stage B.5 target modules source note

## What was inspected
- `analysis/02_validation_GSE115440/scripts/04_module_scoring_GSE115440.R`
- `analysis/02_validation_GSE115440/results/tables/GSE115440_module_gene_overlap.csv`
- `analysis/02_validation_GSE115440/results/tables/GSE115440_module_scores_mean_zscore.csv`

## Authoritative source used by Script 04
Script 04 obtains the standard target modules from **MSigDB via `msigdbr`** for mouse (`species = "Mus musculus"`):
- `HALLMARK_INFLAMMATORY_RESPONSE`
- `HALLMARK_TNFA_SIGNALING_VIA_NFKB`
- `HALLMARK_OXIDATIVE_PHOSPHORYLATION`
- `KEGG_OXIDATIVE_PHOSPHORYLATION`
- `REACTOME_RESPIRATORY_ELECTRON_TRANSPORT`

`MICROGLIA_ACTIVATION_CUSTOM` is defined explicitly inside Script 04 as a fixed custom vector.

## Why `stage_B5_target_modules.csv` is not generated here
This environment currently lacks an R runtime (`Rscript: command not found`), so `msigdbr` cannot be executed to materialize authoritative Hallmark/KEGG/Reactome gene lists into a local file.

To generate `data/modules/stage_B5_target_modules.csv` in an R-enabled environment:

```r
install.packages("msigdbr")
```

Then run a generator that queries `msigdbr::msigdbr(species = "Mus musculus")` and exports long format:
- columns: `module,gene_symbol,module_class`
- includes the five MSigDB modules above + `MICROGLIA_ACTIVATION_CUSTOM` from Script 04.

## Required next input if generation cannot be run locally
If you cannot run `msigdbr` locally, add an official exported file:
- `data/modules/stage_B5_target_modules.csv`
- containing the required modules and columns exactly.
