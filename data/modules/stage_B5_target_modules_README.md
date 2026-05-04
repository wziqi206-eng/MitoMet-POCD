# Stage B.5 target modules source note

## Authoritative source
`analysis/02_validation_GSE115440/scripts/04_module_scoring_GSE115440.R` defines Stage B.5 target modules using:
- `msigdbr::msigdbr(species = "Mus musculus")` for:
  - `HALLMARK_INFLAMMATORY_RESPONSE`
  - `HALLMARK_TNFA_SIGNALING_VIA_NFKB`
  - `HALLMARK_OXIDATIVE_PHOSPHORYLATION`
  - `KEGG_OXIDATIVE_PHOSPHORYLATION`
  - `REACTOME_RESPIRATORY_ELECTRON_TRANSPORT`
- A fixed custom vector for `MICROGLIA_ACTIVATION_CUSTOM` (Script 04 lines 300-305).

## Generator script
Use:
- `data/modules/generate_stage_B5_target_modules.R`

This script writes:
- `data/modules/stage_B5_target_modules.csv`

with columns:
- `module,gene_symbol,module_class`

## Run
```bash
Rscript data/modules/generate_stage_B5_target_modules.R
```

## Package policy
The generator **does not auto-install** packages.
If missing, install manually, e.g.:
```r
install.packages("msigdbr")
install.packages("dplyr")
install.packages("readr")
install.packages("stringr")
```

## Notes
- Do not hand-edit or fabricate module gene lists.
- Re-run the generator to refresh from authoritative sources.
