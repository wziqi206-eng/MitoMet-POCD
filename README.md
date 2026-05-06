# MitoMet-POCD / NeuroMitoMap

A conservative module-level audit of public POCD/PND-related mouse hippocampal transcriptomic datasets.

This repository documents a conservative computational biology project focused on reproducible dataset auditing, annotation-caveat handling, predefined module scoring, random-control sensitivity QC, and transparent cross-dataset comparison.

The source-of-truth report is:

**[technical_report.md](technical_report.md)**

---

## 1. Project status

Current status: **Documentation preparation**

Completed evidence chain:

1. GSE95426 annotation-rescue-aware module coverage QC
2. GSE95426 POCD vs Control module scoring
3. GSE95426 random-control sensitivity QC
4. Conservative cross-dataset module comparison with GSE115440
5. Rank-based concordance and sign-test sensitivity layer
6. Technical report source of truth

This repository should be read as a conservative technical audit, not as a clean shared-pathway claim.

---

## 2. Main finding

The final evidence is mixed.

GSE115440 shows clearer inflammatory module activation in Surgery vs Control.

GSE95426 does **not** show inflammatory activation in POCD vs Control. Its inflammatory Hallmark modules shift downward and the custom microglia module is near-null.

GSE95426 shows strong OXPHOS / respiratory electron transport down-shifts, but these are constrained by a broad random-control downward bias.

The final cross-dataset interpretation is:

> Partial rank-level concordance, weak directional concordance, and caveat-bound biological interpretation.

---

## 3. Datasets used in the primary comparison

| Dataset | Role | Tissue | Main contrast | Use in this repository |
|---|---|---|---|---|
| GSE95426 | Anchor with caveats | Mouse hippocampus | POCD vs Control | Module-level scoring only, PASS_WITH_CAVEAT |
| GSE115440 | Priority comparison dataset | Mouse hippocampus | Surgery vs Control | Independent module-level comparison |

Notes:

- Raw expression values are processed independently within each dataset.
- Cross-dataset comparison is performed only at the module-summary level.
- Surgery+Maresin1 samples in GSE115440 are not part of the primary cross-dataset comparison.
- GSE174412 is not used as hippocampal comparison evidence in this release.

---

## 4. Biological modules

The primary comparison uses six predefined biological target modules:

| Module | Module family |
|---|---|
| HALLMARK_INFLAMMATORY_RESPONSE | inflammatory |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | inflammatory |
| MICROGLIA_ACTIVATION_CUSTOM | inflammatory |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | mitochondrial |
| KEGG_OXIDATIVE_PHOSPHORYLATION | mitochondrial |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | mitochondrial |

Single-gene DEG results are not the main evidence layer.

---

## 5. Key results

### GSE95426

GSE95426 target-module coverage was high enough for predefined module-level scoring under PASS_WITH_CAVEAT.

Key Script 11 result:

- HALLMARK_INFLAMMATORY_RESPONSE: down in POCD
- HALLMARK_TNFA_SIGNALING_VIA_NFKB: down in POCD
- MICROGLIA_ACTIVATION_CUSTOM: near-null
- HALLMARK_OXIDATIVE_PHOSPHORYLATION: down in POCD
- KEGG_OXIDATIVE_PHOSPHORYLATION: down in POCD
- REACTOME_RESPIRATORY_ELECTRON_TRANSPORT: down in POCD

Key Script 11.5 result:

- 45 / 50 random-control modules also shifted downward in POCD.
- This creates a major caution for interpreting GSE95426 pathway-specific down-shifts.
- OXPHOS / respiratory electron transport modules were still more extreme than random controls, but remain hypothesis-generating only.

### GSE115440

GSE115440 shows clearer inflammatory module activation in Surgery vs Control.

Inflammatory modules show positive effect sizes and consistent rescue patterns in the available comparison table.

Mitochondrial/OXPHOS-related modules are weaker or mixed.

### Cross-dataset comparison

The two datasets do not show clean direction-level agreement.

Script 06 result:

| Metric | Value |
|---|---:|
| Spearman rho | 0.886 |
| Bootstrap 95% CI | [0.2, 1] |
| Direction agreement | 1 / 6 |
| Sign test p-value | 0.21875 |

Interpretation:

- Rank-level similarity is present.
- Direction-level agreement is weak.
- The project supports a conservative audit narrative.

---

## 6. Repository map

| Path | Purpose |
|---|---|
| `technical_report.md` | Source-of-truth technical report |
| `AGENTS.md` | Guardrails for future code and interpretation |
| `analysis/01_anchor_GSE95426/` | GSE95426 anchor analysis scripts and outputs |
| `analysis/02_validation_GSE115440/` | GSE115440 comparison analysis scripts and outputs |
| `analysis/03_cross_dataset_concordance/` | Cross-dataset comparison and rank/sign sensitivity |
| `data/metadata/` | Curated metadata and dataset registry |
| `data/modules/` | Predefined biological target modules |
| `docs/` | Project progress notes |
| `report/` | Earlier project notes and supporting documentation |
| `src/` | Earlier preprocessing scripts |

---

## 7. Key scripts

| Script | Purpose |
|---|---|
| `analysis/01_anchor_GSE95426/scripts/10_stage_B5_module_coverage_QC.R` | GSE95426 target-module coverage QC |
| `analysis/01_anchor_GSE95426/scripts/11_module_scoring_GSE95426.R` | GSE95426 module scoring and POCD vs Control contrast |
| `analysis/01_anchor_GSE95426/scripts/11_5_module_sensitivity_QC_GSE95426.R` | GSE95426 random-control sensitivity QC |
| `analysis/03_cross_dataset_concordance/scripts/05_cross_dataset_module_concordance.R` | Conservative cross-dataset module comparison |
| `analysis/03_cross_dataset_concordance/scripts/06_rank_concordance_sign_test.R` | Spearman rank concordance and sign test |

---

## 8. Key output tables

| Output | Purpose |
|---|---|
| `analysis/01_anchor_GSE95426/results/tables/GSE95426_module_pocd_vs_control_results.csv` | GSE95426 module contrast table |
| `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_random_control_summary.csv` | Random-control summary |
| `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_target_vs_random_empirical_QC.csv` | Target-vs-random empirical QC |
| `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_module_concordance_table.csv` | Main cross-dataset comparison table |
| `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_conservative_summary.txt` | Conservative cross-dataset summary |
| `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.csv` | Rank/sign sensitivity result |

---

## 9. Key figures

The `figures/` directory contains a small application-facing figure set copied from committed script outputs.

| Figure | Purpose |
|---|---|
| `figures/01_GSE95426_module_effect_direction_barplot.png` | GSE95426 module-level POCD vs Control effect directions |
| `figures/02_GSE95426_random_control_effect_size_distribution.png` | Random-control effect-size distribution showing broad GSE95426 downward bias |
| `figures/03_GSE95426_target_modules_vs_random_controls.png` | Target modules compared against random-control modules |
| `figures/04_cross_dataset_module_effect_size_heatmap.png` | Cross-dataset module effect-size comparison |
| `figures/05_cross_dataset_directionality_tileplot.png` | Direction-level agreement and disagreement across datasets |

These figures are duplicated for readability only. The original script outputs remain under the corresponding `analysis/` result folders.

---

## 10. Reproducibility notes

This project prioritizes traceability from committed scripts to committed output tables.

Reproducibility expectations:

- Scripts used for the final report are committed.
- Output CSV/TXT files used by the final report are committed.
- Random seeds used in sensitivity analyses are visible in scripts.
- Run logs or session information are retained where available.
- The README summarizes `technical_report.md`; it does not replace it.

Suggested reading order:

1. `README.md`
2. `technical_report.md`
3. `AGENTS.md`
4. Script 11 / 11.5 / 05 / 06 outputs

---

## 11. Interpretation boundaries

This repository does not claim:

- a clean shared POCD/PND pathway across both datasets
- inflammatory activation across both datasets
- a causal disease pathway
- a new biomarker signature
- prediction modeling
- raw cross-platform expression integration

The strongest value of this project is technical:

> It shows how contradictory public transcriptomic evidence can be processed, audited, stress-tested, and reported without forcing a simplified story.

---

## 12. Project value

This project is intended as a public, reproducible repository artifact for computational biology, biomedical informatics, and neuroinformatics documentation use.

It demonstrates:

- public dataset reuse
- transcriptomic workflow organization
- annotation-caveat handling
- predefined module-level analysis
- sensitivity checks with random controls
- conservative cross-dataset comparison
- transparent reporting when the original biological expectation is contradicted by data

---

## Author / Contact

Prepared by the repository maintainer.

For questions or reproducibility issues, please use the GitHub repository issue tracker or contact the maintainer through the GitHub profile associated with this repository.
