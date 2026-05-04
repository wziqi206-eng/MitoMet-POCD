# A Conservative Module-Level Audit of POCD/PND Mouse Hippocampal Transcriptomes

**Subtitle:** Annotation-rescue caveats, random-control sensitivity QC, and partial cross-dataset module concordance

**Project:** MitoMet-POCD / NeuroMitoMap  
**Report status:** v0.1 packaging draft  
**Current phase:** Application-facing technical report  
**Primary evidence layer:** predefined module-level scoring and cross-dataset comparison  
**Data policy:** independent dataset processing only; no raw expression merging across datasets
**Conservative reading rule:** all biological interpretations in this report are descriptive, module-level, and caveat-bound. The report prioritizes transparent handling of discordant findings over a simplified pathway story.

---

## 1. Executive summary

This project reanalyzes two public mouse hippocampal POCD/PND-related transcriptomic datasets using a conservative module-level workflow. The goal is not to assert a clean disease pathway story, but to document what remains supportable after independent processing, annotation caveats, random-control sensitivity checks, and cross-dataset comparison.

The final result is mixed and informative. GSE115440 shows clearer inflammatory module activation in Surgery vs Control. GSE95426 does not show inflammatory activation in POCD vs Control; instead, its inflammatory Hallmark modules shift downward and the custom microglia module is near-null. GSE95426 also shows strong OXPHOS / respiratory electron transport down-shifts, but those signals are constrained by a broad random-control downward bias in the same dataset.

The cross-dataset result is therefore best described as **partial rank-level concordance with weak directional concordance**, not as a clean shared pathway result.

**Interpretation guardrail:** This section does not state that both datasets support inflammatory activation.

---

## 2. Scope and boundaries

This report is a technical audit of module-level public-data reuse. It focuses on:

1. GSE95426 POCD vs Control module contrasts.
2. GSE95426 random-control sensitivity QC.
3. Cross-dataset comparison between GSE95426 and GSE115440 at the module-effect-size and direction level.
4. A final rank/sign concordance sensitivity layer.

This report does not use raw-expression merging across datasets. It does not treat single-gene DEG as the primary evidence layer. It does not use GSE174412. It does not use Surgery+Maresin1 as part of the primary cross-dataset comparison. It does not present causal conclusions.

**Interpretation guardrail:** This section restricts the project to conservative module comparison and avoids overextending the result.

---

## 3. Source-of-truth files used in this report

This report uses only the following project outputs:

### GSE95426 Script 11

- `analysis/01_anchor_GSE95426/results/tables/GSE95426_module_pocd_vs_control_results.csv`
- `analysis/01_anchor_GSE95426/results/tables/GSE95426_module_directionality_summary.csv`
- `analysis/01_anchor_GSE95426/results/tables/GSE95426_random_control_sensitivity_QC.csv`

### GSE95426 Script 11.5

- `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_random_control_summary.csv`
- `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_target_vs_random_empirical_QC.csv`
- `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_module_interpretation_flags.csv`
- `analysis/01_anchor_GSE95426/results/tables/GSE95426_script11_5_interpretation_caveat.txt`

### Cross-dataset Script 05

- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_module_concordance_table.csv`
- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_directionality_summary.csv`
- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_interpretation_flags.csv`
- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_conservative_summary.txt`

### Rank/sign Script 06

- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.csv`
- `analysis/03_cross_dataset_concordance/results/tables/cross_dataset_rank_concordance.txt`

Older planning memos, earlier biological narratives, and external interpretation notes are not used as data sources for this report. They may explain project history, but they do not override the Script 11, 11.5, 05, and 06 outputs.

### Reproducibility note

The report is designed to be traceable from committed repository outputs. Each major interpretation is tied to a generated table or text output from Script 11, Script 11.5, Script 05, or Script 06. The repository should retain script files, result tables, run logs, and session information where available so that an outside reader can audit how each conclusion was produced.

Key reproducibility expectations before application-facing release:

- Script 11, Script 11.5, Script 05, and Script 06 should remain committed.
- Generated CSV and TXT outputs used by this report should remain committed.
- Random seeds used in sensitivity steps should remain visible in scripts.
- Package/session information should remain available through run logs or output text.
- README should later summarize this report, not replace it.

**Interpretation guardrail:** This section makes the report output-driven rather than memo-driven.

---

## 4. Analysis design

The project uses predefined biological modules instead of relying on small-n single-gene DEG lists.

The six biological target modules are:

| Module | Module family |
|---|---|
| HALLMARK_INFLAMMATORY_RESPONSE | inflammatory |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | inflammatory |
| MICROGLIA_ACTIVATION_CUSTOM | inflammatory |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | mitochondrial |
| KEGG_OXIDATIVE_PHOSPHORYLATION | mitochondrial |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | mitochondrial |

For each dataset, expression values were processed independently. Cross-dataset comparison was performed only after module-level summaries had been generated within each dataset.

The central comparison quantities are module-level direction, disease-vs-control delta, Cohen's d, and interpretation flags. P-values, where present, are descriptive only.

**Interpretation guardrail:** This section defines a conservative analysis layer and does not promote single-gene or pathway-specific overreach.

---

## 5. GSE95426 module-level POCD vs Control results

GSE95426 is retained as a PASS_WITH_CAVEAT module-level dataset. In Script 11, all six biological target modules had high detected-gene coverage:

| Module | Detected / total genes | Coverage | Cohen's d | Direction |
|---|---:|---:|---:|---|
| HALLMARK_INFLAMMATORY_RESPONSE | 186 / 201 | 92.54% | -0.776 | down_in_pocd |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 187 / 199 | 93.97% | -0.775 | down_in_pocd |
| MICROGLIA_ACTIVATION_CUSTOM | 22 / 23 | 95.65% | -0.127 | near_null |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 181 / 199 | 90.95% | -1.011 | down_in_pocd |
| KEGG_OXIDATIVE_PHOSPHORYLATION | 113 / 128 | 88.28% | -0.999 | down_in_pocd |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | 139 / 157 | 88.54% | -1.064 | down_in_pocd |

The inflammatory results do not support inflammatory activation in GSE95426. Two inflammatory Hallmark modules are moderately down in POCD, and the custom microglia module is near-null. The mitochondrial/OXPHOS-related modules show stronger negative effect sizes, but these must be interpreted together with Script 11.5 random-control QC.

**Interpretation guardrail:** This section explicitly states that GSE95426 does not support inflammatory activation.

---

## 6. GSE95426 random-control sensitivity QC

Script 11.5 identified a broad downward shift among random-control modules in GSE95426:

| Metric | Value |
|---|---:|
| Random modules | 50 |
| Mean Cohen's d | -0.470 |
| Median Cohen's d | -0.443 |
| Minimum Cohen's d | -0.937 |
| Maximum Cohen's d | -0.075 |
| Down in POCD | 45 / 50 |
| Near-null | 5 / 50 |
| Up in POCD | 0 / 50 |

This creates a major interpretation caveat. GSE95426 module-level contrasts should not be interpreted as clean pathway-specific evidence because many random-control modules also shift downward.

However, the OXPHOS and respiratory electron transport modules remain notable because their absolute effect sizes are more extreme than the random-control distribution:

| Module | Cohen's d | Empirical abs percentile vs random | Stronger than random max |
|---|---:|---:|---|
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | -1.011 | 1.00 | TRUE |
| KEGG_OXIDATIVE_PHOSPHORYLATION | -0.999 | 1.00 | TRUE |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | -1.064 | 1.00 | TRUE |

The correct interpretation is therefore: GSE95426 OXPHOS / respiratory electron transport down-shift is hypothesis-generating and PASS_WITH_CAVEAT, not a clean standalone conclusion.

**Interpretation guardrail:** This section prevents the mitochondrial/OXPHOS result from being overstated.

---

## 7. GSE115440 module evidence carried into cross-dataset comparison

In the Script 05 cross-dataset table, GSE115440 shows clearer inflammatory activation in Surgery vs Control:

| Module | GSE115440 Cohen's d | GSE115440 direction | Rescue consistency |
|---|---:|---|---|
| HALLMARK_INFLAMMATORY_RESPONSE | 1.006 | up | consistent_rescue |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | 1.658 | up | consistent_rescue |
| MICROGLIA_ACTIVATION_CUSTOM | 1.446 | up | consistent_rescue |

The mitochondrial/OXPHOS-related GSE115440 results are weaker or mixed:

| Module | GSE115440 Cohen's d | GSE115440 direction | Rescue consistency |
|---|---:|---|---|
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | 0.425 | up | not_consistent |
| KEGG_OXIDATIVE_PHOSPHORYLATION | 0.141 | up | consistent_rescue |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | -0.020 | down | consistent_rescue |

Thus, GSE115440 is the clearer dataset for inflammatory module activation, while its mitochondrial/OXPHOS module results do not form a strong consistent axis.

**Interpretation guardrail:** This section separates the clearer GSE115440 inflammatory result from the weaker mitochondrial/OXPHOS result.

---

## 8. Cross-dataset module concordance

Script 05 compared the six biological modules between GSE95426 and GSE115440.

The key cross-dataset pattern is discordant:

| Module | GSE95426 d | GSE95426 direction | GSE115440 d | GSE115440 direction | Final interpretation |
|---|---:|---|---:|---|---|
| HALLMARK_INFLAMMATORY_RESPONSE | -0.776 | down | 1.006 | up | discordant_opposite_direction |
| HALLMARK_TNFA_SIGNALING_VIA_NFKB | -0.775 | down | 1.658 | up | discordant_opposite_direction |
| MICROGLIA_ACTIVATION_CUSTOM | -0.127 | near_null | 1.446 | up | partial_or_near_null |
| HALLMARK_OXIDATIVE_PHOSPHORYLATION | -1.011 | down | 0.425 | up | caveated_hypothesis_generating |
| KEGG_OXIDATIVE_PHOSPHORYLATION | -0.999 | down | 0.141 | up | caveated_hypothesis_generating |
| REACTOME_RESPIRATORY_ELECTRON_TRANSPORT | -1.064 | down | -0.020 | down | caveated_hypothesis_generating |

The cross-dataset comparison supports three conservative conclusions:

1. GSE115440 provides the clearer inflammatory activation evidence.
2. GSE95426 does not show inflammatory activation.
3. GSE95426 OXPHOS / respiratory electron transport down-shift is stronger than random controls but remains caveated because of broad random-control downward bias.

**Interpretation guardrail:** This section states partial cross-dataset concordance rather than clean agreement.

---

## 9. Rank-based concordance and sign test

Script 06 added a final lightweight robustness layer using only the six biological target modules.

Results:

| Metric | Value |
|---|---:|
| Number of biological target modules | 6 |
| Spearman rho | 0.886 |
| Bootstrap replicates | 1000 |
| Bootstrap valid replicates | 1000 |
| Bootstrap 95% CI lower bound | 0.200 |
| Bootstrap 95% CI upper bound | 1.000 |
| Direction agreement | 1 / 6 |
| Direction agreement fraction | 0.167 |
| Binomial sign test p-value | 0.21875 |
| Interpretation label | descriptive_partial_concordance_only |

This indicates partial rank-level concordance but weak directional concordance. In practical terms, the relative ordering of module effect sizes has some similarity, but the disease-vs-control directions do not align well across the two datasets.

**Interpretation guardrail:** This section uses Script 06 to strengthen the conservative interpretation, not to rescue an older story.

---

## 10. Core conclusions

The supportable conclusions are:

1. The project successfully builds a reproducible module-level audit workflow for two independently processed public POCD/PND-related mouse hippocampal transcriptomic datasets.
2. GSE115440 shows clear inflammatory module activation in Surgery vs Control.
3. GSE95426 does not show inflammatory module activation in POCD vs Control.
4. GSE95426 shows strong OXPHOS / respiratory electron transport down-shifts, but those signals are constrained by broad random-control downward bias.
5. Cross-dataset evidence is partial: rank-level concordance is present, but direction-level concordance is weak.
6. The project is strongest as a technical report demonstrating reproducible data auditing, annotation-caveat handling, predefined module scoring, random-control sensitivity QC, and transparent reporting of discordant findings.

The project should not be framed as a clean shared disease pathway result.

**Interpretation guardrail:** This section is the final guardrail against overstating the biological interpretation.

---

## 11. Limitations

### Small sample size

Both datasets have small group sizes. This makes single-gene inference unstable and motivates the module-level approach. Even module-level results must remain descriptive and conservative.

### GSE95426 annotation caveat

GSE95426 is retained under PASS_WITH_CAVEAT for predefined module-level scoring. The report does not use it for strong single-gene interpretation.

### GSE95426 broad random-control downward bias

The strongest limitation is the Script 11.5 result showing that 45 / 50 random-control modules also shift downward in POCD. This prevents clean pathway-specific interpretation of GSE95426 module down-shifts.

### Cross-dataset directional discordance

Only 1 / 6 biological target modules have direction agreement between GSE95426 and GSE115440 in Script 06. This limits any strong cross-dataset biological conclusion.

### Dataset design differences

The two datasets differ in study design and technical platform. This report handles them through independent processing and module-level comparison, but this cannot fully remove dataset-specific effects.

### Descriptive statistics

Cohen's d, module deltas, descriptive p-values, Spearman rho, bootstrap intervals, and sign tests are used as descriptive evidence layers. They do not establish causality.

**Interpretation guardrail:** This section makes the weaknesses visible instead of hiding them.

---

## 12. Falsifiers and decision rules

This section states what would weaken or overturn each core conclusion.

### Conclusion 1

**Conclusion:** The project provides a reproducible module-level public-data audit workflow.

**Falsifier:** This would be weakened if the scripts cannot be rerun from documented inputs, if output files are missing, or if the repository lacks sufficient dependency/version information.

### Conclusion 2

**Conclusion:** GSE115440 shows clearer inflammatory module activation.

**Falsifier:** This would be weakened if GSE115440 metadata labels were wrong, if Surgery and Control were reversed, or if rerunning module scoring produced inflammatory module effect sizes near zero or in the opposite direction.

### Conclusion 3

**Conclusion:** GSE95426 does not show inflammatory module activation.

**Falsifier:** This would be weakened if GSE95426 group labels were reversed, if annotation rescue changed inflammatory module coverage or direction materially, or if an independent rerun produced inflammatory modules up in POCD.

### Conclusion 4

**Conclusion:** GSE95426 OXPHOS / respiratory electron transport down-shifts are hypothesis-generating PASS_WITH_CAVEAT signals.

**Falsifier:** This would be weakened if random-control resampling showed that these modules are not more extreme than matched random modules, or if annotation changes removed target-module coverage.

### Conclusion 5

**Conclusion:** Cross-dataset concordance is partial, with rank-level similarity but weak directional agreement.

**Falsifier:** This would be weakened if Script 06 were rerun with corrected inputs and direction agreement became high across the six biological modules, or if Spearman rho became near zero under a corrected module table.

### Conclusion 6

**Conclusion:** The project is best used as a technical audit and application-facing reproducibility artifact.

**Falsifier:** This would be weakened if the repository cannot be understood by an outside reader from the report and README, or if the report overstates biological conclusions beyond the output tables.

**Interpretation guardrail:** This section defines boundaries before application materials are written.

---

## Final report-level interpretation

The strongest final framing is:

> This project is a conservative module-level audit of two public POCD/PND-related mouse hippocampal transcriptomic datasets. It demonstrates independent processing, annotation-caveat handling, predefined module scoring, random-control sensitivity QC, and transparent cross-dataset comparison. The final evidence supports partial rank-level concordance but weak directional concordance, with GSE115440 showing clearer inflammatory activation and GSE95426 showing caveated OXPHOS / respiratory electron transport down-shift under broad random-control downward bias.

This framing is suitable for a GitHub technical report and later application materials as evidence of reproducible computational biology practice, dataset auditing, and disciplined interpretation under contradictory evidence.

