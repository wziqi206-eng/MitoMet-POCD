# MitoMet-POCD

A computational biology project integrating mitochondrial dysfunction, metabolic disruption, oxidative inflammation, and synaptic plasticity in postoperative cognitive dysfunction.

---

## 1. Project Overview

Postoperative cognitive dysfunction (POCD) is a clinically important neurocognitive complication following anesthesia and surgery. Current studies suggest that POCD is not driven by a single pathway, but may involve interacting abnormalities in mitochondrial function, oxidative stress, metabolic regulation, inflammatory signaling, and synaptic plasticity.

This project, **MitoMet-POCD**, aims to build an interpretable computational network model of POCD by integrating literature-curated mechanisms and public transcriptomic datasets.

The long-term goal is to identify candidate molecular drivers and pathway-level signatures that may help guide future experimental validation.

---

## 2. Research Question

Can fragmented mechanisms of postoperative cognitive dysfunction be integrated into a unified computational network to prioritize candidate molecular drivers?

---

## 3. Core Hypothesis

POCD may be driven by a mitochondrial-metabolic-synaptic network involving:

- NAD+ / SIRT1 metabolism
- Fatty acid oxidation
- Mitochondrial respiratory chain dysfunction
- Oxidative stress
- Inflammatory signaling
- Insulin signaling
- TGF-beta / SMAD signaling
- CREB-NR2B-mediated synaptic plasticity

Instead of treating these pathways as independent mechanisms, this project explores whether they form an interconnected computational network.

---

## 4. Biological Rationale

This project is based on four major mechanistic directions:

### 4.1 NAD+ / SIRT1 / Oxidative Stress

NMN pretreatment has been reported to improve isoflurane-induced cognitive impairment by restoring NAD+ related signaling and reducing oxidative stress.

This supports the inclusion of the **NAD_SIRT1** and **OXIDATIVE_STRESS** modules.

### 4.2 PPARα / Fatty Acid Oxidation

Fenofibrate, a PPARα agonist, has been reported to protect against isoflurane-induced cognitive dysfunction by enhancing fatty acid oxidation.

This supports the inclusion of the **FAO** module.

### 4.3 Insulin / TGF-beta / CREB-NR2B Signaling

Insulin and TGF-beta signaling may converge on CREB, which regulates NR2B-related synaptic plasticity in POCD.

This supports the inclusion of the **INSULIN_SIGNALING**, **TGF_BETA**, and **CREB_NR2B_SYNAPSE** modules.

### 4.4 Mitochondrial ROS / Inflammation

Mitochondrial oxidative stress and inflammatory signaling, including MAPK and NF-kB pathways, may represent shared neural injury mechanisms.

This supports the inclusion of the **MITO_RESPIRATORY_CHAIN**, **OXIDATIVE_STRESS**, and **INFLAMMATION** modules.

---

## 5. Core Pathway Modules

The first version of this project includes eight curated gene modules:

| Module | Biological Focus |
|---|---|
| NAD_SIRT1 | NAD+ metabolism and SIRT1-related redox regulation |
| FAO | Fatty acid oxidation and PPARα signaling |
| MITO_RESPIRATORY_CHAIN | Mitochondrial respiratory chain and oxidative phosphorylation |
| OXIDATIVE_STRESS | ROS regulation and antioxidant defense |
| INFLAMMATION | NF-kB, MAPK, and inflammatory cytokine signaling |
| INSULIN_SIGNALING | Insulin receptor and downstream metabolic signaling |
| TGF_BETA | TGF-beta / SMAD signaling |
| CREB_NR2B_SYNAPSE | CREB, NR2B, and synaptic plasticity |

The curated module table is available at:

`references/gene_sets/pocd_modules_v1.csv`

---

## 6. Project Workflow

The project will be developed through the following workflow:

1. Literature curation
2. Gene module construction
3. Public dataset collection
4. Expression matrix preprocessing
5. Differential expression analysis
6. Pathway/module scoring
7. Network and hub gene analysis
8. Interpretable machine learning
9. Technical report and application portfolio

---

## 7. Current Repository Structure

- `README.md`
- `data/`
  - `metadata/`
    - `dataset_registry.csv`
  - `processed/`
  - `raw/`
- `notebooks/`
- `references/`
  - `literature_matrix.csv`
  - `gene_sets/`
    - `pocd_modules_v1.csv`
- `report/`
- `results/`
- `src/`

---

## 8. Current Core Files

| File | Description |
|---|---|
| `references/literature_matrix.csv` | Curated evidence table from key mechanism papers |
| `references/gene_sets/pocd_modules_v1.csv` | First version of POCD-related gene modules |
| `data/metadata/dataset_registry.csv` | Public datasets planned for analysis |

---

## 9. Planned Public Datasets

This project does not treat any single public dataset as the main source of discovery. Instead, datasets are used in different roles within a cross-dataset validation framework.

| Dataset | Role in This Project | Priority |
|---|---|---|
| GSE95426 | Initial pilot dataset for workflow testing and module-level analysis | A |
| GSE178995 | Mitochondrial and respiratory-chain validation dataset | A |
| GSE115440 | Inflammation and NF-kB validation dataset | B |
| GSE95070 | Optional miRNA regulatory-layer dataset | B |
| GSE161340 | Optional aging hippocampus cell-type reference dataset | C |

The goal is not to rediscover differentially expressed genes from a single dataset, but to test whether literature-curated POCD-related modules show coordinated dysregulation across multiple POCD/PND-related datasets.

---

## 10. Planned Analyses

### 10.1 Differential Expression Analysis

The first analysis step will compare POCD-like samples with control samples to identify differentially expressed genes.

Expected outputs:

- DEG table
- Volcano plot
- Heatmap of top differentially expressed genes

### 10.2 Pathway and Module Scoring

Each sample will be scored across the eight curated biological modules.

Expected outputs:

- Module score matrix
- Module score boxplots
- Module correlation heatmap

### 10.3 Network Analysis

Candidate molecular drivers will be prioritized using:

- Literature support
- Differential expression
- Protein-protein interaction network centrality
- Cross-dataset consistency

Expected outputs:

- Hub gene ranking
- POCD mechanism network
- Candidate driver table

### 10.4 Interpretable Machine Learning

Machine learning models will be trained using pathway/module scores as features.

Planned models:

- Logistic Regression
- Random Forest
- XGBoost

Interpretation methods:

- Feature importance
- SHAP analysis

Expected outputs:

- Classification performance table
- SHAP summary plot
- Ranked module importance

---

## 11. Expected Outputs

This project aims to generate:

- A literature-derived POCD mechanism network
- A curated POCD-related gene module library
- Differential expression results from public datasets
- Pathway/module score visualizations
- Hub gene and candidate driver ranking
- Interpretable machine-learning analysis
- A technical report suitable for research portfolio use

---

## 12. Project Status

Current status: **MVP construction phase**

Completed:

- Repository created
- Project folder structure established
- Literature evidence table created
- First version of pathway modules created
- Dataset registry created

Next steps:

1. Download and process GSE95426
2. Organize sample metadata
3. Run first-pass differential expression analysis
4. Generate first volcano plot and heatmap
5. Calculate pathway/module scores

---

## 13. Notes

This project is intended as a computational hypothesis-generating study. It does not claim to prove causal mechanisms. Instead, it aims to prioritize candidate pathways and molecular drivers for future experimental validation.
