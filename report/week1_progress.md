# Week 1 Progress

## Project Name

MitoMet-POCD

## Project Goal

This project aims to build a computational network model of postoperative cognitive dysfunction by integrating mitochondrial dysfunction, metabolic disruption, oxidative inflammation, and synaptic plasticity.

## Completed Work

- Created the GitHub repository: `MitoMet-POCD`
- Established the initial project folder structure
- Created the literature evidence table
- Created the first version of POCD-related pathway modules
- Created the public dataset registry
- Updated the README with project overview, biological rationale, workflow, planned analyses, and expected outputs

## Current Core Files

| File | Purpose |
|---|---|
| `README.md` | Main project overview |
| `references/literature_matrix.csv` | Literature-derived mechanism evidence |
| `references/gene_sets/pocd_modules_v1.csv` | First version of curated POCD-related gene modules |
| `data/metadata/dataset_registry.csv` | Planned public datasets for analysis |

## Current Biological Modules

The current version includes eight modules:

1. NAD+ / SIRT1 metabolism
2. Fatty acid oxidation
3. Mitochondrial respiratory chain
4. Oxidative stress
5. Inflammation / NF-kB / MAPK
6. Insulin signaling
7. TGF-beta / SMAD signaling
8. CREB / NR2B / synaptic plasticity

## Next Tasks

1. Download and inspect GSE95426
2. Create sample metadata for GSE95426
3. Prepare the expression matrix
4. Run first-pass differential expression analysis
5. Generate the first volcano plot and heatmap
6. Calculate module scores for the eight curated gene modules

## Notes

This project is currently in the MVP construction phase. The current goal is not to prove causality, but to build a computational hypothesis-generating framework for prioritizing POCD-related molecular pathways and candidate drivers.
