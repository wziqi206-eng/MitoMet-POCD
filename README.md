# MitoMet-POCD

A computational biology project integrating mitochondrial dysfunction, metabolic disruption, oxidative inflammation, and synaptic plasticity in postoperative cognitive dysfunction.

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

1. **NAD+ / SIRT1 / oxidative stress**  
   NMN pretreatment has been reported to improve isoflurane-induced cognitive impairment by restoring NAD+ related signaling and reducing oxidative stress.

2. **PPARα / fatty acid oxidation**  
   Fenofibrate, a PPARα agonist, has been reported to protect against isoflurane-induced cognitive dysfunction by enhancing fatty acid oxidation.

3. **Insulin / TGF-beta / CREB-NR2B signaling**  
   Insulin and TGF-beta signaling may converge on CREB, which regulates NR2B-related synaptic plasticity in POCD.

4. **Mitochondrial ROS / inflammation**  
   Mitochondrial oxidative stress and inflammatory signaling, including MAPK and NF-kB pathways, may represent shared neural injury mechanisms.

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

```text
references/gene_sets/pocd_modules_v1.csv
