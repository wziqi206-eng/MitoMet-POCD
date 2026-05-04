# AGENTS.md — MitoMet-POCD Repository Rules

This repository is a computational biology Bridge Project for AI+Bio / Computational Biology PhD applications.

## Locked project scope

The project is currently an application-ready GitHub technical report / Bridge Project, not a formal preprint and not a Disease Mechanism Engine.

The main biological narrative is inflammation-first:
- Neuroinflammation, microglial activation, TNF-alpha/NF-kB, and inflammatory response modules are the primary reproducible axis.
- Mitochondrial/OXPHOS transcriptional modules are weak, mixed, context-dependent, and hypothesis-generating unless proven otherwise by strict module-level concordance.

## Dataset rules

GSE95426:
- Tier-1 hippocampal POCD anchor.
- Use only as module-level evidence with annotation-rescue caveats.
- Do not treat as a clean single-gene discovery dataset.

GSE115440:
- Tier-2 priority validation.
- Primary comparison is Surgery vs Control only.
- Surgery+Maresin1 is auxiliary rescue evidence only.

GSE174412:
- Paused / conflicting_unverified.
- Do not use as hippocampal validation unless sample-level GEO tissue metadata is manually reconciled.

## Forbidden actions

Do not:
- Merge raw expression matrices across datasets.
- Treat DEG as primary evidence.
- Use Maresin1 as primary validation.
- Use GSE174412 as validation.
- Claim causal mechanism from small-n transcriptomics.
- Claim AI/ML prediction without a real held-out ML component.
- Make preprint-level claims without explicit user approval.
- Add raw data to git.
- Pop, drop, or modify git stash.
- Commit or push without user approval.
- Rewrite project history.
- Perform broad refactors unrelated to the task.

## Coding rules

- Prefer small, auditable scripts.
- Always check input file existence.
- Always write outputs to task-specific results/tables and results/figures folders.
- Always print sessionInfo() for R scripts when appropriate.
- Do not silently invent missing files or columns.
- If ambiguity exists, stop and report candidate files.

## Interpretation rules

Allowed:
- Report module-level directionality.
- Report effect sizes and concordance.
- Report negative controls transparently.
- Use conservative language.

Forbidden:
- Overclaim mitochondrial/OXPHOS reproducibility.
- Hide negative controls.
- Turn weak/mixed results into mechanism claims.
