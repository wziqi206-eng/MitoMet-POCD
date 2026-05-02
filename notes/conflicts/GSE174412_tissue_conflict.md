# GSE174412 Tissue Conflict Note

## Current status

GSE174412 is temporarily paused as a hippocampal POCD validation dataset.

## Reason

Current secondary evidence suggests that GSE174412 may represent cortex tissue rather than hippocampus. Because the MitoMet-POCD / NeuroMitoMap project is currently framed around hippocampal POCD / perioperative neurocognitive dysfunction signatures, using GSE174412 as a direct hippocampal validation dataset may introduce a tissue mismatch.

## Risk

If GSE174412 is incorrectly described as hippocampus, the project may be criticized for cross-tissue overclaiming.

## Decision

Do not use GSE174412 in the main validation pipeline until GEO sample-level tissue metadata is manually checked.

## Temporary role

- Not Tier 2 validation
- Possible future exploratory / cross-tissue comparison dataset
- Status: conflicting_unverified

## Required next action

Manually inspect GEO sample metadata, including:
- tissue field
- source name
- sample title
- experimental design
- associated publication
