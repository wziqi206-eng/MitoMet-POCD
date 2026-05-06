cat > AGENTS.md << 'EOF'
# AGENTS.md — MitoMet-POCD / NeuroMitoMap

## Project Identity

This repository is for:

**MitoMet-POCD / NeuroMitoMap**

Project positioning:

A conservative module-level public-data audit of small-n mouse hippocampal POCD/PND transcriptomic datasets.

This project is NOT:
- a mechanism proof project
- a disease mechanism engine
- an AI/ML prediction project
- a causal inference project
- a biomarker discovery paper
- a clinical translation claim

This project DOES emphasize:
- reproducibility
- conservative interpretation
- module-level analysis
- QC and annotation caveats
- cross-dataset comparison
- transparent negative-result handling
- application-ready documentation

## Current Stage

The project is in:

**Application Packaging Phase**

Scientific exploration should be treated as mostly frozen.

Do not expand scope.
Do not add exploratory ML.
Do not introduce new biological claims.
Do not reinterpret existing results unless explicitly asked.

## Non-Negotiable Scientific Boundaries

Never write or imply:
- proof of mechanism
- causal pathway
- therapeutic target
- robust biomarker
- clinical relevance
- strong statistical validation
- high-confidence disease mechanism
- AI-driven biological discovery

Avoid wording such as:
- demonstrates mechanism
- reveals causal pathway
- identifies therapeutic target
- strong inflammatory activation
- robust cross-dataset validation
- disease-driving axis
- predictive AI system

Preferred wording:
- conservative audit
- module-level comparison
- descriptive analysis
- partial rank-level concordance
- weak or limited directional concordance
- hypothesis-generating only
- public-data reproducibility project
- small-n caveat
- annotation-caveat handling
- not suitable for strong mechanistic inference

## Allowed Work

You may:
- inspect repository structure
- inspect file contents
- create or edit markdown documentation
- improve README files
- generate conservative summaries
- check consistency across docs
- grep for risky overclaim language
- create small helper scripts for consistency/QC checks
- validate paths and file existence
- explain errors
- propose minimal fixes
- improve logging and comments
- prepare reviewable diffs

## Disallowed Work

You must not:
- fabricate results
- fabricate statistics
- fabricate citations
- fabricate GEO metadata
- fabricate PMIDs
- invent biological interpretation
- silently strengthen claims
- introduce AI/ML framing
- overwrite raw data
- delete files
- rename major directories
- commit or push unless explicitly instructed

## Mandatory Preflight

Before editing files, run:

```bash
git status --short
git log --oneline -5