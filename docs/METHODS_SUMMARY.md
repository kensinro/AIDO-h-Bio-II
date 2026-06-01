# Methods summary

This code package implements the AIDO-h-Biology II workflow.

## Core concepts

1. Biological annotation defines candidate BP/pathway terms.
2. Observation-readiness checks whether each term is sufficiently represented after projection onto measured genes.
3. Endpoint-specific discriminability is quantified after observation-readiness filtering.
4. Random baselines are generated after filtering to compare structured BP/pathway terms with size-aware random controls.
5. Supplementary analyses characterize threshold sensitivity, event-count effects, low-resolution terms, redundancy, and combined-cohort sensitivity.

## Main external-facing terminology

Use:
- observation-readiness
- representation completeness
- endpoint-specific discriminability
- BP/pathway observable
- low-resolution or near-unobservable term

Avoid treating low-resolution terms as biologically unimportant. The classification only describes dataset-level representation adequacy.
