# Reproducibility contract

This document defines the minimum reproducibility contract for the AIDO-h-Biology II code package.

## 1. Separation of code, configuration, data, and outputs

- Version-controlled code lives in `scripts/`.
- Machine-specific paths belong in a local `config/config_paths.json` copied from `config/config_paths_template.json`.
- Raw/external datasets and licensed gene-set resources are not committed to this repository.
- Generated tables, figures, score matrices, and logs are outputs and should not be committed unless intentionally frozen as a release artifact.

## 2. Required provenance

A reproducible analysis run should record:

- repository commit SHA;
- Python version;
- package versions;
- random seed(s);
- input dataset/resource identifiers and versions;
- relevant input checksums when redistribution is not possible;
- configuration used for the run;
- start/end timestamps;
- generated output paths.

## 3. Scientific constants

Scientific thresholds and analysis constants must remain explicit and auditable. In particular, the current workflow documents:

- observation-readiness threshold: `K = 10`;
- discriminability threshold: `D = -log10(p)` with `D >= 1.301` corresponding to `p <= 0.05`;
- survival eligibility requirements and random-baseline settings as defined by the active scripts/manuscript.

Changes to scientific constants must not be introduced as repository-hygiene refactors. They require an explicit scientific rationale and validation.

## 4. Random analyses

Random-baseline scripts must use explicit seeds and report the number of repeats. Final manuscript-facing outputs should be distinguishable from pilot/convergence runs.

## 5. Configuration migration

The repository currently contains historical scripts with machine-specific Windows paths. Migration to configuration-driven execution must preserve the numerical/scientific behavior of the validated analysis. Path/config refactoring should therefore be performed separately from changes to statistical logic.

## 6. Validation gates

Before a repository upgrade is merged to `main`, verify at minimum:

1. Python syntax/compile checks for all active scripts.
2. No accidental raw data or machine-local configuration is committed.
3. Documentation matches actual active script names and execution order.
4. Configuration keys used by scripts match the template.
5. Figure/table documentation maps to generated outputs.
6. Scientific constants have not drifted unintentionally.
7. A smoke run is completed when the required external inputs are available.

A full numerical regression against frozen manuscript-facing outputs is preferred before release when the required datasets are available.
