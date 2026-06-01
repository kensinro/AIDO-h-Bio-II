# AIDO-h-Biology II code

This repository contains the Python analysis code for the manuscript:

**Observation-readiness assessment refines pathway-level molecular interpretation in cancer transcriptomics**

The code implements a representation-first workflow for pathway-level transcriptomic analysis:
observation-readiness assessment, endpoint-specific discriminability analysis, random-baseline calibration,
reviewer-defense analyses, and manuscript/supplementary figure generation.

## Repository structure

```text
scripts/
  01_pathway_observation_readiness_multidb.py
  02_random_baseline_pilot_convergence.py
  03_sizebin_random_baseline_T300.py
  04_reviewer_defense_analyses.py
  05_generate_main_manuscript_figures_2_4.py
  06_generate_supplementary_figure_S2_robust.py
  07_export_supplementary_tables.py

config/
  config_paths_template.json
  PATHS_AND_INPUTS.md

docs/
  CODE_ORDER.md
  OUTPUT_FILES.md
  SUPPLEMENTARY_TABLES.md
  METHODS_SUMMARY.md

archive/
  Original uploaded scripts, pilot scripts, legacy Hallmark-only scripts, and figure-script backups.
```

## Main workflow

1. Run `01_pathway_observation_readiness_multidb.py` to generate all-cancer, multidatabase observation-readiness and discriminability result tables.
2. Run `02_random_baseline_pilot_convergence.py` if random-repeat convergence needs to be checked.
3. Run `03_sizebin_random_baseline_T300.py` to generate the final size-bin matched random-baseline results.
4. Run `04_reviewer_defense_analyses.py` to generate threshold sensitivity, top-D before/after filtering, event-count sensitivity, low-resolution characterization, redundancy, and combined-cohort sensitivity outputs.
5. Run `05_generate_main_manuscript_figures_2_4.py` and `06_generate_supplementary_figure_S2_robust.py` to regenerate manuscript and supplementary figures from the output tables.
6. Optionally run `07_export_supplementary_tables.py` to organize generated result tables into a supplementary-table pack.

## Supplementary tables

Large supplementary tables are not stored directly in this repository because of file-size limits.
The scripts generate the source result tables used to construct the supplementary tables.
Formatted full supplementary tables should be supplied with the manuscript submission or deposited in a data repository.

## Installation

```bash
pip install -r requirements.txt
```

## Notes

The scripts contain Windows-style example paths used in the original analysis. Before running,
edit the input and output paths to match your local folders. See `config/PATHS_AND_INPUTS.md`.
# AIDO-h-Bio-II
