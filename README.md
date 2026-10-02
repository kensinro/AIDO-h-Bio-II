# AIDO-h-Biology II code

This repository contains the Python analysis code for the manuscript:

**Observation-readiness assessment supports pathway-level cancer transcriptomics interpretation through task-specific discriminability analysis**

The code implements a representation-first workflow for pathway-level transcriptomic analysis:
observation-readiness assessment, endpoint-specific discriminability analysis, random-baseline calibration,
reviewer-defense analyses, and manuscript/supplementary figure generation.

## O11 reproducibility package

The reviewer-facing O11 reproducibility package is public.

The complete reviewer-facing supplementary-data release is:

**CIX-26-0144 O11 Reproducibility Package v1.0**  
https://github.com/kensinro/AIDO-h-Bio-II/releases/tag/CIX-26-0144-O11-v1.0

Release assets include:

- complete machine-readable Supplementary Tables **S1–S18** as 17 XLSX assets;
- Supplementary Table **S19** reduced-coverage stability summary;
- Supplementary Table **S20** survival-model upgrade summary;
- Supplementary Table **S21** GSVA/ssGSEA benchmark results;
- Supplementary Table **S22** GSE39582 external-validation results, exact gene mapping used, and exact RFS phenotype used;
- Supplementary Table **S23** BRCA/COAD multivariable Cox results.

The public O11 landing page and repository copies of reviewer-facing revision-result records are available in `docs/`:

- [`docs/O11_REPRODUCIBILITY_PACKAGE.md`](docs/O11_REPRODUCIBILITY_PACKAGE.md)
- [`docs/O11_RELEASE_ASSET_MANIFEST.md`](docs/O11_RELEASE_ASSET_MANIFEST.md)
- [`docs/O11_S19_REDUCED_COVERAGE.csv`](docs/O11_S19_REDUCED_COVERAGE.csv)
- [`docs/O11_S20_SURVIVAL_UPGRADE.csv`](docs/O11_S20_SURVIVAL_UPGRADE.csv)
- [`docs/O11_S21_GSVA_SSGSEA_BENCHMARK.csv`](docs/O11_S21_GSVA_SSGSEA_BENCHMARK.csv)
- [`docs/O11_S22_GSE39582_EXTERNAL_VALIDATION_RESULTS.csv`](docs/O11_S22_GSE39582_EXTERNAL_VALIDATION_RESULTS.csv)
- [`docs/O11_S22_GSE39582_GENE_MAPPING_USED.csv`](docs/O11_S22_GSE39582_GENE_MAPPING_USED.csv)
- [`docs/O11_S22_GSE39582_RFS_PHENOTYPE_USED.csv`](docs/O11_S22_GSE39582_RFS_PHENOTYPE_USED.csv)
- [`docs/O11_S23_MULTIVARIABLE_COX.csv`](docs/O11_S23_MULTIVARIABLE_COX.csv)

The O11 page records the fixed provenance, exact MSigDB release and GMT hashes, analysis controls, frozen target identities, external-validation configuration, and the retained S19 artifact limitation.

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
  O11_REPRODUCIBILITY_PACKAGE.md
  O11_RELEASE_ASSET_MANIFEST.md
  O11_S19_REDUCED_COVERAGE.csv
  O11_S20_SURVIVAL_UPGRADE.csv
  O11_S21_GSVA_SSGSEA_BENCHMARK.csv
  O11_S22_GSE39582_EXTERNAL_VALIDATION_RESULTS.csv
  O11_S22_GSE39582_GENE_MAPPING_USED.csv
  O11_S22_GSE39582_RFS_PHENOTYPE_USED.csv
  O11_S23_MULTIVARIABLE_COX.csv

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

Reviewer-facing Supplementary Tables S1–S23 are publicly delivered through the O11 release linked above. Repository copies of the smaller revision-result records are retained in `docs/` for direct inspection.

## Installation

```bash
pip install -r requirements.txt
```

## Notes

The scripts contain Windows-style example paths used in the original analysis. Before running,
edit the input and output paths to match your local folders. See `config/PATHS_AND_INPUTS.md`.
