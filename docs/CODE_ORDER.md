# Code order

Recommended order for the AIDO-h-Biology II workflow.

## 1. Main all-cancer multidatabase analysis

```bash
python scripts/01_pathway_observation_readiness_multidb.py
```

Generates observation-readiness classes, matched-gene summaries, BP/pathway scores,
survival discriminability summaries, and full Hallmark/GO BP/Reactome result tables.

## 2. Random-baseline pilot

```bash
python scripts/02_random_baseline_pilot_convergence.py
```

Runs repeat-number pilot analyses for T = 50, 100, 200, 300, 400, and 500.

## 3. Final size-bin random baseline

```bash
python scripts/03_sizebin_random_baseline_T300.py
```

Generates the final T=300 size-bin matched random baseline.

## 4. Reviewer-defense analyses

```bash
python scripts/04_reviewer_defense_analyses.py
```

Generates robustness and reviewer-defense analyses, including matched-gene threshold
sensitivity, before/after filtering summaries, event-count sensitivity, low-resolution
characterization, redundancy analysis, and combined-cohort sensitivity.

## 5. Figures

```bash
python scripts/05_generate_main_manuscript_figures_2_4.py
python scripts/06_generate_supplementary_figure_S2_robust.py
```

Regenerates manuscript and supplementary figures from output tables.

## 6. Optional supplementary-table exporter

```bash
python scripts/07_export_supplementary_tables.py ^
  --main-root D:/AIDO-Temp/AIDO-h-Biology-II-AllCancers-MultiDB ^
  --random-root D:/AIDO-Temp/AIDO-h-Biology-II-RandomGlobal-SizeBin-T300 ^
  --reviewer-root D:/AIDO-Temp/AIDO-h-Biology-II-ReviewerAnalyses ^
  --out-dir D:/AIDO-Temp/AIDO-h-Biology-II-SupplementaryTables
```
