# Output files

## Main multidatabase analysis

Expected key outputs from `01_pathway_observation_readiness_multidb.py` include:

- `Run_status_all_cancers_MultiDB.csv`
- `Table_Main_AIDOh_Biology_II_AllCancers_MultiDB_full_results.csv`
- `Table_1_cancer_database_observation_readiness_summary.csv`
- `Table_2_database_overall_observation_readiness_summary.csv`
- `Table_3_before_after_AIDOh_filter_summary.csv`
- `Table_4_low_resolution_and_near_unobservable_terms.csv`
- `Table_5_top_BP_cancer_by_D_survival.csv`
- `Table_full_results_Hallmark.csv`
- `Table_full_results_GO_BP.csv`
- `Table_full_results_Reactome.csv`
- `per_cancer_results/`
- `per_database_results/`
- `bp_scores/`

## Size-bin random baseline

Expected key outputs from `03_sizebin_random_baseline_T300.py` include:

- `Selected_observation_ready_terms_with_size_bins.csv`
- `Table_SizeBinRandom_T300_bin_level_baselines.csv`
- `Table_SizeBinRandom_T300_term_level_mapped_results.csv`
- `Table_SizeBinRandom_T300_database_summary.csv`
- `Table_SizeBinRandom_T300_cancer_database_summary.csv`
- `Table_SizeBinRandom_T300_top_structured_favored_terms.csv`
- `Table_SizeBinRandom_T300_highD_but_not_above_random_q95.csv`
- `Run_status_sizebin_random_T300.csv`

## Supplementary tables

The repository generates source result tables. The optional exporter script can copy,
rename, and split these outputs into a supplementary-table pack. Some formatted
submission-ready supplementary tables may require manual formatting outside the pipeline.
