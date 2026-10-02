# O11 release-asset manifest — CIX-26-0144

This manifest records the reviewer-facing public release state for the Cancer Informatics revision CIX-26-0144.

## Public release

**CIX-26-0144 O11 Reproducibility Package v1.0**  
https://github.com/kensinro/AIDO-h-Bio-II/releases/tag/CIX-26-0144-O11-v1.0

## Public S1–S18 release assets

The release contains 17 uploaded XLSX assets covering Supplementary Tables S1–S18, including a/b/c subdivisions where applicable:

- `Supplementary_Table_S1-S9.xlsx` — 120,897,192 bytes
- `Supplementary_Table_S10_top_structured_favored_terms.xlsx` — 181,851,748 bytes
- `Supplementary_Table_S11_threshold_sensitivity_overall_by_database.xlsx` — 1,616 bytes
- `Supplementary_Table_S12_threshold_sensitivity_by_cancer_database.xlsx` — 37,275 bytes
- `Supplementary_Table_S13a_event_count_correlation_summary.xlsx` — 945 bytes
- `Supplementary_Table_S13b_event_count_sensitivity_by_cancer_database.xlsx` — 7,592 bytes
- `Supplementary_Table_S14a_database_regime_summary.xlsx` — 772 bytes
- `Supplementary_Table_S14b_cancer_database_regime_summary.xlsx` — 13,654 bytes
- `Supplementary_Table_S15a_low_resolution_characterization_summary.xlsx` — 482 bytes
- `Supplementary_Table_S15b_top_highD_low_resolution_examples.xlsx` — 67,287 bytes
- `Supplementary_Table_S15c_low_resolution_full_table.xlsx` — 21,822,633 bytes
- `Supplementary_Table_S16a_redundancy_jaccard_summary_top_terms.xlsx` — 5,596 bytes
- `Supplementary_Table_S16b_redundancy_jaccard_pairwise_top_terms.xlsx` — 32,885,644 bytes
- `Supplementary_Table_S17a_combined_cohort_sensitivity_summary.xlsx` — 1,127 bytes
- `Supplementary_Table_S17b_main_results_excluding_combined_cohorts.xlsx` — 161,931,467 bytes
- `Supplementary_Table_S18a_exact_random_validation_selected_candidates.xlsx` — 159,837 bytes
- `Supplementary_Table_S18b_exact_random_validation_T1000_results.xlsx` — 53,743 bytes

## Public S19–S23 revision-analysis assets

The release also contains the reviewer-facing revision-analysis records:

- `Supplementary_Table_S19_reduced_coverage_stability_summary.csv`
- `Supplementary_Table_S20_survival_model_upgrade_summary.csv`
- `Supplementary_Table_S21_GSVA_ssGSEA_benchmark_results.csv`
- `Supplementary_Table_S22_GSE39582_EXTERNAL_GENE_MAPPING_USED.csv`
- `Supplementary_Table_S22_GSE39582_RFS_EXTERNAL_VALIDATION_RESULTS.csv`
- `Supplementary_Table_S22_GSE39582_RFS_PHENOTYPE_USED.csv`
- `Supplementary_Table_S23_BRCA_COAD_multivariable_Cox_results.csv`

Repository copies of the smaller reviewer-facing records and the O11 landing page are retained in `docs/` for direct inspection.

## Provenance controls

- MSigDB human release: `v2026.1.Hs`
- Hallmark GMT SHA256: `eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596`
- GO Biological Process GMT SHA256: `9be09dd06d6652566eb52eed530d62e6dfecc4365c1e81afd6f0b7f2e86dd4f9`
- Reactome GMT SHA256: `5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341`
- Reviewer-analysis random seed: `20260525`
- Size-bin random calibration: `T=300`
- Exact selected-candidate random validation: `T=1000`
- Reduced-coverage stability: 200 deterministic subsets per target/coverage condition
- Comparator benchmark: R 4.6.0; GSVA 2.7.17
- External validation: GSE39582 / GPL570, frozen RFS endpoint

## Integrity / limitation note

The S19 replicate-level reduced-coverage artifact was not located as a separate retained file during the package audit. The authoritative closure record and structured 200-replicate summary are public; no missing replicate-level file has been reconstructed or fabricated.

## Release status

- Public code repository: **PASS**
- Public S1–S18 release assets: **PASS (17/17 uploaded)**
- Public S19–S23 revision-analysis records: **PASS**
- Public S22 results + exact gene mapping + exact RFS phenotype: **PASS**
- O11 external reviewer delivery: **CLOSED**
