# O11 Reproducibility Package — CIX-26-0144

**Status:** public reviewer-facing reproducibility package  
**Manuscript:** *Observation-readiness assessment supports pathway-level cancer transcriptomics interpretation through task-specific discriminability analysis*  
**Journal / manuscript ID:** Cancer Informatics / CIX-26-0144

## Public release

The complete reviewer-facing supplementary-data release is:

**CIX-26-0144 O11 Reproducibility Package v1.0**  
https://github.com/kensinro/AIDO-h-Bio-II/releases/tag/CIX-26-0144-O11-v1.0

The release contains:

- 17 uploaded `.xlsx` assets covering Supplementary Tables **S1–S18** (including a/b/c subdivisions where applicable);
- Supplementary Table **S19** reduced-coverage stability summary;
- Supplementary Table **S20** survival-model upgrade summary;
- Supplementary Table **S21** GSVA/ssGSEA benchmark results;
- Supplementary Table **S22** GSE39582 external-validation results, exact gene mapping used, and exact RFS phenotype used;
- Supplementary Table **S23** BRCA/COAD multivariable Cox results.

## Public code and reviewer-facing revision outputs

The analysis code is maintained in this public repository. The repository contains the observation-readiness, random-baseline, reviewer-defense, figure-generation, and supplementary-table export workflows together with configuration and reproducibility documentation.

Repository copies of the machine-readable reviewer-facing summaries corresponding to Supplementary Tables **S19–S23** are public in `docs/`:

- `O11_S19_REDUCED_COVERAGE.csv`
- `O11_S20_SURVIVAL_UPGRADE.csv`
- `O11_S21_GSVA_SSGSEA_BENCHMARK.csv`
- `O11_S22_GSE39582_EXTERNAL_VALIDATION_RESULTS.csv`
- `O11_S23_MULTIVARIABLE_COX.csv`

For Supplementary Table **S22**, the exact supporting records used for the external-validation analysis are also public:

- `O11_S22_GSE39582_GENE_MAPPING_USED.csv`
- `O11_S22_GSE39582_RFS_PHENOTYPE_USED.csv`

These records complement the S22 results summary by exposing the mapped frozen target genes and the exact relapse-free-survival phenotype table used in the validation analysis.

## Expression and analysis provenance

- Primary processed expression input: `GE.tsv` for each analysis label.
- Analysis-label accounting: 25 processed labels = 23 primary cancer cohorts + COADREAD + LUNG.
- Endpoint-blind specimen repair: one allowed tumor specimen class per patient; non-tumor/control classes excluded; replicate averaging only within the selected specimen class.
- Exact upstream UCSC Xena dataset version/download dates that cannot be uniquely recovered from retained local filenames are reported as unrecoverable rather than inferred.

## Gene-set provenance

MSigDB human release: **v2026.1.Hs**

- Hallmark GMT SHA256: `eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596`
- GO Biological Process GMT SHA256: `9be09dd06d6652566eb52eed530d62e6dfecc4365c1e81afd6f0b7f2e86dd4f9`
- Reactome GMT SHA256: `5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341`

## Fixed analysis controls

- Primary observation-readiness threshold: K=10 effective matched genes; sensitivity K=5/10/15/20.
- Size-bin random calibration: T=300 per cancer-database-size-bin combination.
- Exact selected-candidate random validation: T=1000.
- Reviewer-analysis random seed: `20260525`.
- Reduced-coverage stability: 200 deterministic subsets per target/coverage condition; five prespecified frozen targets; no survival outcome used for target/gene/subset selection.
- O14 comparator: R 4.6.0; GSVA 2.7.17; same endpoint-blind repaired BRCA/COAD matrices; five frozen targets.
- O15 external validation: GSE39582 / GPL570; relapse-free survival frozen endpoint; no cohort, endpoint, or target replacement.

## Frozen target identities

- BRCA — HALLMARK_ESTROGEN_RESPONSE_LATE: 193 effective matched genes.
- BRCA — GOBP_REGULATION_OF_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_DNA_DAMAGE: 35.
- BRCA — REACTOME_GLYCOGEN_BREAKDOWN_GLYCOGENOLYSIS: 15.
- COAD — GOBP_NEGATIVE_REGULATION_OF_INCLUSION_BODY_ASSEMBLY: 12.
- COAD — REACTOME_TP53_REGULATES_TRANSCRIPTION_OF_CELL_CYCLE_GENES: 45.

## Retained evidence limitation

The S19 replicate-level reduced-coverage artifact was not located as a separate retained file during the final package audit. The public S19 CSV reports the authoritative 200-replicate summary statistics from the frozen closure record; no missing replicate-level file has been reconstructed or fabricated.

## Closure status

- Public code repository: **PASS**
- Public S1–S18 large-output release assets: **PASS (17/17 uploaded)**
- Public S19–S23 revision-analysis release assets: **PASS**
- Public S22 results + mapping + phenotype support records: **PASS**
- O11 external reviewer delivery: **CLOSED**
