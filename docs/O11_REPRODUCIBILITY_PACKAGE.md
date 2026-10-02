# O11 Reproducibility Package — CIX-26-0144

**Status:** public reviewer-facing landing page  
**Manuscript:** *Observation-readiness assessment supports pathway-level cancer transcriptomics interpretation through task-specific discriminability analysis*  
**Journal / manuscript ID:** Cancer Informatics / CIX-26-0144

## Purpose

This page addresses Reviewer 2's reproducibility request by separating the public code repository from large machine-readable supplementary outputs while preserving the provenance required to audit the revised evidence chain.

## Public code

The analysis code is maintained in this public repository. The repository is intentionally code-focused and contains the observation-readiness, random-baseline, reviewer-defense, figure-generation, and supplementary-table export workflows together with configuration and reproducibility documentation.

## Public reviewer-facing revision outputs

The following machine-readable reviewer-facing summaries are public in `docs/`:

- `O11_S19_REDUCED_COVERAGE.csv`
- `O11_S20_SURVIVAL_UPGRADE.csv`
- `O11_S21_GSVA_SSGSEA_BENCHMARK.csv`
- `O11_S22_GSE39582_EXTERNAL_VALIDATION_RESULTS.csv`
- `O11_S23_MULTIVARIABLE_COX.csv`

## Supplementary-output inventory

The assembled O11 package contains the complete S1–S23 evidence chain: S1–S18 original supplementary outputs plus S19–S23 revision-result outputs. Several S1–S18 files are larger than the normal GitHub repository file-size limit and therefore cannot be stored as ordinary Git blobs. The source scripts in this repository regenerate the corresponding analysis outputs, and the full assembled large-output bundle is maintained separately as a revision-delivery artifact.

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

## Large-output delivery status

This public page and the public CSVs close the reviewer-facing provenance and machine-readable summary layer. The full S1–S18 large-output bundle still requires a mechanism that accepts files above ordinary GitHub blob limits, such as GitHub Release assets, or delivery through the journal revision-file mechanism. Public availability is not claimed for those large files until that upload is complete.
