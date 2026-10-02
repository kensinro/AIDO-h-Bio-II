# O11 reproducibility package (CIX-26-0144)

This repository is the **code repository** for the AIDO-h-Bio-II workflow. Large machine-readable supplementary outputs are maintained separately because several files exceed ordinary GitHub file-size limits.

For the Cancer Informatics major revision (CIX-26-0144), an O11 reproducibility package was assembled on 2026-10-02. It contains the full S1-S18 supplementary-output collection, structured reviewer-facing S19-S23 revision-result tables, original GSE39582 mapping/phenotype/validation records, and provenance/analysis-contract records.

## Fixed provenance

- Primary expression input per processed TCGA analysis label: `GE.tsv`.
- Processed labels: 25 analysis labels = 23 primary cohorts + COADREAD + LUNG.
- MSigDB release: `v2026.1.Hs`.
- Hallmark GMT SHA256: `eecaf6dad908334ae885406ec72bdc0646d8917588ed7c219fac92fc5363f596`.
- GO Biological Process GMT SHA256: `9be09dd06d6652566eb52eed530d62e6dfecc4365c1e81afd6f0b7f2e86dd4f9`.
- Reactome GMT SHA256: `5d61f289a2400cddfbb3a3353829fd2284a360bbe50f2093b566c4b7bea93341`.
- Primary observation-readiness threshold: K=10 effective matched genes; sensitivity K=5/10/15/20.
- Size-bin random calibration: T=300 per cancer-database-size-bin combination.
- Exact selected-candidate random validation: T=1000.
- Reviewer-analysis random seed: `20260525`.
- Reduced-coverage analysis: 200 deterministic subsets per target/coverage condition with no outcome-guided target/gene/subset selection.
- Comparator benchmark: R 4.6.0; GSVA 2.7.17.
- External validation: GSE39582 / GPL570; frozen relapse-free-survival endpoint.

Exact upstream UCSC Xena dataset version/download dates could not be uniquely reconstructed from retained local filenames; the revision records this limitation rather than inferring unsupported provenance.

## Package-delivery note

The large-output package is **not stored in this GitHub repository**. Reviewer 2 explicitly allowed complete outputs to be placed in a stable repository **or supplied with the revision**. The manuscript and response letter must therefore identify the final reviewer-accessible delivery mechanism accurately and must not state that the large output files are already hosted in GitHub.

## Retained evidence limitation

During the 2026-10-02 package audit, the S19 replicate-level reduced-coverage artifact was not located as a separate Drive file. The authoritative closure record and a structured summary of the reported 200-replicate results are retained; no missing replicate-level file is reconstructed or fabricated.
