# Bacteriocin protein architecture and MIC analysis v1.1.0

This release substantially expands and strengthens the original literature-derived analysis.

## New in version 1.1.0

- Expanded from 15 MIC observations and 11 proteins to 47 source-defined observations and 22 proteins.
- Preserved 35 exact values, eight bounded intervals, and four right-censored measurements.
- Added validated UniProtKB or NCBI Protein accessions for every analyzed protein.
- Added official IUPred3, InterProScan, and Clustal Omega analyses for all 22 sequences.
- Added protein-cluster bootstrap interval regression with false-discovery-rate control.
- Added family, assay, study-removal, protein-removal, robust-regression, quantile-regression, and protein-median sensitivity analyses.
- Added a nested leave-one-protein-out ridge model and intercept-only predictive baseline.
- Added complete provenance, literature-search, and claim-source tables.
- Added five updated main figures and three supplementary figures.

## Interpretation

Several pooled sequence-feature associations reached the prespecified threshold, but they were not stable after accounting for protein family, assay, and study structure. The multifeature model also failed held-out-protein prediction. The release therefore reports both the candidate pooled signals and the evidence that they are not yet generalizable.

## Archive note

The public release contains reproducibility materials only. The manuscript and supplementary manuscript are intentionally excluded.
