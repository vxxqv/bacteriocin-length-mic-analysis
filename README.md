# Bacteriocin protein architecture and MIC analysis

This repository contains the reproducibility package for a literature-derived analysis of bacteriocin protein architecture and minimum inhibitory concentration.

## Scope

Version 1.1 expands the dataset to 47 source-defined MIC observations from 22 validated proteins, nine primary studies, four protein families, and six assay groups. It preserves exact values, bounded intervals, and right-censored observations without replacing them with artificial point estimates.

The prespecified predictors are total protein length, InterProScan-supported functional-region length, IUPred3 disorder fraction, N-terminal disorder fraction, functional-region charge at pH 7, and functional-region hydropathy.

## Main conclusion

Pooled interval models identified inverse associations for total length, global disorder, N-terminal disorder, and functional-region charge. These associations were attenuated, reversed, or rendered uncertain after family, assay, and study sensitivity analyses. A nested leave-one-protein-out ridge model performed worse than an intercept-only baseline. The current literature therefore does not support a generalizable architecture-potency rule independent of protein family and assay structure.

## Repository contents

- `analysis_package/code` contains the comment-free analysis scripts, the prespecified protocol, and dependency requirements.
- `analysis_package/data` contains validated sequences and compact bioinformatics outputs.
- `analysis_package/tables` contains the expanded MIC dataset, accession table, model results, robustness analyses, literature-search log, and claim-source ledger.
- `analysis_package/figures` contains five main figures and three supplementary figures.
- `analysis_package/bioinformatics` contains official IUPred3, InterProScan, and Clustal Omega outputs and provenance.
- `CITATION.cff` and `.zenodo.json` contain citation and archive metadata.

The manuscript and supplementary manuscript are intentionally excluded from this public repository.

## Reproduction

1. Install Python 3.12.
2. Create and activate a virtual environment.
3. Install the packages in `analysis_package/code/requirements.txt`.
4. Run `analysis_package/code/advanced_analysis.py` from the repository root.

The service-refresh scripts for IUPred3, InterProScan, and Clustal Omega require internet access and depend on the current official service interfaces.

## Bioinformatics tools

- IUPred3 in long-disorder mode with the default 0.5 threshold
- InterProScan 5.78-109.0 through the EMBL-EBI Job Dispatcher
- Clustal Omega 1.2.4 through the EMBL-EBI Job Dispatcher
- Biopython ProtParam for physicochemical sequence descriptors

The associated study cites the publications for all software and public biological databases used.

## Citation

Patni, V., & Atanaskovic, I. (2026). Bacteriocin protein length and minimum inhibitory concentration analysis (Version 1.0.0) [Computer software]. Zenodo. https://doi.org/10.5281/zenodo.22065859
