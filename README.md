# Bacteriocin protein architecture and MIC analysis

This repository contains the data, code, sequences, bioinformatics outputs, and figures used to examine whether bacteriocin protein architecture is associated with minimum inhibitory concentration.

## Study at a glance

- 47 source-defined MIC observations
- 22 sequence-validated bacteriocins
- Nine primary studies, four protein families, and six assay groups
- Exact, bounded, and right-censored MIC values analyzed without forcing ranges into single values
- IUPred3, InterProScan, Clustal Omega, Biopython, interval regression, clustered bootstrapping, robustness tests, and held-out-protein prediction

Four of six prespecified protein features met false-discovery-rate control in pooled models. Greater total length, global disorder, N-terminal disorder, and functional-domain charge were associated with lower MIC. Family, assay, and study analyses showed that these are strong candidates for direct testing, but not yet independent potency predictors.

## Repository map

```text
.
├── analysis_package
│   ├── bioinformatics    Official tool outputs and provenance
│   ├── code              Analysis and service-refresh scripts
│   ├── data              Validated sequences and analysis inputs
│   ├── figures           Main and supplementary figures
│   └── tables            MIC data, accessions, models, and source ledger
├── .zenodo.json          Zenodo release metadata
├── CITATION.cff          GitHub citation metadata
├── LICENSE.txt           MIT license
├── README.md             Project overview
└── RELEASE_NOTES_v1.1.1.md
```

## Reproduce the analysis

Use Python 3.12, install `analysis_package/code/requirements.txt`, then run:

```powershell
python analysis_package/code/advanced_analysis.py
```

## Archive and citation

The archived record and its exportable citation formats are available at:

https://doi.org/10.5281/zenodo.22065859
