# Bacteriocin protein architecture and MIC analysis

This repository contains the data, sequences, bioinformatics outputs, code, tables, and final figures for an analysis of bacteriocin protein architecture and minimum inhibitory concentration.

The dataset includes 47 MIC observations for 22 sequence-validated bacteriocins from nine studies. In the pooled analysis, greater total length, global disorder, N-terminal disorder, and functional-domain charge were associated with lower MIC after false-discovery-rate correction. Additional analyses examined the effects of protein family, assay type, source study, and held-out proteins.

## Files

```text
.
├── analysis_package
│   ├── bioinformatics    Clustal Omega and InterProScan outputs with provenance
│   ├── code              Reproducible analysis and tool-refresh scripts
│   ├── data              Validated sequences and analysis inputs
│   ├── figures           Final numbered figures and editable Figure 1 PowerPoint
│   └── tables            MIC data, accessions, model results, and source ledger
├── .zenodo.json          Archive metadata
├── CITATION.cff          Citation metadata
├── LICENSE.txt
├── README.md
└── RELEASE_NOTES_v1.1.0.md
```

## Run

Install the packages in `analysis_package/code/requirements.txt`, then run:

```powershell
python analysis_package/code/advanced_analysis.py
```

The archived record provides downloadable citation formats:

https://doi.org/10.5281/zenodo.22065858
