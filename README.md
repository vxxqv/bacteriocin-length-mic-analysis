# Bacteriocin activity and protein sequence analysis

This repository holds the data, code, sequence predictions, and figures for our study of bacteriocin activity. The original dataset has 47 MIC records for 22 proteins. An expanded analysis adds a published six-pyocin panel tested against the same 29 isolates, plus separate L-pyocin and S6 isolate panels. MICs, inhibition zones, spot responses, and biofilm readings are analysed as different outcomes.

In the six-pyocin panel, longer proteins and proteins with more predicted disorder tended to inhibit more isolates. Both sequence features tracked each other, and the six-protein screen is exploratory. In matched numeric MICs, the shorter SX1 had a lower value than the closely related SX2 in 14 of 17 non-tied isolate pairs. The code includes every tested feature, pairwise MIC comparison, and exact test result.

## Files

```text
analysis_package/
  bioinformatics/       Original Clustal Omega and InterProScan results
  code/                 Original MIC analysis
  data/                 Original sequences
  expanded_cohorts/     New published panels, code, tool outputs, and figures
  figures/              Earlier publication figures
  tables/               Original MIC results and accession list
```

Start with [the expanded analysis guide](analysis_package/expanded_cohorts/README.md). The latest archive and citation formats are available from [Zenodo](https://doi.org/10.5281/zenodo.22065858). The published source supplements remain credited to their original authors; `LICENSE.txt` covers the code written for this repository.
