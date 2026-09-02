# Advanced analysis protocol

Analysis date: 2026-09-01

## Primary question

Does bacteriocin protein architecture predict source-defined minimum inhibitory concentration after accounting for study, assay, target, and protein-family structure?

## Inclusion criteria for literature expansion

Records must describe a purified, soluble, proteinaceous bacteriocin active against a Gram-negative bacterium, report a quantitative MIC or a concentration endpoint that the source explicitly defines as MIC, provide an amino acid sequence or resolvable protein accession, and provide enough assay detail to classify the method. Small ribosomally synthesized peptides, crude extracts, inhibition-zone-only records, tailocins, phage particles, and records without a convertible concentration are outside the expanded quantitative dataset.

## Analysis units

The observation-level analysis retains one row per bacteriocin-target-assay MIC. A protein-level sensitivity analysis uses the median MIC for each unique protein. Dependence is addressed with study-clustered and protein-clustered uncertainty estimates where estimable.

## Prespecified primary predictors

1. Total sequence length
2. Representative functional-domain length
3. IUPred3 disordered-residue fraction
4. N-terminal 40-residue disorder fraction
5. Functional-domain net charge proxy at pH 7
6. Functional-domain hydropathy

## Exploratory predictors

Whole-protein molecular mass, isoelectric point, net charge proxy at pH 7, hydropathy, aromaticity, instability index, secondary-structure propensities, amino acid composition, longest disordered segment, motif length, domain fraction, N-terminal charge, N-terminal hydropathy, family, import mechanism, assay group, and target genus.

## Model ladder

1. Unadjusted ordinary least squares with HC3 uncertainty
2. Study-clustered and protein-clustered uncertainty when at least three clusters are available
3. Robust Huber regression
4. Median quantile regression
5. Protein-median sensitivity analysis
6. Dilution-assay sensitivity analysis
7. Within-study standardized outcome analysis
8. Family and assay adjusted ordinary least squares when degrees of freedom permit
9. Leave-one-protein-out and leave-one-study-out analyses
10. Protein-level bootstrap confidence intervals and permutation tests
11. Ridge and elastic-net prediction with nested leave-one-protein-out validation, reported only as prediction performance

## Multiplicity and interpretation

The six primary predictors are controlled with the Benjamini-Hochberg procedure. Exploratory predictors are controlled in a separate family. Effect sizes, confidence intervals, stability, and out-of-sample performance take precedence over nominal P values. No model is selected solely because it produces a smaller P value.

## Escalation rule

If the expanded literature search does not yield enough comparable records for stable multivariable inference, the report will broaden mechanistically through family, import, domain, receptor, and assay annotations while retaining the quantitative MIC analysis as exploratory. Missing evidence will not be imputed as an observed result.
