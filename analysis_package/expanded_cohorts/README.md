# Expanded bacteriocin analysis

These files reproduce the published-data analyses in the revised study. The original 47 MIC records, the 2014 L-pyocin spot sheet, the 2026 six-pyocin panel, and the independent S6 receptor cohort are kept as separate datasets. A concentration, an inhibition zone, a spot code, and a biofilm absorbance are different measurements and are never pooled into one response.

The six-pyocin panel has 29 shared isolates and six accession-checked sequences. IUPred3 long mode supplies residue-level disorder scores, Clustal Omega supplies the alignment and percent identity matrix, InterProScan supplies domain signatures, and Biopython calculates simple sequence properties. `thai_sequence_activity_associations.csv` contains all six screened predictors, exact 720-label permutation P values, and Benjamini-Hochberg adjusted values. `thai_sequence_response_pairs.csv` compares sequence identity with response-set overlap; the permutation test shuffles protein labels rather than treating 15 pairs as independent observations.

## Reproduce

Install `requirements.txt` and run these scripts from this directory:

```powershell
python matched_analysis.py
python expanded_analysis.py
python advanced_analysis.py
python final_figures.py
python reproduce_scope_analysis.py
python paired_pyocin_mic.py --output .
python sequence_activity_analysis.py --output .
python summarize_pyocin_interproscan.py --output .
python regional_sequence_identity.py --output .
python sequence_activity_figure.py --output .
```

The archived tool outputs are enough for the analysis above. To request fresh remote predictions, run `run_clustal_omega.py --email YOUR_EMAIL`, `run_pyocin_iupred3.py --output .`, and `run_pyocin_interproscan.py --email YOUR_EMAIL --output .`. The InterProScan and Clustal Omega job IDs and raw IUPred3 scores are kept in `data/`.

The closest pair, SX1 and SX2, shares 75% identity in the full Clustal Omega matrix. The Pfam translocation region annotated in SX1 is 135/140 identical to the aligned SX2 region. The terminal Gene3D DNase-annotated region in SX1 is 19/103 identical among paired alignment columns, with 27 SX2 gap columns. InterProScan did not assign that DNase signature to SX2 at its thresholds. These annotations identify a sequence difference; they do not demonstrate why the MICs differ.

## Sources and limits

The source studies are [Ghequire et al. (2014)](https://doi.org/10.1002/mbo3.210), [Dingemans et al. (2016)](https://doi.org/10.1002/mbo3.339), and [Laosena et al. (2026)](https://doi.org/10.7717/peerj.21850). `Results_Audit.md` records the source checks and endpoint definitions. The six-protein feature screen is exploratory; length and predicted disorder are strongly correlated, so this panel cannot separate their effects. The paired MIC comparisons use only isolates with numeric values for both proteins and keep unquantified NC entries out of the calculation.

Tool references: [IUPred3](https://doi.org/10.1093/nar/gkab408), [InterProScan 5](https://doi.org/10.1093/bioinformatics/btu031), [Clustal Omega](https://doi.org/10.1038/msb.2011.75), and [Biopython](https://doi.org/10.1093/bioinformatics/btp163).
