BACTERIOCIN PROTEIN LENGTH AND MIC ANALYSIS

AUTHORS

Vivaan Patni

Dr Iva Atanaskovic

AFFILIATION

Sapiente Education Program

CORRESPONDING AUTHOR

Vivaan Patni

vivaanpatni2010@gmail.com

CONTENTS

analysis_package/code contains the complete comment-free Python analysis scripts and pinned package requirements.
analysis_package/data contains the 15 source-defined MIC observations, 11 validated accessions, retained protein sequences, IUPred3 summaries, InterProScan annotations, and retained Clustal Omega alignment.
analysis_package/tables contains all processed numerical outputs used in the study.
analysis_package/figures contains the five principal figures and three supplementary figures.
analysis_package/bioinformatics contains service outputs and provenance for IUPred3, InterProScan 5.78-109.0, and Clustal Omega 1.2.4.
CITATION.cff and .zenodo.json provide repository and Zenodo archive metadata.

PRIMARY RESULTS

The primary dataset contains 15 MIC observations from 11 proteins.
Total protein length slope: 0.001354 log10 micrograms per millilitre per amino acid, HC3 P = 0.744.
IUPred3 disordered-region length slope: 0.000158, HC3 P = 0.915.
Representative functional-domain length slope: -0.001752, HC3 P = 0.760.
The representative-domain analysis is exploratory because Pfam PF01024, Gene3D G3DSA:3.30.450.400, and SMART SM00108 represent different protein families.
Median retained full-length pairwise identity: 15.74% across 55 pairs.

REPRODUCTION

1. Install Python 3.12.
2. Open PowerShell in the repository directory.
3. Run python -m venv .venv.
4. Run .\.venv\Scripts\python.exe -m pip install -r analysis_package\code\requirements.txt.
5. Run .\.venv\Scripts\python.exe analysis_package\code\analysis_pipeline.py.
6. Run .\.venv\Scripts\python.exe analysis_package\code\motif_analysis.py.
7. Run .\.venv\Scripts\python.exe analysis_package\code\retained_clustalo.py.

The scripts expect the retained data files to be in the same directory as the scripts. A ready-to-run copy of each required input is included in analysis_package/code.

To refresh IUPred3 predictions from the official service, run .\.venv\Scripts\python.exe analysis_package\code\run_iupred3.py before the analysis pipeline. This step requires internet access and depends on the official server interface.

BIOINFORMATICS TOOLS

IUPred3 was used in long-disorder mode with a 0.5 threshold.
InterProScan 5.78-109.0 supplied Pfam, Gene3D, and SMART domain annotations.
Clustal Omega 1.2.4 supplied the multiple-sequence alignment. Retained pairwise identity and a descriptive neighbor-joining tree were calculated from that alignment with Biopython 1.85.

The associated study cites the publications for IUPred3, InterProScan, the EMBL-EBI Job Dispatcher, Clustal Omega, Biopython, SciPy, statsmodels, Matplotlib, and seaborn.
