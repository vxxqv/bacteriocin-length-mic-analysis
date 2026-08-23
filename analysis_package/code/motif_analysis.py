from pathlib import Path
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import statsmodels.formula.api as smf
from Bio import SeqIO

root = Path(__file__).resolve().parent
table_dir = root / "tables"
figure_dir = root / "figures"
table_dir.mkdir(exist_ok=True)
figure_dir.mkdir(exist_ok=True)

retained = [
    "Pyocin S5",
    "PaeM4",
    "LlpA BW11M1",
    "Pyocin L1",
    "Syringacin M",
    "KpneA",
    "KaerA",
    "KvarIa",
    "KvarM",
    "KpneM",
    "KpneM2"
]

protein_ids = {
    "Pyocin S5": "Pyocin_S5",
    "PaeM4": "PaeM4",
    "LlpA BW11M1": "LlpA_BW11M1",
    "Pyocin L1": "Pyocin_L1",
    "Syringacin M": "Syringacin_M",
    "KpneA": "KpneA",
    "KaerA": "KaerA",
    "KvarIa": "KvarIa",
    "KvarM": "KvarM",
    "KpneM": "KpneM",
    "KpneM2": "KpneM2"
}

records = list(SeqIO.parse(root / "retained_bacteriocins.fasta", "fasta"))

accessions = pd.read_csv(root / "retained_accession_audit.csv")

iupred = pd.read_csv(root / "retained_iupred3_summary.csv")

columns = [
    "protein", "md5", "length", "analysis", "signature_accession", "signature_description",
    "start", "end", "score", "status", "date", "interpro_accession", "interpro_description",
    "go_terms", "pathways"
]
annotations = pd.read_csv(root / "interproscan_retained_results.tsv", sep="\t", header=None, names=columns)

rules = {
    "Pyocin S5": ("Pfam", "PF01024", "pore-forming domain", "single"),
    "PaeM4": ("Gene3D", "G3DSA:3.30.450.400", "ColM catalytic domain", "single"),
    "LlpA BW11M1": ("SMART", "SM00108", "tandem lectin domains", "sum"),
    "Pyocin L1": ("SMART", "SM00108", "tandem lectin domains", "sum"),
    "Syringacin M": ("Gene3D", "G3DSA:3.30.450.400", "ColM catalytic domain", "single"),
    "KpneA": ("Pfam", "PF01024", "pore-forming domain", "single"),
    "KaerA": ("Pfam", "PF01024", "pore-forming domain", "single"),
    "KvarIa": ("Pfam", "PF01024", "pore-forming domain", "single"),
    "KvarM": ("Gene3D", "G3DSA:3.30.450.400", "ColM catalytic domain", "single"),
    "KpneM": ("Gene3D", "G3DSA:3.30.450.400", "ColM catalytic domain", "single"),
    "KpneM2": ("Gene3D", "G3DSA:3.30.450.400", "ColM catalytic domain", "single")
}

motif_rows = []
for bacteriocin in retained:
    method, signature, description, aggregation = rules[bacteriocin]
    protein = protein_ids[bacteriocin]
    hits = annotations.loc[
        annotations["protein"].eq(protein)
        & annotations["analysis"].eq(method)
        & annotations["signature_accession"].eq(signature)
    ].sort_values(["start", "end"])
    spans = [(int(row.start), int(row.end)) for row in hits.itertuples()]
    length = int(sum(end - start + 1 for start, end in spans))
    motif_rows.append({
        "bacteriocin": bacteriocin,
        "protein_id": protein,
        "annotation_method": method,
        "signature_accession": signature,
        "representative_feature": description,
        "aggregation": aggregation,
        "feature_coordinates": "; ".join(f"{start}-{end}" for start, end in spans),
        "motif_length_aa": length
    })

motifs = pd.DataFrame(motif_rows)
motifs.to_csv(table_dir / "representative_motif_lengths.csv", index=False)

features = accessions.merge(iupred, on="bacteriocin", how="left").merge(motifs, on="bacteriocin", how="left")
features.to_csv(table_dir / "retained_protein_features.csv", index=False)

mic = pd.read_csv(table_dir / "clean_primary_mic_dataset.csv")
mic = mic.drop(columns=["motif_length_aa"], errors="ignore")
mic = mic.loc[mic["bacteriocin"].isin(retained)].merge(motifs[["bacteriocin", "motif_length_aa"]], on="bacteriocin", how="left")
mic.to_csv(table_dir / "clean_primary_mic_dataset.csv", index=False)

def regression(data, outcome, analysis):
    fit = smf.ols(f"{outcome} ~ motif_length_aa", data=data).fit()
    robust = fit.get_robustcov_results(cov_type="HC3")
    rho, rho_p = stats.spearmanr(data["motif_length_aa"], data[outcome])
    ci = robust.conf_int()[1]
    return {
        "analysis": analysis,
        "outcome": outcome,
        "predictor": "motif_length_aa",
        "n_records": len(data),
        "n_proteins": data["bacteriocin"].nunique(),
        "slope": robust.params[1],
        "hc3_se": robust.bse[1],
        "ci95_low": ci[0],
        "ci95_high": ci[1],
        "hc3_p": robust.pvalues[1],
        "r_squared": fit.rsquared,
        "spearman_rho": rho,
        "spearman_p": rho_p
    }, fit

median = mic.groupby("bacteriocin", as_index=False).agg({
    "log10_mic_ug_ml": "median",
    "log10_mic_nM": "median",
    "motif_length_aa": "first"
})
direct = mic.loc[~mic["converted_from_nM"]].copy()
dilution = mic.loc[mic["assay_group"].isin(["agar dilution", "microbroth dilution"])].copy()
analyses = [
    (mic, "log10_mic_ug_ml", "primary_all_source_defined_mics"),
    (median, "log10_mic_ug_ml", "one_record_per_protein_median"),
    (direct, "log10_mic_ug_ml", "direct_mass_unit_reports"),
    (dilution, "log10_mic_ug_ml", "dilution_assay_subset"),
    (mic, "log10_mic_nM", "secondary_molar_outcome")
]
rows = []
fits = {}
for data, outcome, analysis in analyses:
    row, fit = regression(data, outcome, analysis)
    rows.append(row)
    fits[analysis] = fit
pd.DataFrame(rows).to_csv(table_dir / "motif_regression_results.csv", index=False)

sns.set_theme(style="ticks", context="paper", font_scale=1.0)
fig, ax = plt.subplots(figsize=(5.7, 4.4))
sns.scatterplot(data=mic, x="motif_length_aa", y="log10_mic_ug_ml", hue="representative_feature" if "representative_feature" in mic.columns else None, s=58, ax=ax)
grid = pd.DataFrame({"motif_length_aa": np.linspace(mic["motif_length_aa"].min(), mic["motif_length_aa"].max(), 100)})
prediction = fits["primary_all_source_defined_mics"].get_prediction(grid).summary_frame()
ax.plot(grid["motif_length_aa"], prediction["mean"], color="black", linewidth=1.2)
ax.fill_between(grid["motif_length_aa"], prediction["mean_ci_lower"], prediction["mean_ci_upper"], color="black", alpha=0.12)
ax.set_xlabel("Representative functional-domain length (aa)")
ax.set_ylabel("log10 MIC (ug/mL)")
legend = ax.get_legend()
if legend is not None:
    legend.remove()
fig.tight_layout()
fig.savefig(figure_dir / "Figure_5_motif_length_vs_MIC.png", bbox_inches="tight", dpi=300)
plt.close(fig)

provenance = {
    "source": "EMBL-EBI InterProScan 5.78-109.0 retained-protein results",
    "rule": "Pfam PF01024 for pore-forming proteins, Gene3D G3DSA:3.30.450.400 for ColM catalytic domains, and summed SMART SM00108 tandem lectin domains",
    "interpretation": "Exploratory cross-family representative functional-domain span"
}
(table_dir / "motif_analysis_provenance.json").write_text(json.dumps(provenance, indent=2), encoding="utf-8")
