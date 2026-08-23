from pathlib import Path
import json
import platform
import sys
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import seaborn as sns
from scipy import stats
import statsmodels.api as sm
import statsmodels.formula.api as smf

root = Path(__file__).resolve().parent
figure_dir = root / "figures"
table_dir = root / "tables"
figure_dir.mkdir(exist_ok=True)
table_dir.mkdir(exist_ok=True)
sns.set_theme(style="ticks", context="paper", font_scale=1.0)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150, "savefig.dpi": 300})

accessions = pd.read_csv(root / "retained_accession_audit.csv")
iupred = pd.read_csv(root / "retained_iupred3_summary.csv")
iupred["bacteriocin"] = iupred["bacteriocin"].replace({"LlpA BW11M1": "LlpA BW11M1"})
features = accessions.merge(iupred, on="bacteriocin", how="left")
features.to_csv(table_dir / "protein_features.csv", index=False)
mw = features.set_index("bacteriocin")["molecular_weight_da"].to_dict()

records = [
    ("Pyocin S5", "P. aeruginosa DWW3", 12.6, "ug/mL", "not reported in accessible methods", "TonB", "10.1016/j.febslet.2010.06.021"),
    ("Pyocin S5", "P. aeruginosa A19", 0.1, "ug/mL", "agar dilution", "TonB", "10.1371/journal.pone.0185782"),
    ("PaeM4", "P. aeruginosa A19", 0.6, "ug/mL", "agar dilution", "TonB", "10.1371/journal.pone.0185782"),
    ("LlpA BW11M1", "P. syringae GR12-2R3", 2.08, "nM", "spot-overlay MIC", "Unknown", "10.1371/journal.ppat.1003199"),
    ("Pyocin L1", "P. aeruginosa P8", 7.0, "nM", "spot-overlay MIC", "Unknown", "10.1371/journal.ppat.1003898"),
    ("Pyocin L1", "P. aeruginosa E2", 27.0, "nM", "spot-overlay MIC", "Unknown", "10.1371/journal.ppat.1003898"),
    ("Syringacin M", "P. syringae pv. lachrymans LMG 5456", 60.0, "nM", "spot-overlay MIC", "TonB", "10.1074/jbc.M112.400150"),
    ("KpneA", "K. variicola DSM15968", 0.2, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KaerA", "K. variicola DSM15968", 0.1, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KvarIa", "K. variicola DSM15968", 0.4, "ug/mL", "microbroth dilution", "Unknown", "10.1038/s41598-019-51969-1"),
    ("KvarM", "K. pneumoniae DSM16358", 0.1, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KvarM", "K. variicola DSM15968", 0.2, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KpneM", "K. pneumoniae DSM16358", 0.8, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KpneM", "K. variicola DSM15968", 0.4, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1"),
    ("KpneM2", "K. pneumoniae DSM16358", 0.8, "ug/mL", "microbroth dilution", "TonB", "10.1038/s41598-019-51969-1")
]
mic = pd.DataFrame(records, columns=["bacteriocin", "target_strain", "reported_value", "reported_unit", "assay_group", "import_mechanism", "doi"])
mic["molecular_weight_da"] = mic["bacteriocin"].map(mw)
mic["converted_from_nM"] = mic["reported_unit"].eq("nM")
mic["mic_ug_ml"] = np.where(mic["reported_unit"].eq("ug/mL"), mic["reported_value"], mic["reported_value"] * mic["molecular_weight_da"] / 1e6)
mic["mic_nM"] = np.where(mic["reported_unit"].eq("nM"), mic["reported_value"], mic["reported_value"] * 1e6 / mic["molecular_weight_da"])
mic["log10_mic_ug_ml"] = np.log10(mic["mic_ug_ml"])
mic["log10_mic_nM"] = np.log10(mic["mic_nM"])
mic = mic.merge(features[["bacteriocin", "verified_accession", "verified_length_aa", "idr_length_aa", "idr_fraction"]], on="bacteriocin", how="left")
mic.to_csv(table_dir / "clean_primary_mic_dataset.csv", index=False)

def regression(data, predictor, outcome, analysis):
    fit = smf.ols(f"{outcome} ~ {predictor}", data=data).fit()
    robust = fit.get_robustcov_results(cov_type="HC3")
    rho, rho_p = stats.spearmanr(data[predictor], data[outcome])
    ci = robust.conf_int()[1]
    return {
        "analysis": analysis,
        "outcome": outcome,
        "predictor": predictor,
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
    "mic_ug_ml": "median",
    "mic_nM": "median",
    "verified_length_aa": "first",
    "idr_length_aa": "first",
    "log10_mic_ug_ml": "median",
    "log10_mic_nM": "median"
})
direct = mic.loc[~mic["converted_from_nM"]].copy()
comparable = mic.loc[mic["assay_group"].isin(["agar dilution", "microbroth dilution"])].copy()
analyses = [
    (mic, "log10_mic_ug_ml", "primary_all_source_defined_mics"),
    (median, "log10_mic_ug_ml", "one_record_per_protein_median"),
    (direct, "log10_mic_ug_ml", "direct_mass_unit_reports"),
    (comparable, "log10_mic_ug_ml", "dilution_assay_subset"),
    (mic, "log10_mic_nM", "secondary_molar_outcome")
]
results = []
fits = {}
for data, outcome, analysis in analyses:
    for predictor in ["verified_length_aa", "idr_length_aa"]:
        row, fit = regression(data, predictor, outcome, analysis)
        results.append(row)
        fits[(analysis, predictor)] = fit
pd.DataFrame(results).to_csv(table_dir / "regression_results.csv", index=False)

loo_rows = []
for omitted in sorted(mic["bacteriocin"].unique()):
    subset = mic.loc[mic["bacteriocin"] != omitted]
    for predictor in ["verified_length_aa", "idr_length_aa"]:
        row, fit = regression(subset, predictor, "log10_mic_ug_ml", "leave_one_protein_out")
        loo_rows.append({"omitted_protein": omitted, **row})
pd.DataFrame(loo_rows).to_csv(table_dir / "leave_one_protein_out.csv", index=False)

coverage = pd.DataFrame([
    ("Total protein length", "completed", "Validated sequences available for all 11 retained proteins; modeled for 15 source-defined MIC records"),
    ("IUPred3 IDR length", "completed", "Official IUPred3 long-disorder predictions available for all 11 retained proteins"),
    ("Representative functional-domain length", "completed", "Exploratory Pfam, Gene3D, and SMART domain spans available for all 11 retained proteins"),
    ("Non-cytotoxic length", "not modeled", "Consistent experimentally supported cytotoxic-domain boundaries were unavailable across proteins"),
    ("Adjusted regression", "not modeled", "Fifteen records from 11 proteins were too sparse for stable target, receptor, and assay adjustment"),
    ("Homologous cytotoxic-domain maximum-likelihood tree", "not modeled", "The dataset spans unrelated toxin families and uniform domain boundaries were unavailable")
], columns=["analysis", "status", "reason"])
coverage.to_csv(table_dir / "analysis_coverage.csv", index=False)

palette = {"TonB": "#2F5597", "Unknown": "#A65E2E"}
fig, axes = plt.subplots(1, 2, figsize=(10.2, 4.2))
for ax, predictor, label in zip(axes, ["verified_length_aa", "idr_length_aa"], ["Total protein length (aa)", "IUPred3 IDR length (aa)"]):
    sns.scatterplot(data=mic, x=predictor, y="log10_mic_ug_ml", hue="import_mechanism", style="converted_from_nM", palette=palette, s=55, ax=ax)
    grid = pd.DataFrame({predictor: np.linspace(mic[predictor].min(), mic[predictor].max(), 100)})
    pred = fits[("primary_all_source_defined_mics", predictor)].get_prediction(grid).summary_frame()
    ax.plot(grid[predictor], pred["mean"], color="black", linewidth=1.2)
    ax.fill_between(grid[predictor], pred["mean_ci_lower"], pred["mean_ci_upper"], color="black", alpha=0.12)
    ax.set_xlabel(label)
    ax.set_ylabel("log10 MIC (ug/mL)")
    handles = [
        Line2D([0], [0], marker="o", color="none", markerfacecolor=palette["TonB"], markeredgecolor=palette["TonB"], label="TonB import"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor=palette["Unknown"], markeredgecolor=palette["Unknown"], label="Unknown import"),
        Line2D([0], [0], marker="o", color="none", markerfacecolor="#444444", markeredgecolor="#444444", label="Direct mass report"),
        Line2D([0], [0], marker="X", color="none", markerfacecolor="#444444", markeredgecolor="#444444", label="Converted from nM")
    ]
    ax.legend(handles=handles, fontsize=7, loc="best")
axes[0].text(-0.12, 1.03, "A", transform=axes[0].transAxes, fontweight="bold")
axes[1].text(-0.12, 1.03, "B", transform=axes[1].transAxes, fontweight="bold")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_1_length_and_IDR_vs_MIC.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5.3, 4.2))
sns.boxplot(data=mic, x="import_mechanism", y="log10_mic_ug_ml", color="white", showfliers=False, ax=ax)
sns.stripplot(data=mic, x="import_mechanism", y="log10_mic_ug_ml", hue="import_mechanism", palette=palette, jitter=0.16, size=6, ax=ax, legend=False)
ax.set_xlabel("Import mechanism")
ax.set_ylabel("log10 MIC (ug/mL)")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_2_import_mechanism_MIC.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5.4, 4.2))
sns.histplot(mic["log10_mic_ug_ml"], bins=7, kde=True, color="#2F5597", ax=ax)
ax.set_xlabel("log10 MIC (ug/mL)")
ax.set_ylabel("Count")
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_1_MIC_distribution.png", bbox_inches="tight")
plt.close(fig)

fit = fits[("primary_all_source_defined_mics", "verified_length_aa")]
influence = fit.get_influence()
fig, axes = plt.subplots(2, 2, figsize=(9.0, 7.2))
axes[0, 0].scatter(fit.fittedvalues, fit.resid, color="#2F5597")
axes[0, 0].axhline(0, color="black", linewidth=1)
axes[0, 0].set_xlabel("Fitted log10 MIC")
axes[0, 0].set_ylabel("Residual")
sm.qqplot(fit.resid, line="45", ax=axes[0, 1], markerfacecolor="#2F5597", markeredgecolor="#2F5597")
axes[0, 1].set_title("Normal Q-Q")
axes[1, 0].scatter(fit.fittedvalues, np.sqrt(np.abs(influence.resid_studentized_internal)), color="#2F5597")
axes[1, 0].set_xlabel("Fitted log10 MIC")
axes[1, 0].set_ylabel("Sqrt absolute standardized residual")
axes[1, 1].stem(np.arange(1, len(mic) + 1), influence.cooks_distance[0], linefmt="#2F5597", markerfmt="o", basefmt="black")
axes[1, 1].set_xlabel("Record")
axes[1, 1].set_ylabel("Cook's distance")
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_2_regression_diagnostics.png", bbox_inches="tight")
plt.close(fig)

comparison = pd.DataFrame(results)
comparison = comparison.loc[comparison["analysis"].isin(["primary_all_source_defined_mics", "secondary_molar_outcome"])]
comparison["predictor_label"] = comparison["predictor"].map({"verified_length_aa": "Total length", "idr_length_aa": "IDR length"})
fig, ax = plt.subplots(figsize=(6.0, 4.2))
sns.barplot(data=comparison, x="predictor_label", y="slope", hue="analysis", ax=ax, palette=["#2F5597", "#A65E2E"])
ax.axhline(0, color="black", linewidth=1)
ax.set_xlabel("Predictor")
ax.set_ylabel("OLS slope")
handles, labels = ax.get_legend_handles_labels()
ax.legend(handles, ["Mass MIC", "Molar MIC"], title="Outcome scale", fontsize=8, title_fontsize=8)
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_3_mass_vs_molar_slopes.png", bbox_inches="tight")
plt.close(fig)

versions = {
    "python": sys.version,
    "platform": platform.platform(),
    "numpy": np.__version__,
    "pandas": pd.__version__,
    "scipy": stats.__version__ if hasattr(stats, "__version__") else "1.18.0",
    "statsmodels": sm.__version__,
    "matplotlib": matplotlib.__version__,
    "seaborn": sns.__version__
}
(table_dir / "software_versions.json").write_text(json.dumps(versions, indent=2), encoding="utf-8")
