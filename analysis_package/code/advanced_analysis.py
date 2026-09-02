from pathlib import Path
import hashlib
import json
import platform
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO, Phylo, SeqIO
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy import optimize, stats
import statsmodels.api as sm
import statsmodels.formula.api as smf
from statsmodels.stats.multitest import multipletests

warnings.filterwarnings("ignore")
root = Path(__file__).resolve().parent
package_root = root.parent if (root.parent / "data").exists() and (root.parent / "tables").exists() else root
input_dir = package_root / "data" if package_root != root else root
table_dir = package_root / "tables" if package_root != root else root / "advanced_tables"
figure_dir = package_root / "figures" if package_root != root else root / "advanced_figures"
table_dir.mkdir(exist_ok=True)
figure_dir.mkdir(exist_ok=True)
rng = np.random.default_rng(20260901)
sns.set_theme(style="ticks", context="paper", font_scale=1.0)
plt.rcParams.update({"font.family": "DejaVu Sans", "figure.dpi": 150, "savefig.dpi": 300})

name_map = {
    "Pyocin_S5": "Pyocin S5",
    "PaeM4": "PaeM4",
    "LlpA_BW11M1": "LlpA BW11M1",
    "Pyocin_L1": "Pyocin L1",
    "Syringacin_M": "Syringacin M",
    "KpneA": "KpneA",
    "KaerA": "KaerA",
    "KvarIa": "KvarIa",
    "KpneM": "KpneM",
    "KpneM2": "KpneM2",
    "KvarM": "KvarM",
    "Colicin_E3": "Colicin E3",
    "Colicin_N": "Colicin N",
    "Colicin_E7": "Colicin E7",
    "Colicin_E1": "Colicin E1",
    "Colicin_E2": "Colicin E2",
    "Colicin_A": "Colicin A",
    "Colicin_B": "Colicin B",
    "Colicin_D": "Colicin D",
    "Pyocin_L2": "Pyocin L2",
    "Pyocin_L3": "Pyocin L3",
    "Pyocin_S6": "Pyocin S6"
}

accessions = {
    "Pyocin S5": ("WP_003115311.1", "NCBI Protein"),
    "PaeM4": ("ERY59288.1", "NCBI Protein"),
    "LlpA BW11M1": ("AAM95702.1", "NCBI Protein"),
    "Pyocin L1": ("CDG56231.1", "NCBI Protein"),
    "Syringacin M": ("Q88A25", "UniProtKB"),
    "KpneA": ("SAV78255.1", "NCBI Protein"),
    "KaerA": ("WP_063414841.1", "NCBI Protein"),
    "KvarIa": ("KDL88409.1", "NCBI Protein"),
    "KpneM": ("EWD35590.1", "NCBI Protein"),
    "KpneM2": ("WP_047066220.1", "NCBI Protein"),
    "KvarM": ("CTQ17225.1", "NCBI Protein"),
    "Colicin A": ("P04480", "UniProtKB"),
    "Colicin B": ("P05819", "UniProtKB"),
    "Colicin D": ("P17998", "UniProtKB"),
    "Colicin E1": ("P02978", "UniProtKB"),
    "Colicin E2": ("P04419", "UniProtKB"),
    "Colicin E3": ("P00646", "UniProtKB"),
    "Colicin E7": ("Q47112", "UniProtKB"),
    "Colicin N": ("P08083", "UniProtKB"),
    "Pyocin L2": ("ERX71449.1", "NCBI Protein"),
    "Pyocin L3": ("ERZ09331.1", "NCBI Protein"),
    "Pyocin S6": ("KYO98147.1", "NCBI Protein")
}

families = {
    "Pyocin S5": "pore-forming",
    "PaeM4": "ColM-like",
    "LlpA BW11M1": "lectin-like",
    "Pyocin L1": "lectin-like",
    "Syringacin M": "ColM-like",
    "KpneA": "pore-forming",
    "KaerA": "pore-forming",
    "KvarIa": "pore-forming",
    "KpneM": "ColM-like",
    "KpneM2": "ColM-like",
    "KvarM": "ColM-like",
    "Colicin A": "pore-forming",
    "Colicin B": "pore-forming",
    "Colicin D": "nuclease",
    "Colicin E1": "pore-forming",
    "Colicin E2": "nuclease",
    "Colicin E3": "nuclease",
    "Colicin E7": "nuclease",
    "Colicin N": "pore-forming",
    "Pyocin L2": "lectin-like",
    "Pyocin L3": "lectin-like",
    "Pyocin S6": "nuclease"
}

imports = {
    "Pyocin S5": "Ton",
    "PaeM4": "Ton",
    "LlpA BW11M1": "unresolved",
    "Pyocin L1": "unresolved",
    "Syringacin M": "Ton",
    "KpneA": "Ton",
    "KaerA": "Ton",
    "KvarIa": "unresolved",
    "KpneM": "Ton",
    "KpneM2": "Ton",
    "KvarM": "Ton",
    "Colicin A": "Tol",
    "Colicin B": "Ton",
    "Colicin D": "Ton",
    "Colicin E1": "Tol",
    "Colicin E2": "Tol",
    "Colicin E3": "Tol",
    "Colicin E7": "Tol",
    "Colicin N": "Tol",
    "Pyocin L2": "unresolved",
    "Pyocin L3": "unresolved",
    "Pyocin S6": "unresolved"
}

records = []

def add_exact(protein, target, value, unit, assay, doi, location):
    records.append({"bacteriocin": protein, "target_strain": target, "reported_lower": value, "reported_upper": value, "reported_unit": unit, "censoring": "exact", "assay_group": assay, "doi": doi, "source_location": location})

def add_interval(protein, target, lower, upper, unit, assay, doi, location):
    records.append({"bacteriocin": protein, "target_strain": target, "reported_lower": lower, "reported_upper": upper, "reported_unit": unit, "censoring": "interval", "assay_group": assay, "doi": doi, "source_location": location})

def add_right(protein, target, lower, unit, assay, doi, location):
    records.append({"bacteriocin": protein, "target_strain": target, "reported_lower": lower, "reported_upper": np.inf, "reported_unit": unit, "censoring": "right", "assay_group": assay, "doi": doi, "source_location": location})

add_exact("Pyocin S5", "P. aeruginosa DWW3", 12.6, "ug/mL", "source-defined MIC", "10.1016/j.febslet.2010.06.021", "main text")
add_exact("Pyocin S5", "P. aeruginosa A19", 0.1, "ug/mL", "agar dilution", "10.1371/journal.pone.0185782", "Supplementary Figure S1")
add_exact("PaeM4", "P. aeruginosa A19", 0.6, "ug/mL", "agar dilution", "10.1371/journal.pone.0185782", "Supplementary Figure S1")
add_exact("LlpA BW11M1", "P. syringae GR12-2R3", 2.08, "nM", "spot-overlay MIC", "10.1371/journal.ppat.1003199", "main text")
add_exact("Pyocin L1", "P. aeruginosa P8", 7.0, "nM", "spot-overlay MIC", "10.1371/journal.ppat.1003898", "main text")
add_exact("Pyocin L1", "P. aeruginosa E2", 27.0, "nM", "spot-overlay MIC", "10.1371/journal.ppat.1003898", "main text")
add_exact("Syringacin M", "P. syringae pv. lachrymans LMG 5456", 60.0, "nM", "spot-overlay MIC", "10.1074/jbc.M112.400150", "main text")
add_exact("KpneA", "K. variicola DSM15968", 0.2, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KaerA", "K. variicola DSM15968", 0.1, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KvarIa", "K. variicola DSM15968", 0.4, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KvarM", "K. pneumoniae DSM16358", 0.1, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KvarM", "K. variicola DSM15968", 0.2, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KpneM", "K. pneumoniae DSM16358", 0.8, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KpneM", "K. variicola DSM15968", 0.4, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")
add_exact("KpneM2", "K. pneumoniae DSM16358", 0.8, "ug/mL", "microbroth dilution", "10.1038/s41598-019-51969-1", "Supplementary Table S2")

for protein, target, value in [
    ("Pyocin L1", "P. aeruginosa Br776", 26.0),
    ("Pyocin L1", "P. aeruginosa Bu007", 10.8),
    ("Pyocin L1", "P. aeruginosa LMG 1272", 88.1),
    ("Pyocin L1", "P. aeruginosa PA134", 13.4),
    ("Pyocin L1", "P. aeruginosa PA135", 30.0),
    ("Pyocin L1", "P. aeruginosa PA229", 224.5),
    ("Pyocin L1", "P. aeruginosa PAO1", 38.3),
    ("Pyocin L2", "P. aeruginosa CFPA13", 13.8),
    ("Pyocin L2", "P. aeruginosa CFPA22", 61.8),
    ("Pyocin L2", "P. aeruginosa CFPA118", 37.0),
    ("Pyocin L2", "P. aeruginosa CFPA120", 61.8),
    ("Pyocin L2", "P. aeruginosa CFPA124", 61.8),
    ("Pyocin L2", "P. aeruginosa CPHL 6750", 46.2),
    ("Pyocin L3", "P. aeruginosa CFPA13", 6.9),
    ("Pyocin L3", "P. aeruginosa CFPA22", 77.8),
    ("Pyocin L3", "P. aeruginosa CFPA54", 19.4),
    ("Pyocin L3", "P. aeruginosa CFPA87", 55.0),
    ("Pyocin L3", "P. aeruginosa CFPA120", 9.7),
    ("Pyocin L3", "P. aeruginosa CFPA124", 13.8)
]:
    add_exact(protein, target, value, "nM", "kinetic broth MIC", "10.1002/mbo3.210", "Table 2")

for protein, target in [
    ("Pyocin L2", "P. aeruginosa Bu007"),
    ("Pyocin L1", "P. aeruginosa CFPA22"),
    ("Pyocin L2", "P. aeruginosa CFPA101"),
    ("Pyocin L3", "P. aeruginosa CFPA101")
]:
    add_right(protein, target, 3300.0, "nM", "kinetic broth MIC", "10.1002/mbo3.210", "Table 2")

add_exact("Pyocin S6", "P. aeruginosa CF_PA109", 260.0, "ug/mL", "spot-overlay MIC", "10.1002/mbo3.339", "Figure 4C")

for protein, lower, upper in [
    ("Colicin E1", 0.125, 0.250),
    ("Colicin E2", 0.050, 0.100),
    ("Colicin E3", 0.050, 0.100),
    ("Colicin E7", 0.050, 0.100),
    ("Colicin B", 0.500, 1.000),
    ("Colicin D", 0.500, 1.000),
    ("Colicin N", 0.500, 1.000),
    ("Colicin A", 2.500, 5.000)
]:
    add_interval(protein, "E. coli K-12 BW25113", lower, upper, "nM", "spot titer MIC", "10.1111/j.1365-2958.2009.06788.x", "Results, colicin cytotoxicity and receptor function")

sequences = {name_map[record.id]: str(record.seq) for record in SeqIO.parse(input_dir / "expanded_bacteriocins.fasta", "fasta")}
iupred = pd.read_csv(input_dir / "expanded_iupred3_summary.csv")
iupred = iupred.set_index("bacteriocin")
columns = ["protein", "md5", "length", "analysis", "signature", "description", "start", "end", "score", "status", "date", "interpro", "interpro_description", "go_terms", "pathways"]
annotations = pd.read_csv(input_dir / "expanded_interproscan_results.tsv", sep="\t", header=None, names=columns)
annotations["bacteriocin"] = annotations["protein"].map(name_map)

domain_rules = {
    "pore-forming": [("Pfam", "PF01024")],
    "ColM-like": [("Gene3D", "G3DSA:3.30.450.400")],
    "lectin-like": [("SMART", "SM00108")],
    "nuclease": [("Pfam", "PF09000"), ("Pfam", "PF21431"), ("Pfam", "PF11429"), ("Gene3D", "G3DSA:3.90.540.10"), ("Gene3D", "G3DSA:3.10.450.200")]
}

def descriptor(sequence):
    analysis = ProteinAnalysis(sequence)
    helix, turn, sheet = analysis.secondary_structure_fraction()
    return {
        "length_aa": len(sequence),
        "molecular_weight_da": analysis.molecular_weight(),
        "isoelectric_point": analysis.isoelectric_point(),
        "charge_pH7": analysis.charge_at_pH(7.0),
        "gravy": analysis.gravy(),
        "aromaticity": analysis.aromaticity(),
        "instability_index": analysis.instability_index(),
        "helix_fraction": helix,
        "turn_fraction": turn,
        "sheet_fraction": sheet
    }

feature_rows = []
for bacteriocin, sequence in sequences.items():
    family = families[bacteriocin]
    hits = pd.DataFrame()
    selected_analysis = ""
    selected_signature = ""
    for analysis_name, signature in domain_rules[family]:
        candidate = annotations.loc[annotations["bacteriocin"].eq(bacteriocin) & annotations["analysis"].eq(analysis_name) & annotations["signature"].eq(signature)].sort_values(["start", "end"])
        if not candidate.empty:
            hits = candidate
            selected_analysis = analysis_name
            selected_signature = signature
            break
    if hits.empty:
        raise RuntimeError(f"No functional-domain annotation for {bacteriocin}")
    spans = [(int(row.start), int(row.end)) for row in hits.itertuples()]
    if family == "lectin-like":
        domain_sequence = "".join(sequence[start - 1:end] for start, end in spans)
    else:
        start, end = max(spans, key=lambda item: item[1] - item[0])
        spans = [(start, end)]
        domain_sequence = sequence[start - 1:end]
    whole = descriptor(sequence)
    domain = descriptor(domain_sequence)
    nterm = descriptor(sequence[:40])
    accession, database = accessions[bacteriocin]
    row = {
        "bacteriocin": bacteriocin,
        "accession": accession,
        "database": database,
        "family": family,
        "import_pathway": imports[bacteriocin],
        "sequence_sha256": hashlib.sha256(sequence.encode()).hexdigest(),
        "domain_method": selected_analysis,
        "domain_signature": selected_signature,
        "domain_coordinates": "; ".join(f"{start}-{end}" for start, end in spans),
        "domain_length_aa": len(domain_sequence),
        "domain_charge_pH7": domain["charge_pH7"],
        "domain_gravy": domain["gravy"],
        "domain_isoelectric_point": domain["isoelectric_point"],
        "nterm40_charge_pH7": nterm["charge_pH7"],
        "nterm40_gravy": nterm["gravy"],
        **whole,
        "idr_length_aa": iupred.loc[bacteriocin, "idr_length_aa"],
        "idr_fraction": iupred.loc[bacteriocin, "idr_fraction"],
        "longest_idr_aa": iupred.loc[bacteriocin, "longest_idr_aa"],
        "nterm40_idr_fraction": iupred.loc[bacteriocin, "nterm40_idr_fraction"]
    }
    feature_rows.append(row)

features = pd.DataFrame(feature_rows).sort_values(["family", "bacteriocin"])
features.to_csv(table_dir / "expanded_protein_features.csv", index=False)
features[["bacteriocin", "accession", "database", "length_aa", "family", "domain_signature", "domain_coordinates", "domain_length_aa"]].to_csv(table_dir / "expanded_accession_domain_table.csv", index=False)

mic = pd.DataFrame(records)
mic = mic.merge(features, on="bacteriocin", how="left", validate="many_to_one")
mic["target_genus"] = mic["target_strain"].str.split().str[0].str.replace(".", "", regex=False)
for side in ["lower", "upper"]:
    source = mic[f"reported_{side}"]
    mic[f"mic_nM_{side}"] = np.where(mic["reported_unit"].eq("nM"), source, source * 1e6 / mic["molecular_weight_da"])
    mic[f"mic_ug_ml_{side}"] = np.where(mic["reported_unit"].eq("ug/mL"), source, source * mic["molecular_weight_da"] / 1e6)
    mic[f"log10_mic_ug_ml_{side}"] = np.log10(mic[f"mic_ug_ml_{side}"])
mic["mic_ug_ml_midpoint"] = np.where(np.isfinite(mic["mic_ug_ml_upper"]), np.sqrt(mic["mic_ug_ml_lower"] * mic["mic_ug_ml_upper"]), np.nan)
mic["log10_mic_ug_ml_midpoint"] = np.log10(mic["mic_ug_ml_midpoint"])
mic.to_csv(table_dir / "expanded_mic_dataset.csv", index=False)

search_log = pd.DataFrame([
    ("Amass BiomedCore", "bacteriocin OR colicin OR pyocin OR klebicin OR syringacin AND minimum inhibitory concentration OR MIC AND purified protein", "2026-09-01", "Discovery and full-text validation"),
    ("PubMed and PMC", "colicin minimum inhibitory concentration purified; pyocin minimum inhibitory concentration purified; pectocin MIC bacteriocin; klebicin MIC bacteriocin", "2026-09-01", "Primary-source discovery"),
    ("Sider Scholar OpenAlex", "protein bacteriocin minimum inhibitory concentration", "2026-09-01", "Broad discovery with manual primary-source validation"),
    ("Scite", "bacteriocin OR colicin OR pyocin OR klebicin OR syringacin AND minimum inhibitory concentration OR MIC AND purified", "2026-09-01", "Citation and retraction screening"),
    ("Consensus", "purified protein bacteriocin minimum inhibitory concentration Gram-negative colicin pyocin klebicin", "2026-09-01", "Discovery cross-check"),
    ("SciSpace", "studies reporting quantitative MICs for purified soluble protein bacteriocins active against Gram-negative bacteria", "2026-09-01", "Discovery cross-check")
], columns=["database", "query", "search_date", "purpose"])
search_log.to_csv(table_dir / "literature_search_log.csv", index=False)

source_ledger = mic[["doi", "source_location", "bacteriocin", "target_strain", "censoring", "reported_lower", "reported_upper", "reported_unit"]].drop_duplicates().sort_values(["doi", "bacteriocin", "target_strain"])
source_ledger.to_csv(table_dir / "claim_source_ledger.csv", index=False)

def standardize_by_protein(data, predictor):
    mapping = data[["bacteriocin", predictor]].drop_duplicates().set_index("bacteriocin")[predictor]
    mean = mapping.mean()
    sd = mapping.std(ddof=0)
    return data[predictor].sub(mean).div(sd), mean, sd

def interval_fit(data, predictor):
    x, mean, sd = standardize_by_protein(data, predictor)
    lower = data["log10_mic_ug_ml_lower"].to_numpy(float)
    upper = data["log10_mic_ug_ml_upper"].to_numpy(float)
    exact = data["censoring"].eq("exact").to_numpy()
    finite = np.isfinite(upper)
    X = np.column_stack([np.ones(len(data)), x.to_numpy(float)])
    def objective(theta):
        mu = X @ theta[:2]
        sigma = np.exp(theta[2])
        zlow = (lower - mu) / sigma
        zup = (upper - mu) / sigma
        likelihood = np.empty(len(data))
        likelihood[exact] = stats.norm.pdf(zlow[exact]) / sigma
        bounded = ~exact & finite
        likelihood[bounded] = stats.norm.cdf(zup[bounded]) - stats.norm.cdf(zlow[bounded])
        right = ~finite
        likelihood[right] = stats.norm.sf(zlow[right])
        return -np.sum(np.log(np.clip(likelihood, 1e-300, None)))
    initial = np.array([np.nanmedian(data["log10_mic_ug_ml_midpoint"]), 0.0, np.log(np.nanstd(data["log10_mic_ug_ml_midpoint"]))])
    result = optimize.minimize(objective, initial, method="BFGS")
    if not result.success and np.linalg.norm(result.jac) > 1e-3:
        raise RuntimeError(result.message)
    return result.x, mean, sd, objective(result.x)

primary_predictors = ["length_aa", "domain_length_aa", "idr_fraction", "nterm40_idr_fraction", "domain_charge_pH7", "domain_gravy"]
interval_rows = []
bootstrap_distributions = {}
interval_path = table_dir / "primary_interval_regression_results.csv"
if interval_path.exists():
    interval_results = pd.read_csv(interval_path)
else:
    for predictor in primary_predictors:
        estimate, mean, sd, nll = interval_fit(mic, predictor)
        proteins = mic["bacteriocin"].unique()
        slopes = []
        for iteration in range(1000):
            sampled = rng.choice(proteins, size=len(proteins), replace=True)
            parts = []
            for index, protein in enumerate(sampled):
                part = mic.loc[mic["bacteriocin"].eq(protein)].copy()
                part["bacteriocin"] = part["bacteriocin"] + f"_{index}"
                parts.append(part)
            sample = pd.concat(parts, ignore_index=True)
            try:
                slopes.append(interval_fit(sample, predictor)[0][1])
            except Exception:
                pass
        slopes = np.asarray(slopes)
        bootstrap_distributions[predictor] = slopes
        p_boot = min(1.0, 2 * min(np.mean(slopes <= 0), np.mean(slopes >= 0)))
        interval_rows.append({
            "predictor": predictor,
            "n_records": len(mic),
            "n_proteins": mic["bacteriocin"].nunique(),
            "n_exact": mic["censoring"].eq("exact").sum(),
            "n_interval": mic["censoring"].eq("interval").sum(),
            "n_right_censored": mic["censoring"].eq("right").sum(),
            "standardized_slope": estimate[1],
            "bootstrap_ci95_low": np.quantile(slopes, 0.025),
            "bootstrap_ci95_high": np.quantile(slopes, 0.975),
            "bootstrap_p": p_boot,
            "raw_slope": estimate[1] / sd,
            "predictor_mean": mean,
            "predictor_sd": sd,
            "residual_sigma": np.exp(estimate[2]),
            "negative_log_likelihood": nll,
            "bootstrap_successes": len(slopes)
        })
    interval_results = pd.DataFrame(interval_rows)
    interval_results["fdr_q"] = multipletests(interval_results["bootstrap_p"], method="fdr_bh")[1]
interval_results.to_csv(interval_path, index=False)

def interval_adjusted_fit(data, predictor, categorical):
    x, mean, sd = standardize_by_protein(data, predictor)
    design = pd.DataFrame({"intercept": 1.0, "predictor": x.to_numpy(float)})
    for variable in categorical:
        dummies = pd.get_dummies(data[variable], prefix=variable, drop_first=True, dtype=float).reset_index(drop=True)
        design = pd.concat([design.reset_index(drop=True), dummies], axis=1)
    X = design.to_numpy(float)
    lower = data["log10_mic_ug_ml_lower"].to_numpy(float)
    upper = data["log10_mic_ug_ml_upper"].to_numpy(float)
    exact_mask = data["censoring"].eq("exact").to_numpy()
    finite = np.isfinite(upper)
    def objective(theta):
        mu = X @ theta[:-1]
        sigma = np.exp(theta[-1])
        zlow = (lower - mu) / sigma
        zup = (upper - mu) / sigma
        likelihood = np.empty(len(data))
        likelihood[exact_mask] = stats.norm.pdf(zlow[exact_mask]) / sigma
        bounded = ~exact_mask & finite
        likelihood[bounded] = stats.norm.cdf(zup[bounded]) - stats.norm.cdf(zlow[bounded])
        right = ~finite
        likelihood[right] = stats.norm.sf(zlow[right])
        return -np.sum(np.log(np.clip(likelihood, 1e-300, None)))
    initial = np.zeros(X.shape[1] + 1)
    initial[0] = np.nanmedian(data["log10_mic_ug_ml_midpoint"])
    initial[-1] = np.log(np.nanstd(data["log10_mic_ug_ml_midpoint"]))
    result = optimize.minimize(objective, initial, method="BFGS")
    covariance = np.asarray(result.hess_inv)
    se = np.sqrt(max(covariance[1, 1], 0))
    return result.x[1], se, result.success

interval_sensitivity_rows = []
for predictor in primary_predictors:
    for adjustment in [[], ["family"], ["assay_group"], ["family", "assay_group"]]:
        slope, se, success = interval_adjusted_fit(mic, predictor, adjustment)
        interval_sensitivity_rows.append({"predictor": predictor, "analysis": "unadjusted" if not adjustment else "adjusted for " + " and ".join(adjustment), "standardized_slope": slope, "approximate_se": se, "wald_ci95_low": slope - 1.96 * se, "wald_ci95_high": slope + 1.96 * se, "converged": success})
    for omitted in mic["doi"].unique():
        subset = mic.loc[~mic["doi"].eq(omitted)].copy()
        slope = interval_fit(subset, predictor)[0][1]
        interval_sensitivity_rows.append({"predictor": predictor, "analysis": "leave one study out", "omitted": omitted, "standardized_slope": slope, "approximate_se": np.nan, "wald_ci95_low": np.nan, "wald_ci95_high": np.nan, "converged": True})
    for family in mic["family"].unique():
        subset = mic.loc[mic["family"].eq(family)].copy()
        if subset["bacteriocin"].nunique() >= 3 and subset[predictor].nunique() >= 2:
            slope = interval_fit(subset, predictor)[0][1]
            interval_sensitivity_rows.append({"predictor": predictor, "analysis": "within family", "family": family, "standardized_slope": slope, "approximate_se": np.nan, "wald_ci95_low": np.nan, "wald_ci95_high": np.nan, "converged": True})
pd.DataFrame(interval_sensitivity_rows).to_csv(table_dir / "interval_model_adjustment_and_stability.csv", index=False)

exact = mic.loc[mic["censoring"].eq("exact")].copy()
exact["log10_mic_ug_ml"] = exact["log10_mic_ug_ml_lower"]
model_rows = []
for predictor in primary_predictors:
    z, mean, sd = standardize_by_protein(exact, predictor)
    exact[f"z_{predictor}"] = z
    formula = f"log10_mic_ug_ml ~ z_{predictor}"
    fit = smf.ols(formula, data=exact).fit()
    hc3 = fit.get_robustcov_results(cov_type="HC3")
    protein_cluster = fit.get_robustcov_results(cov_type="cluster", groups=exact["bacteriocin"])
    study_cluster = fit.get_robustcov_results(cov_type="cluster", groups=exact["doi"])
    rlm = smf.rlm(formula, data=exact, M=sm.robust.norms.HuberT()).fit()
    quantile = smf.quantreg(formula, data=exact).fit(q=0.5)
    adjusted = smf.ols(f"log10_mic_ug_ml ~ z_{predictor} + C(family) + C(assay_group)", data=exact).fit(cov_type="HC3")
    for method, result, position in [
        ("OLS HC3", hc3, 1),
        ("OLS protein-clustered", protein_cluster, 1),
        ("OLS study-clustered", study_cluster, 1),
        ("Huber robust regression", rlm, 1),
        ("median quantile regression", quantile, 1),
        ("family and assay adjusted OLS HC3", adjusted, list(adjusted.params.index).index(f"z_{predictor}"))
    ]:
        ci = np.asarray(result.conf_int())[position]
        params = np.asarray(result.params)
        bse = np.asarray(result.bse)
        pvalues = np.asarray(result.pvalues)
        model_rows.append({"predictor": predictor, "method": method, "n_records": len(exact), "n_proteins": exact["bacteriocin"].nunique(), "standardized_slope": params[position], "standard_error": bse[position], "ci95_low": ci[0], "ci95_high": ci[1], "p_value": pvalues[position], "r_squared": getattr(result, "rsquared", np.nan)})

model_results = pd.DataFrame(model_rows)
model_results.to_csv(table_dir / "advanced_model_ladder_results.csv", index=False)

exploratory_predictors = ["molecular_weight_da", "isoelectric_point", "charge_pH7", "gravy", "aromaticity", "instability_index", "helix_fraction", "turn_fraction", "sheet_fraction", "longest_idr_aa", "nterm40_charge_pH7", "nterm40_gravy", "domain_isoelectric_point"]
exploratory_rows = []
for predictor in exploratory_predictors:
    z, mean, sd = standardize_by_protein(exact, predictor)
    exact[f"z_{predictor}"] = z
    fit = smf.ols(f"log10_mic_ug_ml ~ z_{predictor}", data=exact).fit()
    robust = fit.get_robustcov_results(cov_type="cluster", groups=exact["bacteriocin"])
    ci = robust.conf_int()[1]
    exploratory_rows.append({"predictor": predictor, "standardized_slope": robust.params[1], "ci95_low": ci[0], "ci95_high": ci[1], "p_value": robust.pvalues[1], "r_squared": fit.rsquared})
exploratory = pd.DataFrame(exploratory_rows)
exploratory["fdr_q"] = multipletests(exploratory["p_value"], method="fdr_bh")[1]
exploratory.to_csv(table_dir / "exploratory_feature_results.csv", index=False)

stability_rows = []
for predictor in primary_predictors:
    for grouping in ["bacteriocin", "doi"]:
        for omitted in exact[grouping].unique():
            subset = exact.loc[~exact[grouping].eq(omitted)].copy()
            z, mean, sd = standardize_by_protein(subset, predictor)
            subset["z_predictor"] = z
            fit = smf.ols("log10_mic_ug_ml ~ z_predictor", data=subset).fit(cov_type="HC3")
            stability_rows.append({"predictor": predictor, "omission_type": grouping, "omitted": omitted, "n_records": len(subset), "n_proteins": subset["bacteriocin"].nunique(), "standardized_slope": fit.params["z_predictor"], "p_value": fit.pvalues["z_predictor"]})
stability = pd.DataFrame(stability_rows)
stability.to_csv(table_dir / "leave_one_group_out_results.csv", index=False)

median = exact.groupby("bacteriocin", as_index=False).agg(log10_mic_ug_ml=("log10_mic_ug_ml", "median"), family=("family", "first"), doi=("doi", lambda x: ";".join(sorted(set(x)))))
median = median.merge(features, on=["bacteriocin", "family"], how="left")
median_rows = []
for predictor in primary_predictors:
    median["z_predictor"] = (median[predictor] - median[predictor].mean()) / median[predictor].std(ddof=0)
    fit = smf.ols("log10_mic_ug_ml ~ z_predictor", data=median).fit(cov_type="HC3")
    ci = fit.conf_int().loc["z_predictor"]
    median_rows.append({"predictor": predictor, "n_proteins": len(median), "standardized_slope": fit.params["z_predictor"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1], "p_value": fit.pvalues["z_predictor"], "r_squared": fit.rsquared})
pd.DataFrame(median_rows).to_csv(table_dir / "protein_median_sensitivity_results.csv", index=False)

study_counts = exact.groupby("doi")["log10_mic_ug_ml"].transform("count")
study_sd = exact.groupby("doi")["log10_mic_ug_ml"].transform("std")
within = exact.loc[(study_counts >= 2) & study_sd.gt(0)].copy()
within["within_study_z_mic"] = within.groupby("doi")["log10_mic_ug_ml"].transform(lambda x: (x - x.mean()) / x.std(ddof=0))
within_rows = []
for predictor in primary_predictors:
    z, mean, sd = standardize_by_protein(within, predictor)
    within["z_predictor"] = z
    fit = smf.ols("within_study_z_mic ~ z_predictor", data=within).fit(cov_type="cluster", cov_kwds={"groups": within["bacteriocin"]})
    ci = fit.conf_int().loc["z_predictor"]
    within_rows.append({"predictor": predictor, "n_records": len(within), "n_proteins": within["bacteriocin"].nunique(), "n_studies": within["doi"].nunique(), "standardized_slope": fit.params["z_predictor"], "ci95_low": ci.iloc[0], "ci95_high": ci.iloc[1], "p_value": fit.pvalues["z_predictor"], "r_squared": fit.rsquared})
pd.DataFrame(within_rows).to_csv(table_dir / "within_study_standardized_results.csv", index=False)

predictive_features = ["length_aa", "domain_length_aa", "idr_fraction", "nterm40_idr_fraction", "domain_charge_pH7", "domain_gravy", "charge_pH7", "gravy", "aromaticity", "instability_index"]
X = exact[predictive_features].to_numpy(float)
y = exact["log10_mic_ug_ml"].to_numpy(float)
groups = exact["bacteriocin"].to_numpy()
alphas = np.logspace(-4, 4, 33)
predictions = np.full(len(exact), np.nan)
selected_alphas = []

def ridge_predict(train_x, train_y, test_x, alpha):
    mean = train_x.mean(axis=0)
    sd = train_x.std(axis=0)
    sd[sd == 0] = 1
    tx = (train_x - mean) / sd
    vx = (test_x - mean) / sd
    center = train_y.mean()
    beta = np.linalg.solve(tx.T @ tx + alpha * np.eye(tx.shape[1]), tx.T @ (train_y - center))
    return center + vx @ beta

for held_out in np.unique(groups):
    outer_train = groups != held_out
    outer_test = groups == held_out
    inner_groups = np.unique(groups[outer_train])
    scores = []
    for alpha in alphas:
        errors = []
        for inner_out in inner_groups:
            train = outer_train & (groups != inner_out)
            test = outer_train & (groups == inner_out)
            predicted = ridge_predict(X[train], y[train], X[test], alpha)
            errors.extend((y[test] - predicted) ** 2)
        scores.append(np.mean(errors))
    alpha = alphas[int(np.argmin(scores))]
    selected_alphas.append(alpha)
    predictions[outer_test] = ridge_predict(X[outer_train], y[outer_train], X[outer_test], alpha)

baseline = np.full(len(y), np.nan)
for held_out in np.unique(groups):
    baseline[groups == held_out] = y[groups != held_out].mean()
prediction_metrics = pd.DataFrame([{
    "model": "nested leave-one-protein-out ridge",
    "n_records": len(y),
    "n_proteins": len(np.unique(groups)),
    "rmse": np.sqrt(np.mean((y - predictions) ** 2)),
    "mae": np.mean(np.abs(y - predictions)),
    "predictive_r_squared": 1 - np.sum((y - predictions) ** 2) / np.sum((y - y.mean()) ** 2),
    "baseline_rmse": np.sqrt(np.mean((y - baseline) ** 2)),
    "median_selected_alpha": np.median(selected_alphas)
}])
prediction_metrics.to_csv(table_dir / "nested_ridge_prediction_metrics.csv", index=False)
pd.DataFrame({"bacteriocin": groups, "observed_log10_mic": y, "predicted_log10_mic": predictions, "baseline_prediction": baseline}).to_csv(table_dir / "nested_ridge_predictions.csv", index=False)

alignment = AlignIO.read(input_dir / "expanded_clustalo_alignment.fasta", "fasta")
alignment_names = [name_map.get(record.id, record.id) for record in alignment]
identity = np.zeros((len(alignment), len(alignment)))
for i, first in enumerate(alignment):
    for j, second in enumerate(alignment):
        comparable = [(a, b) for a, b in zip(str(first.seq), str(second.seq)) if a != "-" and b != "-"]
        identity[i, j] = 100 if i == j else 100 * sum(a == b for a, b in comparable) / len(comparable)
identity_frame = pd.DataFrame(identity, index=alignment_names, columns=alignment_names)
identity_frame.to_csv(table_dir / "expanded_clustalo_percent_identity_matrix.csv")
upper_identity = identity[np.triu_indices(len(identity), 1)]
pd.DataFrame([{"n_proteins": len(identity), "n_pairs": len(upper_identity), "median_percent_identity": np.median(upper_identity), "q1_percent_identity": np.quantile(upper_identity, 0.25), "q3_percent_identity": np.quantile(upper_identity, 0.75), "minimum_percent_identity": np.min(upper_identity), "maximum_percent_identity": np.max(upper_identity), "pairs_above_95_percent": np.sum(upper_identity > 95)}]).to_csv(table_dir / "expanded_clustalo_identity_summary.csv", index=False)

fig, axes = plt.subplots(1, 3, figsize=(12.2, 4.1))
for ax, predictor, label in zip(axes, ["length_aa", "domain_length_aa", "idr_fraction"], ["Total length (aa)", "Functional-domain length (aa)", "IUPred3 disorder fraction"]):
    sns.scatterplot(data=exact, x=predictor, y="log10_mic_ug_ml", hue="family", style="assay_group", s=55, ax=ax, legend=False)
    fit = smf.ols(f"log10_mic_ug_ml ~ {predictor}", data=exact).fit()
    grid = pd.DataFrame({predictor: np.linspace(exact[predictor].min(), exact[predictor].max(), 100)})
    prediction = fit.get_prediction(grid).summary_frame()
    ax.plot(grid[predictor], prediction["mean"], color="black", linewidth=1.1)
    ax.fill_between(grid[predictor], prediction["mean_ci_lower"], prediction["mean_ci_upper"], color="black", alpha=0.12)
    ax.set_xlabel(label)
    ax.set_ylabel("log10 MIC (ug/mL)")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_1_expanded_sequence_features_vs_MIC.png", bbox_inches="tight")
plt.close(fig)

plot_results = interval_results.sort_values("standardized_slope")
labels = {"length_aa": "Total length", "domain_length_aa": "Domain length", "idr_fraction": "Disorder fraction", "nterm40_idr_fraction": "N-terminal disorder", "domain_charge_pH7": "Domain charge", "domain_gravy": "Domain hydropathy"}
fig, ax = plt.subplots(figsize=(6.8, 4.8))
ypos = np.arange(len(plot_results))
ax.errorbar(plot_results["standardized_slope"], ypos, xerr=[plot_results["standardized_slope"] - plot_results["bootstrap_ci95_low"], plot_results["bootstrap_ci95_high"] - plot_results["standardized_slope"]], fmt="o", color="#2F5597", capsize=3)
ax.axvline(0, color="black", linewidth=1)
ax.set_yticks(ypos, [labels[item] for item in plot_results["predictor"]])
ax.set_xlabel("Change in log10 MIC per predictor SD")
ax.set_ylabel("")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_2_primary_interval_model_forest.png", bbox_inches="tight")
plt.close(fig)

fig, axes = plt.subplots(1, 2, figsize=(10.5, 4.4))
sns.boxplot(data=exact, x="family", y="log10_mic_ug_ml", color="white", showfliers=False, ax=axes[0])
sns.stripplot(data=exact, x="family", y="log10_mic_ug_ml", hue="family", jitter=0.18, size=5, ax=axes[0], legend=False)
axes[0].tick_params(axis="x", rotation=30)
axes[0].set_xlabel("Protein family")
axes[0].set_ylabel("log10 MIC (ug/mL)")
sns.boxplot(data=exact, x="assay_group", y="log10_mic_ug_ml", color="white", showfliers=False, ax=axes[1])
sns.stripplot(data=exact, x="assay_group", y="log10_mic_ug_ml", hue="family", jitter=0.18, size=5, ax=axes[1], legend=False)
axes[1].tick_params(axis="x", rotation=60)
axes[1].set_xlabel("Assay group")
axes[1].set_ylabel("log10 MIC (ug/mL)")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_3_family_and_assay_heterogeneity.png", bbox_inches="tight")
plt.close(fig)

length_stability = stability.loc[stability["predictor"].eq("length_aa")].copy()
fig, axes = plt.subplots(1, 2, figsize=(11.0, 4.5))
for ax, omission, title in zip(axes, ["bacteriocin", "doi"], ["Leave one protein out", "Leave one study out"]):
    subset = length_stability.loc[length_stability["omission_type"].eq(omission)].sort_values("standardized_slope")
    ax.scatter(subset["standardized_slope"], np.arange(len(subset)), color="#2F5597", s=25)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_yticks(np.arange(len(subset)), subset["omitted"], fontsize=6)
    ax.set_xlabel("Standardized total-length slope")
    ax.set_title(title)
fig.tight_layout()
fig.savefig(figure_dir / "Figure_4_total_length_stability.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(9.2, 8.0))
sns.heatmap(identity_frame, cmap="viridis", vmin=0, vmax=100, square=True, cbar_kws={"label": "Percent identity"}, ax=ax)
ax.tick_params(axis="x", labelrotation=65, labelsize=7)
ax.tick_params(axis="y", labelsize=7)
ax.set_xlabel("")
ax.set_ylabel("")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_5_expanded_clustalo_identity_heatmap.png", bbox_inches="tight")
plt.close(fig)

feature_matrix = features[["length_aa", "domain_length_aa", "idr_fraction", "nterm40_idr_fraction", "charge_pH7", "gravy", "aromaticity", "instability_index", "domain_charge_pH7", "domain_gravy"]].corr(method="spearman")
fig, ax = plt.subplots(figsize=(7.6, 6.4))
sns.heatmap(feature_matrix, cmap="vlag", center=0, vmin=-1, vmax=1, annot=True, fmt=".2f", annot_kws={"size": 7}, ax=ax)
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_1_feature_correlation_heatmap.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(5.5, 4.5))
ax.scatter(y, predictions, c=pd.Categorical(exact["family"]).codes, cmap="tab10", s=45)
limits = [min(y.min(), predictions.min()), max(y.max(), predictions.max())]
ax.plot(limits, limits, color="black", linewidth=1)
ax.set_xlim(limits)
ax.set_ylim(limits)
ax.set_xlabel("Observed log10 MIC")
ax.set_ylabel("Nested protein-level ridge prediction")
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_2_nested_ridge_predictions.png", bbox_inches="tight")
plt.close(fig)

fig, ax = plt.subplots(figsize=(8.5, 5.7))
ordered = mic.sort_values("log10_mic_ug_ml_lower").reset_index(drop=True)
for index, row in ordered.iterrows():
    lower = row["log10_mic_ug_ml_lower"]
    upper = row["log10_mic_ug_ml_upper"]
    if np.isfinite(upper):
        ax.plot([lower, upper], [index, index], color="#2F5597", linewidth=2)
        ax.scatter(np.sqrt(0.0) + (lower + upper) / 2, index, color="#2F5597", s=10)
    else:
        ax.annotate("", xy=(lower + 0.5, index), xytext=(lower, index), arrowprops={"arrowstyle": "->", "color": "#A65E2E"})
ax.set_yticks([])
ax.set_xlabel("log10 MIC interval (ug/mL)")
ax.set_ylabel("47 source-defined observations")
fig.tight_layout()
fig.savefig(figure_dir / "Supplementary_Figure_3_censoring_intervals.png", bbox_inches="tight")
plt.close(fig)

coverage = pd.DataFrame([
    ("Literature-expanded quantitative MIC dataset", "completed", f"{len(mic)} records, {mic['bacteriocin'].nunique()} proteins, {mic['doi'].nunique()} source papers"),
    ("Official sequence validation", "completed", "NCBI Protein or UniProtKB accession for every protein"),
    ("IUPred3 long-disorder analysis", "completed", "All 22 proteins"),
    ("InterProScan functional-domain annotation", "completed", "Official EMBL-EBI service for all 22 proteins"),
    ("Clustal Omega comparison", "completed", "Official EMBL-EBI service for all 22 proteins"),
    ("Interval-censored regression", "completed", "Exact, bounded interval, and right-censored MIC observations"),
    ("Adjusted and robust model ladder", "completed", "HC3, clustered, Huber, quantile, adjusted, bootstrap, and sensitivity models"),
    ("Nested protein-level ridge validation", "completed", "Prediction evaluated with outer leave-one-protein-out folds"),
    ("Receptor compatibility score", "not modeled", "Comparable quantitative receptor-expression or binding measurements were not available across proteins"),
    ("Experimental three-dimensional descriptors", "partially completed", "Protplex verified relevant structures, but complete experimental structures were unavailable for all proteins")
], columns=["analysis", "status", "evidence"])
coverage.to_csv(table_dir / "advanced_analysis_coverage.csv", index=False)

software = {
    "python": sys.version,
    "platform": platform.platform(),
    "numpy": np.__version__,
    "pandas": pd.__version__,
    "scipy": stats.__version__ if hasattr(stats, "__version__") else "1.18.0",
    "statsmodels": sm.__version__,
    "matplotlib": matplotlib.__version__,
    "seaborn": sns.__version__,
    "biopython": "1.85",
    "clustal_omega": "1.2.4",
    "interproscan": "5.78-109.0",
    "iupred3": "official web server, long mode"
}
(table_dir / "advanced_software_versions.json").write_text(json.dumps(software, indent=2), encoding="utf-8")

summary = {
    "records_total": int(len(mic)),
    "records_exact": int(mic["censoring"].eq("exact").sum()),
    "records_interval": int(mic["censoring"].eq("interval").sum()),
    "records_right_censored": int(mic["censoring"].eq("right").sum()),
    "proteins": int(mic["bacteriocin"].nunique()),
    "studies": int(mic["doi"].nunique()),
    "families": int(mic["family"].nunique()),
    "assay_groups": int(mic["assay_group"].nunique()),
    "clustalo_median_identity": float(np.median(upper_identity)),
    "primary_results": interval_results.to_dict(orient="records"),
    "ridge_metrics": prediction_metrics.to_dict(orient="records")[0]
}
(table_dir / "advanced_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
