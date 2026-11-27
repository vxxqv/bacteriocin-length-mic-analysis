from pathlib import Path
from itertools import permutations
from collections import Counter
import csv
import hashlib
import io
import json
import shutil
import zipfile
import numpy as np
import pandas as pd
from Bio import SeqIO
from scipy.stats import spearmanr, friedmanchisquare, binomtest, fisher_exact
from statsmodels.stats.multitest import multipletests
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.gridspec import GridSpec

root = Path(__file__).resolve().parent
out = root
raw = (out / 'data' / 'Laosena_2026_Data_S1.xlsx').read_bytes()
wb = pd.ExcelFile(io.BytesIO(raw))
raw_sus = pd.read_excel(wb, "Pyocin susceptibility", header=None)
raw_form = pd.read_excel(wb, "Biofilm formation A550", header=None)
raw_bio = pd.read_excel(wb, "Biofilm eradication A405", header=None)
names = ["AP41", "SX1", "SX2", "S5", "G", "M1"]
acc = dict(zip(names, ["BAA02196.1", "UYS84543.1", "UYS84545.1", "WP_003115311.1", "WP_170835785.1", "ERZ09841.1"]))
seqs = {r.id: str(r.seq) for r in SeqIO.parse(out / "data" / "Laosena_2026_reference_pyocins.fasta", "fasta")}
lengths = {n: len(seqs[a]) for n, a in acc.items()}
assert len(seqs) == 6
assert raw_sus.shape[1] == 22

sus_rows = []
mic_rows = []
for _, row in raw_sus.iloc[3:].iterrows():
    isolate = row.iloc[0]
    if not isinstance(isolate, str) or not isolate.startswith("PATH-"):
        continue
    for j, name in enumerate(names):
        vals = [float(row.iloc[1 + j * 3 + k]) for k in range(3)]
        sus_rows.append({"isolate": isolate, "pyocin": name, "length_aa": lengths[name], "zone_1_mm": vals[0], "zone_2_mm": vals[1], "zone_3_mm": vals[2], "zone_mean_mm": float(np.mean(vals)), "all_positive": int(all(v > 0 for v in vals)), "any_positive": int(any(v > 0 for v in vals))})
    for j, name in enumerate(names[:3]):
        value = row.iloc[19 + j]
        mic_rows.append({"isolate": isolate, "pyocin": name, "mic_raw": str(value), "mic_ug_ml": float(value) if isinstance(value, (int, float, np.number)) and not pd.isna(value) else np.nan})
sus = pd.DataFrame(sus_rows)
mic = pd.DataFrame(mic_rows)
assert sus.isolate.nunique() == 29 and len(sus) == 174
assert len(mic) == 87
sus.to_csv(out / "thai_isolate_pyocin_zones.csv", index=False)
mic.to_csv(out / "thai_isolate_pyocin_mic.csv", index=False)
seq_tab = pd.DataFrame([{"pyocin": n, "accession": acc[n], "length_aa": lengths[n], "sequence_sha256": hashlib.sha256(seqs[acc[n]].encode()).hexdigest()} for n in names])
seq_tab.to_csv(out / "thai_pyocin_accessions.csv", index=False)
coverage = sus.groupby("pyocin", sort=False).agg(length_aa=("length_aa", "first"), all_positive=("all_positive", "sum"), any_positive=("any_positive", "sum"), mean_zone_mm=("zone_mean_mm", "mean")).reset_index()
coverage["pyocin"] = pd.Categorical(coverage.pyocin, names, ordered=True)
coverage = coverage.sort_values("pyocin")
coverage.to_csv(out / "thai_coverage.csv", index=False)
zone_matrix = sus.pivot(index="isolate", columns="pyocin", values="zone_mean_mm")[names]
binary_matrix = sus.pivot(index="isolate", columns="pyocin", values="all_positive")[names]
discordant = int((sus.any_positive != sus.all_positive).sum())
friedman_zone = friedmanchisquare(*[zone_matrix[n] for n in names])
rho_cov = spearmanr(coverage.length_aa, coverage.all_positive).statistic
rho_zone = spearmanr(coverage.length_aa, coverage.mean_zone_mm).statistic
perm_cov = np.array([spearmanr(coverage.length_aa, p).statistic for p in permutations(coverage.all_positive)])
perm_zone = np.array([spearmanr(coverage.length_aa, p).statistic for p in permutations(coverage.mean_zone_mm)])
coverage_perm_p = float(np.mean(abs(perm_cov) >= abs(rho_cov) - 1e-12))
zone_perm_p = float(np.mean(abs(perm_zone) >= abs(rho_zone) - 1e-12))
loo = []
for name in names:
    sub = coverage[coverage.pyocin != name]
    loo.append({"omitted_pyocin": name, "rho_length_coverage": float(spearmanr(sub.length_aa, sub.all_positive).statistic)})
pd.DataFrame(loo).to_csv(out / "thai_length_leave_one_out.csv", index=False)
union_rows = []
from itertools import combinations
for k in range(1, 7):
    for group in combinations(names, k):
        union_rows.append({"pyocins": "+".join(group), "n_pyocins": k, "covered_isolates": int(binary_matrix[list(group)].any(axis=1).sum())})
pd.DataFrame(union_rows).to_csv(out / "thai_response_set_unions.csv", index=False)
pair_rows = []
for name in names[1:]:
    a = binary_matrix.AP41.astype(int)
    b = binary_matrix[name].astype(int)
    ab = int(((a == 1) & (b == 0)).sum())
    ba = int(((a == 0) & (b == 1)).sum())
    p = binomtest(min(ab, ba), ab + ba, .5).pvalue if ab + ba else 1.0
    pair_rows.append({"comparison": f"AP41 vs {name}", "AP41_only": ab, "other_only": ba, "p": p})
pair = pd.DataFrame(pair_rows)
pair["holm_p"] = multipletests(pair.p, method="holm")[1]
pair.to_csv(out / "thai_paired_coverage.csv", index=False)
mic_matrix = mic.pivot(index="isolate", columns="pyocin", values="mic_ug_ml")[names[:3]]
mic_complete = mic_matrix.dropna()
friedman_mic = friedmanchisquare(*[np.log10(mic_complete[n]) for n in names[:3]]) if len(mic_complete) >= 3 else None
mic_summary = mic.groupby("pyocin", sort=False).agg(n_numeric=("mic_ug_ml", "count"), median_mic_ug_ml=("mic_ug_ml", "median")).reset_index()
mic_summary.to_csv(out / "thai_mic_summary.csv", index=False)

form_rows = []
for _, row in raw_form.iloc[2:].iterrows():
    isolate = row.iloc[0]
    if isinstance(isolate, str) and isolate.startswith("PATH-"):
        vals = [float(row.iloc[k]) for k in range(1, 5)]
        form_rows.append({"isolate": isolate, "mean_A550": float(np.mean(vals)), "sd_A550": float(np.std(vals, ddof=1))})
formation = pd.DataFrame(form_rows)
assert len(formation) == 29
formation.to_csv(out / "thai_biofilm_formation.csv", index=False)
formation_breadth = formation.merge(sus.groupby("isolate").all_positive.sum().rename("positive_pyocins"), on="isolate")
biofilm_breadth_rho, biofilm_breadth_p = spearmanr(formation_breadth.mean_A550, formation_breadth.positive_pyocins)
formation_breadth.to_csv(out / "thai_biofilm_formation_breadth.csv", index=False)

bio_rows = []
current_isolate = None
for _, row in raw_bio.iloc[3:].iterrows():
    if isinstance(row.iloc[0], str) and row.iloc[0].startswith("PATH-"):
        current_isolate = row.iloc[0]
    concentration = row.iloc[1]
    if current_isolate is None or not isinstance(concentration, (int, float, np.number)) or pd.isna(concentration):
        continue
    for j, treatment in enumerate(["AP41", "SX1", "SX2", "Amikacin", "Gentamicin"]):
        vals = [float(row.iloc[2 + 4 * j + k]) for k in range(4)]
        bio_rows.append({"isolate": current_isolate, "concentration_ug_ml": float(concentration), "treatment": treatment, "mean_A405": float(np.mean(vals)), "sd_A405": float(np.std(vals, ddof=1))})
bio = pd.DataFrame(bio_rows)
assert bio.isolate.nunique() == 3
base_a405 = bio[bio.concentration_ug_ml == 0].set_index(["isolate", "treatment"]).mean_A405
bio["relative_A405"] = [r.mean_A405 / base_a405.loc[(r.isolate, r.treatment)] for r in bio.itertuples()]
bio.to_csv(out / "thai_biofilm_dose_response.csv", index=False)

old = pd.read_csv(out / "spot_isolate_profiles.csv")
old["paeruginosa"] = old.Species.eq("Pseudomonas aeruginosa")
for rule in ["published", "strict"]:
    positive = ["+"] if rule == "strict" else ["+", "T"]
    old["any_" + rule] = old[["PyoL1", "PyoL2", "PyoL3"]].isin(positive).any(axis=1).astype(int)
    old["count_" + rule] = old[["PyoL1", "PyoL2", "PyoL3"]].isin(positive).sum(axis=1)
species = old.groupby("paeruginosa").agg(n=("Isolate", "size"), any_published=("any_published", "sum"), any_strict=("any_strict", "sum")).reset_index()
species.to_csv(out / "ltype_species_boundary.csv", index=False)
abr = old[old.paeruginosa & old.ABR.astype(str).str.fullmatch("[0-5]")].copy()
abr["ABR_num"] = abr.ABR.astype(int)
abr_rows = []
for source_name, sub in [("all_annotated", abr), ("Pirnay_only", abr[abr.Source == "JP. Pirnay"])]:
    for rule in ["published", "strict"]:
        rho, p = spearmanr(sub.ABR_num, sub["count_" + rule])
        abr_rows.append({"subset": source_name, "rule": rule, "n": len(sub), "rho": rho, "p_unadjusted": p, "positive_any": int(sub["any_" + rule].sum())})
abr_results = pd.DataFrame(abr_rows)
abr_results.to_csv(out / "ltype_antibiotic_resistance.csv", index=False)
abr.groupby(["Source", "ABR_num"]).agg(n=("Isolate", "size"), any_published=("any_published", "sum"), any_strict=("any_strict", "sum")).reset_index().to_csv(out / "ltype_antibiotic_resistance_strata.csv", index=False)

s6 = pd.DataFrame([
    ("FpvAI", 0, 2), ("FpvAI + FpvB", 23, 11),
    ("FpvAIIb", 20, 0), ("FpvAIIb + FpvB", 16, 7),
    ("FpvAIII", 0, 1), ("FpvAIII + FpvB", 29, 1)
], columns=["receptor_category", "resistant", "sensitive"])
s6["n"] = s6.resistant + s6.sensitive
s6["sensitive_fraction"] = s6.sensitive / s6.n
assert s6.n.sum() == 110 and s6.sensitive.sum() == 22
s6.to_csv(out / "s6_receptor_external_cohort.csv", index=False)

summary = {
    "source": "Laosena et al. 2026, DOI 10.7717/peerj.21850, supplemental Data S1 and Table S2",
    "raw_xlsx_sha256": hashlib.sha256(raw).hexdigest(),
    "thai_n_isolates": int(sus.isolate.nunique()),
    "thai_six_pyocins": names,
    "thai_zone_replicate_rows": len(sus),
    "thai_zone_discordant_triplicates": discordant,
    "thai_friedman_zone_p": float(friedman_zone.pvalue),
    "thai_friedman_zone_statistic": float(friedman_zone.statistic),
    "thai_length_coverage_spearman": float(rho_cov),
    "thai_length_coverage_exact_permutation_p": coverage_perm_p,
    "thai_length_zone_spearman": float(rho_zone),
    "thai_length_zone_exact_permutation_p": zone_perm_p,
    "thai_length_coverage_loo_rho_min": float(min(r["rho_length_coverage"] for r in loo)),
    "thai_length_coverage_loo_rho_max": float(max(r["rho_length_coverage"] for r in loo)),
    "thai_union_AP41_SX1": int(binary_matrix[["AP41", "SX1"]].any(axis=1).sum()),
    "thai_mic_numeric_counts": {r.pyocin: int(r.n_numeric) for r in mic_summary.itertuples()},
    "thai_mic_complete_cases": len(mic_complete),
    "thai_friedman_mic_p": float(friedman_mic.pvalue) if friedman_mic is not None else None,
    "thai_biofilm_formation_n": len(formation),
    "thai_biofilm_formation_breadth_rho": float(biofilm_breadth_rho),
    "thai_biofilm_formation_breadth_p": float(biofilm_breadth_p),
    "thai_biofilm_dose_response_strains": int(bio.isolate.nunique()),
    "ltype_all_species_n": len(old),
    "ltype_paeruginosa_n": int(old.paeruginosa.sum()),
    "ltype_outgroup_n": int((~old.paeruginosa).sum()),
    "ltype_abr_annotated_n": len(abr),
    "ltype_abr_pirnay_n": int((abr.Source == "JP. Pirnay").sum()),
    "s6_n": int(s6.n.sum()),
    "s6_sensitive_n": int(s6.sensitive.sum())
}
(out / "scope_expansion_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12, "axes.titleweight": "bold", "axes.spines.top": False, "axes.spines.right": False})
fig = plt.figure(figsize=(18, 11), constrained_layout=True)
gs = GridSpec(2, 3, figure=fig, width_ratios=[1.5, 1, 1], height_ratios=[1, 1])
ax = fig.add_subplot(gs[:, 0])
ordered = zone_matrix.sort_values(["AP41", "SX1", "SX2"], ascending=False)
cmap = LinearSegmentedColormap.from_list("zones", ["#f8f9fc", "#dbe9ef", "#55b8ba", "#167a9d", "#193e68"])
im = ax.imshow(ordered.values, cmap=cmap, vmin=0, vmax=max(20, np.ceil(ordered.values.max())))
ax.set_xticks(range(6), names, rotation=45, ha="right")
ax.set_yticks(range(29), ordered.index, fontsize=7)
ax.set_title("A  Six proteins, one 29-isolate panel", loc="left", pad=10)
ax.set_xlabel("Pyocin")
fig.colorbar(im, ax=ax, fraction=.035, pad=.02, label="Zone diameter (mm)")

ax = fig.add_subplot(gs[0, 1])
for r in coverage.itertuples():
    ax.scatter(r.length_aa, r.all_positive, s=130, color="#1e819b" if r.pyocin in names[:3] else "#c77e4c", edgecolor="white", linewidth=1.1)
    ax.annotate(str(r.pyocin), (r.length_aa, r.all_positive), xytext=(7, 3), textcoords="offset points", fontsize=9)
ax.set_xlim(250, 820)
ax.set_ylim(0, 31)
ax.set_xlabel("Reference protein length (aa)")
ax.set_ylabel("Isolates with inhibition (of 29)")
ax.set_title("B  Length and response breadth", loc="left", pad=10)
ax.text(.02, .02, "Six protein-level observations", transform=ax.transAxes, fontsize=8, color="#555b69")

ax = fig.add_subplot(gs[0, 2])
for i, name in enumerate(names[:3]):
    vals = mic_matrix[name].dropna().values
    jitter = np.random.default_rng(75 + i).normal(0, .06, len(vals))
    ax.scatter(np.full(len(vals), i) + jitter, vals, s=38, alpha=.82, color=["#194f73", "#2c9ba8", "#c67c3e"][i], edgecolor="white", linewidth=.4)
    ax.hlines(np.median(vals), i-.22, i+.22, color="#202839", linewidth=2)
ax.set_yscale("log")
ax.set_xticks(range(3), names[:3])
ax.set_ylabel("Published spot-dilution MIC (µg/mL)")
ax.set_title("C  MICs for three selected pyocins", loc="left", pad=10)
ax.text(.02, .02, "NC retained as unquantified", transform=ax.transAxes, fontsize=8, color="#555b69")

ax = fig.add_subplot(gs[1, 1])
bio100 = bio[bio.concentration_ug_ml == 100].pivot(index="isolate", columns="treatment", values="relative_A405")[["AP41", "SX1", "SX2", "Amikacin", "Gentamicin"]]
biocmap = LinearSegmentedColormap.from_list("bio", ["#20486c", "#64b7b7", "#f4dfaa", "#bd634e"])
im2 = ax.imshow(bio100.values, cmap=biocmap, vmin=0, vmax=1.2, aspect="auto")
ax.set_xticks(range(5), bio100.columns)
ax.set_yticks(range(3), bio100.index)
ax.set_title("D  Biofilm signal at 100 µg/mL", loc="left", pad=10)
for i in range(3):
    for j in range(5):
        v = bio100.iloc[i, j]
        ax.text(j, i, f"{v:.2f}", ha="center", va="center", color="white" if v < .4 or v > .9 else "#1b2433", fontsize=9)
fig.colorbar(im2, ax=ax, fraction=.035, pad=.02, label="Fraction of untreated A405")

ax = fig.add_subplot(gs[1, 2])
plot_s6 = s6[s6.n >= 20].iloc[::-1]
ax.barh(plot_s6.receptor_category, plot_s6.sensitive_fraction * 100, color=["#65b9b1", "#287c9d", "#91c7c3", "#bd7a53"])
for i, row in enumerate(plot_s6.itertuples()):
    ax.text(row.sensitive_fraction * 100 + 1, i, f"{row.sensitive}/{row.n}", va="center", fontsize=9)
ax.set_xlim(0, 50)
ax.set_xlabel("Sensitive isolates (%)")
ax.set_title("E  Independent S6 receptor panel", loc="left", pad=10)
fig.savefig(out / "Vivaan_Final_Figure_3.png", dpi=300, facecolor="white")
fig.savefig(out / "Vivaan_Final_Figure_3.pdf", facecolor="white")
fig.savefig(out / "Vivaan_Final_Figure_3.svg", facecolor="white")
plt.close(fig)
