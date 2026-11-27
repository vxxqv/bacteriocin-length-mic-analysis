from pathlib import Path
import hashlib
import json
import math
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap, LinearSegmentedColormap, LogNorm, Normalize, TwoSlopeNorm
from matplotlib.patches import Circle, Ellipse, FancyArrowPatch, FancyBboxPatch, Rectangle, PathPatch
from matplotlib.path import Path as MplPath
from matplotlib.colorbar import ColorbarBase


ROOT = Path(__file__).resolve().parent
DATA = pd.read_csv(ROOT / "data" / "expanded_mic_dataset.csv")
FEATURES = pd.read_csv(ROOT / "data" / "expanded_protein_features.csv").set_index("bacteriocin")
IDENTITY = pd.read_csv(ROOT / "data" / "expanded_clustalo_percent_identity_matrix.csv", index_col=0)
IUPRED_PROFILE = pd.read_csv(ROOT / "data" / "iupred3_l_type_profiles.csv")
if len(IUPRED_PROFILE) != 802 or IUPRED_PROFILE.groupby("protein").size().to_dict() != {"L1": 256, "L2": 256, "L3": 290}:
    raise ValueError("Unexpected L-type IUPred3 profile size")
SOURCE = ROOT / "data" / "Ghequire_2014_Table_S1.xlsx"
EXPECTED_HASH = "1dc915c57866cb464d7dd6ecd158ab6e823562cde44d0402de429d47f2b2ca3f"
actual_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
if actual_hash != EXPECTED_HASH:
    raise ValueError("Published supplementary workbook hash changed")

raw = pd.read_excel(SOURCE, header=2)
if len(raw) != 242 or raw["Isolate"].nunique() != 242:
    raise ValueError("Unexpected supplementary table size or isolate identifiers")
raw["Species"] = raw["Species"].ffill()
for name in ["PyoL1", "PyoL2", "PyoL3"]:
    if not set(raw[name].unique()).issubset({"+", "-", "T"}):
        raise ValueError((name, sorted(set(raw[name].unique()))))
pa = raw.loc[raw["Species"].eq("Pseudomonas aeruginosa")].copy()
if len(pa) != 150 or pa["Isolate"].nunique() != 150:
    raise ValueError("Pseudomonas aeruginosa panel is not 150 distinct isolates")
raw.to_csv(ROOT / "spot_isolate_profiles.csv", index=False)


def exact_p(b, c):
    n = b + c
    return min(1.0, 2.0 * sum(math.comb(n, k) for k in range(min(b, c) + 1)) / 2 ** n) if n else 1.0


pair_rows = []
summary = {"source_doi": "10.1002/mbo3.210", "source_supplement": "Table S1", "panel_isolates": 150, "source_all_species_isolates": 242}
names = ["PyoL1", "PyoL2", "PyoL3"]
for definition, positive_codes in [("published_counting", {"+", "T"}), ("strict_plus_only", {"+"})]:
    indicators = {name: pa[name].isin(positive_codes).to_numpy(dtype=bool) for name in names}
    summary[definition] = {name: int(indicators[name].sum()) for name in names}
    summary[definition]["any_pyocin"] = int(np.logical_or.reduce(list(indicators.values())).sum())
    for i, first in enumerate(names):
        for second in names[i + 1:]:
            a = indicators[first]
            b = indicators[second]
            both = int(np.logical_and(a, b).sum())
            first_only = int(np.logical_and(a, ~b).sum())
            second_only = int(np.logical_and(~a, b).sum())
            neither = int(np.logical_and(~a, ~b).sum())
            pair_rows.append({"definition": definition, "first": first, "second": second,
                              "both": both, "first_only": first_only, "second_only": second_only,
                              "neither": neither, "union": both + first_only + second_only,
                              "jaccard": both / (both + first_only + second_only),
                              "paired_risk_difference_second_minus_first": (second_only - first_only) / len(pa),
                              "exact_mcnemar_p": exact_p(first_only, second_only)})
pairs = pd.DataFrame(pair_rows)
for definition, part in pairs.groupby("definition"):
    indices = part.sort_values("exact_mcnemar_p").index.tolist()
    running = 0.0
    for rank, index in enumerate(indices):
        running = max(running, min(1.0, (len(indices) - rank) * float(pairs.loc[index, "exact_mcnemar_p"])))
        pairs.loc[index, "holm_adjusted_p"] = running
pairs.to_csv(ROOT / "spot_pairwise_stats.csv", index=False)
patterns = pa[names].astype(str).agg("".join, axis=1).value_counts().rename_axis("raw_pattern").reset_index(name="isolates")
patterns.to_csv(ROOT / "spot_response_patterns.csv", index=False)
summary["raw_codes"] = {name: {code: int((pa[name] == code).sum()) for code in ["+", "T", "-"]} for name in names}
summary["paired_results"] = pairs.to_dict(orient="records")
summary["all_negative_published_counting"] = int((pa[names] == "-").all(axis=1).sum())
(ROOT / "spot_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
(ROOT / "source_manifest.json").write_text(json.dumps({"source_doi": "10.1002/mbo3.210", "pmcid": "PMC4263511", "source_file": "Ghequire_2014_Table_S1.xlsx", "source_url": "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC4263511/supplementaryFiles", "sha256": actual_hash, "sheet": "Sheet1", "header_row_one_based": 3, "raw_codes": ["+", "T", "-"], "analysis_note": "The published positive totals include + and T. Strict-plus-only is reported as a sensitivity analysis because the supplementary spreadsheet does not define T."}, indent=2), encoding="utf-8")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8, "text.color": "#142439", "axes.labelcolor": "#142439",
                     "svg.fonttype": "none", "pdf.fonttype": 42, "savefig.pad_inches": 0.08})
ink = "#142439"
muted = "#596b7d"
teal = "#087f82"
orange = "#d57239"
soft = "#edf1f3"
mic_cmap = LinearSegmentedColormap.from_list("mic", ["#0b526e", "#2795a1", "#86c9aa", "#f2d37a", "#e98552", "#993f53"])
fold_cmap = LinearSegmentedColormap.from_list("fold", ["#e3f6f0", "#86c9bc", "#238f91", "#315c87", "#433f76"])
contrast_cmap = LinearSegmentedColormap.from_list("contrast", ["#205b9a", "#a9c5dc", "#f8faf9", "#f1b28b", "#bd572c"])


def save(fig, number):
    (ROOT / "supporting_figures").mkdir(exist_ok=True)
    stem = ROOT / "supporting_figures" / f"Vivaan_Reanalysis_Figure_{number}"
    for suffix in ["png", "svg", "pdf"]:
        fig.savefig(stem.with_suffix("." + suffix), dpi=360, facecolor="white", bbox_inches="tight")
    plt.close(fig)


def grid(ax, rows, cols, values, labels, cmap, norm, text_sizes=None, dark_low=True):
    masked = np.ma.masked_invalid(values)
    plot_cmap = cmap.copy()
    plot_cmap.set_bad(soft)
    ax.imshow(masked, cmap=plot_cmap, norm=norm, aspect="auto", interpolation="nearest")
    ax.set_xticks(range(len(cols)), cols)
    ax.set_yticks(range(len(rows)), rows)
    ax.tick_params(axis="both", length=0, pad=7)
    ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False)
    ax.set_xticks(np.arange(-0.5, len(cols), 1), minor=True)
    ax.set_yticks(np.arange(-0.5, len(rows), 1), minor=True)
    ax.grid(which="minor", color="white", linewidth=2.5)
    ax.tick_params(which="minor", bottom=False, left=False)
    for spine in ax.spines.values():
        spine.set_visible(False)
    for i in range(len(rows)):
        for j in range(len(cols)):
            value = values[i, j]
            text_value = labels[i][j]
            if not np.isfinite(value):
                color = "#718092"
            else:
                fraction = norm(value)
                color = "white" if (dark_low and fraction < 0.18) or fraction > 0.77 else ink
            size = text_sizes[i][j] if text_sizes is not None else 8
            ax.text(j, i, text_value, ha="center", va="center", color=color, fontsize=size, fontweight="bold", linespacing=1.15)


def mic_row(protein, strain):
    part = DATA.loc[DATA.doi.eq("10.1002/mbo3.210") & DATA.bacteriocin.eq("Pyocin " + protein) & DATA.target_strain.eq("P. aeruginosa " + strain)]
    if len(part) > 1:
        raise ValueError((protein, strain, len(part)))
    return None if part.empty else part.iloc[0]


strains = ["Bu007", "CFPA13", "CFPA22", "CFPA101", "CFPA120", "CFPA124"]
proteins = ["L1", "L2", "L3"]
values = np.full((3, 6), np.nan)
labels = [["NR" for _ in strains] for _ in proteins]
censored = []
for i, protein in enumerate(proteins):
    for j, strain in enumerate(strains):
        row = mic_row(protein, strain)
        if row is None:
            continue
        values[i, j] = float(row.mic_nM_lower)
        if row.censoring == "right":
            labels[i][j] = ">3300"
            censored.append((i, j))
        elif row.censoring == "exact":
            labels[i][j] = f"{float(row.mic_nM_lower):g}"
        else:
            raise ValueError((protein, strain, row.censoring))

fig = plt.figure(figsize=(8.5, 8.0), facecolor="white")
fig.text(0.105, 0.971, "The MIC ranking reverses", fontsize=17, fontweight="bold")
fig.text(0.105, 0.933, "Two equal-length pyocins, two strains, and opposite outcomes under one study protocol.", fontsize=9, color=muted)
fig.text(0.105, 0.899, "L1 and L2: 256 amino acids each  |  86.3% sequence identity", fontsize=8.8, color="#315d79", fontweight="bold")
for x, background, stripe, strain, number, winner in [(0.105, "#edf3fa", "#265c9b", "Bu007", ">305×", "L1 has the lower MIC"),
                                                      (0.525, "#fff0e7", "#c46634", "CFPA22", ">53×", "L2 has the lower MIC")]:
    fig.add_artist(FancyBboxPatch((x, 0.745), 0.37, 0.132, boxstyle="round,pad=0.004,rounding_size=0.018", transform=fig.transFigure, facecolor=background, edgecolor="none"))
    fig.add_artist(Rectangle((x + 0.009, 0.757), 0.006, 0.108, transform=fig.transFigure, facecolor=stripe, edgecolor="none"))
    fig.text(x + 0.035, 0.846, strain, fontsize=9.6, fontweight="bold", color=ink)
    fig.text(x + 0.035, 0.79, number, fontsize=22, fontweight="bold", color=stripe)
    fig.text(x + 0.19, 0.794, winner, fontsize=8.6, color=ink)
fig.text(0.105, 0.705, "A   Published MIC landscape", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.155, 0.435, 0.79, 0.225])
mic_norm = LogNorm(vmin=6.9, vmax=3300)
grid(ax, ["L1  256 aa", "L2  256 aa", "L3  290 aa"], strains, values, labels, mic_cmap, mic_norm)
for i, j in censored:
    ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, facecolor="none", hatch="///", edgecolor="#4d3147", linewidth=0.5))
fig.text(0.155, 0.402, "Low MIC / stronger inhibition", fontsize=8, color=muted)
fig.text(0.945, 0.402, "High MIC", fontsize=8, color=muted, ha="right")
cb = fig.add_axes([0.155, 0.382, 0.79, 0.014])
ColorbarBase(cb, cmap=mic_cmap, norm=mic_norm, orientation="horizontal", ticks=[10, 30, 100, 300, 1000, 3300])
cb.set_xticklabels(["10", "30", "100", "300", "1000", "3300"], fontsize=7)
fig.text(0.105, 0.323, "B   Lower MIC within matched strain pairs", fontsize=11, fontweight="bold")
comparison_rows = [("L1 : L2", "L1", "L2"), ("L2 : L3", "L2", "L3")]
contrast = np.full((2, 6), np.nan)
contrast_labels = [["NR" for _ in strains] for _ in comparison_rows]
uncertain = []
for i, (_, first, second) in enumerate(comparison_rows):
    for j, strain in enumerate(strains):
        a = mic_row(first, strain)
        b = mic_row(second, strain)
        if a is None or b is None:
            continue
        if a.censoring == "right" and b.censoring == "right":
            contrast_labels[i][j] = "Both\n>3300"
            uncertain.append((i, j))
            continue
        ratio = float(a.mic_nM_lower) / float(b.mic_nM_lower)
        contrast[i, j] = np.log2(ratio)
        if a.censoring == "right":
            contrast_labels[i][j] = f"{second}\n>{math.floor(ratio)}x"
        elif b.censoring == "right":
            contrast_labels[i][j] = f"{first}\n>{math.floor(1 / ratio)}x"
        else:
            winner = second if ratio >= 1 else first
            contrast_labels[i][j] = f"{winner}\n{max(ratio, 1 / ratio):.2g}x"
ax = fig.add_axes([0.155, 0.135, 0.79, 0.15])
contrast_norm = TwoSlopeNorm(vmin=-3, vcenter=0, vmax=3)
grid(ax, [x[0] for x in comparison_rows], strains, contrast, contrast_labels, contrast_cmap, contrast_norm)
ax.tick_params(axis="x", top=False, labeltop=False, bottom=True, labelbottom=True)
for i, j in uncertain:
    ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, facecolor="none", hatch="xx", edgecolor="#8e9ba6", linewidth=0.4))
fig.text(0.105, 0.080, "Blue: first protein lower   |   Orange: second protein lower   |   NR: no paired MIC", fontsize=8.2, color=muted)
fig.text(0.105, 0.048, "Contrast color saturates at eightfold; labels preserve larger bounds and exact ratios.", fontsize=8.2, color=muted)
save(fig, 1)


def block_matrix(targets, columns, family, doi):
    result = np.full((len(targets), len(columns)), np.nan)
    cell_text = [["NR" for _ in columns] for _ in targets]
    for i, target in enumerate(targets):
        part = DATA.loc[DATA.doi.eq(doi) & DATA.target_strain.eq(target) & DATA.family.eq(family) & DATA.censoring.eq("exact")]
        if not len(part):
            raise ValueError((doi, target, family))
        minimum = float(part.mic_nM_lower.min())
        for j, protein in enumerate(columns):
            row = part.loc[part.bacteriocin.eq(protein)]
            if len(row) > 1:
                raise ValueError((target, protein))
            if row.empty:
                continue
            mic = float(row.iloc[0].mic_nM_lower)
            fold = mic / minimum
            result[i, j] = math.log2(fold)
            cell_text[i][j] = f"{mic:.3g} nM\n{fold:.2g}x"
    return result, cell_text


fig = plt.figure(figsize=(8.0, 7.0), facecolor="white")
fig.text(0.12, 0.967, "Inside the matched comparisons", fontsize=16, fontweight="bold")
fig.text(0.12, 0.929, "Each colored row compares proteins on one target under one study protocol.", fontsize=9, color=muted)
fig.text(0.12, 0.875, "A   How much of the MIC collection can be compared?", fontsize=11, fontweight="bold")
stage_data = [("47", "all MICs", "22 proteins  |  9 studies"),
              ("29", "matched MICs", "11 blocks  |  17 proteins  |  3 studies"),
              ("16", "exact matched", "7 blocks  |  8 proteins  |  2 studies")]
for index, (count, title, detail) in enumerate(stage_data):
    x = 0.12 + index * 0.31
    color = ["#eaf0f3", "#d7eae8", "#bce3dc"][index]
    fig.add_artist(FancyBboxPatch((x, 0.754), 0.245, 0.097, boxstyle="round,pad=0.009,rounding_size=0.012", transform=fig.transFigure, facecolor=color, edgecolor="none"))
    fig.text(x + 0.018, 0.81, count, fontsize=17, fontweight="bold", color=ink)
    fig.text(x + 0.084, 0.812, title, fontsize=9, fontweight="bold", color=ink)
    fig.text(x + 0.018, 0.775, detail, fontsize=6.9, color=muted)
    if index < 2:
        fig.add_artist(FancyArrowPatch((x + 0.255, 0.802), (x + 0.296, 0.802), arrowstyle="-|>", mutation_scale=12, linewidth=1.4, color="#7c9aa8", transform=fig.transFigure))
fig.text(0.15, 0.722, "B   L-type pyocins", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.15, 0.269, 0.31, 0.385])
targets = ["P. aeruginosa " + x for x in ["CFPA13", "CFPA22", "CFPA120", "CFPA124"]]
matrix, cell_text = block_matrix(targets, ["Pyocin L2", "Pyocin L3"], "lectin-like", "10.1002/mbo3.210")
fold_norm = Normalize(vmin=0, vmax=3)
grid(ax, [x.split()[-1] for x in targets], ["L2\n256 aa", "L3\n290 aa"], matrix, cell_text, fold_cmap, fold_norm, dark_low=False)
fig.text(0.592, 0.722, "C   ColM-like klebicins", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.59, 0.462, 0.32, 0.192])
targets = ["K. pneumoniae DSM16358", "K. variicola DSM15968"]
matrix, cell_text = block_matrix(targets, ["KpneM2", "KvarM", "KpneM"], "ColM-like", "10.1038/s41598-019-51969-1")
grid(ax, ["K. pneumoniae", "K. variicola"], ["KpneM2\n272 aa", "KvarM\n276 aa", "KpneM\n278 aa"], matrix, cell_text, fold_cmap, fold_norm, dark_low=False)
fig.text(0.592, 0.42, "D   Pore-forming klebicins", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.59, 0.269, 0.32, 0.083])
matrix, cell_text = block_matrix(["K. variicola DSM15968"], ["KpneA", "KaerA", "KvarIa"], "pore-forming", "10.1038/s41598-019-51969-1")
grid(ax, ["K. variicola"], ["KpneA\n377 aa", "KaerA\n379 aa", "KvarIa\n414 aa"], matrix, cell_text, fold_cmap, fold_norm, dark_low=False)
fig.text(0.15, 0.202, "Within each row: MIC relative to the lowest MIC", fontsize=9, color=muted)
cb = fig.add_axes([0.15, 0.165, 0.66, 0.025])
ColorbarBase(cb, cmap=fold_cmap, norm=fold_norm, orientation="horizontal", ticks=[0, 1, 2, 3])
cb.set_xticklabels(["1x  lowest", "2x", "4x", "8x or more"], fontsize=8)
fig.text(0.15, 0.075, "Values inside cells are exact MICs in nM and fold differences within a matched row. Gray means unreported.", fontsize=8.3, color=muted)
save(fig, 2)


code_map = {"-": 0, "T": 1, "+": 2}
pattern_order = patterns.sort_values(["isolates", "raw_pattern"], ascending=[False, True]).reset_index(drop=True)
code_matrix = np.array([[code_map[code] for code in pattern] for pattern in pattern_order.raw_pattern], dtype=int)
fig = plt.figure(figsize=(8.6, 7.1), facecolor="white")
fig.text(0.11, 0.967, "The 150-isolate response fingerprint", fontsize=16, fontweight="bold")
fig.text(0.11, 0.929, "Fourteen raw response signatures summarize one shared P. aeruginosa spot-assay panel.", fontsize=9, color=muted)
fig.text(0.11, 0.875, "A   Raw signatures and their sizes", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.15, 0.165, 0.25, 0.655])
code_cmap = ListedColormap(["#edf1f3", "#f0b45e", "#087f82"])
ax.imshow(code_matrix, cmap=code_cmap, vmin=-0.5, vmax=2.5, aspect="auto", interpolation="nearest")
ax.set_xticks(range(3), ["L1", "L2", "L3"])
ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=8, labelsize=10)
ax.set_yticks(range(len(pattern_order)), [str(x) for x in pattern_order.isolates])
ax.tick_params(axis="y", length=0, pad=7, labelsize=8)
ax.set_ylabel("Number of isolates", labelpad=8)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.set_xticks(np.arange(-0.5, 3, 1), minor=True)
ax.set_yticks(np.arange(-0.5, len(pattern_order), 1), minor=True)
ax.grid(which="minor", color="white", linewidth=2)
ax.tick_params(which="minor", bottom=False, left=False)
for row, pattern in enumerate(pattern_order.raw_pattern):
    for column, code in enumerate(pattern):
        ax.text(column, row, code, ha="center", va="center", fontsize=8.5, fontweight="bold", color="white" if code == "+" else ink)
ax_bar = fig.add_axes([0.425, 0.165, 0.105, 0.655])
bar_colors = LinearSegmentedColormap.from_list("signature_abundance", ["#badad6", "#218d90", "#154e75"])
for row, count in enumerate(pattern_order.isolates):
    ax_bar.barh(row, count, height=0.73, color=bar_colors(math.sqrt(count / 52)), edgecolor="none")
ax_bar.set_xlim(0, 52)
ax_bar.set_ylim(len(pattern_order) - 0.5, -0.5)
ax_bar.set_xticks([])
ax_bar.set_yticks([])
for spine in ax_bar.spines.values():
    spine.set_visible(False)
fig.text(0.427, 0.832, "Frequency", fontsize=7.5, color=muted)
for x, label, color in [(0.15, "-  no inhibition", "#edf1f3"), (0.275, "T  reported separately", "#f0b45e"), (0.458, "+  positive symbol", "#087f82")]:
    fig.add_artist(Rectangle((x, 0.082), 0.018, 0.018, transform=fig.transFigure, facecolor=color, edgecolor="#afbbc4" if color == "#edf1f3" else "none"))
    fig.text(x + 0.024, 0.084, label, fontsize=7.7, color=muted)
fig.text(0.57, 0.875, "B   Coverage depends on the code rule", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.59, 0.505, 0.32, 0.305])
clear = np.array([summary["raw_codes"][x]["+"] for x in names])
unclear = np.array([summary["raw_codes"][x]["T"] for x in names])
ax.barh([0, 1, 2], clear, color=teal, height=0.55, label="+")
ax.barh([0, 1, 2], unclear, left=clear, color="#f0b45e", height=0.55, label="T")
ax.set_yticks([0, 1, 2], ["L1", "L2", "L3"])
ax.set_ylim(2.5, -0.5)
ax.set_xlim(0, 80)
ax.set_xticks([0, 20, 40, 60, 80])
ax.grid(axis="x", color="#e2e8eb", linewidth=0.8)
ax.set_axisbelow(True)
ax.tick_params(axis="y", length=0)
for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
for y, (a, b) in enumerate(zip(clear, unclear)):
    ax.text(a + b + 1.4, y, f"{a + b}", va="center", fontsize=9, fontweight="bold")
fig.text(0.57, 0.442, "C   Shared susceptible isolates", fontsize=11, fontweight="bold")
jaccard = np.full((3, 3), np.nan)
for _, row in pairs.loc[pairs.definition.eq("published_counting")].iterrows():
    i = names.index(row["first"])
    j = names.index(row["second"])
    jaccard[i, j] = jaccard[j, i] = float(row.jaccard)
ax = fig.add_axes([0.60, 0.145, 0.285, 0.255])
overlap_cmap = LinearSegmentedColormap.from_list("overlap", ["#f1f3f4", "#b3d7d0", "#439e9c", "#245d83"])
grid(ax, ["L1", "L2", "L3"], ["L1", "L2", "L3"], jaccard,
     [["same" if i == j else f"{100 * jaccard[i, j]:.1f}%" for j in range(3)] for i in range(3)], overlap_cmap, Normalize(vmin=0, vmax=0.4), dark_low=False)
fig.text(0.11, 0.035, "Published totals include both + and T. The meaning of T is not defined in the source sheet, so strict + results are reported separately.", fontsize=7.6, color=muted)
save(fig, 3)


members = ["Pyocin L1", "Pyocin L2", "Pyocin L3"]
identity_values = IDENTITY.loc[members, members].to_numpy(dtype=float)
fig = plt.figure(figsize=(8.6, 7.4), facecolor="white")
fig.text(0.09, 0.963, "A sequence map beyond total length", fontsize=16, fontweight="bold")
fig.text(0.09, 0.926, "Alignment, predicted regions, and residue-level disorder are shown for the same three pyocins.", fontsize=8.8, color=muted)
fig.text(0.09, 0.859, "A   Clustal Omega identity", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.11, 0.542, 0.36, 0.265])
identity_cmap = LinearSegmentedColormap.from_list("identity", ["#eef2f4", "#b8d7d1", "#57aaa7", "#0e6680", "#173f68"])
grid(ax, ["L1", "L2", "L3"], ["L1", "L2", "L3"], identity_values,
     [[f"{identity_values[i, j]:.1f}%" for j in range(3)] for i in range(3)], identity_cmap, Normalize(vmin=0, vmax=100), dark_low=False)
cb = fig.add_axes([0.11, 0.508, 0.36, 0.016])
ColorbarBase(cb, cmap=identity_cmap, norm=Normalize(vmin=0, vmax=100), orientation="horizontal", ticks=[0, 25, 50, 75, 100])
cb.set_xticklabels(["0", "25", "50", "75", "100%"], fontsize=8)
fig.text(0.55, 0.859, "B   Protein length and SMART hits", fontsize=11, fontweight="bold")
ax = fig.add_axes([0.58, 0.566, 0.36, 0.225])
for y, member in zip([2, 1, 0], members):
    feature = FEATURES.loc[member]
    length = int(feature.length_aa)
    ax.add_patch(FancyBboxPatch((0, y - 0.13), length, 0.26, boxstyle="round,pad=0.02,rounding_size=0.12", facecolor="#dbe5e9", edgecolor="none"))
    for item in str(feature.domain_coordinates).split(";"):
        first, last = [int(x.strip()) for x in item.split("-")]
        ax.add_patch(FancyBboxPatch((first, y - 0.13), last - first + 1, 0.26, boxstyle="round,pad=0.02,rounding_size=0.08", facecolor="#d57239", edgecolor="none"))
    ax.text(length + 7, y, f"{length} aa", va="center", fontsize=8.2, color=ink)
ax.set_xlim(0, 340)
ax.set_ylim(-0.6, 2.6)
ax.set_yticks([2, 1, 0], ["L1", "L2", "L3"])
ax.set_xticks([0, 100, 200, 300])
ax.set_xlabel("Amino acid position", labelpad=4)
ax.tick_params(axis="y", length=0, pad=9)
ax.grid(axis="x", color="#e7ecef", linewidth=0.8)
ax.set_axisbelow(True)
for spine in ["top", "right", "left"]:
    ax.spines[spine].set_visible(False)
fig.add_artist(Rectangle((0.58, 0.478), 0.018, 0.018, transform=fig.transFigure, facecolor="#d57239", edgecolor="none"))
fig.text(0.605, 0.480, "SMART SM00108 predicted region", fontsize=8.4, color=muted)
fig.text(0.09, 0.438, "C   IUPred3 long-mode disorder profile", fontsize=11, fontweight="bold")
profile_colors = {"L1": "#275d9b", "L2": "#c96939", "L3": "#0d8684"}
for index, name in enumerate(["L1", "L2", "L3"]):
    part = IUPRED_PROFILE.loc[IUPRED_PROFILE.protein.eq(name)]
    ax = fig.add_axes([0.15, 0.335 - index * 0.095, 0.79, 0.077])
    ax.fill_between(part.position, 0, part.score, color=profile_colors[name], alpha=0.24, linewidth=0)
    ax.plot(part.position, part.score, color=profile_colors[name], linewidth=1.5)
    ax.axhline(0.5, color="#7f8d99", linestyle=(0, (3, 3)), linewidth=0.8)
    ax.set_xlim(0, 300)
    ax.set_ylim(0, 0.65)
    ax.set_yticks([0, 0.5], ["0", "0.5"])
    ax.tick_params(axis="y", length=0, labelsize=7, colors=muted)
    ax.tick_params(axis="x", length=0, labelsize=7, colors=muted, labelbottom=index == 2)
    if index == 2:
        ax.set_xticks([0, 50, 100, 150, 200, 250, 300])
        ax.set_xlabel("Amino acid position", fontsize=8.5, labelpad=2)
    else:
        ax.set_xticks([])
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.spines["bottom"].set_color("#dbe4e8")
    fig.text(0.09, 0.36 - index * 0.095, name, fontsize=9.3, fontweight="bold", color=profile_colors[name])
fig.text(0.09, 0.046, "Only four L1 residues exceed the 0.5 prediction threshold; L2 and L3 have none. These are predictions, not measured structure.", fontsize=8.25, color=muted)
save(fig, 5)


def mechanism_scene(ax, stage):
    ax.set_xlim(0, 100)
    ax.set_ylim(0, 100)
    ax.set_aspect("equal")
    ax.axis("off")
    ax.add_patch(FancyBboxPatch((1, 1), 98, 98, boxstyle="round,pad=0.015,rounding_size=5", facecolor="#f9fbfb", edgecolor="#dbe6e8", linewidth=1.2))
    ax.text(8, 88, ["1   Recognition", "2   BAM inhibition", "3   Envelope response"][stage], fontsize=10.5, fontweight="bold", color=ink)
    ax.text(8, 79, "Outside target cell", fontsize=7.8, color=muted)
    ax.add_patch(Rectangle((7, 42), 86, 9, facecolor="#72b8b6", edgecolor="none"))
    ax.add_patch(Rectangle((7, 34), 86, 7, facecolor="#b7ded7", edgecolor="none"))
    ax.plot([7, 93], [51, 51], color="#3d8f9c", linewidth=1.2)
    ax.plot([7, 93], [34, 34], color="#3d8f9c", linewidth=1.2)
    for xpos in range(12, 94, 8):
        ax.add_patch(Circle((xpos, 47), 1.25, facecolor="#b8e4de", edgecolor="none"))
        ax.add_patch(Circle((xpos, 37), 1.25, facecolor="#e0f2ed", edgecolor="none"))
    ax.text(8, 27, "Periplasm", fontsize=7.8, color=muted)
    ax.add_patch(FancyBboxPatch((44, 29), 15, 30, boxstyle="round,pad=0.02,rounding_size=4", facecolor="#244b70", edgecolor="#18354e", linewidth=1.2))
    ax.add_patch(FancyBboxPatch((49, 36), 5, 17, boxstyle="round,pad=0.02,rounding_size=2", facecolor="#f9fbfb", edgecolor="none"))
    ax.plot([46, 44, 47, 48], [58, 63, 67, 64], color="#244b70", linewidth=2.5, solid_capstyle="round")
    ax.text(72, 56, "BamA", fontsize=8.5, color="#244b70", fontweight="bold")
    if stage == 0:
        ax.add_patch(Ellipse((22, 69), 21, 14, angle=-15, facecolor="#806fbd", edgecolor="#5a5295", linewidth=1.2))
        ax.add_patch(Ellipse((36, 67), 21, 14, angle=-15, facecolor="#5189bb", edgecolor="#39709e", linewidth=1.2))
        ax.plot([40, 46], [64, 60], color=orange, linewidth=4, solid_capstyle="round")
        ax.add_patch(FancyArrowPatch((42, 66), (48, 64), arrowstyle="-|>", mutation_scale=12, color="#607b91", linewidth=1.5))
        ax.text(9, 15, "L-type pyocin meets\nsurface-exposed BamA", fontsize=8.5, color=ink)
    elif stage == 1:
        ax.add_patch(Ellipse((27, 67), 21, 14, angle=-10, facecolor="#806fbd", edgecolor="#5a5295", linewidth=1.2))
        ax.add_patch(Ellipse((41, 66), 21, 14, angle=-10, facecolor="#5189bb", edgecolor="#39709e", linewidth=1.2))
        ax.plot([49, 51, 51], [62, 55, 44], color=orange, linewidth=4, solid_capstyle="round")
        ax.add_patch(Circle((51, 44), 2.2, facecolor=orange, edgecolor="none"))
        ax.text(7, 15, "C-terminal peptide docks;\nprotein body stays outside", fontsize=8.5, color=ink)
    else:
        ax.add_patch(Ellipse((25, 68), 21, 14, angle=-10, facecolor="#806fbd", edgecolor="#5a5295", linewidth=1.2))
        ax.add_patch(Ellipse((39, 67), 21, 14, angle=-10, facecolor="#5189bb", edgecolor="#39709e", linewidth=1.2))
        ax.plot([47, 51, 51], [62, 55, 44], color=orange, linewidth=4, solid_capstyle="round")
        ax.add_patch(Ellipse((76, 23), 15, 9, facecolor="#d9eacf", edgecolor="#6b9c78", linewidth=1.2))
        ax.text(76, 23, "OMP", ha="center", va="center", fontsize=6.8, fontweight="bold", color="#366e57")
        ax.add_patch(FancyArrowPatch((76, 29), (76, 40), arrowstyle="-|>", mutation_scale=13, color="#5b9a74", linewidth=2.2))
        ax.plot([68, 84], [37, 47], color="#bd543f", linewidth=3)
        ax.plot([68, 84], [47, 37], color="#bd543f", linewidth=3)
        ax.text(7, 15, "OMP insertion stalls;\nenvelope stress follows", fontsize=8.5, color=ink)


fig = plt.figure(figsize=(9.0, 4.5), facecolor="white")
fig.text(0.07, 0.955, "How an L-type pyocin can stop outer-membrane assembly", fontsize=15, fontweight="bold")
fig.text(0.07, 0.905, "A literature-based mechanism, shown separately from the MIC and spot-assay measurements.", fontsize=8.8, color=muted)
for index, x in enumerate([0.07, 0.375, 0.68]):
    ax = fig.add_axes([x, 0.19, 0.255, 0.65])
    mechanism_scene(ax, index)
for x in [0.337, 0.642]:
    fig.add_artist(FancyArrowPatch((x, 0.51), (x + 0.028, 0.51), arrowstyle="-|>", mutation_scale=18, linewidth=1.8, color="#7191a0", transform=fig.transFigure))
fig.text(0.07, 0.075, "Mechanism supported by Munder et al. (2026): BamA binding, C-terminal peptide engagement, reduced OMP assembly, and envelope stress.", fontsize=8.2, color=muted)
save(fig, 7)

fig = plt.figure(figsize=(10.6, 7.0), facecolor="white")
fig.text(0.07, 0.962, "One code, two different coverage conclusions", fontsize=16, fontweight="bold")
fig.text(0.07, 0.919, "Each tile counts the same 150 isolates. Only the definition of a positive spot response changes.", fontsize=9, color=muted)
fig.add_artist(FancyBboxPatch((0.07, 0.859), 0.86, 0.044, boxstyle="round,pad=0.003,rounding_size=0.012", transform=fig.transFigure, facecolor="#e8f3f1", edgecolor="none"))
fig.text(0.088, 0.870, "L3 minus L2:  +34 isolates with + or T", fontsize=11, fontweight="bold", color="#17616e")
fig.text(0.665, 0.870, "→  +1 with + only", fontsize=11, fontweight="bold", color="#b15b31")
pair_labels = [("PyoL1", "PyoL2"), ("PyoL1", "PyoL3"), ("PyoL2", "PyoL3")]
cell_labels = [["both", "first only"], ["second only", "neither"]]
cell_colors = [["#248c87", "#346ea1"], ["#d47745", "#84939d"]]
for column, (first, second) in enumerate(pair_labels):
    x = 0.125 + column * 0.285
    fig.text(x + 0.10, 0.817, f"{first[-2:]} vs {second[-2:]}", ha="center", fontsize=11, fontweight="bold")
    for row, definition in enumerate(["published_counting", "strict_plus_only"]):
        record = pairs.loc[pairs.definition.eq(definition) & pairs.first.eq(first) & pairs.second.eq(second)].iloc[0]
        y = 0.535 if row == 0 else 0.175
        ax = fig.add_axes([x, y, 0.20, 0.235])
        counts = [[int(record.both), int(record.first_only)], [int(record.second_only), int(record.neither)]]
        for i in range(2):
            for j in range(2):
                base = np.array(matplotlib.colors.to_rgb(cell_colors[i][j]))
                strength = 0.25 + 0.67 * math.sqrt(counts[i][j] / 150)
                shade = 1 - (1 - base) * strength
                ax.add_patch(FancyBboxPatch((j + 0.015, i + 0.015), 0.97, 0.97, boxstyle="round,pad=0.015,rounding_size=0.06", facecolor=shade, edgecolor="white", linewidth=2))
                ax.text(j + 0.5, i + 0.39, str(counts[i][j]), ha="center", va="center", fontsize=16, fontweight="bold", color=ink)
                ax.text(j + 0.5, i + 0.72, cell_labels[i][j], ha="center", va="center", fontsize=7.5, color=ink)
        ax.set_xlim(0, 2)
        ax.set_ylim(2, 0)
        ax.set_xticks([0.5, 1.5], [f"{second[-2:]} +", f"{second[-2:]} -"])
        ax.tick_params(axis="x", top=True, labeltop=True, bottom=False, labelbottom=False, length=0, pad=3, labelsize=8)
        ax.set_yticks([0.5, 1.5], [f"{first[-2:]} +", f"{first[-2:]} -"])
        ax.tick_params(axis="y", length=0, pad=5, labelsize=8)
        for spine in ax.spines.values():
            spine.set_visible(False)
        excess = int(record.second_only - record.first_only)
        fig.text(x + 0.10, y - 0.038, f"{second[-2:]} minus {first[-2:]}: {excess:+d}/150", ha="center", fontsize=8.2, color=ink)
        p_text = "< 0.001" if record.holm_adjusted_p < 0.001 else f"= {record.holm_adjusted_p:.3g}"
        fig.text(x + 0.10, y - 0.071, f"Holm-adjusted P {p_text}", ha="center", fontsize=8.2, color=muted)
fig.text(0.028, 0.655, "Published\n+ or T", fontsize=9, color=ink, fontweight="bold", va="center")
fig.text(0.028, 0.295, "Strict\n+ only", fontsize=9, color=ink, fontweight="bold", va="center")
fig.text(0.07, 0.065, "Blue and orange cells are the discordant isolates used by the paired test. Positive counts are spot responses, not MICs.", fontsize=8.2, color=muted)
save(fig, 4)

pattern_counts = {str(row.raw_pattern): int(row.isolates) for row in patterns.itertuples()}
if len(pattern_counts) != 14 or sum(pattern_counts.values()) != 150:
    raise ValueError("Unexpected raw response pattern distribution")
priority = {"+": 0, "T": 1, "-": 2}
positions = []
group_totals = []
for stage in range(3):
    ordered_patterns = sorted(pattern_counts, key=lambda pattern: (priority[pattern[stage]], priority[pattern[(stage + 1) % 3]], priority[pattern[(stage + 2) % 3]], pattern))
    cursor = 0.80
    stage_positions = {}
    stage_groups = {}
    for code in ["+", "T", "-"]:
        members = [pattern for pattern in ordered_patterns if pattern[stage] == code]
        high = cursor
        for pattern in members:
            height = pattern_counts[pattern] / 150 * 0.62
            stage_positions[pattern] = (cursor, cursor - height)
            cursor -= height
        stage_groups[code] = (high, cursor, sum(pattern_counts[pattern] for pattern in members))
        cursor -= 0.020
    positions.append(stage_positions)
    group_totals.append(stage_groups)

fig = plt.figure(figsize=(11.2, 6.8), facecolor="white")
ax = fig.add_axes([0, 0, 1, 1])
ax.set_xlim(0, 1)
ax.set_ylim(0, 1)
ax.axis("off")
fig.text(0.055, 0.952, "One isolate panel, fourteen response signatures", fontsize=17, fontweight="bold")
fig.text(0.055, 0.907, "Flow widths represent isolates. Every path retains its original +, T, or - code at each pyocin.", fontsize=9, color=muted)
stage_x = [0.095, 0.345, 0.595]
stage_colors = {"+": "#087f82", "T": "#e5a446", "-": "#d7e1e5"}
for stage, x in enumerate(stage_x):
    fig.text(x + 0.018, 0.844, f"L{stage + 1}", fontsize=13, fontweight="bold", ha="center")
    totals = group_totals[stage]
    fig.text(x + 0.018, 0.819, f"+ {totals['+'][2]}   T {totals['T'][2]}   - {totals['-'][2]}", fontsize=7.5, color=muted, ha="center")
for stage in [0, 1]:
    start_x = stage_x[stage] + 0.036
    end_x = stage_x[stage + 1]
    width = end_x - start_x
    for pattern in sorted(pattern_counts, key=lambda item: -pattern_counts[item]):
        top0, bottom0 = positions[stage][pattern]
        top1, bottom1 = positions[stage + 1][pattern]
        vertices = [(start_x, top0), (start_x + width * 0.46, top0), (end_x - width * 0.46, top1), (end_x, top1),
                    (end_x, bottom1), (end_x - width * 0.46, bottom1), (start_x + width * 0.46, bottom0), (start_x, bottom0), (start_x, top0)]
        codes = [MplPath.MOVETO, MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.LINETO,
                 MplPath.CURVE4, MplPath.CURVE4, MplPath.CURVE4, MplPath.CLOSEPOLY]
        ax.add_patch(PathPatch(MplPath(vertices, codes), facecolor=stage_colors[pattern[2]], alpha=0.70, edgecolor="white", linewidth=0.65))
for stage, x in enumerate(stage_x):
    for code in ["+", "T", "-"]:
        high, low, count = group_totals[stage][code]
        ax.add_patch(FancyBboxPatch((x, low), 0.036, high - low, boxstyle="round,pad=0.001,rounding_size=0.005", facecolor=stage_colors[code], edgecolor="white", linewidth=1.2))
        if high - low > 0.035:
            ax.text(x + 0.018, (high + low) / 2, code, ha="center", va="center", fontsize=10, fontweight="bold", color="white" if code == "+" else ink)
fig.text(0.735, 0.844, "Observed signatures", fontsize=11.5, fontweight="bold")
fig.text(0.735, 0.813, "Raw L1 / L2 / L3 codes", fontsize=8, color=muted)
for index, (pattern, count) in enumerate(sorted(pattern_counts.items(), key=lambda item: (-item[1], item[0]))):
    y = 0.773 - index * 0.044
    ax.text(0.737, y, pattern, fontsize=8.5, color=ink, va="center", family="DejaVu Sans Mono")
    ax.add_patch(FancyBboxPatch((0.795, y - 0.010), 0.145 * count / 52, 0.020, boxstyle="round,pad=0.001,rounding_size=0.006", facecolor=stage_colors[pattern[2]], edgecolor="none"))
    ax.text(0.953, y, str(count), fontsize=8.5, color=ink, va="center", ha="right", fontweight="bold")
for x, code, label in [(0.055, "+", "positive symbol"), (0.25, "T", "undefined source code"), (0.48, "-", "no reported inhibition")]:
    ax.add_patch(Rectangle((x, 0.088), 0.017, 0.017, facecolor=stage_colors[code], edgecolor="#aebdc5" if code == "-" else "none"))
    ax.text(x + 0.024, 0.090, f"{code}  {label}", fontsize=8.2, color=muted, va="bottom")
fig.text(0.055, 0.037, "Most common: 52 isolates had no response to any pyocin. The source counted T as positive, but did not define it in Table S1.", fontsize=8.4, color=muted)
save(fig, 6)

print(json.dumps({"panel_isolates": len(pa), "plus_or_t": summary["published_counting"], "strict_plus": summary["strict_plus_only"]}, indent=2))
