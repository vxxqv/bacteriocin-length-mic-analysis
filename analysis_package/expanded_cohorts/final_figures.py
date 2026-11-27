from pathlib import Path
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, Normalize
from matplotlib.patches import Rectangle, Circle, FancyArrowPatch, FancyBboxPatch

root = Path(__file__).resolve().parent
navy = "#17243a"
blue = "#216b99"
teal = "#008f87"
orange = "#d76b3a"
purple = "#7a60a7"
red = "#ba4b55"
grey = "#697789"
light = "#eef2f4"
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 12, "axes.labelsize": 10, "xtick.labelsize": 9, "ytick.labelsize": 9, "svg.fonttype": "none", "pdf.fonttype": 42, "axes.spines.top": False, "axes.spines.right": False})

def label(ax, letter, title):
    ax.text(0, 1.10, letter, transform=ax.transAxes, fontweight="bold", fontsize=15, color=navy, va="bottom")
    ax.text(0.065, 1.105, title, transform=ax.transAxes, fontweight="bold", fontsize=12, color=navy, va="bottom")

def save(fig, number):
    for suffix in ["png", "pdf", "svg"]:
        fig.savefig(root / f"Vivaan_Final_Figure_{number}.{suffix}", dpi=300, bbox_inches="tight", facecolor="white")
    plt.close(fig)

mic = pd.read_csv(root / "data" / "expanded_mic_dataset.csv")
mic = mic.loc[mic.doi.eq("10.1002/mbo3.210") & mic.bacteriocin.isin(["Pyocin L1", "Pyocin L2", "Pyocin L3"])].copy()
mic["short_strain"] = mic.target_strain.str.split().str[-1]
strains = ["Bu007", "CFPA13", "CFPA22", "CFPA120", "CFPA124", "CFPA101"]
proteins = ["Pyocin L1", "Pyocin L2", "Pyocin L3"]
lookup = {(r.short_strain, r.bacteriocin): r for r in mic.itertuples()}
fig = plt.figure(figsize=(14.5, 10.5), layout="constrained")
gs = fig.add_gridspec(2, 2, width_ratios=[1.24, 1], height_ratios=[1.1, 0.9], wspace=0.14, hspace=0.28)
ax = fig.add_subplot(gs[0, 0])
label(ax, "A", "Reported MICs within one study")
cmap = LinearSegmentedColormap.from_list("mic", ["#e0f4ec", "#77c7ba", "#237f9e", "#282f62"])
norm = Normalize(vmin=np.log10(5), vmax=np.log10(3300))
for i, strain in enumerate(strains):
    for j, protein in enumerate(proteins):
        record = lookup.get((strain, protein))
        if record is None:
            color, txt, hatch = light, "NR", None
        else:
            value = float(record.mic_nM_lower)
            color = cmap(norm(np.log10(value)))
            txt = ">3300" if record.censoring == "right" else f"{value:g}"
            hatch = "////" if record.censoring == "right" else None
        ax.add_patch(Rectangle((j - 0.47, i - 0.45), 0.94, 0.9, facecolor=color, edgecolor="white", lw=2, hatch=hatch))
        ax.text(j, i, txt, ha="center", va="center", fontsize=10, fontweight="bold", color="white" if record is not None and float(record.mic_nM_lower) >= 60 else navy, bbox={"facecolor": "#283461", "edgecolor": "none", "pad": 1.4} if record is not None and record.censoring == "right" else None)
ax.set_xlim(-0.5, 2.5)
ax.set_ylim(len(strains) - 0.5, -0.5)
ax.set_xticks(range(3), ["L1", "L2", "L3"])
ax.xaxis.tick_top()
ax.set_yticks(range(len(strains)), strains)
ax.tick_params(length=0)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.text(0, -0.10, "nM; hatched cells are right-censored; NR = not reported", transform=ax.transAxes, color=grey, fontsize=9, va="top")

ax = fig.add_subplot(gs[0, 1])
label(ax, "B", "L3 versus L2, exact paired MICs")
ratio_rows = [("CFPA120", 9.7 / 61.8), ("CFPA124", 13.8 / 61.8), ("CFPA13", 6.9 / 13.8), ("CFPA22", 77.8 / 61.8)]
for i, (name, ratio) in enumerate(ratio_rows):
    ax.plot([0, np.log2(ratio)], [i, i], lw=4, color=teal if ratio < 1 else orange, solid_capstyle="round")
    ax.scatter(np.log2(ratio), i, s=100, color=teal if ratio < 1 else orange, zorder=3)
    ax.text(np.log2(ratio) + (0.13 if ratio > 1 else -0.13), i, f"{ratio:.2f}×", ha="left" if ratio > 1 else "right", va="center", fontweight="bold", color=navy)
ax.axvline(0, color=navy, lw=1.2)
ax.set_yticks(range(4), [x[0] for x in ratio_rows])
ax.set_ylim(3.6, -0.6)
ax.set_xlim(-3.7, 1.6)
ax.set_xticks([-3, -2, -1, 0, 1], ["⅛", "¼", "½", "1", "2"])
ax.set_xlabel("L3 MIC / L2 MIC (log2 scale)")
ax.grid(axis="x", alpha=0.15)

ss = pd.read_csv(root / "advanced_mic_decomposition.csv")
ax = fig.add_subplot(gs[1, 0])
label(ax, "C", "Variation across four matched strains")
parts = [("Target strain", ss.iloc[1].share_of_total_percent, blue), ("Protein", ss.iloc[0].share_of_total_percent, teal), ("Interaction + assay error", ss.iloc[2].share_of_total_percent, purple)]
left = 0
for name, value, color in parts:
    ax.barh(0, value, left=left, color=color, height=0.45)
    ax.text(left + value / 2, 0, f"{value:.1f}%", color="white", ha="center", va="center", fontweight="bold", fontsize=11)
    left += value
ax.set_xlim(0, 100)
ax.set_ylim(-0.75, 0.75)
ax.set_yticks([])
ax.set_xlabel("Share of observed log10 MIC sum of squares")
ax.legend([Rectangle((0, 0), 1, 1, facecolor=c) for _, _, c in parts], [n for n, _, _ in parts], loc="lower center", bbox_to_anchor=(0.5, -0.48), ncol=3, frameon=False, fontsize=9)
ax.text(0, 0.69, "One MIC per cell: interaction and measurement error cannot be separated.", transform=ax.transAxes, color=grey, fontsize=9)

blocks = pd.read_csv(root / "matched_block_length_slopes.csv")
ax = fig.add_subplot(gs[1, 1])
label(ax, "D", "Matched length slopes disagree")
slopes = blocks.slope_log10_nM_per_aa.to_numpy(float)
names = [f"{row.family}: {row.target_strain.split()[-1]}" for row in blocks.itertuples()]
order = np.argsort(slopes)
colors = [teal if slopes[k] < 0 else orange for k in order]
ax.barh(range(len(order)), slopes[order], color=colors, height=0.66)
ax.axvline(0, color=navy, lw=1)
ax.set_yticks(range(len(order)), [names[k] for k in order])
ax.set_xlabel("Within-block slope, log10 nM per residue")
ax.grid(axis="x", alpha=0.15)
save(fig, 1)

profile = pd.read_csv(root / "spot_isolate_profiles.csv")
panel = profile.loc[profile.Species.eq("Pseudomonas aeruginosa")].copy()
pattern = pd.read_csv(root / "spot_response_patterns.csv").sort_values("isolates", ascending=False)
boot = pd.read_csv(root / "advanced_paired_bootstrap.csv")
model = pd.read_csv(root / "advanced_conditional_logit.csv")
union = pd.read_csv(root / "advanced_code_bounds.csv")
fig = plt.figure(figsize=(14.5, 12.2), layout="constrained")
gs = fig.add_gridspec(2, 2, width_ratios=[1.0, 1.05], height_ratios=[1.17, 1], wspace=0.17, hspace=0.25)
ax = fig.add_subplot(gs[0, 0])
label(ax, "A", "Fourteen observed response signatures")
code_colors = {"+": teal, "T": orange, "-": "#cbd4dc"}
for i, row in enumerate(pattern.itertuples()):
    ax.barh(i, row.isolates, height=0.74, color=blue if "T" not in row.raw_pattern else orange, alpha=0.85)
    ax.text(row.isolates + 0.4, i, str(row.isolates), va="center", color=navy, fontsize=9)
ax.set_xlim(0, 61)
ax.set_ylim(len(pattern) - 0.45, -0.55)
ax.set_xticks([])
ax.set_yticks(range(len(pattern)), pattern.raw_pattern.tolist(), fontfamily="DejaVu Sans Mono", fontsize=10)
for spine in ax.spines.values():
    spine.set_visible(False)
ax.text(0, -0.06, "+ explicit response     T undefined     – no response", transform=ax.transAxes, fontsize=9, color=grey, va="top")

ax = fig.add_subplot(gs[0, 1])
label(ax, "B", "Coverage depends on the T code")
strict = panel[["PyoL1", "PyoL2", "PyoL3"]].eq("+").sum().to_numpy()
published = panel[["PyoL1", "PyoL2", "PyoL3"]].isin(["+", "T"]).sum().to_numpy()
strict = np.append(strict, panel[["PyoL1", "PyoL2", "PyoL3"]].eq("+").any(axis=1).sum())
published = np.append(published, panel[["PyoL1", "PyoL2", "PyoL3"]].isin(["+", "T"]).any(axis=1).sum())
yy = np.arange(4)
ax.barh(yy, published, height=0.66, color="#e8a47c", label="Additional if T counts")
ax.barh(yy, strict, height=0.66, color=teal, label="Explicit +")
for i, (low, high) in enumerate(zip(strict, published)):
    ax.text(high + 1.5, i, f"{low}–{high}", va="center", fontweight="bold", color=navy)
ax.set_yticks(yy, ["L1", "L2", "L3", "Any of three"])
ax.set_ylim(3.5, -0.5)
ax.set_xlim(0, 118)
ax.set_xlabel("Isolates out of 150")
ax.grid(axis="x", alpha=0.15)

ax = fig.add_subplot(gs[1, 0])
label(ax, "C", "Paired L3 minus L2 coverage")
subset = boot.loc[(boot.first.eq("PyoL2")) & (boot.second.eq("PyoL3"))]
for i, definition in enumerate(["published_counting", "strict_plus_only"]):
    row = subset.loc[subset.definition.eq(definition)].iloc[0]
    value = 100 * row.risk_difference_second_minus_first
    lo = 100 * row.isolate_bootstrap_low
    hi = 100 * row.isolate_bootstrap_high
    ax.plot([lo, hi], [i, i], lw=4, color=blue if i == 0 else orange, solid_capstyle="round")
    ax.scatter(value, i, s=110, color=blue if i == 0 else orange, zorder=3)
    ax.text(hi + 1.5, i, f"{value:+.1f} points", va="center", color=navy, fontweight="bold")
ax.axvline(0, color=navy, lw=1)
ax.set_yticks([0, 1], ["+ or T", "+ only"])
ax.set_xlim(-14, 45)
ax.set_ylim(1.7, -0.7)
ax.set_xlabel("Difference in percentage points; paired isolate bootstrap interval")
ax.grid(axis="x", alpha=0.15)

ax = fig.add_subplot(gs[1, 1])
label(ax, "D", "Within-isolate model, L3 versus L1")
rows = model.loc[model.protein.eq("PyoL3")]
for i, definition in enumerate(["published_counting", "strict_plus_only"]):
    row = rows.loc[rows.definition.eq(definition)].iloc[0]
    color = blue if i == 0 else orange
    ax.plot([row.ci95_low, row.ci95_high], [i, i], color=color, lw=4, solid_capstyle="round")
    ax.scatter(row.odds_ratio, i, s=110, color=color, zorder=3)
    ax.text(row.ci95_high + 0.15, i, f"OR {row.odds_ratio:.2f}", va="center", fontweight="bold", color=navy)
ax.axvline(1, color=navy, lw=1)
ax.set_yticks([0, 1], ["+ or T", "+ only"])
ax.set_ylim(1.7, -0.7)
ax.set_xlim(0, 8.4)
ax.set_xlabel("Conditional odds ratio, 95% Wald interval")
ax.grid(axis="x", alpha=0.15)
save(fig, 2)

identity = pd.read_csv(root / "data" / "expanded_clustalo_percent_identity_matrix.csv", index_col=0)
features = pd.read_csv(root / "data" / "expanded_protein_features.csv").set_index("bacteriocin")
disorder = pd.read_csv(root / "data" / "iupred3_l_type_profiles.csv")
fig = plt.figure(figsize=(14.5, 10.0), layout="constrained")
gs = fig.add_gridspec(2, 2, height_ratios=[0.95, 1.18], width_ratios=[0.7, 1.3], hspace=0.22, wspace=0.16)
ax = fig.add_subplot(gs[0, 0])
label(ax, "A", "Sequence identity")
short = ["L1", "L2", "L3"]
ids = ["Pyocin L1", "Pyocin L2", "Pyocin L3"]
matrix = identity.loc[ids, ids].to_numpy(float)
cmap2 = LinearSegmentedColormap.from_list("identity", ["#ecf0f3", "#9ed4d0", "#178a91", "#1c4266"])
ax.imshow(matrix, cmap=cmap2, vmin=0, vmax=100)
ax.set_xticks(range(3), short)
ax.set_yticks(range(3), short)
for i in range(3):
    for j in range(3):
        ax.text(j, i, f"{matrix[i, j]:.1f}%", ha="center", va="center", color="white" if matrix[i, j] > 60 else navy, fontweight="bold", fontsize=11)
ax.tick_params(length=0)

ax = fig.add_subplot(gs[0, 1])
label(ax, "B", "Length, SMART regions, and IUPred3")
for i, protein in enumerate(ids):
    row = features.loc[protein]
    length = int(row.length_aa)
    y = 2.6 - i * 1.2
    ax.add_patch(Rectangle((0, y - 0.13), length, 0.26, color="#d4dce3", ec="none"))
    if pd.notna(row.domain_coordinates):
        for seg in str(row.domain_coordinates).replace(";", ",").split(","):
            if "-" in seg:
                try:
                    start, end = [int(x.strip()) for x in seg.split("-", 1)]
                    ax.add_patch(Rectangle((start, y - 0.13), end - start + 1, 0.26, color=orange, ec="white", lw=0.5))
                except ValueError:
                    pass
    ax.text(-10, y, short[i], ha="right", va="center", fontweight="bold", color=navy)
    ax.text(length + 8, y, f"{length} aa", va="center", color=navy, fontsize=9)
    part = disorder.loc[disorder.protein.eq(short[i])]
    ax.plot(part.position, y - 0.34 + part.score * 0.33, color=purple, lw=1.3)
ax.set_xlim(-25, 340)
ax.set_ylim(-0.5, 3.3)
ax.set_xlabel("Residue position")
ax.set_yticks([])
ax.grid(axis="x", alpha=0.15)

ax = fig.add_subplot(gs[1, :])
label(ax, "C", "L-type pyocin mechanism reported in the literature")
ax.set_xlim(0, 14)
ax.set_ylim(0, 5)
ax.axis("off")
stage_x = [0.35, 5.05, 9.75]
stage_names = ["Surface recognition", "BamA engagement", "Envelope disruption"]
for k, x in enumerate(stage_x):
    ax.add_patch(FancyBboxPatch((x, 0.48), 3.85, 3.68, boxstyle="round,pad=0.12,rounding_size=0.13", ec="#d8e0e5", fc="#f8fafb", lw=1.2))
    ax.text(x + 0.17, 3.88, stage_names[k], color=navy, fontsize=11, fontweight="bold", va="top")
    ax.add_patch(Rectangle((x + 0.24, 1.55), 3.37, 0.46, facecolor="#75b6b3", edgecolor="#3d888c"))
    ax.add_patch(FancyBboxPatch((x + 1.84, 1.22), 0.52, 1.12, boxstyle="round,pad=0.02", fc=blue, ec=navy, lw=1.2))
    ax.add_patch(FancyBboxPatch((x + 2.02, 1.33), 0.16, 0.9, boxstyle="round,pad=0.01", fc="white", ec="none"))
    ax.text(x + 2.10, 1.02, "BamA", ha="center", va="top", fontsize=9, color=navy)
    ax.add_patch(Circle((x + 1.05, 2.66), 0.24, fc=purple, ec="none"))
    ax.add_patch(Circle((x + 1.50, 2.66), 0.24, fc=purple, ec="none"))
    ax.add_patch(FancyArrowPatch((x + 1.65, 2.58), (x + 2.10, 2.29 if k == 0 else 1.79), arrowstyle="-", color=orange, lw=5, mutation_scale=10))
    ax.text(x + 0.98, 3.10, "Pyocin", ha="center", fontsize=9, color=navy)
    if k == 2:
        ax.add_patch(FancyArrowPatch((x + 2.68, 1.35), (x + 2.68, 2.07), arrowstyle="-|>", mutation_scale=15, color=red, lw=2))
        ax.text(x + 2.68, 2.36, "OMP insertion\nreduced", ha="center", fontsize=9, color=red)
for x in [4.38, 9.08]:
    ax.add_patch(FancyArrowPatch((x, 2.27), (x + 0.55, 2.27), arrowstyle="-|>", mutation_scale=18, color=grey, lw=2))
save(fig, 3)
