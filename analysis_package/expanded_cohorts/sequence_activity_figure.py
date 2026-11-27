from pathlib import Path
import argparse
import json
import re
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
names = ["AP41", "SX1", "SX2", "S5", "G", "M1"]
features = pd.read_csv(root / "thai_sequence_activity_features.csv").set_index("pyocin").loc[names]
pairs = pd.read_csv(root / "thai_sequence_response_pairs.csv")
disorder = np.full((6, int(features.length_aa.max())), np.nan)
for i, accession in enumerate(features.accession):
    scores = json.loads((root / "data" / "thai_pyocin_iupred3_raw" / f"{accession}.json").read_text(encoding="utf-8"))["iupred2"]
    disorder[i, :len(scores)] = scores
identity = np.eye(6) * 100
overlap = np.eye(6)
for row in pairs.itertuples():
    i = names.index(row.pyocin_1)
    j = names.index(row.pyocin_2)
    identity[i, j] = identity[j, i] = row.clustalo_identity_percent
    overlap[i, j] = overlap[j, i] = row.response_jaccard

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10, "axes.titlesize": 14, "axes.labelsize": 11, "xtick.labelsize": 9, "ytick.labelsize": 10, "savefig.facecolor": "white"})
fig = plt.figure(figsize=(15, 10.2), facecolor="white")
grid = fig.add_gridspec(2, 2, height_ratios=[1, 1.05], width_ratios=[1.22, 1], left=.075, right=.96, bottom=.095, top=.92, wspace=.24, hspace=.34)
ax_a = fig.add_subplot(grid[0, 0])
ax_b = fig.add_subplot(grid[0, 1])
ax_c = fig.add_subplot(grid[1, 0])
ax_d = fig.add_subplot(grid[1, 1])
palette = LinearSegmentedColormap.from_list("disorder", ["#f1f5f9", "#75c4c9", "#007d91", "#1d365f"])
heat = ax_a.imshow(np.ma.masked_invalid(disorder), aspect="auto", interpolation="nearest", cmap=palette, vmin=0, vmax=1)
ax_a.set_yticks(range(6), [f"{name}  ({features.loc[name, 'length_aa']} aa)" for name in names])
ax_a.set_xlabel("Residue position")
ax_a.set_title("A  Predicted disorder along each protein", loc="left", fontweight="bold", color="#17243a")
ax_a.set_xticks([0, 200, 400, 600, 776])
ax_a.spines[:].set_visible(False)
ax_a.tick_params(length=0)
cb = fig.colorbar(heat, ax=ax_a, fraction=.025, pad=.025)
cb.set_label("IUPred3 score")

ax_b.scatter(features.disorder_fraction, features.covered_isolates, s=185, color="#008a9b", edgecolor="white", linewidth=1.5, zorder=3)
offsets = {"AP41": (.008, -1.7), "SX1": (-.01, 1.5), "SX2": (.006, .4), "S5": (.006, .4), "G": (.006, .4), "M1": (.006, .4)}
for name in names:
    x = features.loc[name, "disorder_fraction"]
    y = features.loc[name, "covered_isolates"]
    dx, dy = offsets[name]
    ax_b.text(x + dx, y + dy, name, fontsize=10, va="center", ha="right" if name == "SX1" else "left", color="#17243a")
ax_b.set_xlim(.025, .43)
ax_b.set_ylim(2, 31)
ax_b.set_xlabel("Predicted disordered residue fraction")
ax_b.set_ylabel("Responding isolates out of 29")
ax_b.set_title("B  Disorder and response breadth", loc="left", fontweight="bold", color="#17243a")
ax_b.grid(color="#e3e8ee", linewidth=.8)
ax_b.set_axisbelow(True)
for spine in ["top", "right"]:
    ax_b.spines[spine].set_visible(False)

seq_cmap = LinearSegmentedColormap.from_list("similarity", ["#f3f7f9", "#c8e8ea", "#1d9eaa", "#173a66"])
im_c = ax_c.imshow(identity, vmin=0, vmax=100, cmap=seq_cmap, aspect="equal")
ax_c.set_title("C  Clustal Omega identity (%)", loc="left", fontweight="bold", color="#17243a")
ax_d.set_title("D  Shared response sets (%)", loc="left", fontweight="bold", color="#17243a")
im_d = ax_d.imshow(overlap * 100, vmin=0, vmax=100, cmap=seq_cmap, aspect="equal")
for ax, matrix in [(ax_c, identity), (ax_d, overlap * 100)]:
    ax.set_xticks(range(6), names, rotation=45, ha="right")
    ax.set_yticks(range(6), names)
    ax.tick_params(length=0)
    ax.spines[:].set_visible(False)
    for i in range(6):
        for j in range(6):
            value = matrix[i, j]
            ax.text(j, i, f"{value:.0f}", ha="center", va="center", fontsize=8.8, color="white" if value > 52 else "#203149")
fig.text(.075, .035, "Six proteins from one published 29-isolate panel. Disorder is predicted, and response overlap is descriptive.", fontsize=9.5, color="#5a6879")
for suffix in ["png", "pdf", "svg"]:
    fig.savefig(root / f"Vivaan_Sequence_Analysis_Figure.{suffix}", dpi=300)
plt.close(fig)
