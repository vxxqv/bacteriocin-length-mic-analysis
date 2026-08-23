from pathlib import Path
import re
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import seaborn as sns
from Bio import AlignIO, Phylo
from Bio.Align import MultipleSeqAlignment
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

root = Path(__file__).resolve().parent
table_dir = root / "tables"
figure_dir = root / "figures"
table_dir.mkdir(exist_ok=True)
figure_dir.mkdir(exist_ok=True)

retained = [
    "Pyocin_S5",
    "PaeM4",
    "LlpA_BW11M1",
    "Pyocin_L1",
    "Syringacin_M",
    "KpneA",
    "KaerA",
    "KvarIa",
    "KvarM",
    "KpneM",
    "KpneM2"
]

alignment = AlignIO.read(root / "clustalo_retained_alignment.fasta", "fasta")
selected = list(alignment)
keep_columns = [index for index in range(alignment.get_alignment_length()) if any(record.seq[index] != "-" for record in selected)]
trimmed_records = []
for record in selected:
    sequence = "".join(record.seq[index] for index in keep_columns)
    trimmed_records.append(SeqRecord(Seq(sequence), id=record.id, description=""))
trimmed = MultipleSeqAlignment(trimmed_records)
AlignIO.write(trimmed, root / "clustalo_retained_alignment.fasta", "fasta")
AlignIO.write(trimmed, root / "clustalo_retained_alignment.clustal", "clustal")

names = [record.id for record in trimmed]
matrix = np.zeros((len(names), len(names)))
for i, first in enumerate(trimmed):
    for j, second in enumerate(trimmed):
        comparable = [(a, b) for a, b in zip(str(first.seq), str(second.seq)) if a != "-" and b != "-"]
        matrix[i, j] = 100 if i == j else 100 * sum(a == b for a, b in comparable) / len(comparable)
pim = pd.DataFrame(matrix, index=names, columns=names)
pim.to_csv(table_dir / "clustalo_retained_percent_identity_matrix.csv")

calculator = DistanceCalculator("identity")
distance_matrix = calculator.get_distance(trimmed)
tree = DistanceTreeConstructor().nj(distance_matrix)
for clade in tree.get_nonterminals():
    clade.name = None
Phylo.write(tree, root / "clustalo_retained_nj_tree.newick", "newick")

sns.set_theme(style="ticks", context="paper", font_scale=1.0)
fig, ax = plt.subplots(figsize=(7.5, 6.3))
sns.heatmap(pim, cmap="viridis", vmin=0, vmax=100, square=True, cbar_kws={"label": "Percent identity"}, ax=ax)
ax.set_xlabel("")
ax.set_ylabel("")
ax.tick_params(axis="x", labelrotation=55, labelsize=8)
ax.tick_params(axis="y", labelsize=8)
fig.tight_layout()
fig.savefig(figure_dir / "Figure_3_clustalo_identity_heatmap.png", bbox_inches="tight", dpi=300)
plt.close(fig)

fig, ax = plt.subplots(figsize=(7.5, 5.2))
Phylo.draw(tree, axes=ax, do_show=False, show_confidence=False)
ax.set_xlabel("Neighbor-joining distance from retained Clustal Omega alignment")
ax.set_ylabel("")
fig.tight_layout()
fig.savefig(figure_dir / "Figure_4_clustalo_retained_tree.png", bbox_inches="tight", dpi=300)
plt.close(fig)

upper = matrix[np.triu_indices(len(names), 1)]
summary = pd.DataFrame([{
    "n_proteins": len(names),
    "n_pairs": len(upper),
    "median_percent_identity": np.median(upper),
    "q1_percent_identity": np.quantile(upper, 0.25),
    "q3_percent_identity": np.quantile(upper, 0.75),
    "minimum_percent_identity": np.min(upper),
    "maximum_percent_identity": np.max(upper)
}])
summary.to_csv(table_dir / "clustalo_retained_identity_summary.csv", index=False)
