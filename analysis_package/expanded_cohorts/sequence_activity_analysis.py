from pathlib import Path
from itertools import permutations
import argparse
import json
import re
import numpy as np
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
names = ["AP41", "SX1", "SX2", "S5", "G", "M1"]
accessions = pd.read_csv(root / "thai_pyocin_accessions.csv").set_index("pyocin")
disorder = pd.read_csv(root / "thai_pyocin_iupred3.csv").set_index("accession")
coverage = pd.read_csv(root / "thai_coverage.csv").set_index("pyocin")
fasta = {record.id: str(record.seq) for record in SeqIO.parse(root / "data" / "Laosena_2026_reference_pyocins.fasta", "fasta")}
features = []
for name in names:
    accession = accessions.loc[name, "accession"]
    sequence = fasta[accession]
    protein = ProteinAnalysis(sequence)
    features.append({"pyocin": name, "accession": accession, "length_aa": len(sequence), "disorder_fraction": disorder.loc[accession, "predicted_disorder_fraction"], "nterm40_disorder_fraction": disorder.loc[accession, "nterm40_disorder_fraction"], "charge_pH7": protein.charge_at_pH(7.0), "gravy": protein.gravy(), "aromaticity": protein.aromaticity(), "covered_isolates": int(coverage.loc[name, "all_positive"])})
feature_table = pd.DataFrame(features)
feature_table.to_csv(root / "thai_sequence_activity_features.csv", index=False)

def exact_spearman(x, y):
    observed = float(spearmanr(x, y).statistic)
    null = np.array([spearmanr(x, labels).statistic for labels in permutations(y)], dtype=float)
    return observed, float(np.mean(np.abs(null) >= abs(observed) - 1e-12))

predictors = ["length_aa", "disorder_fraction", "nterm40_disorder_fraction", "charge_pH7", "gravy", "aromaticity"]
association_rows = []
for predictor in predictors:
    rho, p = exact_spearman(feature_table[predictor].to_numpy(), feature_table.covered_isolates.to_numpy())
    association_rows.append({"predictor": predictor, "spearman_rho": rho, "exact_permutation_p": p})
association = pd.DataFrame(association_rows)
association["bh_q_six_features"] = multipletests(association.exact_permutation_p, method="fdr_bh")[1]
association.to_csv(root / "thai_sequence_activity_associations.csv", index=False)

matrix_text = (root / "data" / "thai_pyocin_clustalo_percent_identity.txt").read_text(encoding="utf-8")
identity = np.zeros((6, 6), dtype=float)
for line in matrix_text.splitlines():
    match = re.match(r"\s*([1-6]):\s+([^\s]+)\s+(.+)", line)
    if match:
        identity[int(match.group(1)) - 1] = [float(value) for value in match.group(3).split()]
if np.any(np.diag(identity) != 100) or not np.allclose(identity, identity.T):
    raise ValueError("Clustal Omega percent identity matrix did not parse correctly")
zones = pd.read_csv(root / "thai_isolate_pyocin_zones.csv")
responses = zones.pivot(index="isolate", columns="pyocin", values="all_positive")[names].astype(bool)
jaccard = np.eye(6)
for i in range(6):
    for j in range(i + 1, 6):
        overlap = (responses.iloc[:, i] & responses.iloc[:, j]).sum()
        union = (responses.iloc[:, i] | responses.iloc[:, j]).sum()
        jaccard[i, j] = jaccard[j, i] = overlap / union if union else 0.0
upper = np.triu_indices(6, 1)
observed = float(spearmanr(identity[upper], jaccard[upper]).statistic)
null = np.array([spearmanr(identity[upper], jaccard[np.ix_(perm, perm)][upper]).statistic for perm in permutations(range(6))], dtype=float)
mantel_p = float(np.mean(np.abs(null) >= abs(observed) - 1e-12))
pair_rows = []
for i, j in zip(*upper):
    pair_rows.append({"pyocin_1": names[i], "pyocin_2": names[j], "clustalo_identity_percent": identity[i, j], "response_jaccard": jaccard[i, j]})
pd.DataFrame(pair_rows).to_csv(root / "thai_sequence_response_pairs.csv", index=False)
summary = {"n_proteins": 6, "n_response_pairs": 15, "sequence_identity_response_overlap_spearman": observed, "exact_label_permutation_p": mantel_p, "n_permutations": 720, "length_disorder_fraction_spearman": float(spearmanr(feature_table.length_aa, feature_table.disorder_fraction).statistic), "feature_screen_q_method": "Benjamini-Hochberg across six sequence features"}
(root / "thai_sequence_activity_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(feature_table.to_string(index=False))
print(association.to_string(index=False))
print(summary)
