from pathlib import Path
import json
import math
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / "data" / "expanded_mic_dataset.csv"
FEATURES = ROOT / "data" / "expanded_protein_features.csv"
IDENTITY = ROOT / "data" / "expanded_clustalo_percent_identity_matrix.csv"
BLOCK = ["doi", "target_strain", "assay_group", "family"]

data = pd.read_csv(SOURCE)
features = pd.read_csv(FEATURES)
identity = pd.read_csv(IDENTITY, index_col=0)
data["block_id"] = data[BLOCK].astype(str).agg(" | ".join, axis=1)
data["log10_mic_nM"] = np.log10(data["mic_nM_lower"])
groups = data.groupby(BLOCK, sort=True)
matched = pd.concat([part for _, part in groups if part.bacteriocin.nunique() > 1], ignore_index=True)
exact = matched.loc[matched.censoring.eq("exact")].copy()
exact_groups = exact.groupby(BLOCK, sort=True)
exact_matched = pd.concat([part for _, part in exact_groups if part.bacteriocin.nunique() > 1], ignore_index=True)

coverage = []
for key, part in groups:
    if part.bacteriocin.nunique() < 2:
        continue
    coverage.append(dict(zip(BLOCK, key)) | {
        "records": len(part),
        "proteins": part.bacteriocin.nunique(),
        "exact": int(part.censoring.eq("exact").sum()),
        "interval": int(part.censoring.eq("interval").sum()),
        "right_censored": int(part.censoring.eq("right").sum()),
        "distinct_lengths": part.length_aa.nunique(),
        "protein_names": "; ".join(sorted(part.bacteriocin.unique()))
    })
coverage = pd.DataFrame(coverage)
coverage.to_csv(ROOT / "matched_block_coverage.csv", index=False)

rows = []
for key, part in exact_groups:
    if part.bacteriocin.nunique() < 2:
        continue
    x = part.length_aa.to_numpy(dtype=float)
    y = part.log10_mic_nM.to_numpy(dtype=float)
    slope = np.nan if np.var(x) == 0 else np.cov(x, y, bias=True)[0, 1] / np.var(x)
    rows.append(dict(zip(BLOCK, key)) | {
        "proteins": part.bacteriocin.nunique(),
        "length_range_aa": int(x.max() - x.min()),
        "slope_log10_nM_per_aa": slope,
        "log10_mic_span": float(y.max() - y.min()),
        "proteins_ordered_by_length": "; ".join(part.sort_values("length_aa").bacteriocin)
    })
slopes = pd.DataFrame(rows)
slopes.to_csv(ROOT / "matched_block_length_slopes.csv", index=False)

pair_rows = []
for strain, part in data.loc[data.doi.eq("10.1002/mbo3.210")].groupby("target_strain"):
    for a, b in [("Pyocin L1", "Pyocin L2"), ("Pyocin L2", "Pyocin L3")]:
        left = part.loc[part.bacteriocin.eq(a)]
        right = part.loc[part.bacteriocin.eq(b)]
        if len(left) != 1 or len(right) != 1:
            continue
        left = left.iloc[0]
        right = right.iloc[0]
        pair_rows.append({
            "target_strain": strain, "protein_a": a, "protein_b": b,
            "length_a_aa": int(left.length_aa), "length_b_aa": int(right.length_aa),
            "mic_a_nM_lower": left.mic_nM_lower, "mic_b_nM_lower": right.mic_nM_lower,
            "censoring_a": left.censoring, "censoring_b": right.censoring,
            "ratio_a_over_b_lower_or_exact": left.mic_nM_lower / right.mic_nM_lower if right.censoring == "exact" else np.nan,
            "ratio_b_over_a_lower_or_exact": right.mic_nM_lower / left.mic_nM_lower if left.censoring == "exact" else np.nan
        })
pairs = pd.DataFrame(pair_rows)
pairs.to_csv(ROOT / "matched_pyocin_pairs.csv", index=False)

span_rows = []
for (doi, protein), part in data.groupby(["doi", "bacteriocin"]):
    if doi != "10.1002/mbo3.210":
        continue
    exact_part = part.loc[part.censoring.eq("exact")]
    if len(exact_part) < 2:
        continue
    minimum = exact_part.mic_nM_lower.min()
    maximum = exact_part.mic_nM_lower.max()
    right = part.loc[part.censoring.eq("right")]
    if len(right) and right.mic_nM_lower.max() > maximum:
        maximum = right.mic_nM_lower.max()
        bound = "strict_lower_bound"
    else:
        bound = "exact_observed_span"
    span_rows.append({"doi": doi, "bacteriocin": protein, "n_target_strains": part.target_strain.nunique(),
                      "minimum_exact_mic_nM": minimum, "maximum_or_lower_bound_nM": maximum,
                      "fold_span": maximum / minimum, "span_type": bound})
spans = pd.DataFrame(span_rows).sort_values("fold_span", ascending=False)
spans.to_csv(ROOT / "within_protein_target_spans.csv", index=False)

key_studies = data.groupby("doi").agg(records=("bacteriocin", "size"), proteins=("bacteriocin", "nunique"), targets=("target_strain", "nunique")).reset_index()
key_studies.to_csv(ROOT / "study_coverage.csv", index=False)

exact_l23 = pairs.loc[pairs.protein_a.eq("Pyocin L2") & pairs.protein_b.eq("Pyocin L3") & pairs.censoring_a.eq("exact") & pairs.censoring_b.eq("exact")].copy()
exact_l23["fold"] = exact_l23.mic_a_nM_lower / exact_l23.mic_b_nM_lower
rank = slopes.loc[slopes.length_range_aa.gt(0)].copy()

summary = {
    "all_records": int(len(data)), "all_proteins": int(data.bacteriocin.nunique()),
    "all_studies": int(data.doi.nunique()), "matched_blocks": int(len(coverage)),
    "matched_records": int(len(matched)), "matched_proteins": int(matched.bacteriocin.nunique()),
    "matched_studies": int(matched.doi.nunique()), "exact_matched_blocks": int(len(slopes)),
    "exact_matched_records": int(len(exact_matched)),
    "length_varying_exact_blocks": int(rank.shape[0]),
    "positive_length_slopes": int(rank.slope_log10_nM_per_aa.gt(0).sum()),
    "negative_length_slopes": int(rank.slope_log10_nM_per_aa.lt(0).sum()),
    "l1_l2_sequence_identity_percent": float(identity.loc["Pyocin L1", "Pyocin L2"]),
    "l1_l2_length_aa": int(features.set_index("bacteriocin").loc["Pyocin L1", "length_aa"]),
    "bu007_l2_to_l1_mic_lower_bound_fold": float(3300 / 10.8),
    "cfpa22_l1_to_l2_mic_lower_bound_fold": float(3300 / 61.8),
    "l2_l3_exact_strains": int(len(exact_l23)),
    "l2_l3_fold_min": float(exact_l23.fold.min()),
    "l2_l3_fold_max": float(exact_l23.fold.max())
}
(ROOT / "matched_analysis_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
print(json.dumps(summary, indent=2))
