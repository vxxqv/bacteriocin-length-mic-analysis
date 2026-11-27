from pathlib import Path
from itertools import combinations
import json
import math
import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, entropy, chi2
from statsmodels.discrete.conditional_models import ConditionalLogit

root = Path(__file__).resolve().parent
mic = pd.read_csv(root / "data" / "expanded_mic_dataset.csv")
profiles = pd.read_csv(root / "spot_isolate_profiles.csv")
panel = profiles.loc[profiles.Species.eq("Pseudomonas aeruginosa")].copy()
if len(panel) != 150 or panel.Isolate.nunique() != 150:
    raise ValueError("The paired spot panel must contain 150 distinct isolates")
names = ["PyoL1", "PyoL2", "PyoL3"]
if not set(panel[names].stack()).issubset({"+", "T", "-"}):
    raise ValueError("Unexpected spot response code")

mic_part = mic.loc[mic.doi.eq("10.1002/mbo3.210") & mic.bacteriocin.isin(["Pyocin L2", "Pyocin L3"]) & mic.censoring.eq("exact")]
matched_strains = [strain for strain, part in mic_part.groupby("target_strain") if set(part.bacteriocin) == {"Pyocin L2", "Pyocin L3"}]
if len(matched_strains) != 4:
    raise ValueError("Expected four exact L2/L3 matched strains")
wide = mic_part.loc[mic_part.target_strain.isin(matched_strains)].pivot(index="target_strain", columns="bacteriocin", values="mic_nM_lower")
wide = wide[["Pyocin L2", "Pyocin L3"]].sort_index()
log_values = np.log10(wide.to_numpy(dtype=float))
grand = log_values.mean()
protein_mean = log_values.mean(axis=0)
strain_mean = log_values.mean(axis=1)
protein_component = np.tile(protein_mean - grand, (len(wide), 1))
strain_component = np.repeat((strain_mean - grand)[:, None], 2, axis=1)
interaction_component = log_values - grand - protein_component - strain_component
total_ss = float(np.sum((log_values - grand) ** 2))
ss_rows = [
    {"component": "protein", "sum_squares": float(np.sum(protein_component ** 2))},
    {"component": "strain", "sum_squares": float(np.sum(strain_component ** 2))},
    {"component": "protein_by_strain_plus_error", "sum_squares": float(np.sum(interaction_component ** 2))},
]
for row in ss_rows:
    row["share_of_total_percent"] = 100 * row["sum_squares"] / total_ss
pd.DataFrame(ss_rows).to_csv(root / "advanced_mic_decomposition.csv", index=False)
ratio = wide["Pyocin L3"] / wide["Pyocin L2"]
ratio_rows = []
for strain in wide.index:
    remaining = ratio.drop(strain)
    ratio_rows.append({"held_out_strain": strain, "remaining_strains": len(remaining), "geometric_mean_l3_to_l2": float(np.exp(np.log(remaining).mean())), "l3_lower_count": int((remaining < 1).sum()), "l2_lower_count": int((remaining > 1).sum())})
pd.DataFrame(ratio_rows).to_csv(root / "advanced_leave_one_strain_out.csv", index=False)

indicators = {
    "published_counting": panel[names].isin(["+", "T"]).to_numpy(dtype=int),
    "strict_plus_only": panel[names].eq("+").to_numpy(dtype=int),
}
model_rows = []
bootstrap_rows = []
association_rows = []
rng = np.random.default_rng(20261102)
bootstrap_indices = rng.integers(0, len(panel), size=(20000, len(panel)))
design = np.tile(np.array([[0, 0], [1, 0], [0, 1]], dtype=float), (len(panel), 1))
groups = np.repeat(np.arange(len(panel)), 3)
for definition, matrix in indicators.items():
    informative = (matrix.sum(axis=1) > 0) & (matrix.sum(axis=1) < 3)
    mask = np.repeat(informative, 3)
    fit = ConditionalLogit(matrix.reshape(-1)[mask], design[mask], groups=groups[mask]).fit(disp=False, maxiter=1000)
    interval = np.asarray(fit.conf_int())
    for index, protein in enumerate(names[1:]):
        model_rows.append({"definition": definition, "protein": protein, "reference": names[0], "informative_isolates": int(informative.sum()), "test": "coefficient_wald", "log_odds_coefficient": float(fit.params[index]), "likelihood_ratio_chisq": None, "odds_ratio": float(np.exp(fit.params[index])), "ci95_low": float(np.exp(interval[index, 0])), "ci95_high": float(np.exp(interval[index, 1])), "p_value": float(fit.pvalues[index])})
    null_log_likelihood = float(fit.model.loglike(np.zeros(2)))
    model_rows.append({"definition": definition, "protein": "joint_model", "reference": names[0], "informative_isolates": int(informative.sum()), "test": "likelihood_ratio_2df", "log_odds_coefficient": None, "likelihood_ratio_chisq": float(2 * (fit.llf - null_log_likelihood)), "odds_ratio": None, "ci95_low": None, "ci95_high": None, "p_value": float(chi2.sf(2 * (fit.llf - null_log_likelihood), 2))})
    for first_index, second_index in combinations(range(3), 2):
        first = matrix[:, first_index]
        second = matrix[:, second_index]
        observed_diff = float((second - first).mean())
        boot_diff = (second[bootstrap_indices] - first[bootstrap_indices]).mean(axis=1)
        low, high = np.quantile(boot_diff, [0.025, 0.975])
        bootstrap_rows.append({"definition": definition, "first": names[first_index], "second": names[second_index], "risk_difference_second_minus_first": observed_diff, "isolate_bootstrap_low": float(low), "isolate_bootstrap_high": float(high), "resamples": len(boot_diff), "seed": 20261102})
        both = int(np.sum((first == 1) & (second == 1)))
        first_only = int(np.sum((first == 1) & (second == 0)))
        second_only = int(np.sum((first == 0) & (second == 1)))
        neither = int(np.sum((first == 0) & (second == 0)))
        odds_ratio, p_value = fisher_exact([[both, first_only], [second_only, neither]], alternative="two-sided")
        denominator = math.sqrt((both + first_only) * (second_only + neither) * (both + second_only) * (first_only + neither))
        association_rows.append({"definition": definition, "first": names[first_index], "second": names[second_index], "observed_overlap": both, "independence_expected_overlap": float((both + first_only) * (both + second_only) / len(panel)), "phi": float((both * neither - first_only * second_only) / denominator), "fisher_odds_ratio": float(odds_ratio), "fisher_p": float(p_value)})
model_frame = pd.DataFrame(model_rows)
model_frame.to_csv(root / "advanced_conditional_logit.csv", index=False)
pd.DataFrame(bootstrap_rows).to_csv(root / "advanced_paired_bootstrap.csv", index=False)
association_frame = pd.DataFrame(association_rows)
for definition, part in association_frame.groupby("definition"):
    ordered = part.sort_values("fisher_p").index.tolist()
    running = 0.0
    for rank, index in enumerate(ordered):
        running = max(running, min(1.0, (len(ordered) - rank) * association_frame.loc[index, "fisher_p"]))
        association_frame.loc[index, "holm_fisher_p"] = running
association_frame.to_csv(root / "advanced_overlap_association.csv", index=False)

combination_rows = []
bound_rows = []
for size in range(1, 4):
    for indices in combinations(range(3), size):
        members = [names[index] for index in indices]
        lower = int(indicators["strict_plus_only"][:, indices].any(axis=1).sum())
        upper = int(indicators["published_counting"][:, indices].any(axis=1).sum())
        bound_rows.append({"proteins": "+".join(members), "number_of_proteins": size, "strict_plus_lower": lower, "plus_or_t_upper": upper, "uncertain_isolates": upper - lower})
        for definition, matrix in indicators.items():
            count = int(matrix[:, indices].any(axis=1).sum())
            best_single = max(int(matrix[:, index].sum()) for index in indices)
            combination_rows.append({"definition": definition, "proteins": "+".join(members), "number_of_proteins": size, "covered_isolates": count, "coverage_percent": 100 * count / len(panel), "gain_over_best_member": count - best_single})
pd.DataFrame(combination_rows).to_csv(root / "advanced_coverage_combinations.csv", index=False)
pd.DataFrame(bound_rows).to_csv(root / "advanced_code_bounds.csv", index=False)
pair_bound_rows = []
for first_index, second_index in combinations(range(3), 2):
    first_low = int(indicators["strict_plus_only"][:, first_index].sum())
    first_high = int(indicators["published_counting"][:, first_index].sum())
    second_low = int(indicators["strict_plus_only"][:, second_index].sum())
    second_high = int(indicators["published_counting"][:, second_index].sum())
    pair_bound_rows.append({"first": names[first_index], "second": names[second_index], "minimum_second_minus_first": second_low - first_high, "maximum_second_minus_first": second_high - first_low, "denominator": len(panel)})
pd.DataFrame(pair_bound_rows).to_csv(root / "advanced_pairwise_code_bounds.csv", index=False)
pattern_counts = panel[names].astype(str).agg("".join, axis=1).value_counts().sort_index()
diversity = {"observed_raw_signatures": len(pattern_counts), "possible_raw_signatures": 27, "shannon_bits": float(entropy(pattern_counts.to_numpy(), base=2)), "effective_number_of_signatures": float(2 ** entropy(pattern_counts.to_numpy(), base=2))}
(root / "advanced_diversity.json").write_text(json.dumps(diversity, indent=2), encoding="utf-8")
report = {"mic_sum_squares": ss_rows, "mic_l3_to_l2_ratios": {strain: float(value) for strain, value in ratio.items()}, "conditional_logit": model_rows, "diversity": diversity}
(root / "advanced_summary.json").write_text(json.dumps(report, indent=2, allow_nan=False), encoding="utf-8")
print(json.dumps({"mic_sum_squares": ss_rows, "conditional_logit": model_rows, "diversity": diversity}, indent=2, allow_nan=False))
