from pathlib import Path
from itertools import combinations
import argparse
import numpy as np
import pandas as pd
from scipy.stats import binomtest
from statsmodels.stats.multitest import multipletests

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
mic = pd.read_csv(root / "thai_isolate_pyocin_mic.csv").pivot(index="isolate", columns="pyocin", values="mic_ug_ml")
zones = pd.read_csv(root / "thai_isolate_pyocin_zones.csv").pivot(index="isolate", columns="pyocin", values="all_positive")
rows = []
for first, second in combinations(["AP41", "SX1", "SX2"], 2):
    paired = mic[[first, second]].dropna()
    ratios = paired[second] / paired[first]
    higher = int((ratios > 1).sum())
    lower = int((ratios < 1).sum())
    ties = int((ratios == 1).sum())
    sign_result = binomtest(higher, higher + lower, 0.5) if higher + lower else None
    sign_interval = sign_result.proportion_ci(confidence_level=0.95, method="exact") if sign_result else None
    response_first_only = int(((zones[first] == 1) & (zones[second] == 0)).sum())
    response_second_only = int(((zones[first] == 0) & (zones[second] == 1)).sum())
    rows.append({"first": first, "second": second, "n_paired_numeric_mic": len(paired), "second_mic_higher": higher, "second_mic_lower": lower, "mic_ties": ties, "second_mic_higher_fraction_non_ties": higher / (higher + lower) if higher + lower else np.nan, "second_mic_higher_exact_ci_low": float(sign_interval.low) if sign_interval else np.nan, "second_mic_higher_exact_ci_high": float(sign_interval.high) if sign_interval else np.nan, "median_second_over_first_mic": float(np.median(ratios)), "mic_exact_sign_p": float(sign_result.pvalue) if sign_result else 1.0, "first_only_zone_response": response_first_only, "second_only_zone_response": response_second_only, "response_exact_mcnemar_p": float(binomtest(response_first_only, response_first_only + response_second_only, 0.5).pvalue) if response_first_only + response_second_only else 1.0})
results = pd.DataFrame(rows)
results["mic_holm_p_three_pairs"] = multipletests(results.mic_exact_sign_p, method="holm")[1]
results["response_holm_p_three_pairs"] = multipletests(results.response_exact_mcnemar_p, method="holm")[1]
results.to_csv(root / "thai_paired_mic_and_response.csv", index=False)
print(results.to_string(index=False))
