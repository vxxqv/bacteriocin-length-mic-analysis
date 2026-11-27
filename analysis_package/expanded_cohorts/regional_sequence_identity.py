from pathlib import Path
import argparse
import pandas as pd
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
hits = pd.read_csv(root / "thai_pyocin_interproscan_hits.csv")
alignment = {record.id: str(record.seq) for record in SeqIO.parse(root / "data" / "thai_pyocin_clustalo_alignment.fasta", "fasta")}
sx1 = alignment["UYS84543.1"]
sx2 = alignment["UYS84545.1"]
if len(sx1) != len(sx2):
    raise ValueError("Aligned sequences have different lengths")
regions = [("translocation", "Pfam", "PF06958"), ("terminal DNase annotation", "Gene3D", "G3DSA:3.90.540.10")]
rows = []
for label, analysis, signature in regions:
    match = hits[(hits.pyocin == "SX1") & (hits.analysis == analysis) & (hits.signature == signature)]
    if len(match) != 1:
        raise ValueError((analysis, signature, len(match)))
    start, end = int(match.iloc[0].start), int(match.iloc[0].end)
    position = 0
    paired = 0
    identical = 0
    gaps_in_sx2 = 0
    for first, second in zip(sx1, sx2):
        if first != "-":
            position += 1
        if start <= position <= end and first != "-":
            if second == "-":
                gaps_in_sx2 += 1
            else:
                paired += 1
                identical += first == second
    if paired + gaps_in_sx2 != end - start + 1:
        raise ValueError(label)
    rows.append({"sx1_region": label, "analysis": analysis, "signature": signature, "sx1_start": start, "sx1_end": end, "region_length_aa": end - start + 1, "paired_alignment_columns": paired, "identical_paired_columns": identical, "sx2_gap_columns": gaps_in_sx2, "identity_percent_among_paired_columns": 100 * identical / paired})
results = pd.DataFrame(rows)
results.to_csv(root / "thai_sx1_sx2_regional_identity.csv", index=False)
print(results.to_string(index=False))
