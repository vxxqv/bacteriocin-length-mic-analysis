from pathlib import Path
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
names = ["accession", "sequence_md5", "sequence_length_aa", "analysis", "signature", "signature_description", "start", "end", "score", "status", "date", "interpro_accession", "interpro_description", "go_terms", "pathways"]
raw = pd.read_csv(root / "data" / "thai_pyocin_interproscan.tsv", sep="\t", header=None, names=names, dtype=str, keep_default_na=False)
accessions = pd.read_csv(root / "thai_pyocin_accessions.csv").set_index("accession")
if set(raw.accession) != set(accessions.index):
    raise ValueError("InterProScan accessions do not match the six reference sequences")
raw["pyocin"] = raw.accession.map(accessions.pyocin)
raw["start"] = raw.start.astype(int)
raw["end"] = raw.end.astype(int)
raw["sequence_length_aa"] = raw.sequence_length_aa.astype(int)
if any(raw.end > raw.sequence_length_aa) or any(raw.start < 1):
    raise ValueError("InterProScan coordinates are out of range")
columns = ["pyocin", "accession", "analysis", "signature", "signature_description", "start", "end", "interpro_accession", "interpro_description", "score"]
raw[columns].to_csv(root / "thai_pyocin_interproscan_hits.csv", index=False)
summary = raw[raw.interpro_accession.ne("-")].groupby(["pyocin", "accession"], as_index=False).agg(n_interpro_terms=("interpro_accession", "nunique"), interpro_terms=("interpro_accession", lambda values: ";".join(sorted(set(values)))))
summary = accessions.reset_index()[["pyocin", "accession", "length_aa"]].merge(summary, on=["pyocin", "accession"], how="left")
summary["n_interpro_terms"] = summary.n_interpro_terms.fillna(0).astype(int)
summary["interpro_terms"] = summary.interpro_terms.fillna("")
summary.to_csv(root / "thai_pyocin_interproscan_summary.csv", index=False)
print(summary.to_string(index=False))
print(raw.groupby(["pyocin", "analysis"]).size().to_string())
