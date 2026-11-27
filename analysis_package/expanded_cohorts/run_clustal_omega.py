from pathlib import Path
import argparse
import json
import time
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--email", required=True)
args = parser.parse_args()
out = Path(__file__).resolve().parent
base = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
fasta = (out / "data" / "Laosena_2026_reference_pyocins.fasta").read_text(encoding="utf-8")
response = requests.post(base + "/run", data={"email": args.email, "title": "Pyocin six-protein alignment", "stype": "protein", "outfmt": "clustal_num", "sequence": fasta}, timeout=60)
response.raise_for_status()
job_id = response.text.strip()
for _ in range(90):
    status = requests.get(base + "/status/" + job_id, timeout=30).text.strip()
    if status == "FINISHED":
        break
    if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
        raise RuntimeError(status)
    time.sleep(5)
else:
    raise TimeoutError(job_id)
types = requests.get(base + "/resulttypes/" + job_id, timeout=30).text
result_names = {"aln-clustal_num": "thai_pyocin_clustalo_alignment.clustal", "fa": "thai_pyocin_clustalo_alignment.fasta", "pim": "thai_pyocin_clustalo_percent_identity.txt"}
retrieved = {}
for kind, name in result_names.items():
    if kind in types:
        result = requests.get(base + "/result/" + job_id + "/" + kind, timeout=60)
        result.raise_for_status()
        (out / "data" / name).write_bytes(result.content)
        retrieved[kind] = name
(out / "data" / "thai_pyocin_clustal_provenance.json").write_text(json.dumps({"job_id": job_id, "service": base, "outputs": retrieved}, indent=2), encoding="utf-8")
print(job_id, retrieved)
