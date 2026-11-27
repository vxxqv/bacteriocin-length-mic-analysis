from pathlib import Path
import argparse
import json
import time
import requests

parser = argparse.ArgumentParser()
parser.add_argument("--email", required=True)
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
root.mkdir(parents=True, exist_ok=True)
sequence = (root / "data" / "Laosena_2026_reference_pyocins.fasta").read_text(encoding="utf-8")
base = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
job_file = root / "data" / "thai_pyocin_interproscan_job.json"
if job_file.exists():
    job = json.loads(job_file.read_text(encoding="utf-8"))
    job_id = job["job_id"]
else:
    response = requests.post(base + "/run", data={"email": args.email, "title": "Six published pyocin sequences", "goterms": "true", "pathways": "false", "stype": "p", "sequence": sequence}, timeout=180)
    response.raise_for_status()
    job_id = response.text.strip()
    job_file.write_text(json.dumps({"job_id": job_id, "service": base, "source": "Laosena et al. 2026 Table S2"}, indent=2), encoding="utf-8")
for _ in range(180):
    status_response = requests.get(base + "/status/" + job_id, timeout=60)
    status_response.raise_for_status()
    status = status_response.text.strip()
    print(job_id, status, flush=True)
    if status == "FINISHED":
        break
    if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
        raise RuntimeError(status)
    time.sleep(10)
else:
    raise TimeoutError(job_id)
result_types = requests.get(base + "/resulttypes/" + job_id, timeout=60)
result_types.raise_for_status()
retrieved = {}
for result_type, filename in {"tsv": "thai_pyocin_interproscan.tsv", "json": "thai_pyocin_interproscan.json", "gff": "thai_pyocin_interproscan.gff3"}.items():
    if result_type in result_types.text:
        response = requests.get(base + "/result/" + job_id + "/" + result_type, timeout=300)
        response.raise_for_status()
        (root / "data" / filename).write_bytes(response.content)
        retrieved[result_type] = filename
job_file.write_text(json.dumps({"job_id": job_id, "service": base, "source": "Laosena et al. 2026 Table S2", "outputs": retrieved}, indent=2), encoding="utf-8")
print(retrieved)
