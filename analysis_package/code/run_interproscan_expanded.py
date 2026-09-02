from pathlib import Path
import json
import time
import requests

root = Path(__file__).resolve().parent
base = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5"
sequence = (root / "expanded_bacteriocins.fasta").read_text(encoding="utf-8")
payload = {
    "email": "vivaanpatni2010@gmail.com",
    "title": "Expanded bacteriocin domain annotation",
    "goterms": "true",
    "pathways": "true",
    "sequence": sequence
}
response = requests.post(f"{base}/run", data=payload, timeout=180)
response.raise_for_status()
job_id = response.text.strip()
while True:
    status = requests.get(f"{base}/status/{job_id}", timeout=60).text.strip()
    if status == "FINISHED":
        break
    if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
        raise RuntimeError(status)
    time.sleep(10)
types = requests.get(f"{base}/resulttypes/{job_id}", timeout=60).text
outputs = {
    "tsv": "expanded_interproscan_results.tsv",
    "xml": "expanded_interproscan_results.xml",
    "gff": "expanded_interproscan_results.gff3",
    "json": "expanded_interproscan_results.json",
    "sequence": "expanded_interproscan_input.fasta",
    "submission": "expanded_interproscan_submission.xml"
}
retrieved = {}
for result_type, filename in outputs.items():
    if result_type in types:
        item = requests.get(f"{base}/result/{job_id}/{result_type}", timeout=300)
        item.raise_for_status()
        (root / filename).write_bytes(item.content)
        retrieved[result_type] = filename
(root / "expanded_interproscan_provenance.json").write_text(json.dumps({"job_id": job_id, "service": "EMBL-EBI Job Dispatcher InterProScan 5", "retrieval_date": "2026-09-01", "outputs": retrieved}, indent=2), encoding="utf-8")
print(job_id)
