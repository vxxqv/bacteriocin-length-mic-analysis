from pathlib import Path
import json
import time
import requests

root = Path(__file__).resolve().parent
base = "https://www.ebi.ac.uk/Tools/services/rest/clustalo"
sequence = (root / "expanded_bacteriocins.fasta").read_text(encoding="utf-8")
payload = {
    "email": "vivaanpatni2010@gmail.com",
    "title": "Expanded bacteriocin full-length alignment",
    "guidetreeout": "true",
    "dismatout": "true",
    "mbed": "false",
    "mbediteration": "false",
    "stype": "protein",
    "outfmt": "clustal_num",
    "order": "aligned",
    "sequence": sequence
}
response = requests.post(f"{base}/run", data=payload, timeout=180)
response.raise_for_status()
job_id = response.text.strip()
(root / "expanded_clustalo_job_id.txt").write_text(job_id, encoding="utf-8")
while True:
    status = requests.get(f"{base}/status/{job_id}", timeout=60).text.strip()
    if status == "FINISHED":
        break
    if status in {"ERROR", "FAILURE", "NOT_FOUND"}:
        details = requests.get(f"{base}/result/{job_id}/out", timeout=60).text
        raise RuntimeError(f"{status}: {details}")
    time.sleep(5)
types = requests.get(f"{base}/resulttypes/{job_id}", timeout=60).text
outputs = {
    "aln-clustal_num": "expanded_clustalo_alignment.clustal",
    "fa": "expanded_clustalo_alignment.fasta",
    "tree": "expanded_clustalo_guide_tree.dnd",
    "phylotree": "expanded_clustalo_tree.newick",
    "pim": "expanded_clustalo_percent_identity.txt",
    "matrix": "expanded_clustalo_distance_matrix.txt",
    "submission": "expanded_clustalo_submission.xml"
}
retrieved = {}
for result_type, filename in outputs.items():
    if result_type in types:
        item = requests.get(f"{base}/result/{job_id}/{result_type}", timeout=180)
        item.raise_for_status()
        (root / filename).write_bytes(item.content)
        retrieved[result_type] = filename
(root / "expanded_clustalo_provenance.json").write_text(json.dumps({"job_id": job_id, "service": "EMBL-EBI Job Dispatcher Clustal Omega 1.2.4", "retrieval_date": "2026-09-01", "outputs": retrieved}, indent=2), encoding="utf-8")
print(job_id)
