from pathlib import Path
import json
import re
import time
import pandas as pd
import requests
from Bio import SeqIO

root = Path(__file__).resolve().parent
result_dir = root / "iupred3_raw"
result_dir.mkdir(exist_ok=True)
base = "https://iupred3.elte.hu"
rows = []

def segments(scores, threshold=0.5):
    found = []
    start = None
    for index, score in enumerate(scores, start=1):
        if score >= threshold and start is None:
            start = index
        if score < threshold and start is not None:
            found.append((start, index - 1))
            start = None
    if start is not None:
        found.append((start, len(scores)))
    return found

for record in SeqIO.parse(root / "retained_bacteriocins.fasta", "fasta"):
    output = result_dir / f"{record.id}.json"
    if output.exists():
        data = json.loads(output.read_text(encoding="utf-8"))
    else:
        session = requests.Session()
        home = session.get(f"{base}/", timeout=60)
        home.raise_for_status()
        token = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', home.text).group(1)
        payload = {
            "csrfmiddlewaretoken": token,
            "accession": "",
            "inp_seq": str(record.seq),
            "context": "long",
            "smoothing": "savgol"
        }
        result = session.post(f"{base}/plot", data=payload, headers={"Referer": f"{base}/"}, timeout=180)
        result.raise_for_status()
        raw_path = re.search(r'href="(/raw_json%3F[0-9]+)"', result.text).group(1)
        raw = session.get(f"{base}{raw_path}", timeout=60)
        raw.raise_for_status()
        data = raw.json()
        output.write_text(json.dumps(data, indent=2), encoding="utf-8")
        time.sleep(0.5)
    scores = data["iupred2"]
    found = segments(scores)
    rows.append({
        "bacteriocin": record.id.replace("_", " "),
        "sequence_length_aa": len(scores),
        "idr_length_aa": sum(score >= 0.5 for score in scores),
        "idr_fraction": sum(score >= 0.5 for score in scores) / len(scores),
        "longest_idr_aa": max((end - start + 1 for start, end in found), default=0),
        "idr_segments": "; ".join(f"{start}-{end}" for start, end in found)
    })

pd.DataFrame(rows).to_csv(root / "retained_iupred3_summary.csv", index=False)
