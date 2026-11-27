from pathlib import Path
import argparse
import json
import re
import time
import pandas as pd
import requests
from Bio import SeqIO

parser = argparse.ArgumentParser()
parser.add_argument("--output", type=Path, required=True)
args = parser.parse_args()
root = args.output.resolve()
raw_dir = root / "data" / "thai_pyocin_iupred3_raw"
raw_dir.mkdir(parents=True, exist_ok=True)
base = "https://iupred3.elte.hu"
rows = []
for record in SeqIO.parse(root / "data" / "Laosena_2026_reference_pyocins.fasta", "fasta"):
    path = raw_dir / f"{record.id}.json"
    if path.exists():
        data = json.loads(path.read_text(encoding="utf-8"))
    else:
        session = requests.Session()
        page = session.get(base + "/", timeout=60)
        page.raise_for_status()
        token = re.search(r'name="csrfmiddlewaretoken" value="([^"]+)"', page.text)
        if token is None:
            raise RuntimeError("IUPred3 form token was not found")
        response = session.post(base + "/plot", data={"csrfmiddlewaretoken": token.group(1), "accession": "", "inp_seq": str(record.seq), "context": "long", "smoothing": "savgol"}, headers={"Referer": base + "/"}, timeout=180)
        response.raise_for_status()
        raw_path = re.search(r'href="(/raw_json%3F[0-9]+)"', response.text)
        if raw_path is None:
            raise RuntimeError("IUPred3 raw result link was not found")
        result = session.get(base + raw_path.group(1), timeout=60)
        result.raise_for_status()
        data = result.json()
        path.write_text(json.dumps(data), encoding="utf-8")
        time.sleep(1)
    scores = data["iupred2"]
    if len(scores) != len(record.seq):
        raise ValueError(record.id)
    rows.append({"accession": record.id, "length_aa": len(record.seq), "predicted_disordered_residues": sum(score >= 0.5 for score in scores), "predicted_disorder_fraction": sum(score >= 0.5 for score in scores) / len(scores), "nterm40_disorder_fraction": sum(score >= 0.5 for score in scores[:40]) / min(len(scores), 40)})
    print(record.id, flush=True)
pd.DataFrame(rows).to_csv(root / "thai_pyocin_iupred3.csv", index=False)
