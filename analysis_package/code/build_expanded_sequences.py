from pathlib import Path
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

root = Path(__file__).resolve().parent
existing = {record.id: record for record in SeqIO.parse(root.parent / "outputs" / "analysis_package" / "data" / "retained_bacteriocins.fasta", "fasta")}
colicin_names = {
    "P04480": "Colicin_A",
    "P05819": "Colicin_B",
    "P17998": "Colicin_D",
    "P02978": "Colicin_E1",
    "P04419": "Colicin_E2",
    "P00646": "Colicin_E3",
    "Q47112": "Colicin_E7",
    "P08083": "Colicin_N"
}
pyocin_names = {
    "ERX71449.1": "Pyocin_L2",
    "ERZ09331.1": "Pyocin_L3",
    "KYO98147.1": "Pyocin_S6"
}
records = list(existing.values())
for record in SeqIO.parse(root / "expanded_colicins_uniprot.fasta", "fasta"):
    accession = record.id.split("|")[1]
    records.append(SeqRecord(record.seq, id=colicin_names[accession], description=f"accession={accession}"))
for record in SeqIO.parse(root / "expanded_pyocins_ncbi.fasta", "fasta"):
    accession = record.id
    records.append(SeqRecord(record.seq, id=pyocin_names[accession], description=f"accession={accession}"))
SeqIO.write(records, root / "expanded_bacteriocins.fasta", "fasta")
