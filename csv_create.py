import pandas as pd
from Bio import SeqIO

thermo_file = "thermo_cleaned.fasta"
meso_file = "meso_cleaned.fasta"

output_csv = "final_dataset.csv"

amino_acids = list("ACDEFGHIKLMNPQRSTVWY")

data = []

def extract_features(record, label):
    seq = str(record.seq).upper()
    length = len(seq)

    features = {}
    features["id"] = record.id
    features["length"] = length

    for aa in amino_acids:
        features[f"freq_{aa}"] = seq.count(aa) / length

    features["label"] = label

    return features

for record in SeqIO.parse(thermo_file, "fasta"):
    data.append(extract_features(record, 1))

for record in SeqIO.parse(meso_file, "fasta"):
    data.append(extract_features(record, 0))

df = pd.DataFrame(data)

df.to_csv(output_csv, index=False)

print("CSV file created:", output_csv)
print("Total sequences:", len(df))
