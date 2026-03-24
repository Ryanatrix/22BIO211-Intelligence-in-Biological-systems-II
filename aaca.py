from Bio import SeqIO
from collections import Counter
import numpy as np
import pandas as pd

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

def compute_avg_composition(fasta_file):
    vectors = []

    for record in SeqIO.parse(fasta_file, "fasta"):
        seq = str(record.seq).replace("X", "")
        length = len(seq)
        if length == 0:
            continue

        count = Counter(seq)
        vec = [count[aa] / length for aa in amino_acids]
        vectors.append(vec)

    return np.mean(vectors, axis=0)


thermo = compute_avg_composition("thermo_raw.fasta")
meso = compute_avg_composition("meso_raw.fasta")

df = pd.DataFrame({
    "AA": list(amino_acids),
    "Thermo": thermo,
    "Meso": meso,
    "Diff": thermo - meso
})

print(df)
df.to_csv("composition_analysis.csv", index=False)