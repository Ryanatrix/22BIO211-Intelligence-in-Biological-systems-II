from Bio import SeqIO
from collections import Counter
import numpy as np
from sklearn.metrics.pairwise import cosine_similarity

min_len = 100
max_len = 300
max_x = 0.05
limit = 500
threshold = 0.95

thermo_in = "thermo_raw.fasta"
meso_in = "meso_raw.fasta"

thermo_out = "thermo_curated.fasta"
meso_out = "meso_curated.fasta"

combined_out = "final_dataset.fasta"

amino_acids = "ACDEFGHIKLMNPQRSTVWY"

def seq_to_vec(seq):
    seq = seq.replace("X", "")
    length = len(seq)
    count = Counter(seq)

    return np.array([
        count[aa] / length if length > 0 else 0
        for aa in amino_acids
    ])


def curate_cosine(input_file, output_file):

    total = 0
    short_removed = 0
    long_removed = 0
    x_removed = 0
    sim_removed = 0

    kept_vecs = []
    kept_records = []

    for record in SeqIO.parse(input_file, "fasta"):
        total += 1
        seq = str(record.seq).upper()

        if len(seq) < min_len:
            short_removed += 1
            continue
        if len(seq) > max_len:
            long_removed += 1
            continue

        if seq.count("X") / len(seq) > max_x:
            x_removed += 1
            continue

        vec = seq_to_vec(seq)

        is_similar = False
        for old_vec in kept_vecs:
            sim = cosine_similarity([vec], [old_vec])[0][0]
            if sim >= threshold:
                sim_removed += 1
                is_similar = True
                break

        if is_similar:
            continue

        kept_vecs.append(vec)
        kept_records.append(record)

        if len(kept_records) == limit:
            break

    SeqIO.write(kept_records, output_file, "fasta")

    print("\n--------------------------------------")
    print("File:", input_file)
    print("Total:", total)
    print(f"Removed (<{min_len}):", short_removed)
    print(f"Removed (>{max_len}):", long_removed)
    print("Removed (X):", x_removed)
    print("Removed (similar):", sim_removed)
    print("Final saved:", len(kept_records))
    print("--------------------------------------")

    return kept_records

if __name__ == "__main__":

    print("\nProcessing Thermophiles...")
    thermo_records = curate_cosine(thermo_in, thermo_out)

    print("\nProcessing Mesophiles...")
    meso_records = curate_cosine(meso_in, meso_out)

    for i, rec in enumerate(thermo_records):
        rec.id = f"thermo_{i}"
        rec.description = ""

    for i, rec in enumerate(meso_records):
        rec.id = f"meso_{i}"
        rec.description = ""

    combined = thermo_records + meso_records

    SeqIO.write(combined, combined_out, "fasta")

    print("\n✅ Combined dataset created with labels!")
    print("Total combined sequences:", len(combined))
    print("Saved as:", combined_out)