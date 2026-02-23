import torch
import pandas as pd
from Bio import SeqIO
from transformers import BertTokenizer, BertModel
from tqdm import tqdm

print("Loading ProtBERT model...")

tokenizer = BertTokenizer.from_pretrained("Rostlab/prot_bert", do_lower_case=False)
model = BertModel.from_pretrained("Rostlab/prot_bert")

model.eval()
device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
model = model.to(device)

def get_protein_embedding(sequence):

    sequence = " ".join(list(sequence))
    sequence = "[CLS] " + sequence + " [SEP]"

    encoded = tokenizer(sequence,
                        return_tensors="pt",
                        padding=True,
                        truncation=True,
                        max_length=1024)

    encoded = {k: v.to(device) for k, v in encoded.items()}

    with torch.no_grad():
        output = model(**encoded)

    embedding = output.last_hidden_state[0,0,:].cpu().numpy()

    return embedding

def process_fasta(fasta_file, label):

    rows = []

    for record in tqdm(list(SeqIO.parse(fasta_file, "fasta"))):
        seq = str(record.seq)

        emb = get_protein_embedding(seq)

        row = {"id": record.id, "label": label}

        for i, val in enumerate(emb):
            row[f"bert_{i}"] = float(val)

        rows.append(row)

    return pd.DataFrame(rows)

print("\nProcessing Thermophiles...")
thermo_df = process_fasta("thermo_cleaned.fasta", 1)

print("\nProcessing Mesophiles...")
meso_df = process_fasta("meso_cleaned.fasta", 0)

final_df = pd.concat([thermo_df, meso_df], ignore_index=True)

final_df.to_csv("bert_features.csv", index=False)

print("\nBERT feature extraction completed!")
print("Saved file: bert_features.csv")
print("Total proteins:", len(final_df))