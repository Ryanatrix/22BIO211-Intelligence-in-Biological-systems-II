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

def get_embeddings_batch(sequences, batch_size=8):

    embeddings = []

    for i in tqdm(range(0, len(sequences), batch_size), desc="Embedding batches"):

        batch_seqs = sequences[i:i + batch_size]

        batch_seqs = [" ".join(list(seq)) for seq in batch_seqs]

        encoded = tokenizer(batch_seqs,
                            return_tensors="pt",
                            padding=True,
                            truncation=True,
                            max_length=1024)

        encoded = {k: v.to(device) for k, v in encoded.items()}

        with torch.no_grad():
            output = model(**encoded)

        batch_embeddings = output.last_hidden_state.mean(dim=1).cpu().numpy()

        embeddings.extend(batch_embeddings)

    return embeddings

def process_fasta(fasta_file):

    records = list(SeqIO.parse(fasta_file, "fasta"))

    sequences = []
    ids = []
    labels = []

    for record in records:
        seq = str(record.seq)
        seq_id = record.id.lower()

        if "thermo" in seq_id:
            label = 1
        elif "meso" in seq_id:
            label = 0
        else:
            continue

        sequences.append(seq)
        ids.append(record.id)
        labels.append(label)

    print(f"\nTotal valid sequences: {len(sequences)}")

    embeddings = get_embeddings_batch(sequences)

    data = []

    for i in range(len(embeddings)):
        row = {"id": ids[i], "label": labels[i]}

        for j, val in enumerate(embeddings[i]):
            row[f"bert_{j}"] = float(val)

        data.append(row)

    return pd.DataFrame(data)

print("\nProcessing Combined Dataset...")

df = process_fasta("final_dataset.fasta")

df.to_csv("bert_features2.csv", index=False)

print("\n✅ BERT feature extraction completed!")
print("Saved file: bert_features.csv")
print("Total proteins:", len(df))