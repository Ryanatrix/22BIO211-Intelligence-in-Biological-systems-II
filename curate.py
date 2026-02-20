from Bio import SeqIO
from Bio.Align import PairwiseAligner

min_length = 50
max_x_ratio = 0.05
similarity_threshold = 0.90
required_count = 500

thermo_input = "thermo2.fasta"
meso_input = "meso2.fasta"

thermo_output = "thermo_cleaned.fasta"
meso_output = "meso_cleaned.fasta"

global_align = PairwiseAligner()
global_align.mode = "global"
global_align.match_score = 2
global_align.mismatch_score = -1
global_align.open_gap_score = -2
global_align.extend_gap_score = -1

local_align = PairwiseAligner()
local_align.mode = "local"
local_align.match_score = 2
local_align.mismatch_score = -1
local_align.open_gap_score = -2
local_align.extend_gap_score = -1


def compute_similarity(seq1, seq2):
    score = global_align.score(seq1, seq2)
    max_score = min(len(seq1), len(seq2)) * global_align.match_score
    return score / max_score


def clean_fasta(input_file, output_file):

    kept_sequences = []
    kept_records = []

    for record in SeqIO.parse(input_file, "fasta"):

        seq = str(record.seq).upper()

        if len(seq) < min_length:
            continue

        if seq.count("X") / len(seq) > max_x_ratio:
            continue

        if seq in kept_sequences:
            continue

        redundant = False
        for old_seq in kept_sequences:
            if compute_similarity(seq, old_seq) >= similarity_threshold:
                redundant = True
                break

        if redundant:
            continue

        kept_sequences.append(seq)
        kept_records.append(record)

        if len(kept_records) == required_count:
            break

    SeqIO.write(kept_records, output_file, "fasta")

    print(f"{output_file} saved with {len(kept_records)} sequences.")

    if len(kept_sequences) >= 2:
        g_score = global_align.score(kept_sequences[0], kept_sequences[1])
        l_score = local_align.score(kept_sequences[0], kept_sequences[1])

        print("Example alignment scores:")
        print("Global:", g_score)
        print("Local :", l_score)
        print()

clean_fasta(thermo_input, thermo_output)
clean_fasta(meso_input, meso_output)
