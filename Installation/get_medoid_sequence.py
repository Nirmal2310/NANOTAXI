import os

n_threads = 1

os.environ["OMP_NUM_THREADS"] = str(n_threads)
os.environ["MKL_NUM_THREADS"] = str(n_threads)
os.environ["OPENBLAS_NUM_THREADS"] = str(n_threads)
os.environ["VECLIB_MAXIMUM_THREADS"] = str(n_threads)
os.environ["NUMEXPR_NUM_THREADS"] = str(n_threads)

import sys
from Bio import SeqIO
import Levenshtein
import numpy as np
import re

def calculate_medoid(input_fasta, output_fasta):
    sequences = list(SeqIO.parse(input_fasta, "fasta"))
    
    num_seqs = len(sequences)
    
    dist_sums = np.zeros(num_seqs)
    
    for i in range(num_seqs):
        seq_i = str(sequences[i].seq)
        for j in range(i + 1, num_seqs):
            seq_j = str(sequences[j].seq)
            d = Levenshtein.distance(seq_i, seq_j) / max(len(seq_i), len(seq_j))
            dist_sums[i] += d
            dist_sums[j] += d
    
    medoid_idx = np.argmin(dist_sums)
    medoid_seq = sequences[medoid_idx]
    
    medoid_seq.id = re.sub(r".fasta", "", output_fasta)
    medoid_seq.description = f"Medoid_of_{num_seqs}_sequences"
    
    with open(output_fasta, "w") as out_f:
        SeqIO.write(medoid_seq, out_f, "fasta")

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python get_medoid_sequence.py <input.fasta> <output.fasta>")
    else:
        calculate_medoid(sys.argv[1], sys.argv[2])