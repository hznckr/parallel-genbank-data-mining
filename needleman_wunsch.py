from mpi4py import MPI
from Bio import SeqIO
import numpy as np
import time
import os

def needleman_wunsch(seq1, seq2):
    len_seq1 = len(seq1) + 1
    len_seq2 = len(seq2) + 1
    score_matrix = np.zeros((len_seq1, len_seq2))

    for i in range(len_seq1):
        score_matrix[i][0] = -i
    for j in range(len_seq2):
        score_matrix[0][j] = -j

    for i in range(1, len_seq1):
        for j in range(1, len_seq2):
            match = score_matrix[i-1][j-1] + (1 if seq1[i-1] == seq2[j-1] else -1)
            delete = score_matrix[i-1][j] - 1
            insert = score_matrix[i][j-1] - 1
            score_matrix[i][j] = max(match, delete, insert)

    return score_matrix[len_seq1 - 1][len_seq2 - 1]

def optimized_parallel_alignment(file_paths):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Rank 0: sadece bir kere dosyaları okur
    if rank == 0:
        sequences = []
        names = []

        for path in file_paths:
            record = SeqIO.read(path, "fasta")
            sequences.append(str(record.seq))
            names.append(os.path.splitext(os.path.basename(path))[0])

    else:
        sequences = None
        names = None

    # Tüm proseslere sekansları ve isimleri gönder
    sequences = comm.bcast(sequences, root=0)
    names = comm.bcast(names, root=0)

    # Tüm karşılaştırmaları belirle
    all_jobs = [(i, j) for i in range(len(sequences)) for j in range(i+1, len(sequences))]

    # İşleri eşit olarak paylaştır
    local_jobs = np.array_split(all_jobs, size)[rank]

    # Her worker kendi işini yapar
    local_results = []
    local_start = time.time()
    for i, j in local_jobs:
        score = needleman_wunsch(sequences[i], sequences[j])
        local_results.append((names[i], names[j], score))
    local_end = time.time()

    # Sonuçları rank 0’a gönder
    all_results = comm.gather(local_results, root=0)

    if rank == 0:
        flat_results = [item for sublist in all_results for item in sublist]
        for seq1, seq2, score in flat_results:
            print(f"{seq1} vs {seq2}: Score = {score}")
        print(f"[Rank 0] Toplam süre: {time.time() - local_start:.4f} saniye")

if __name__ == "__main__":
    folder_path = r"C:\Users\hazan\Desktop\COX1_all"
    fasta_files = [os.path.join(folder_path, f) for f in os.listdir(folder_path) if f.endswith(".fasta")]
    optimized_parallel_alignment(fasta_files)

