# 🧬 Parallel Gene Sequence Alignment (Needleman-Wunsch & MPI)

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![MPI](https://img.shields.io/badge/Parallel-MPI-red?style=for-the-badge)
![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-success?style=for-the-badge)

This project implements a high-performance parallel version of the **Needleman-Wunsch** global alignment algorithm. It is designed to process and compare **COX1 genetic sequences** across multiple species by distributing the computational workload across multiple CPU cores using `mpi4py`.

## 🚀 Key Features
- **Parallel Optimization:** Uses Message Passing Interface (MPI) to distribute sequence comparison pairs across worker nodes, significantly reducing total alignment time.
- **Needleman-Wunsch Algorithm:** A robust dynamic programming implementation for optimal global alignment and similarity scoring.
- **Biopython Integration:** Seamlessly parses FASTA records from GenBank for large-scale data processing.
- **Dynamic Load Balancing:** Utilizes `np.array_split` to ensure an even distribution of alignment tasks among available processors.

## 🛠 Tech Stack
- **Language:** Python
- **Parallel Computing:** `mpi4py` (MPI for Python)
- **Bioinformatics:** `Biopython` (SeqIO)
- **Data Handling:** `NumPy`

## 📊 How It Works
1. **Rank 0 (Master):** Reads all `.fasta` files from the specified directory.
2. **Broadcast:** Sequences are broadcasted to all worker processes using `comm.bcast`.
3. **Parallel Alignment:** Each process calculates scores for a specific subset of sequence pairs.
4. **Gather:** Results are collected and summarized by the master process.

## 🏃 Usage
To run the project in parallel (e.g., using 4 cores):
```bash
mpiexec -n 4 python parallel_alignment.py
