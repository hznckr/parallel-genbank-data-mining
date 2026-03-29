# 🧬 Parallel Gene Sequence Alignment (Needleman-Wunsch & MPI)

![Python](https://img.shields.io/badge/Python-3.9+-3776AB?style=for-the-badge&logo=python&logoColor=white)
![MPI](https://img.shields.io/badge/Parallel-MPI-red?style=for-the-badge)
![Bioinformatics](https://img.shields.io/badge/Domain-Bioinformatics-success?style=for-the-badge)

This repository contains a high-performance implementation of the **Needleman-Wunsch algorithm** for global sequence alignment. Using **MPI (Message Passing Interface)**, the system parallelizes the comparison of multiple genetic sequences, specifically focusing on the **COX1 (Cytochrome c oxidase subunit I)** gene across different species.

## 🚀 Key Features
- **Parallel Computing with MPI:** Distributes alignment tasks across multiple CPU cores using `mpi4py`, significantly reducing the time complexity for large-scale comparisons.
- **Global Alignment Algorithm:** A robust dynamic programming implementation of the Needleman-Wunsch algorithm for optimal similarity scoring.
- **Bioinformatics Integration:** Seamlessly parses FASTA records using `Biopython` for systematic genomic analysis.
- **Scalable Architecture:** Designed to scale linearly with the number of available processors.

## 🛠 Tech Stack
- **Programming Language:** Python
- **Parallelism:** `mpi4py` (MPI for Python)
- **Library Support:** `Biopython` (SeqIO), `NumPy`
- **Dataset:** 16 distinct COX1 gene sequences retrieved from GenBank (FASTA format).

## 📊 How It Works
1. **Master Node (Rank 0):** Loads 16 COX1 sequence files from the `/data` directory.
2. **Broadcast:** Distributes sequence data to all worker nodes via `comm.bcast`.
3. **Parallel Tasking:** Worker nodes independently calculate alignment scores for assigned pairs.
4. **Gather:** Rank 0 collects all scores and outputs the final evolutionary distance matrix/summary.

## 📂 Project Structure
- `parallel_alignment.py`: Core parallelized Python script.
- `/data`: Contains 16 sample `.fasta` files for immediate testing.
- `requirements.txt`: Project dependencies.

## 🏃 Installation & Usage
### 1. Requirements
Ensure you have an MPI implementation (like MS-MPI or MPICH) installed on your system.
```bash
pip install mpi4py biopython numpy
