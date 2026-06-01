# Gene Variant Clustering

## About the project

This project processes `.fastq` sequencing files, filters reads by sequence length, computes pairwise distances using the Needleman-Wunsch alignment algorithm, and groups similar gene variants into clusters.

For each input sample, the program:

- reads all FASTQ records from a given directory
- keeps only sequences whose length is between 291 and 301 bases
- writes the filtered sequences to the `sequences` directory
- builds a distance matrix and saves it to the `matrices` directory
- runs Python clustering to generate cluster summaries in the `clustering` directory

## Installation and Run Instructions

### Requirements

- C++17 compiler with OpenMP support
- Python 3
- Python packages:

```bash
pip install pyspoa scikit-learn-extra numpy
```

### Build the C++ program

```bash
g++ -std=c++17 -fopenmp -O2 main.cpp algorithm.cpp parse.cpp clustering.cpp -o main
```

### Run the program

On Linux or WSL:

```bash
./main
```

On Windows:

```bash
main.exe
```

The program automatically creates the output folders `matrices`, `sequences`, and `clustering` if they do not already exist.