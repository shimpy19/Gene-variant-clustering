compile (requires C++17):
g++ -std=c++17 -fopenmp -O2 main.cpp algorithm.cpp parse.cpp clustering.cpp -o main

run:
./main.exe

Python dependencies (for clustering):
pip install pyspoa scikit-learn-extra numpy

Notes:
- The program now processes all `.fastq` files in `Bioinformatika - jeleni-2/fastq`.
- For each sample it writes `sequences_<sample>.txt` and `distance_matrix_<sample>.txt`, then runs the clustering script.