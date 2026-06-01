# Gene Variant Clustering

## About the project

This project processes `.fastq` sequencing files (deer samples) to produce clusters of similar gene variants. The pipeline is split into two main phases:

- **Phase 1 — per-sample clustering**
	- Parse FASTQ files and extract reads (`parseFastq`). Sequences are filtered by length (the code expects ~291–301 bp samples).
	- Compute an all-vs-all alignment score per sample using a Needleman–Wunsch implementation (scoring: match=+1, mismatch=-1, gap=-2).
	- Convert alignment scores to a distance measure using `distance = max(0, 296 - score)` and write the distance matrix to `matrices/`.
	- Call the Python clustering script (`clustering.py`) to produce cluster summaries and cluster consensus representatives in `clustering/`.

- **Phase 2 — representative grouping (post-clustering)**
	- Read cluster summaries and keep representatives from clusters whose `size >= MINIMUM_CLUSTER_SIZE` (default 20).
	- Compare representatives pairwise with the same Needleman–Wunsch routine; if two representatives have `distance < DISTANCE_THRESHOLD` (default 10) they are merged using a union-find structure.
	- For each final group the code selects a medoid representative (the member with the smallest total distance to other members) and writes grouped output to `clustering/representative_groups.txt`.

## Methods and key parameters

- Alignment: Needleman–Wunsch (global alignment), parameters in `main.cpp` and `algorithm.cpp`.
- Distance conversion: `distance = max(0, 296 - alignment_score)` (project-specific reference value).
- Group thresholds: `MINIMUM_CLUSTER_SIZE` (default 20) and `DISTANCE_THRESHOLD` (default 10) are defined in `main.cpp`.

## Installation and Run Instructions

### Requirements

- C++17 compiler with OpenMP support
- Python 3 (for `clustering.py` and its dependencies)

Install required Python packages inside the same Python environment that will be used when `main` calls `clustering.py` (recommended: WSL virtualenv `venv_linux`):

```bash
source venv_linux/bin/activate
pip install numpy==1.26.4 pyspoa scikit-learn-extra==0.3.0
```

### Build the C++ program

```bash
g++ -std=c++17 -fopenmp -O2 main.cpp algorithm.cpp parse.cpp clustering.cpp -o main
```

### Run the program (recommended inside WSL)

```bash
source venv_linux/bin/activate   # ensure the Python environment is active in WSL
./main.exe                       # or ./main if you built without .exe
```

- Output locations: `matrices/`, `sequences/`, `clustering/`.

**Runtime warning:** Clustering the full deer dataset is compute-intensive — expect approximately **30 minutes** on a typical laptop/desktop. Monitor CPU and disk usage, and run on a machine with sufficient resources.

### Notes

- If you run the compiled binary from Windows rather than WSL, the `system()` call will launch the Windows Python interpreter; make sure `numpy` and the other Python deps are installed for that interpreter as well.
- To adjust thresholds, edit `MINIMUM_CLUSTER_SIZE` and `DISTANCE_THRESHOLD` in `main.cpp` and recompile.