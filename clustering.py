import numpy as np
from sklearn_extra.cluster import KMedoids
import argparse
from pathlib import Path
from spoa import poa

# clustering function that uses KMedoids from scikit-learn-extra to cluster the sequences
# based on the distance matrix 
def clustering(distance_matrix_file, sequences_file):
    distance_matrix = np.loadtxt(distance_matrix_file)

    if distance_matrix.ndim == 1:
        distance_matrix = np.atleast_2d(distance_matrix)

    n_samples = distance_matrix.shape[0]
    if n_samples < 2:
        print("Not enough sequences for clustering (need at least 2).")
        return

    inertia = []
    kmedoids = None
    labels = None
    max_k = min(6, n_samples)

    # try different values of k from 2 to max_k and check the inertia to find the optimal number of clusters
    for num in range(2, max_k + 1):
        print(f"Testing k={num}")
        kmedoids = KMedoids(n_clusters=num, metric='precomputed', method='pam', init='k-medoids++', random_state=42)
        labels = kmedoids.fit_predict(distance_matrix)
        inertia.append(kmedoids.inertia_)

        if len(inertia) >= 2:
            previous_inertia = inertia[-2]
            current_inertia = inertia[-1]
            improvement = 0.0 if previous_inertia == 0 else (previous_inertia - current_inertia) / previous_inertia

            if improvement <= 0.1:
                break

            print("Inertia is too high, try more clusters")

    sequences = {}
    with open(sequences_file) as f:
        for line in f:
            idx, seq = line.strip().split("\t", 1)
            sequences[int(idx)] = seq

    if kmedoids is not None and labels is not None:
        clusters = {}
        for i, label in enumerate(labels):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append(sequences[i])

        # create output directory for clustering results if it doesn't exist
        clustering_dir = Path("clustering")
        clustering_dir.mkdir(parents=True, exist_ok=True)

        sequences_path = Path(sequences_file)
        sample_name = sequences_path.stem
        if sample_name.startswith("sequences_"):
            sample_name = sample_name[len("sequences_"):]

        clustering_output_path = clustering_dir / f"clustering_{sample_name}.txt"
        with open(clustering_output_path, "w") as clustering_output_file:
            for cluster_id, cluster_sequences in clusters.items():
                if len(cluster_sequences) > 1:
                    # perform partial order alignment (POA) to get the consensus sequence for the cluster
                    consensus, msa = poa(cluster_sequences)
                    line = f"Cluster {cluster_id} | size: {len(cluster_sequences)} | consensus: {consensus}"
                else:
                    line = f"Cluster {cluster_id} | size: 1 | consensus: {cluster_sequences[0]}"

                print(line)
                clustering_output_file.write(line + "\n")

        print(f"Saved clustering results to: {clustering_output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--distance_matrix', type=str, required=True, help='Path to the distance matrix txt file.')
    parser.add_argument('--sequences', type=str, required=True, help='Path to the sequences txt file.')
    args = parser.parse_args()
    
    clustering(args.distance_matrix, args.sequences)