import numpy as np
from sklearn_extra.cluster import KMedoids
import argparse
from spoa import poa

def clustering(distance_matrix_file, sequences_file):
    distance_matrix = np.loadtxt(distance_matrix_file)

    inertia = []
    cluster_num = 0
    for num in range(2,7):
        kmedoids = KMedoids(n_clusters=num, metric='precomputed', method='pam', init='k-medoids++', random_state=42)
        labels = kmedoids.fit_predict(distance_matrix)
        inertia.append(kmedoids.inertia_)
        if num > 2 and (inertia[num-2]-inertia[num-1])/inertia[num-2] <= 0.1:
            cluster_num = num
            break
        else:
            print("Inertia is too high, try more clusters")

    sequences = {}
    with open(sequences_file) as f:
        for line in f:
            idx, seq = line.strip().split("\t", 1)
            sequences[int(idx)] = seq

    if kmedoids and labels is not None:
        clusters = {}
        for i, label in enumerate(labels):
            if label not in clusters:
                clusters[label] = []
            clusters[label].append(sequences[i])

        for cluster_id, cluster_sequences in clusters.items():
        #this print is added to check the cluster size and consensus sequence you can remove it when you check the results
            if len(cluster_sequences) > 1:
                consensus, msa = poa(cluster_sequences)
                print(f"Cluster {cluster_id} (size: {len(cluster_sequences)}) consensus: {consensus}")
            else:
                print(f"Cluster {cluster_id} (size: 1) consensus: {cluster_sequences[0]}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('--distance_matrix', type=str, required=True, help='Path to the distance matrix txt file.')
    parser.add_argument('--sequences', type=str, required=True, help='Path to the sequences txt file.')
    args = parser.parse_args()
    
    clustering(args.distance_matrix, args.sequences)