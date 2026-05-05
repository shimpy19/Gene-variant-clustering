import numpy as np
from sklearn_extra.cluster import KMedoids

filename = "distance_matrix_J29.txt"
distance_matrix = np.loadtxt(filename)

kmedoids = KMedoids(n_clusters=3, metric='precomputed', method='pam', init='k-medoids++', random_state=42)

labels = kmedoids.fit_predict(distance_matrix)

sequences = {}
with open("sequences_J29.txt") as f:
    for line in f:
        idx, seq = line.strip().split("\t", 1)
        sequences[int(idx)] = seq

for idx in kmedoids.medoid_indices_:
    cluster_id = labels[idx]
    count = np.sum(labels == cluster_id)

    print(sequences[idx])
    print(str(len(sequences[idx])) + "bp")
    print("Cluster size:", count)
    print()