import numpy as np
from sklearn_extra.cluster import KMedoids

filename = "distance_matrix_J29.txt"
distance_matrix = np.loadtxt(filename)

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
        print("inertia is too high, try more clusters")

sequences = {}
with open("sequences_J29.txt") as f:
    for line in f:
        idx, seq = line.strip().split("\t", 1)
        sequences[int(idx)] = seq

for idx in kmedoids.medoid_indices_:
    cluster_id = labels[idx]
    count = np.sum(labels == cluster_id)
    print("Cluster size:", count)