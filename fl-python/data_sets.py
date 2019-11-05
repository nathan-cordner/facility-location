import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
from sortedcontainers import *
from timeit import default_timer as timer

import csv

# Import that came from the point clusters code
# %matplotlib inline
import seaborn as sns; sns.set()  # for plot styling
from sklearn.datasets.samples_generator import make_blobs
from sklearn.datasets import make_moons
from sklearn.datasets import make_circles

# Methods to generate starting data

def compute_distances(clients):
    """
        Compute the distance between each facility i and client j
        Return a sorted list of where each entry has the form
        [distance, facility_index, client_index]
    """    
    print("Starting to compute distances")
    start = timer()
    client_distances = []
    counter = [0 for x in range(len(clients))]
    biggest_edge = [0 for x in range(len(clients))]
    c = [[0 for j in range(len(clients))] for i in range(len(clients))]
    for i in range(len(clients)):
        for j in range(i, len(clients)):
            current_distance = np.linalg.norm(clients[i]-clients[j])
            client_distances += [[current_distance, i, j]]
            counter[i] += current_distance
            c[i][j] = current_distance
            if biggest_edge[i] < current_distance:
                biggest_edge[i] = current_distance
            # Distance is symmetric
            if not i == j:
               client_distances += [[current_distance, j, i]]
               counter[j] += current_distance
               c[j][i] = current_distance
               if biggest_edge[j] < current_distance:
                   biggest_edge[j] = current_distance
    client_distances = SortedList(client_distances)
    end = timer()
    results = [len(clients) * biggest_edge[i] - counter[i] for i in range(len(clients))]
    print("Finished computing distances. Time:", end-start)

    return client_distances, c

def plot_initial_points(clients):
    """
        Plot facilities and clients on the same graph
    """

    client_coords = [[],[]]

    for j in clients:
        client_coords[0] += [j[0]]
        client_coords[1] += [j[1]]

    # Clients are red circles
    plt.plot(client_coords[0], client_coords[1], 'ro')
    plt.show()

def make_data_set(num_samples, X):
    """
        Prepare data set X for primal-dual algorithm
    """
    clients = [np.array((X[i][0], X[i][1])) for i in range(num_samples)]
    distances, c = compute_distances(clients)

    # plot_initial_points(clients)    
    return clients, distances, c

def cluster_data(num_clusters, num_samples, var=0.40):
    X, y_true = make_blobs(n_samples=num_samples, centers=num_clusters, cluster_std=var, random_state=0)
    return make_data_set(num_samples, X)

def crescent_moons_data(num_samples):
    X, y_true = make_moons(num_samples, noise=.05, random_state=0)
    return make_data_set(num_samples, X)

def noisy_circles_data(num_samples):
    X, y_true = make_circles(num_samples, factor=.5, noise=.05, random_state=0)
    return make_data_set(num_samples, X)

def anisotropic_data(num_samples):
    X, y = make_blobs(n_samples=num_samples, random_state=170)
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    return make_data_set(num_samples, X_aniso)

def varied_variance_data(num_samples):
    X, y = make_blobs(n_samples=num_samples, cluster_std=[0.4, 1.2, 0.8], random_state=170)
    return make_data_set(num_samples, X)

def read_data_set(file_name):
    """
        Prepare data set X for primal-dual algorithm
    """

    # read dataset from file
    f = open(file_name, "r")
    lines = f.readlines()
    points_list = []
    for line in lines:
        points = line.split(",")
        points_list += [[float(points[0]), float(points[1][:-1])]]
    f.close()

    return make_data_set(len(points_list), points_list)

def write_sorted_distances(distances, file_name):
    """
        Write sorted dataset to file
    """
    num_samples = len(distances)

    # Write dataset to file
    f = open(file_name, "w+")
    for i in range(num_samples - 1):
        f.write(str(distances[i][0]) + "," + str(distances[i][1]) + "," + str(distances[i][2]) + "\n")
    f.write(str(distances[i][0]) + "," + str(distances[i][1]) + "," + str(distances[i][2]))
    f.close()

