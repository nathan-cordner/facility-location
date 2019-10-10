import numpy as np
import random
import matplotlib.pyplot as plt
from sortedcontainers import *
from timeit import default_timer as timer

import sys

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
    for i in range(len(clients)):
        for j in range(i, len(clients)):
            current_distance = np.linalg.norm(clients[i]-clients[j])
            client_distances += [[current_distance, i, j]]
            # Distance is symmetric
            if not i == j:
               client_distances += [[current_distance, j, i]]
    client_distances = SortedList(client_distances)
    end = timer()
    print("Finished computing distances. Time:", end-start)
    return client_distances

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

def random_data(num_clients, size):
    """
        Returns a list of n unique random points in the
        Euclidean plane. Coordinates are limited to
        being between -size and size.
    """
    points = []
    while len(points) < num_clients:
        x_coord = random.randint(-size,size)
        y_coord = random.randint(-size,size)
        new_point = np.array((x_coord, y_coord))
        unique_point = True
        for i in points:
            if np.array_equal(new_point,i):
                unique_point = False
        if unique_point:
            points += [new_point]
    distances = compute_distances(points)

    plot_initial_points(points)
    return points, distances

def make_data_set(num_samples, X, file_name):
    """
        Prepare data set X for primal-dual algorithm
    """
    clients = [np.array((X[i][0], X[i][1])) for i in range(num_samples)]

    # Write dataset to file
    f = open(file_name, "w+")
    for i in range(num_samples - 1):
        f.write(str(clients[i][0]) + "," + str(clients[i][1]) + "\n")
    f.write(str(clients[num_samples - 1][0]) + "," + str(clients[num_samples - 1][1]))
    f.close()

    # Display dataset
    # plot_initial_points(clients)

    # Maybe I will do this in C++
    # distances = compute_distances(clients) 
    # return clients, distances, X

def cluster_data(num_clusters, num_samples, var=0.40):
    X, y_true = make_blobs(n_samples=num_samples, centers=num_clusters, cluster_std=var, random_state=0)
    name = str(num_samples) + "_cluster.txt"
    return make_data_set(num_samples, X, name)

def crescent_moons_data(num_samples):
    X, y_true = make_moons(num_samples, noise=.05, random_state=0)
    name = str(num_samples) + "_moon.txt"
    return make_data_set(num_samples, X, name)

def noisy_circles_data(num_samples):
    X, y_true = make_circles(num_samples, factor=.5, noise=.05, random_state=0)
    name = str(num_samples) + "_circle.txt"
    return make_data_set(num_samples, X, name)

def anisotropic_data(num_samples):
    X, y = make_blobs(n_samples=num_samples, random_state=170)
    transformation = [[0.6, -0.6], [-0.4, 0.8]]
    X_aniso = np.dot(X, transformation)
    name = str(num_samples) + "_aniso.txt"
    return make_data_set(num_samples, X_aniso, name)

def varied_variance_data(num_samples):
    X, y = make_blobs(n_samples=num_samples, cluster_std=[0.4, 1.2, 0.8], random_state=170)
    name = str(num_samples) + "_variedvar.txt"
    return make_data_set(num_samples, X, name)

def main():
    data_name = sys.argv[1]
    num_points = int(sys.argv[2])
    if data_name == "cluster":
        cluster_data(4, num_points)
    elif data_name == "moon":
        crescent_moons_data(num_points)
    elif data_name == "circle":
        noisy_circles_data(num_points)
    elif data_name == "aniso":
        anisotropic_data(num_points)
    elif data_name == "variedvar":
        varied_variance_data(num_points)

if __name__ == "__main__": main()

