import data_sets as ds
import primal_dual_facility_location as pd
import sys

def main():
    num_args = len(sys.argv)
    if num_args == 3:
        file_name = sys.argv[1]
        cost = int(sys.argv[2])
        clients, distances, c = ds.read_data_set(file_name)
    else:
        data_name = sys.argv[1]
        num_points = int(sys.argv[2])
        cost = int(sys.argv[3])
        if data_name == "cluster":
            clients, distances, c = ds.cluster_data(4, num_points)
        elif data_name == "moon":
            clients, distances, c = ds.crescent_moons_data(num_points)
        elif data_name == "circle":
            clients, distances, c = ds.noisy_circles_data(num_points)
        elif data_name == "aniso":
            clients, distances, c = ds.anisotropic_data(num_points)
        elif data_name == "variedvar":
            clients, distances, c = ds.varied_variance_data(num_points)
    pd.facility_location(clients, clients, distances, c, cost)

if __name__ == "__main__": main()
