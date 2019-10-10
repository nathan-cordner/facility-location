import matplotlib.pyplot as plt
import sys

import seaborn as sns; sns.set()  # for plot styling

def read_data_set(file_name):
    """
        Prepare data set X for primal-dual algorithm
    """

    # read dataset from file
    f = open(file_name, "r")
    lines = f.readlines()
    facility = True
    clusters = []
    current_entry = [[],[]]
    for line in lines:
        if line == "\n":
            clusters += [current_entry]
            current_entry = [[],[]]
            facility = True
        else:
            points = line.split(",")
            if facility:
                current_entry[0] = [float(points[0]), float(points[1][:-1])]
                facility = False
            else:
                current_entry[1] += [[float(points[0]), float(points[1][:-1])]]
    f.close()
    # print(clusters)

    # Display dataset
    plot_final_points(clusters)

def plot_final_points(open_facilities):
    """
        Here the clients are grouped by open facility
    """

    ax = plt.gca()
    facility_centers = []
    for i in open_facilities:
        color = next(ax._get_lines.prop_cycler)['color']
        # iterate over clients assigned to this facility
        for j in i[1]:
            plt.plot([j[0]], [j[1]], marker='o', color=color)
        facility_centers += [[i[0][0], i[0][1]]]
    # plot the facilities
    for i in facility_centers:
        plt.plot([i[0]], [i[1]], marker='s', color='black')
    plt.show()

def main():
    read_data_set(sys.argv[1])

if __name__ == "__main__": main()
