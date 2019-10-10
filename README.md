# facility-location
Implementations of the Jain and Vazirani uncapacitated facility location algorithm 
 * We assume a uniform facility cost given by some parameter "lambda" >= 0
 * We also assume that facilities and clients are given by the same list of points
 * We return and plot a clustering result of clients assigned to their nearest open facility

Data input:
 * We allow for 2-dimensional Euclidean data
 * CSV files must list one pair of coordinates x,y on each row, with no header row
 * CSV files must have extension ".txt" and be in the same directory
 * We also have some synthetic data sets that we can generate with the following keywords:
   * cluster: 4 circular clusters 
   * moon: 2 crescent moon shapes
   * circle: 2 concentric circle shapes
   * aniso: 2 skewed circular clusters
   * variedvar: 3 circular clusters of varying sizes and densities
 * Synthetic data sets need a value "num_points" >= 0 specified

C++ Implementation:
 * To compile: g++ facility_location.C -o pd_alg -std=c++11
 * To run (custom data set): bash run_pd_alg.sh <file_name> <lambda>
 * To run (synthetic data set): bash run_pd_alg.sh <keyword> <num_points> <lambda>
