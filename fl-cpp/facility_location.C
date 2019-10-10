// a C++ implementation of the primal-dual clustering algorithm

#include "algorithm"
#include "array"
#include "chrono"
#include "cmath"
#include "fstream"
#include "iostream"
#include "string"
#include "vector"
using namespace std;

vector<array<double, 2> > parseCSV(const char* filename) {
    // Reads a text file, inserts data points into a vector
    ifstream file(filename);
    string value_1;
    string value_2;
    vector<array<double, 2> > data_points;

    while (file.good()) {
        getline(file, value_1, ',');
        getline(file, value_2, '\n');
        if (!value_1.empty()) {
            array<double, 2> data_pt = {stod(value_1), stod(value_2)};
            data_points.push_back(data_pt);
        }
    }

    return data_points;
}

vector<array<double, 3> > parseCSVSorted(const char* filename) {
    // Reads a text file, inserts data points into a vector
    ifstream file(filename);
    string value_1;
    string value_2;
    string value_3;
    vector<array<double, 3> > data_points;

    while (file.good()) {
        getline(file, value_1, ',');
        getline(file, value_2, ',');
        getline(file, value_3, '\n');
        // This seems to be necessary here...
        if (!value_1.empty()) {
            array<double, 3> data_pt = {stod(value_1), stod(value_2), stod(value_3)};
            data_points.push_back(data_pt);
        }
    }

    return data_points;    
}

vector<vector<double> > getDistanceGrid(vector<array<double, 3> > &distanceList, int numClients) {
    // put distances from sorted list into grid form

    vector<vector<double> > distanceGrid;
    vector<double> current_vector;
    current_vector.resize(numClients, 0);
    distanceGrid.resize(numClients, current_vector);

    int x; int y; double distance;    
    // cout << "Distance List Size" << distanceList.size() << endl;

    for (int i = 0; i < distanceList.size(); i++) {
        x = (int) distanceList.at(i)[1];
        y = (int) distanceList.at(i)[2];
        distance = distanceList.at(i)[0];

        distanceGrid.at(x).at(y) = distance;
    }

    return distanceGrid;
}

// -------- Methods for computing and sorting distances ---------

// http://www.cplusplus.com/forum/beginner/178293/
double euclideanDistance(double x1, double y1, double x2, double y2) {
	double x = x1 - x2;
	double y = y1 - y2;
	double dist;

	dist = pow(x, 2) + pow(y, 2);
	dist = sqrt(dist);                  

	return dist;
}

// Used for sorting distance vector
bool compareArrays(array<double, 3> a1, array<double, 3> a2) {
    return (a1[0] < a2[0]);
}

// I hope this will create a "reverse" sort, to be able to use the vector pop() operation
bool compareArraysTwo(array<double, 2> a1, array<double, 2> a2) {
    return (a1[0] > a2[0]);
}

vector<array<double, 3> > computeDistances(vector<array<double, 2> > clients){
    /*
     * Compute the distance between each facility i and client j
     * Return a sorted list of where each entry has the form
     * {distance, facility_index, client_index}
     */    
    cout << "Starting to compute distances" << endl;
    vector<array<double, 3> > client_distances;
    vector<array<double, 2> >::const_iterator i;
    vector<array<double, 2> >::const_iterator j;
    double current_distance;
    double x1; double x2; double y1; double y2;
    array<double, 2> entry_1;
    array<double, 2> entry_2;
    double facility_index = 0;
    double client_index = 0;

    for(i = clients.begin(); i != clients.end(); ++i) {
        for (j = i; j != clients.end(); ++j) {
            entry_1 = *i; entry_2 = *j;
            x1 = entry_1[0]; x2 = entry_1[1];
            y1 = entry_2[0]; y2 = entry_2[1];
            current_distance = euclideanDistance(x1, x2, y1, y2);
            array<double, 3> current_entry = {current_distance, facility_index, client_index};
            client_distances.push_back(current_entry);
            // Distance is symmetric
            if (facility_index != client_index) {
               current_entry = {current_distance, client_index, facility_index};
               client_distances.push_back(current_entry);
            }
            client_index++;
        }
        facility_index++;
        client_index = facility_index;
    }
    sort(client_distances.begin(), client_distances.end(), compareArrays);
    cout << "Finished computing distances" << endl;
    return client_distances;
}

// https://stackoverflow.com/questions/12774207/fastest-way-to-check-if-a-file-exist-using-standard-c-c11-c
bool file_exists (const std::string& name) {
    ifstream f(name.c_str());
    return f.good();
}

vector<array<double, 3> > retrieveDistances (const char* filename, vector<array<double, 2> > &data_points) {

    vector<array<double, 3> > distances;
    string sorted_filename(filename);
    sorted_filename.insert(sorted_filename.size()-4,"_sorted");
    
    // Check to see if data set has already been computed
    if (file_exists(sorted_filename)) {
        // File already exists
        distances = parseCSVSorted(sorted_filename.c_str());
    }
    else {
        // compute distances and store results
        auto start = chrono::high_resolution_clock::now();
        distances = computeDistances(data_points);
        auto finish = chrono::high_resolution_clock::now();
        chrono::duration<double> elapsed = finish - start;
        cout << "Elapsed time: " << elapsed.count() << " s\n";

        // Write distances to file
        ofstream output_file;
        output_file.open(sorted_filename.c_str());
        array<double, 3> assignment;

        vector<array<double, 3> >::const_iterator results;

        for (results = distances.begin(); results < distances.end(); results++) {
            assignment = *results;
            output_file << assignment[0] << "," << assignment[1] << "," << assignment[2] << "\n";
        }
        output_file.close();
    }

    return distances;

}


// -------- Methods for the primal dual algorithm ----------

bool openClient (vector<int> client) {
    // Need to avoid looking up elements in empty vectors
    if (!client.empty()) {
    //    if (client.back() == -1) {
        if (find(client.begin(), client.end(), -1) != client.end()) { 
           return false;
        }
    }
    return true;
}

int update_facilities(int i, vector<int> &update_clients, vector<vector<int> > &open_facilities, vector<array<double, 2> > &facility_pay_schedule, vector<array<double, 2> > &next_paid_facility, vector<vector<int> > &clients_copy, vector<vector<double> > &client_w, int num_open_clients, double current_time, vector<double> &client_contribute_times, vector<double> &facility_contributions, vector<double> &facility_opening_cost, vector<array<double, 2> > &client_v) {
    /*
     *  Facility i is now paid for, so we now prevent all clients that are contributing
     *  to facility i from contributing to any other facilities.
     *
     *  We only update clients specified in update_clients
     */

    // cout << "update facilities" << endl;

    // iterate over clients contributing to facility i
    vector<int> clients_to_add;
    vector<int>::const_iterator client_index;
    vector<int> assigned_facilities;
    vector<int>::const_iterator facilities_index;
    vector<int> assigned_clients;
    vector<int>::iterator client_index_2;
    array<double, 2> old_value;
    array<double, 2> new_value;
    double num_clients_facility_f;
    double update_client_pay;
    double current_client_pay;
    double remaining_time;
    int j; int f; int h;


    if (!openClient(update_clients)) {
        return num_open_clients;
    }

    for (client_index = update_clients.begin(); client_index < update_clients.end(); client_index++) {
        j = *client_index;
        // cout << "client index is: " << j << endl;
        if (j == -1) {
            break;
        }

        assigned_facilities = clients_copy.at(j);
        // iterate over facilities that client j contributes to
        if (openClient(assigned_facilities)) {            
            for (facilities_index = assigned_facilities.begin(); facilities_index < assigned_facilities.end(); facilities_index++) {
                f = *facilities_index;
                //assigned_clients = open_facilities.at(f);
                if (openClient(open_facilities.at(f)) && f != i) {
                    if (facility_pay_schedule.at(f)[0] != -1) {
                        if (client_w.at(f).at(j) != -1) {


                            num_clients_facility_f = open_facilities.at(f).size();
                            update_client_pay = num_clients_facility_f * (current_time - client_contribute_times.at(f));
                            old_value = {facility_pay_schedule.at(f)[0], (double) f};
                            facility_contributions.at(f) += update_client_pay;
                            current_client_pay = facility_contributions.at(f);
                            facility_pay_schedule.at(f)[0] -= facility_contributions.at(f);
                            client_contribute_times.at(f) = current_time;
                            if (open_facilities.at(f).size() != 0) {
                                if(find(open_facilities.at(f).begin(), open_facilities.at(f).end(), j) != open_facilities.at(f).end()) {
                                    if (open_facilities.at(f).size() > 1) {
                                        remaining_time = (facility_opening_cost.at(f) - current_client_pay) / (open_facilities.at(f).size() - 1);
                                    }
                                    else {
                                        remaining_time = facility_opening_cost.at(f) - current_time + 1;
                                    }
                                }
                                else {
                                    remaining_time = (facility_opening_cost.at(f) - current_client_pay) / (open_facilities.at(f).size());  
                                }
                            }
                            else {
                                remaining_time = facility_opening_cost.at(f) - current_time + 1;
                            }
                            facility_pay_schedule.at(f)[0] = current_time + remaining_time;

                            new_value = {facility_pay_schedule.at(f)[0], (double) f};
                            // TODO: make this replace & sort more efficient.
                            replace(next_paid_facility.begin(), next_paid_facility.end(), old_value, new_value);
                            sort(next_paid_facility.begin(), next_paid_facility.end(), compareArraysTwo);

                        }
                    
                        // remove all copies of client j from facility f
                        // TODO: make this a log(n) operation
                        open_facilities.at(f).erase(remove(open_facilities.at(f).begin(), open_facilities.at(f).end(), j), open_facilities.at(f).end()); 

                    }

                    // set contribution of client j to facility f equal to 0
                    // client_w.at(f).at(j) = 0;
                }

            }

            clients_copy.at(j).push_back(-1);
            client_v.at(j) = {current_time, (double) i};
            num_open_clients -= 1;
            // cout << "Current Number of Open Clients: " << num_open_clients << endl;
        }
    }

    // Remove facility i from the list of facilities waiting to be paid for
    if (facility_pay_schedule.at(i)[0] != -1) {
        old_value = {facility_pay_schedule.at(i)[0], (double) i};
        // TODO: make this a log(n) operation
        next_paid_facility.erase(remove(next_paid_facility.begin(), next_paid_facility.end(), old_value), next_paid_facility.end());
        facility_pay_schedule.at(i)[0] = -1;
    }

    return num_open_clients;
}


int service_tight_edge(vector<vector<int> > &open_facilities, vector<array<double, 2> > &facilities, vector<vector<int> > &clients_copy, array<double, 3> distance, vector<vector<double> > &client_w, vector<double> &facility_opening_cost, vector<array<double, 2> > &facility_pay_schedule, vector<array<double, 2> > &next_paid_facility, int num_open_clients, vector<double> &client_contribute_times, vector<double> &facility_contributions, vector<array<double, 2> > &client_v){
    /*
     *  Subroutine for the facility location algorithm for the event that
     *  an edge i,j goes tight
     */

    // cout << "service tight edge" << endl;

    vector<int> update_clients;
    array<double, 2> old_value;
    array<double, 2> new_value;
    double update_client_pay;
    double current_client_pay;
    int num_clients_facility_i;

    int i = (int) distance[1];
    int j = (int) distance[2];

    // set the time for when client j starts to contribute
    if (openClient(clients_copy.at(j))) {
        client_w.at(i).at(j) = distance[0];
    }

    // check to see if client j neighbors facility i
    // Assign and remove it if it does, provided the facility is paid for
    if (facility_pay_schedule.at(i)[0] == -1) {
        open_facilities.at(i).push_back(j);
        update_clients.push_back(j);
        num_open_clients = update_facilities(i, update_clients, open_facilities, facility_pay_schedule, next_paid_facility, clients_copy, client_w, num_open_clients, distance[0], client_contribute_times, facility_contributions, facility_opening_cost, client_v);
    }

    // Otherwise, client j contributes to the cost of facility i
    // Update contributions of other clients too
    else {
        if (openClient(clients_copy.at(j))) {
            // Add facility i to client j
            clients_copy.at(j).push_back(i);
        }
        // If facility i currently has no contributing clients
        if (!openClient(open_facilities.at(i))) {
            // Add client j to facility i. Find and replace the -1
            replace(open_facilities.at(i).begin(), open_facilities.at(i).end(), -1, j);
            // Track time last contribution update to facility i
            client_contribute_times.at(i) = distance[0];
            // Initial facility contribution is already set to 0

            // Update when facility i will be paid for
            old_value = {facility_pay_schedule.at(i)[0], (double) i};
            facility_pay_schedule.at(i)[0] = facility_opening_cost.at(i);
            new_value = {facility_pay_schedule.at(i)[0], (double) i};

            // TODO: make this replace & sort more efficient.
            replace(next_paid_facility.begin(), next_paid_facility.end(), old_value, new_value);
            sort(next_paid_facility.begin(), next_paid_facility.end(), compareArraysTwo);

        }

        // Facility i already has some clients
        else {
            // Update client contribution to facility i
            num_clients_facility_i = open_facilities.at(i).size();
            update_client_pay = num_clients_facility_i * (distance[0] - client_contribute_times.at(i));
            facility_contributions.at(i) += update_client_pay;
            current_client_pay = facility_contributions.at(i);

            // Add client j to facility i
            if (openClient(clients_copy.at(j))) {
                open_facilities.at(i).push_back(j);
                num_clients_facility_i++;
            }
            client_contribute_times.at(i) = distance[0];
            // sort(client_contribute_times.at(i).begin(), client_contribute_times.at(i).end());
               
            // Check to see if facility is paid for
            if (current_client_pay >= facility_opening_cost.at(i)) {
                // Clients connected to facility i are removed from contributing to other facilities
                num_open_clients = update_facilities(i, open_facilities.at(i), open_facilities, facility_pay_schedule, next_paid_facility, clients_copy, client_w, num_open_clients, distance[0], client_contribute_times, facility_contributions, facility_opening_cost, client_v);
            }
            else {               
                // Update when facility i will be paid for
                old_value = {facility_pay_schedule.at(i)[0], (double) i};
                double remaining_time;
                if (open_facilities.at(i).size() != 0) {
                    remaining_time = (facility_opening_cost.at(i) - current_client_pay) / num_clients_facility_i;
                }
                else {
                    remaining_time = facility_opening_cost.at(i) - distance[0] + 1;
                }

                facility_pay_schedule.at(i)[0] = distance[0] + remaining_time;
                new_value = {facility_pay_schedule.at(i)[0], (double) i};

                // TODO: make this replace & sort more efficient.
                replace(next_paid_facility.begin(), next_paid_facility.end(), old_value, new_value);
                sort(next_paid_facility.begin(), next_paid_facility.end(), compareArraysTwo);
            }
        }
    }

    return num_open_clients;
}

vector<vector<int> > output_cluster_results(vector<array<double, 2> > &clients, vector<vector<double> > &client_w, vector<vector<double> > &distanceGrid, vector<array<double, 2> > &facility_pay_schedule, vector<vector<int> > &open_facilities, vector<array<double, 2> > &client_v, double current_time) {

    // I'm going to check contributions based off of current_time - client_w (e.g. end time minus start time)

    vector<int> assigned_clients;
    vector<int> temp_open_facilities;
    for (int i = 0; i < facility_pay_schedule.size(); i++) {
        if (facility_pay_schedule.at(i)[0] == -1)
            temp_open_facilities.push_back(i);
    }
    vector<vector<int> > opened_facilities;
    int i;

    while (!temp_open_facilities.empty()) {
        // get and remove a facility
        i = temp_open_facilities.back();
        temp_open_facilities.pop_back();

        // first coordinate is facility, rest are indices of assigned clients
        vector<int> facility_assignment;
        facility_assignment.push_back(i);
        
        for (int j = 0; j < client_v.size(); j++) {
      
            //int current_client = open_facilities.at(i).at(j);
            // double current_contribution = current_time - client_w.at(i).at(current_client);
            int current_client = j;

            // need to use time between when client starts contributing to facility
            // and when client has its contributions capped
            int witness = (int) client_v.at(current_client)[1];
            double witness_cost = (client_v.at(current_client)[0] - client_w.at(witness).at(current_client)) + distanceGrid.at(witness).at(current_client);
            double current_cost = (client_v.at(current_client)[0] - client_w.at(i).at(current_client)) + distanceGrid.at(i).at(current_client);

            if (current_cost <= witness_cost && client_w.at(i).at(current_client) != -1) {
                assigned_clients.push_back(current_client);
                facility_assignment.push_back(current_client);

                vector<int> facilities_to_remove;
                // track other potential open facilities that are counting on contribution from j
                for (int h = 0; h < temp_open_facilities.size(); h++) {
                    int current_facility = temp_open_facilities.at(h);
                    vector<int> current_client_list = open_facilities.at(current_facility);

                    if (client_w.at(i).at(current_client) > 0 && client_w.at(current_facility).at(current_client) > 0) {
                        // blacklist the current facility
                        facilities_to_remove.push_back(current_facility);
                    }
                }
                // erase the offending facilities
                for (int h = 0; h < facilities_to_remove.size(); h++) {
                    temp_open_facilities.erase(remove(temp_open_facilities.begin(), temp_open_facilities.end(), facilities_to_remove.at(h)), temp_open_facilities.end());
                }
            }
        }

        if (facility_assignment.size() > 1) {
            opened_facilities.push_back(facility_assignment);
        }
    }

    for (int j = 0; j < clients.size(); j++) {
        if (find(assigned_clients.begin(), assigned_clients.end(), j) == assigned_clients.end()) {
            // cout << "Assigning client " << j << endl;
            int current_facility;
            double current_distance;
            double min_distance = -1;
            int min_facility_index;

            for (int h = 0; h < opened_facilities.size(); h++) {
                current_facility = opened_facilities.at(h).at(0);
                current_distance = distanceGrid.at(current_facility).at(j);

                if (min_distance == -1) {
                    min_distance = current_distance;
                    min_facility_index = h;
                }
                else if (current_distance < min_distance) {
                    min_distance = current_distance;
                    min_facility_index = h;
                }
            }
            opened_facilities.at(min_facility_index).push_back(j);
        }
    }

    return opened_facilities;

}

int facility_location(vector<array<double, 2> > clients, vector<array<double, 3> > client_distances, double opening_cost, string file_name){
    /*
     * Implementation of algorithm outline found on page 183 of
     * "The Design of Approximation Algorithms" by Williamson and Shmoys,
     * and page 282 of "Approximation Algorithms for Metric Facility
     * Location and k-Median Problems Using the Primal-Dual Schema and
     * Lagrangian Relaxation" by Jain and Vazirani
     *
     * Takes as input a lists representing the facilites and clients,
     * a (sorted) list of distances, and an integer for the facility opening cost.
     * Executes the primal-dual uncapacitated facility location algorithm
     * and finds a good approximate assignment for clients to facilities.
     */

    array<double, 3> distance;
    array<double, 2> next_facility;
    double next_edge;

    cout << "Starting primal-dual algorithm" << endl;
    auto start_2 = chrono::high_resolution_clock::now();

    // note that facilities and clients are the same here
    int num_client_distances = client_distances.size();
    int num_clients = clients.size();

    // obtain the grid of distances c[i][j] to use for clustering assignments
    vector<vector<double> > distanceGrid = getDistanceGrid(client_distances, num_clients);

    // Variable w -- how much a client contributes to a particular facility
    // client_w[i][j] stores the time when client j starts to contribute to facility i
    vector<vector<double> > client_w;
    vector<double> current_vector;
    current_vector.resize(num_clients, -1);
    client_w.resize(num_clients, current_vector);

    // Record time and facility when clients are declared connected
    vector<array<double, 2> > client_v;
    client_v.resize(num_clients, {-1, -1});

    // We assume all facilitiy costs are equal
    vector<double> facility_opening_cost;
    facility_opening_cost.resize(num_clients, opening_cost);
    
    // Variable S -- gets a copy of the list of clients
    // Each client stores a list of facilities that it's assigned to
    vector<vector<int> > clients_copy;
    clients_copy.resize(num_clients);
    int num_open_clients = clients_copy.size();    

    // Variable T -- starts empty, but will accept facilities
    vector<vector<int> > open_facilities;
    vector<int> initial_list_open_facilities;
    initial_list_open_facilities.push_back(-1);
    open_facilities.resize(num_clients, initial_list_open_facilities);

    // Want to also track when clients start to contribute to a given facility
    vector<double> client_contribute_times;
    client_contribute_times.resize(num_clients, 0);

    // We will also track how much is currently being contributed to a given facility
    vector<double> facility_contributions;
    facility_contributions.resize(num_clients, 0);

    // Initial time each facility will be paid for is one more than the max edge length
    // Second coordinate keeps track of the facility index
    // double maximum_distance = client_distances.at(num_client_distances - 1)[0];
    vector<array<double, 2> > facility_pay_schedule;
    facility_pay_schedule.resize(num_clients);
    for (int i = 0; i < num_clients; i++){
        array<double, 2> initial_facility_pay = {opening_cost + 1, (double) i};
        facility_pay_schedule.at(i) = initial_facility_pay;
    } 
    vector<array<double, 2> > next_paid_facility = facility_pay_schedule;

    // Keep track of the index of the edge we service next
    int current_edge_index = 0;

    // Keep track of how much time has passed
    double current_time = 0;    

    // Here we flag "closed" clients and facilities with the value -1
    while (num_open_clients > 0) {
        
        if (current_edge_index == num_client_distances) {
            cout << "All edges have been evaluated. Continuing algorithm" << endl;
            while (!next_paid_facility.empty()) {
                next_facility = next_paid_facility.back();
                next_paid_facility.pop_back();
                if (next_facility[0] >= 0) {
                    current_time = next_facility[0];
                }
                num_open_clients = update_facilities(next_facility[1], open_facilities.at(next_facility[1]), open_facilities, facility_pay_schedule, next_paid_facility, clients_copy, client_w, num_open_clients, current_time, client_contribute_times, facility_contributions, facility_opening_cost, client_v);
                if (num_open_clients <= 0) {
                    break;
                }
            }
            break;
        }
        // cout << "Looking at edge " << current_edge_index << endl;

        distance = client_distances.at(current_edge_index);
        next_edge = distance[0];

        // Decide which event happens next
        if (!next_paid_facility.empty()) {
            next_facility = next_paid_facility.back();
            next_paid_facility.pop_back();
        }
        else {
            // Don't go here!
            next_facility = {next_edge + 1, 0};
        }
        
        if (next_edge <= next_facility[0]) {
            // An edge goes tight
            // Save some time by checking this first
            // if (openClient(clients_copy.at(distance[2]))) {
            num_open_clients = service_tight_edge(open_facilities, clients, clients_copy, distance, client_w, facility_opening_cost, facility_pay_schedule, next_paid_facility, num_open_clients, client_contribute_times, facility_contributions, client_v);
            // }
            current_time = distance[0];
            current_edge_index += 1;
        }
        else {
            // A facility is now paid for
            num_open_clients = update_facilities(next_facility[1], open_facilities.at(next_facility[1]), open_facilities, facility_pay_schedule, next_paid_facility, clients_copy, client_w, num_open_clients, current_time, client_contribute_times, facility_contributions, facility_opening_cost, client_v);
            current_time = next_facility[0];
        }
        // cout << "Current number of open clients: " << num_open_clients << endl;

    }

    auto finish_2 = chrono::high_resolution_clock::now();
    chrono::duration<double> elapsed_2 = finish_2 - start_2;
    cout << "Finished primal-dual algorithm" << endl;
    cout << "Elapsed time: " << elapsed_2.count() << " s\n";

    // Get results of pruning step
    vector<vector<int> > opened_facilities = output_cluster_results(clients, client_w, distanceGrid, facility_pay_schedule, open_facilities, client_v, current_time);

    vector<vector<int> >::const_iterator results;
    vector<int> assignment;
    vector<int>::const_iterator results_2; // list of clients assigned to facility
    int client_assigned_to_facility;
    int num_facilities_opened = opened_facilities.size();

    ofstream output_file;
    output_file.open("pd_result_"+file_name);

    // output algorithm results
    for (results = opened_facilities.begin(); results < opened_facilities.end(); results++) {
        assignment = *results;
        output_file << clients.at(*assignment.begin())[0] << "," << clients.at(*assignment.begin())[1] << endl;
        // both facility and assigned clients are written here
        for (results_2 = assignment.begin(); results_2 < assignment.end(); results_2++) {
            client_assigned_to_facility = *results_2;
            output_file << clients.at(client_assigned_to_facility)[0] << "," << clients.at(client_assigned_to_facility)[1] << endl;
        }
        output_file << endl;

    }

    output_file.close();
    // output_file.close();
    cout << "Number of facilities: " << num_facilities_opened << endl;

    return 0;
}

// ---------------------------------------------------------

int main(int argc, char** argv) {

    vector<array<double, 2> > data_points = parseCSV(argv[1]);
    vector<array<double, 3> > distances = retrieveDistances(argv[1], data_points);
    double lambda = atof(argv[2]);

    // Run primal-dual algorithm
    facility_location(data_points, distances, lambda, argv[1]);
    /*
    while (true){
        cout << "Enter lambda: ";
        cin >> lambda;
        facility_location(data_points, distances, lambda, argv[1]);
    }
    */

    return 0;
}
