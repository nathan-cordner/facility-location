import numpy as np
import matplotlib.pyplot as plt
from sortedcontainers import *
from timeit import default_timer as timer

import seaborn as sns; sns.set()  # for plot styling
from matplotlib.gridspec import GridSpec

def facility_location(clients, facilities, distances_list, c, opening_cost):
    """
        Implementation of algorithm outline found on page 183 of
        "The Design of Approximation Algorithms" by Williamson and Shmoys,
        and page 282 of "Approximation Algorithms for Metric Facility
        Location and k-Median Problems Using the Primal-Dual Schema and
        Lagrangian Relaxation" by Jain and Vazirani

        Takes as input a lists representing the facilites and clients,
        a (sorted) list of distances, and an integer for the facility opening cost.
        Executes the primal-dual uncapacitated facility location algorithm
        and finds a good approximate assignment for clients to facilities.
    """
    print("Starting primal-dual algorithm")

    start = timer()
    
    # Variable c_ij -- the distance between facility i and client j
    # Assume this is already in sorted order
    client_distances = distances_list
    num_client_distances = len(client_distances)

    # Variable w -- how much a client contributes to a particular facility
    # client_w[i][j] stores the time when client j starts to contribute to facility i
    client_w = [[-1 for j in range(len(clients))] for i in range(len(facilities))]

    # Record contribution when clients are declared connected
    client_v = [-1 for j in range(len(clients))]

    # We assume all facilitiy costs are equal
    facility_opening_cost = [opening_cost] * len(facilities)

    # Variable S -- gets a copy of the list of clients
    # Each client stores a list of facilities that it's assigned to
    # Coordinate 2 says if the client is currently open
    clients_copy = [[x, []] for x in clients]
    num_open_clients = len(clients_copy)
    
    # Variable T -- starts empty, but will accept facilities
    open_facilities = [-1] * len(facilities)

    # Initial time each facility will be paid for is one more than lambda
    # Second coordinate keeps track of the facility index
    
    maximum_distance_plus_one = opening_cost + 1
    facility_pay_schedule = [[maximum_distance_plus_one, x] for x in range(len(facilities))]    
    next_paid_facility = SortedList(facility_pay_schedule)

    # Keep track of the index of the edge we service next
    current_edge_index = 0

    # Keep track of how much time has passed
    current_time = 0

    # Want to track how many times we update facility cost
    counter = [0]
       
    # Here we flag "closed" facilities with the value -1
    while num_open_clients > 0:
        
        if current_edge_index == num_client_distances:
            print("All edges have been evaluated. Continuing algorithm")
            # What if the algorithm needs more time?
            while not next_paid_facility == []:
                next_facility = next_paid_facility[0]
                if next_facility[0] >= 0:
                    current_time = next_facility[0]
                num_open_clients = update_facilities(next_facility[1], open_facilities[next_facility[1]][1], open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, current_time)
                if num_open_clients <= 0:
                    break
            break

        distance = client_distances[current_edge_index]
        next_edge = distance[0]

        # Decide which event happens next
        if not next_paid_facility == []:
            next_facility = next_paid_facility[0]
        else:
            # Don't go here!
            next_facility = [next_edge + 1, 0]
        if next_edge <= next_facility[0]:
            # An edge goes tight
            # if not clients_copy[distance[2]] == -1:
            num_open_clients = service_tight_edge(open_facilities, facilities, clients_copy, distance, client_w, facility_opening_cost, facility_pay_schedule, next_paid_facility, num_open_clients, counter, client_v)
            current_time = distance[0]
            current_edge_index += 1
        else:
            # A facility is now paid for
            if next_facility[0] >= 0:
                current_time = next_facility[0]
            num_open_clients = update_facilities(next_facility[1], open_facilities[next_facility[1]][1], open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, current_time)

    end = timer()
    print("Finished primal-dual algorithm. Time:", end - start)
    print("current_time:",current_time)

    # Some facilities might be paid for now, but haven't officially opened.
    # We fix that here.

    for i in range(len(facilities)):
        if not facility_pay_schedule[i] == -1:
            # Update client contribution to facility i
            # Only update pay from clients that are not yet connected
            
            num_clients_facility_i = len(open_facilities[i][1])
            update_client_pay = num_clients_facility_i * (current_time - open_facilities[i][2])
            open_facilities[i][3] += update_client_pay
            current_client_pay = open_facilities[i][3]

            if current_client_pay >= facility_opening_cost[i]:
                # Clients connected to facility i are removed from contributing to other facilities
                num_open_clients = update_facilities(i, open_facilities[i][1], open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, current_time)


    pruning_primal_obj, pruning_opened_facilities, extra_cost = cluster_assignment_via_pruning(clients, facilities, client_w, c, facility_pay_schedule, client_v, opening_cost)

    print("num open facilities:", len(pruning_opened_facilities))

    # plot_final_points(clients, client_assignments)
    for x in pruning_opened_facilities:
        x[0] = facilities[x[0]]

    plot_final_clusters(pruning_opened_facilities)

def service_tight_edge(open_facilities, facilities, clients_copy, distance, client_w, facility_opening_cost, facility_pay_schedule, next_paid_facility, num_open_clients, counter, client_v):
    """
        Subroutine for the facility location algorithm for the event that
        an edge i,j goes tight
    """

    i = distance[1]
    j = distance[2]
    # print("Checking connection", i, j, "at time", distance[0])

    # Set the time for when client j starts to contribute
    if not clients_copy[j] == -1:
        client_w[i][j] = distance[0]

    # check to see if client j neighbors facility i
    # Assign and remove it if it does, provided the facility is paid for
    if facility_pay_schedule[i] == -1:
        if not clients_copy[j] == -1:
            open_facilities[i][1].add(j)
        num_open_clients = update_facilities(i, [j], open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, distance[0])

    # Otherwise, client j contributes to the cost of facility i
    # Update contributions of other clients too
    else:
        if not clients_copy[j] == -1:
            clients_copy[j][1] += [i]
        # If facility i currently has no contributing clients
        if open_facilities[i] == -1:
            # Second coordinate -- track clients contributing to facility i
            # Third coordinate tracks times of last cost update
            # Fourth coordinate tracks running contribution of clients so far
            open_facilities[i] = [facilities[i], SortedList(), distance[0], 0]
            open_facilities[i][1].add(j)
            # Update when facility i will be paid for.
            next_paid_facility.remove([facility_pay_schedule[i][0], i])
            facility_pay_schedule[i][0] = facility_opening_cost[i]
            next_paid_facility.add([facility_pay_schedule[i][0], i])

        # Facility i already has some clients
        else:
            # Update client contribution to facility i
            num_clients_facility_i = len(open_facilities[i][1])

            update_client_pay = num_clients_facility_i * (distance[0] - open_facilities[i][2])
            open_facilities[i][3] += update_client_pay
            # open_facilities[i][3] += new_client_pay
            current_client_pay = open_facilities[i][3]
            # print("Facility",i,"now will be paid in",facility_pay_schedule[i][0] - current_client_pay)

            # Add client j to facility i
            if not clients_copy[j] == -1:
                open_facilities[i][1].add(j)
            open_facilities[i][2] = distance[0]

            if current_client_pay >= facility_opening_cost[i]:
                # Clients connected to facility i are removed from contributing to other facilities
                num_open_clients = update_facilities(i, open_facilities[i][1], open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, distance[0])

            else:               
                # Update when facility i will be paid for
                counter[0] += 1

                next_paid_facility.remove([facility_pay_schedule[i][0], i])
                if not len(open_facilities[i][1]) == 0:
                    remaining_time = (facility_opening_cost[i] - current_client_pay) / len(open_facilities[i][1])
                else:
                    # No currently contributing clients, so it should not get paid for anytime soon
                    remaining_time = facility_opening_cost[i] - distance[0] + 1
                # print("new pay schedule is", facility_pay_schedule[i][0])
                facility_pay_schedule[i][0] = distance[0] + remaining_time
                next_paid_facility.add([facility_pay_schedule[i][0], i])

    return num_open_clients

def update_facilities(i, update_clients, open_facilities, facility_pay_schedule, facility_opening_cost, next_paid_facility, clients_copy, num_open_clients, client_w, client_v, current_time):
    """
        Facility i is now paid for, so we declare its clients to be connected to it.

        We only update clients specified in update_clients
    """
    # iterate over clients contributing to facility i
    clients_to_add = []
    for j in update_clients:
        if not clients_copy[j] == -1:
            # declare client j to be connected
            # print("Client",j,"connects to facility", i,"at time", current_time)

            # Remove client j from contributing to anymore facilities
            # Update when other facilities are going to be paid for
            for f in clients_copy[j][1]:
                if not open_facilities[f] == -1 and not f == i:

                    if not facility_pay_schedule[f] == -1:
                        if not client_w[f][j] == -1:
                            next_paid_facility.remove([facility_pay_schedule[f][0], f])
                            # print("Updating", f, "old val", facility_pay_schedule[f][0])
                            # Record and update contribution of client j to facility f
                            num_clients_facility_f = len(open_facilities[f][1])
                            update_client_pay = num_clients_facility_f * (current_time - open_facilities[f][2])
                            open_facilities[f][3] += update_client_pay
                            current_client_pay = open_facilities[f][3]
                            facility_pay_schedule[f][0] = facility_opening_cost[f] - current_client_pay
                            open_facilities[f][2] = current_time
                            if not len(open_facilities[f][1]) == 0:
                                if j in open_facilities[f][1]:
                                    if len(open_facilities[f][1]) > 1:
                                        remaining_time = (facility_opening_cost[f] - current_client_pay) / (len(open_facilities[f][1]) - 1)
                                    else:
                                        remaining_time = facility_opening_cost[f] - current_time + 1
                                else:
                                    remaining_time = (facility_opening_cost[f] - current_client_pay) / len(open_facilities[f][1])
                            else:
                                # No currently contributing clients, so it should not get paid for anytime soon
                                remaining_time = facility_opening_cost[f] - current_time + 1
                            # print("new pay schedule is", facility_pay_schedule[f][0])
                            facility_pay_schedule[f][0] = current_time + remaining_time

                            next_paid_facility.add([facility_pay_schedule[f][0], f])
                            # print("new val", facility_pay_schedule[f][0])

                    # remove all copies of client j from facility f
                    open_facilities[f][1].remove(j)

            if num_open_clients == len(facility_opening_cost):
                print("Time of first facility opening", current_time)
            client_v[j] = [current_time, i] # record time and facility when client j is connected
            clients_copy[j] = -1
            num_open_clients -= 1
   
    # Remove facility i from the list of facilities waiting to be paid for
    if not facility_pay_schedule[i] == -1:
        next_paid_facility.remove([facility_pay_schedule[i][0], i])
        facility_pay_schedule[i] = -1
        # print("Facility",i,"is now open")

    return num_open_clients

def cluster_assignment_via_pruning(clients, facilities, client_w, c, facility_pay_schedule, client_v, opening_cost):
    """
        Implementing the pruning procedure
    """
    primal_obj = 0
    assigned_clients = []
    # Candidate facilities to remain open
    temp_open_facilities = [i for i in range(len(facility_pay_schedule)) if facility_pay_schedule[i] == -1]
    opened_facilities = []

    while temp_open_facilities != []:
        i = temp_open_facilities[0]
        temp_open_facilities = temp_open_facilities[1:] # slice off current facility
        facility_assignment = [i, []]
        for j in range(len(client_v)):
            # How do I tell if client j ever contributed to facility i? Should I check if client_w[i][j] != -1?
            witness = client_v[j][1]
            # Hopefully client_w is up-to-date
            witness_cost = (client_v[j][0] - client_w[witness][j]) + c[witness][j]
            current_cost = (client_v[j][0] - client_w[i][j]) + c[i][j]
            if current_cost <= witness_cost and client_w[i][j] != -1:
                primal_obj += c[i][j]
                assigned_clients += [j]
                facility_assignment[1] += [clients[j]]
                temp_open_facilities = [h for h in temp_open_facilities if not (client_w[i][j] > 0 and client_w[h][j] > 0)]
        if not len(facility_assignment[1]) == 0:
            opened_facilities += [facility_assignment]

    # Continuing Method 2: assigned clients not yet connected to their nearest open facility
    extra_cost = 0
    for j in range(len(client_v)):
        if not j in assigned_clients:
            nearest_facility_list = [[c[opened_facilities[i][0]][j], i, opened_facilities[i][0]] for i in range(len(opened_facilities))]
            nearest_facility_sublist = min(nearest_facility_list)
            nearest_facility_index = nearest_facility_sublist[1]
            primal_obj += c[nearest_facility_sublist[2]][j]
            extra_cost += c[nearest_facility_sublist[2]][j]
            opened_facilities[nearest_facility_index][1] += [clients[j]]

    primal_obj += len(opened_facilities) * opening_cost
    print("Pruning Assigned Clients", len(assigned_clients))
    print("Pruning Open Facilities", len(opened_facilities))
    return primal_obj, opened_facilities, extra_cost

def plot_final_points(clients, client_assignments):
    """
        Here the clients are grouped by open facility
    """

    facility_centers = []
    for j in range(len(client_assignments)):
        plt.plot([clients[j][0]], [clients[j][1]], marker='o', color='red')
        # iterate over clients assigned to this facility
        if j == 0:
            for k in client_assignments[j]:
                facility_centers += [[clients[k[0]][0], clients[k[0]][1]]]
    # plot the facilities
    for i in facility_centers:
        plt.plot([i[0]], [i[1]], marker='s', color='black')
    plt.show()

def plot_final_clusters(open_facilities):
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



