import igraph as ig
import csv
import collections
import json
import random
import pickle
import sys

# Question 6
try:
    # Some workspace variables take long time to generate, so we store them to disk
    # and use the pre-generated object if available
    with open("g_gcc.pickle", 'rb') as file_object:
        g_gcc = pickle.load(file_object)
except:
    print("Generating Graph and its GCC...")
    with open("g_gcc.pickle", 'wb') as file_object:
        weighted_edges = collections.OrderedDict()  # key: (source node ID, destination node ID), value (sum of weights, number of edges)
        with open("san_francisco-censustracts-2017-4-All-MonthlyAggregate.csv", mode='r') as fgraph:
            reader = csv.DictReader(fgraph, delimiter=',')
            for index, line in enumerate(reader):
                if int(line['month']) != 12:
                    continue
                _sid, _did = int(line['sourceid']), int(line['dstid'])
                edge = (_sid, _did) if _sid < _did else (_did, _sid)
                if edge in weighted_edges:  # consider multiple edges between the same two nodes
                    total_weight, num_edges = weighted_edges[edge]
                    weighted_edges[edge] = (float(line['mean_travel_time']) + total_weight, num_edges + 1)
                else:
                    weighted_edges[edge] = (float(line['mean_travel_time']), 1)

        # merge duplicated edges by averaging their weights
        for edge in weighted_edges.keys():
            total_weight, num_edges = weighted_edges[edge]
            weighted_edges[edge] = total_weight / num_edges

        g = ig.Graph([e for e in weighted_edges.keys()], 
            edge_attrs=dict(weight=[w for w in weighted_edges.values()]))
        # set vertex indices attribute to keep track of indices in later graph manipulations
        for index, vertex in enumerate(g.vs):
            vertex['index'] = index
        print("The graph has {0} vertices and {1} edges.".format(g.vcount(), g.ecount()))
        g_gcc = g.clusters().giant()

        pickle.dump(g_gcc, file_object)
        print("Graph and its GCC generated...")

print("The Giant Connected Component has {0} vertices and {1} edges.".format(g_gcc.vcount(), g_gcc.ecount()))
g_gcc_indices_lookup = {vertex['index']:i for i, vertex in enumerate(g_gcc.vs)}

# Each node is represented by a Geolocation object, 
# with a street address as Display Name and a coordinate (a tuple of latitude and longitude) as Location
class Geolocation():
    def __init__(self, display_name, location):
        self.name= display_name
        self.location = location

geolocations = {}  # key: node ID, value: Geolocation object
with open("san_francisco_censustracts.json") as fgeo:
    geodata = json.load(fgeo)
    for feature in geodata["features"]:
        feature_id = feature["properties"]["MOVEMENT_ID"]
        display_name = feature["properties"]["DISPLAY_NAME"]
        coordinates = [point for l1 in feature["geometry"]["coordinates"] for l2 in l1 for point in l2]
        # average all coordinates of this node
        coordinate_sum = (0, 0)
        for coord in coordinates:
            coordinate_sum = (coordinate_sum[0]+coord[0], coordinate_sum[1]+coord[1])
        coordinate_avg = (coordinate_sum[0]/len(coordinates), coordinate_sum[1]/len(coordinates))

        geolocations[feature_id] = Geolocation(display_name, coordinate_avg)

# Question 7
def pretty_print_loc(loc):
    return (round(loc[0], 3), round(loc[1], 3))

weighted_edges_gcc = collections.OrderedDict()
for index, edge in enumerate(ig.EdgeSeq(g_gcc)):
    weighted_edges_gcc[index] = edge["weight"]
    
g_mst = g_gcc.spanning_tree(weights=list(weighted_edges_gcc.values()))
print("The Minimum Spanning Tree has {0} vertices and {1} edges.".format(g_mst.vcount(), g_mst.ecount()))
g_edge_seq = ig.EdgeSeq(g_mst)
for index, edge in enumerate(g_edge_seq):
    if index % 150 == 1:  # random sampling of edges
        print("One edge in the MST: {0}, with weight {1}".format(edge.tuple, edge["weight"]))
        print("\tSource node street address: {0}; location: {1}\n\tTarget node street address: {2}; location: {3}".format(
            geolocations[str(edge.tuple[0])].name, pretty_print_loc(geolocations[str(edge.tuple[0])].location), 
            geolocations[str(edge.tuple[1])].name, pretty_print_loc(geolocations[str(edge.tuple[1])].location)))
        # print("\"{0}\",\"{1}\",\"{2}\",\"{3}\",\"{4}\",\"{5}\"".format(edge.tuple, edge["weight"],
            # geolocations[str(edge.tuple[0])].name, pretty_print_loc(geolocations[str(edge.tuple[0])].location), 
            # geolocations[str(edge.tuple[1])].name, pretty_print_loc(geolocations[str(edge.tuple[1])].location)))

g_mst_indices_lookup = {vertex['index']:i for i, vertex in enumerate(g_mst.vs)}

# Question 8
try:
    with open("graph_triangles.pickle", 'rb') as file_object:
        graph_triangles = pickle.load(file_object)
except:
    print("Generating GCC graph cliques...")
    with open("graph_triangles.pickle", 'wb') as file_object:
        graph_triangles = list(g_gcc.cliques(min=3, max=3))
        pickle.dump(graph_triangles, file_object)
    print("GCC cliques generated.")

# a random sampling of 1000 triangles
random_indices, count = [], 0
while count < 1000:
    val = random.randint(0, len(graph_triangles))
    if val not in random_indices:
        random_indices.append(val)
        count += 1
print("Random sampling indices generated.")

def get_edge_weight(graph, edge_vertices):
    # Given edge_vertices as a tuple of two integers, return the edge's weight in the undirected graph
    weights = (graph.es.select(_source=edge_vertices[0], _target=edge_vertices[1])["weight"] +
               graph.es.select(_source=edge_vertices[1], _target=edge_vertices[0])["weight"])
    return sum(weights)
    
count_satisfied = 0
for cnt, index in enumerate(random_indices):
    if cnt % 10 == 1:
        print("Checking triangle #{0}".format(cnt))
    triangle = graph_triangles[index]
    triangle_weights = [get_edge_weight(g_gcc, (triangle[0], triangle[1])), 
                        get_edge_weight(g_gcc, (triangle[0], triangle[2])),
                        get_edge_weight(g_gcc, (triangle[1], triangle[2]))]
    triangle_weights.sort()
    if triangle_weights[0] + triangle_weights[1] > triangle_weights[2]:  # triangle inequality
        count_satisfied += 1

print("Percentage of triangles satisfying the triangle inequality: {0}".format(count_satisfied / 1000))

# Question 9
def find_cyclic_walk(graph, start_id, path):
    # Given a graph in adjacency matrix (2D list) and ID of starting node,
    # return the path of a cyclic walk
    path += [start_id]
        
    # if a cycle is found, simply return the path
    if len(path) > 1 and path[0] == start_id:
        return path
    else:
        _node_candidates = []  # adjacent nodes to starting node
        for next_id in range(len(graph)):
            if graph[start_id][next_id] > 0:
                _node_candidates += [next_id]
        # print("\t_node_candidates = {0}".format(_node_candidates))
        for next_id in _node_candidates:
            # copy the graph and path so far, and remove the traversed edge
            _graph_prune = [[item for item in row] for row in graph]
            _graph_prune[start_id][next_id] -= 1
            _graph_prune[next_id][start_id] -= 1
            _path_append = [node for node in path]

            path_res = find_cyclic_walk(_graph_prune, next_id, _path_append)
            if path_res != None:
                return path_res
        print("WARN: no cyclic walk found at node {0}".format(start_id))
        return None
    
def find_eulerian_walk(graph, start_id):
    # Given a graph in adjacency matrix and ID of starting node,
    # return the path of a walk that traverses every edge exactly once
    if graph[start_id] == [0] * len(graph):
        return [start_id]  # empty walk

    path = find_cyclic_walk(graph, start_id, [])
    print("Initial Eulerian Path found with length {0}, path = {1}".format(len(path), path))
    # remove edges used in the first cycle
    i_node = 0
    while i_node < len(path) - 1:
        graph[path[i_node]][path[i_node+1]] -= 1
        graph[path[i_node+1]][path[i_node]] -= 1
        i_node += 1

    # find nodes on the path that have existing edges incident to it,
    # find cycles from these nodes and add the sub-path into the whole path
    # repeat this process until the path does not change anymore
    while True:
        path_new = []
        for index, node_id in enumerate(path):
            if graph[node_id] == [0] * len(graph):  # no incident edge
                path_new += [node_id]
            else:
                subpath = find_cyclic_walk(graph, node_id, [])
                print("One more Eulerian Path found with length {0}, path = {1}".format(len(subpath), subpath))
                path_new += subpath
                # remove used edges in this new subpath
                i_node = 0
                while i_node < len(subpath) - 1:
                    graph[subpath[i_node]][subpath[i_node+1]] -= 1
                    graph[subpath[i_node+1]][subpath[i_node]] -= 1
                    i_node += 1
        if path_new == path:
            break
        else:
            path = path_new
    return path
    
g_mst_matrix = list(g_mst.get_adjacency())
g_mst_indices = [vertex['index'] for vertex in g_mst.vs]

g_mst_matrix = list(map(lambda row: list(map(lambda val: val*2, row)), g_mst_matrix))

try:
    with open("eulerian_path.pickle", 'rb') as file_object:
        eulerian_path = pickle.load(file_object)
except:
    with open("eulerian_path.pickle", 'wb') as file_object:
        print("Finding Eulerian Walk of multi-graph...")
        eulerian_path = find_eulerian_walk(g_mst_matrix, 0)
        eulerian_path = list(map(lambda i_node: g_mst_indices[i_node], eulerian_path))  # map vertices to their canonical indices
        print("Eulerian Walk found.")

        pickle.dump(eulerian_path, file_object)

print("Approximate Traveling Salesman Eulerian Path: \n{0}\nLength: {1}".format(eulerian_path, len(eulerian_path)))

try:
    with open("tour_path_adjacent.pickle", 'rb') as file_object_1:
        with open("tour_cost.pickle", 'rb') as file_object_2:
            tour_path_adjacent = pickle.load(file_object_1)
            tour_cost = pickle.load(file_object_2)
except:
    with open("tour_path_adjacent.pickle", 'wb') as file_object_1:
        with open("tour_cost.pickle", 'wb') as file_object_2:
            # To get the embedded tour from the Eulerian cycle, sequentially extract unique nodes from the path
            tour_path = []
            for node in eulerian_path:
                if node not in tour_path:
                    tour_path.append(node)
            # print("Embedded tour has length: {0}".format(len(tour_path)))

            # Calculate Approximate TSP Cost and the tour path
            tour_path_adjacent, tour_cost = [], 0
            for index in range(len(tour_path)-1):
                if index % 50 == 1:
                    print("Calculating tour cost at node {0}...".format(index))

                edge_cost = get_edge_weight(g_gcc, (
                    g_gcc_indices_lookup[tour_path[index]], g_gcc_indices_lookup[tour_path[index+1]]))
                if edge_cost == 0:  # no edge exists between these two nodes
                    _shortest_path = g_mst.get_shortest_paths(
                        g_mst_indices_lookup[tour_path[index]], to=g_mst_indices_lookup[tour_path[index+1]], output='vpath')[0]
                    for _i in range(len(_shortest_path)-1):
                        tour_path_adjacent.append(g_mst_indices[_shortest_path[_i]])
                        edge_cost += get_edge_weight(g_gcc, (
                            g_gcc_indices_lookup[g_mst_indices[_shortest_path[_i]]], 
                            g_gcc_indices_lookup[g_mst_indices[_shortest_path[_i+1]]]))
                else:
                    tour_path_adjacent.append(tour_path[index])
                tour_cost += edge_cost
            tour_path_adjacent.append(tour_path[len(tour_path)-1])

            pickle.dump(tour_path_adjacent, file_object_1)
            pickle.dump(tour_cost, file_object_2)

print("Traveling Salesman Tour: \n{0}\nCost: {1}".format(tour_path_adjacent, tour_cost))

# Generate a list of nodes represented as a tuple of latitude and longitude
try:
    with open("tour_path_adjacent_coord.pickle", 'rb') as file_object:
        tour_path_adjacent_coord = pickle.load(file_object)
except:
    with open("tour_path_adjacent_coord.pickle", 'wb') as file_object:
        tour_path_adjacent_coord = []
        for node_id in tour_path_adjacent:
            tour_path_adjacent_coord.append(geolocations[str(node_id)].location)
        pickle.dump(tour_path_adjacent_coord, file_object)

# print("Traveling Salesman Tour with Coordinates: \n{0}".format(tour_path_adjacent_coord))
with open("tour_path_coordinates.csv", 'w') as fout:
    for coordinate in tour_path_adjacent_coord:
        fout.write("{0},{1}\n".format(coordinate[0], coordinate[1]))

# Calculate the total cost of the MST, which is the lower bound of the optimal TSP
mst_cost = sum([edge["weight"] for edge in g_mst.es])
print("Minimum Spanning Tree Cost: {0}".format(mst_cost))
