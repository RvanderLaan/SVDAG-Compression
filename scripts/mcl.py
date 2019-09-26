import markov_clustering as mc
import networkx as nx
import matplotlib

file_in = './bin/edges10.txt'
edges = []
with open(file_in) as f:
    for line in f:
        n1, n2, w = line.split(' ')
        edges.append((int(n1), int(n2), float(w)))
        # edges.append((int(n2), int(n1), float(w)))

g = nx.Graph()
g.add_weighted_edges_from(edges)

# matrix = nx.to_numpy_matrix(g)
matrix = nx.to_scipy_sparse_matrix(g)

result = mc.run_mcl(matrix)           # run MCL with default parameters
clusters = mc.get_clusters(result)    # get clusters

mc.draw_graph(matrix, clusters, node_size=50, with_labels=False, edge_color="silver")
matplotlib.pyplot.show()
