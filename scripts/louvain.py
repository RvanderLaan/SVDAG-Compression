import community
from infomap import infomap
import networkx as nx
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors

color_keys = list(mcolors.TABLEAU_COLORS.keys())

file_in = './bin/edges12.txt'
edges = []
with open(file_in) as f:
    for line in f:
        n1, n2, w = line.split(' ')
        edges.append((int(n1), int(n2), float(w)))
        # edges.append((int(n2), int(n1), float(w)))

G = nx.Graph()
G.add_weighted_edges_from(edges)
# G.add_edge(1, 2, weight=1)

# nx.layou

#first compute the best partition
# louvain
partition = community.best_partition(G, resolution=0.2)

# dendogram = community.generate_dendrogram(G, resolution=0.2)
# partition = community.partition_at_level(dendogram, 1)

#drawing
size = float(len(set(partition.values())))
pos = nx.kamada_kawai_layout(G)

# infomap
infomapWrapper = infomap.Infomap("--two-level")

print("Building Infomap network from a NetworkX graph...")
for e in G.edges():
    infomapWrapper.addLink(*e)

print("Find communities with Infomap...")
infomapWrapper.run();

tree = infomapWrapper.tree

print("Found %d top modules with codelength: %f" % (tree.numTopModules(), tree.codelength()))

partition = {}
for node in tree.leafIter():
    partition[node] = node.moduleIndex()


count = 0.
for com in set(partition.values()) :
    count = count + 1.
    list_nodes = [nodes for nodes in partition.keys()
                                if partition[nodes] == com]
    nx.draw_networkx_nodes(G, pos, list_nodes, node_size = 60,
                                node_color = mcolors.TABLEAU_COLORS[color_keys[int(count) % len(color_keys)]])


nx.draw_networkx_edges(G, pos, alpha=0.5)
plt.show()
