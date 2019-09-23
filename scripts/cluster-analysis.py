import numpy as np
import matplotlib.pyplot as plt

filepath = '/home/remi/Documents/Projects/SymVox/bin/edges11.txt'

num_edges = 0

edges = {}
nodes_with_single_edge = {}


with open(filepath) as f:
  p1, p2, pw = None, None, 0
  for line in f:
    num_edges += 1
    n1, n2, nw = line.split(' ')

    if n1 not in edges:
      edges[n1] = [n2]
    else:
      edges[n1].append(n2)
    
    if p1 is not None and p1 is not n1:
      nodes_with_single_edge[n1] = n2

    p1, p2, pw, = n1, n2, nw

num_node_pairs = 0
num_node_triplets = 0

for n in nodes_with_single_edge:
  n2 = nodes_with_single_edge[n]
  if n2 not in edges:
    num_node_pairs += 1

for n in edges:
  n_edges = edges[n]
  if len(n_edges) is 2:
    n2, n3 = n_edges

    num_triplet_edges = 2
    if n2 in edges:
      num_triplet_edges += len(edges[n2])
    if n3 in edges:
      num_triplet_edges += len(edges[n3]) 

    if num_triplet_edges is 3:
      num_node_triplets += 1


subgraph_sizes = []
num_subgraphs = 0
while len(edges) > 0:
  n = next(iter(edges))
  queue = edges.pop(n)
  num_subgraphs += 1

  subgraph_nodes = {n}
  subgraph_nodes.update(queue)

  while len(queue) > 0:
    n2 = queue.pop()
    if n2 in edges:
      subgraph_nodes.update(edges[n2])
      queue.extend(edges.pop(n2))
  subgraph_sizes.append(len(subgraph_nodes))


print('Nodes: {} edges: {} nodes w/ single edge: {} pairs: {} triplets: {}'.format(len(edges), num_edges, len(nodes_with_single_edge), num_node_pairs, num_node_triplets))
print('Subgraphs: {}'.format(num_subgraphs))

print("Smallest subgraph: {} \t Largest subgraph: {}".format(min(subgraph_sizes), max(subgraph_sizes)))

hist, bins = np.histogram(subgraph_sizes, bins=32)

width = 0.7 * (bins[1] - bins[0])
center = (bins[:-1] + bins[1:]) / 2

fig, ax = plt.subplots()
ax.bar(center, hist, align='center', width=width)
ax.set_yscale('log')
plt.show()
