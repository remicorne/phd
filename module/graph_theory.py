import networkx as nx
import matplotlib.pyplot as plt
import numpy as np

#https://networkx.org/documentation/stable/tutorial.html
#constructors
# G=nx.Graph()
# G=nx.DiGraph() #directed 
G=nx.MultiGraph() #multiple edges between nodes
# G=nx.MultiDiGraph()

edge_list = [(1,2),(2,4),(4,1)]
G=nx.from_edgelist(edge_list)
plt.show


G=nx.Graph()
G.add_edges_from(edge_list)


G.add_edge(1,2) #will create nodes 1 and 2 
G.add_edge(2,3,weight=0.7) #weighting will be correlations 
G.add_edge(1, "d") #can pass any object! 

nx.draw(G, with_labels=True, font_weight='bold')
plt.show

#adjency matricies - can they be done with weighting ?





########## GENERATE GRAPHS

#multiple graphs for each correlation type i.e. DA-DA

#single graph multiple edge types corrisponding to each compound correlation DA-DA, DA-GLU


########## COMPARE GRAPHS

#graph metrics : average edge weight, average degree, clustering coefficient, degree distribution, average shortest path length (stats where possible)

#network alignment : GRAAL or NetAlign




