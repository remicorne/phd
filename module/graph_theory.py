import networkx as nx
import matplotlib.pyplot as plt

#https://networkx.org/documentation/stable/tutorial.html
#constructors
G=nx.Graph()
# G=nx.DiGraph() #directed 
# G=nx.MultiGraph() #multiple edges between nodes
# G=nx.MultiDiGraph()

G.add_edge(1,2) #will create nodes 1 and 2 
G.add_edge(2,3,weight=0.7) #weighting will be correlations 
G.add_edge(1, "d") #can pass any object! 

nx.draw(G, with_labels=True, font_weight='bold')
plt.show