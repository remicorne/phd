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

import os
import networkx as nx
import matplotlib.pyplot as plt

ROOT = os.getcwd()  # This gives terminal location (terminal working dir)
INPUT_DIR = f"{ROOT}/input"
OUTPUT_DIR = f"{ROOT}/output"
CACHE_DIR = f"{INPUT_DIR}/cache"

print (INPUT_DIR)


filename = "TCB2_data_HPLC.csv"  # TCB2 #using current working directory plus file name
HT_filename = "TCB2_data_HT.csv"

from main import *
from module.correlogram import buildExperimentCorrMatricies
filename,experiment,correlogram_type,to_correlate,p_value_threshold,n_minimum,columns,corr_method,from_scratch=corrSelector( # generic prompter for creating corr matricies, probably need to add pearson/spearman
    filename,
    experiment='dose_response',
    correlogram_type='compound',
    to_correlate='GLU',
    p_value_threshold=0.05,
    n_minimum=5,
    columns=["OF","PL","CC", "IC","M", "SJ","SL1", "SL6", "SR6", "SR1", "AC", "V",  
                "Am", "dH", "vH", "NAc", "VM", "DM","VL", "DL", "MD",  "VPL",  "VPR", 
                "DG", "Y",  "SC","SN", "VTA", "DR","MR", "CE"],
    corr_method = 'pearson',
    from_scratch=True,
    # hierarchical_clustering=None #need to firgue out correlogram plotting and how it will intergrate 
)

# df_to_corr , correlation_matrix , T_F_significance_matrix , treatment , correlated
matricies=buildExperimentCorrMatricies(
    filename,
    experiment,
    correlogram_type, #compound / ratio
    to_correlate, # whihc compound/s / ratio/s
    p_value_threshold,
    n_minimum,
    columns, #ordered to plot
    corr_method, # 'pearson' 'spearman' 'kendall'
    from_scratch,  # Used in decorator
)


for matrix in matricies:
    print (f" treatment {matrix[3]} correlating {matrix[4]} ")

    df_to_corr=matrix[0]
    correlation_matrix=matrix[1]
    T_F_significance_matrix=matrix[2]
    treatment=matrix[3]
    correlated=matrix[4]

# TODO need to have consistent placement of nodes for visual comparison 
#AllenSDK MouseConnectivityApi

# TODO metrics from graphs
# 1. Node-Level Metrics:
# Degree Centrality: Number of edges connected to a node. High degree nodes might be more influential.
# Betweenness Centrality: Nodes that act as bridges between other nodes in the graph.
# Closeness Centrality: How close a node is to all other nodes in the graph.
# 2. Graph-Level Metrics:
# Average Degree: Average number of edges connected to nodes in the graph.
# Graph Density: Ratio of actual edges to possible edges in the graph.
# Average Shortest Path Length: Average shortest path length between nodes in the graph.
# 3. Community Structure:
# Modularity: Measure of the strength of the division of a network into communities.
# Community Detection: Algorithms to identify groups of densely connected nodes.
# 4. Network Comparison:
# Graph Edit Distance: Quantifies dissimilarity between graphs based on operations needed to transform one graph into another.
# Jaccard Index: Measure of similarity between graphs based on shared edges.


    # Create a new graph
    G = nx.Graph()

    # Iterate through the correlation matrix and significance matrix
    for i, col in enumerate(correlation_matrix.columns):
        for j, row in enumerate(correlation_matrix.index):
            if i < j:
                correlation = correlation_matrix.iloc[j, i]  # As correlation_matrix is symmetric
                significance = T_F_significance_matrix.iloc[j, i]

                # Check for significance and add edge if significant
                if significance:
                    if correlation > 0:
                        edge_color = 'red'  # Positive correlation
                    else:
                        edge_color = 'blue'  # Negative correlation

                    # Add edge to the graph with edge weight and color
                    G.add_edge(col, row, weight=abs(correlation), color=edge_color)

    # Draw the graph
    pos = nx.spring_layout(G)  # You can use other layout algorithms as well
    edges = G.edges()
    weights = [G[u][v]['weight'] for u, v in edges]
    edge_colors = [G[u][v]['color'] for u, v in edges]

    # Draw nodes and edges
    nx.draw_networkx_nodes(G, pos, node_size=300, alpha=0.7)
    nx.draw_networkx_edges(G, pos, edgelist=edges, width=weights, edge_color=edge_colors)

    # Add labels to nodes
    node_labels = {node: node for node in G.nodes()}  # Label nodes with their names
    nx.draw_networkx_labels(G, pos, labels=node_labels, font_size=8)

    # Set title for the graph
    plt.title(f"{correlated} correlations in {treatment}")
    # Display the graph
    plt.show()




#single graph multiple edge types corrisponding to each compound correlation DA-DA, DA-GLU


########## COMPARE GRAPHS

#graph metrics : average edge weight, average degree, clustering coefficient, degree distribution, average shortest path length (stats where possible)

#network alignment : GRAAL or NetAlign




