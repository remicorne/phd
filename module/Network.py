import numpy as np
import networkx as nx
from module.Matrix import Matrix
from dataclasses import dataclass


@dataclass
class Network:
    matrix: Matrix

    def __post_init__(self):
        self.G = nx.MultiDiGraph() if self.matrix.is_square else nx.Graph()
        # directed edge -  to_correlate[0] --> to_correlate[1]
        self.G.clear()

        self.G.add_nodes_from(
            self.matrix.corr_masked.columns.tolist()
        )  # adds every BR as a node
        self.edge_labels = {}
        for (row, col), correlation in self.matrix.corr_masked.stack().dropna().items():
            # Add edge to the graph with edge weight and color
            # Avoid self sorrelation
            if not(row == col and not self.is_bidirectional):
                self.G.add_edge(
                    row,
                    col,
                    weight=abs(correlation),
                    color="red" if correlation > 0 else "blue",
                )
                
                self.edge_labels[(row, col)] = f"{correlation:.2f}"


            ##### Draw the graph
        # pos = nx.spring_layout(G, seed=42)  # using a seed for consistency need allensdk working

        angles = np.linspace(
            0, 2 * np.pi, len(self.matrix.corr_masked.columns), endpoint=False
        )
        self.pos = {
            col: (np.cos(angles[i]), np.sin(angles[i]))
            for i, col in enumerate(self.matrix.corr_masked.columns)
        }

    @property
    def is_bidirectional(self):
        return self.matrix.is_square
    
    def get_title(self):
        title = self.matrix.get_title()
        return title.replace('-', '->') if self.is_bidirectional else title

    def plot_ax(self, ax):
        nx.draw_networkx_nodes(
            self.G,
            self.pos,
            node_size=1100,
            alpha=0.95,
            node_color="white",
            edgecolors="black",
            ax=ax,
        )
        nx.draw_networkx_edges(
            self.G,
            self.pos,
            width=list(nx.get_edge_attributes(self.G, "weight").values()),
            edge_color=list(nx.get_edge_attributes(self.G, "color").values()),
            ax=ax,
            node_size=1100,
            **(
                {"arrowstyle": "->", "arrowsize": 20}
                if self.is_bidirectional
                else {}
            ),
        )
        # Add labels to nodes
        node_labels = {
            node: node for node in self.G.nodes()
        }  # Label nodes with their names
        nx.draw_networkx_labels(
            self.G, self.pos, labels=node_labels, font_size=18, ax=ax
        )    
        # nx.draw_networkx_edge_labels(
        #     self.G, self.pos, edge_labels=self.edge_labels, font_size=18, ax=ax
        # )

        # Set title for the graph
        ax.set_frame_on(False)
        ax.set_title(self.get_title(), fontsize=28, pad=-10, y=1)
        self.matrix.is_square
