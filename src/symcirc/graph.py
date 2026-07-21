from collections import defaultdict
from random import sample

class Graph:
    def __init__(self, name = "Circuit Graph"):
        self.adj = defaultdict(dict)   # node -> list of components
        self.nodes = set()
        self.components = dict()
        self.identities = []
        self.name = name

    @property
    def seed_node(self):
        return sample(self.nodes, 1)[0]

    def add_component(self, component):
        nodes = component.nodes()

        # connect every node pair
        for n in nodes:
            self.nodes.add(n)

        # store component for each node
        for n in nodes:
            self.adj[n][component.name] = component

        self.components[component.name] = component

    def neighbors(self, node):
        """Return nodes connected to given node."""
        nbrs = set()
        for comp in self.adj[node].values():
            for n in comp.nodes():
                if n != node:
                    nbrs.add(n)
        return nbrs

    def components_between(self, node1, node2):
        """Return all components connecting node1 and node2."""
        adj1 = self.adj[node1]
        adj2 = self.adj[node2]

        # iterate over the smaller adjacency dict
        if len(adj1) > len(adj2):
            node1, node2 = node2, node1
            adj1 = adj2

        result = []
        for comp in adj1.values():
            nodes = comp.nodes()
            if node2 in nodes:
                result.append(comp)

        return result

    def incident(self, node):
        """Return components connected to node."""
        return self.adj[node].values()

    def collapse_nodes(self, keep, remove):
        for comp in self.adj[remove].values():
            for attr in ["node1", "node2", "node3", "node4"]:
                if hasattr(comp, attr):
                    if getattr(comp, attr) == remove:
                        setattr(comp, attr, keep)

            self.adj[keep][comp.name] = comp

        del self.adj[remove]

        if remove in self.nodes:
            self.nodes.remove(remove)

    def remove_component(self, component):
        """Remove a component (edge) from the graph."""

        if component.name not in self.components:
            return

        # remove from adjacency lists
        for node in component.nodes():
            if node in self.adj:
                if component.name in self.adj[node]:
                    del self.adj[node][component.name]

        # remove from component set
        del self.components[component.name]

    def parallel_edges(self):
        """
        Return all sets of components (edges) between nodes as a list of lists.
        Each inner list contains all components connecting a specific pair of nodes.
        """
        from collections import defaultdict

        # key: frozenset({node1, node2}) -> list of components
        parallel_dict = defaultdict(list)

        for comp in self.components.values():
            nodes = list(comp.nodes())
            # consider all node pairs in the component
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    key = frozenset([nodes[i], nodes[j]])
                    parallel_dict[key].append(comp)

        # return just the lists of parallel components
        parallel_groups = list(parallel_dict.values())
        return parallel_groups

    def visualize(self):
        import networkx as nx
        import matplotlib.pyplot as plt

        G = nx.MultiGraph()

        # add nodes
        for node in self.nodes:
            G.add_node(node)

        # add edges
        for comp in self.components.values():
            nodes = list(comp.nodes())
            if len(nodes) == 2:
                n1, n2 = nodes
                G.add_edge(n1, n2, label=str(comp.name))
            else:
                for i in range(len(nodes)):
                    for j in range(i + 1, len(nodes)):
                        G.add_edge(nodes[i], nodes[j], label=str(comp.name))

        # layout with minimal spacing
        pos = nx.spring_layout(G, seed=42, k=0.8, iterations=100)
        #pos = nx.shell_layout(G)
        #pos = nx.kamada_kawai_layout(G, scale=2)

        nx.draw_networkx_nodes(G, pos, node_size=1200, node_color="lightblue")
        nx.draw_networkx_labels(G, pos)

        # draw edges with alternating curvature
        edges = list(G.edges(keys=True))
        for u, v, k in edges:
            sign = 1 if k % 2 == 0 else -1  # alternate up/down
            idx = (k + 1) // 2  # magnitude of bulge
            nx.draw_networkx_edges(
                G,
                pos,
                edgelist=[(u, v)],
                connectionstyle=f"arc3,rad={sign * 0.15 * idx}",
                arrows=True
            )

        # draw edge labels correctly along curves
        ax = plt.gca()
        for u, v, k, data in G.edges(keys=True, data=True):
            label = data["label"]
            x1, y1 = pos[u]
            x2, y2 = pos[v]

            # compute label offset matching the curvature
            sign = 1 if k % 2 == 0 else -1
            idx = (k + 1) // 2
            offset = 0.08 * idx * sign

            x = (x1 + x2) / 2
            y = (y1 + y2) / 2 + offset

            ax.text(
                x,
                y,
                label,
                fontsize=9,
                horizontalalignment="center",
                verticalalignment="center",
                bbox=dict(facecolor="white", edgecolor="none", alpha=0.5),
            )

        plt.title(self.name)
        plt.axis("off")
        plt.show()
