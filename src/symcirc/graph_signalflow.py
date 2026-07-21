from collections import defaultdict
from random import sample
from symcirc import utils


class Branch:
    def __init__(self, name, src, dst, sym_gain, num_gain):
        self.name = name
        self.src = src
        self.dst = dst
        self.sym_gain = sym_gain
        self.num_gain = num_gain

    @property
    def label(self):
        return str(self.sym_gain)
    def nodes(self):
        return self.src, self.dst


from collections import defaultdict
import sympy


class SignalFlowGraph:

    def __init__(self, name):
        self.name = name
        self.vertices = set()
        self.edges = defaultdict(list)

    def add_edge(self, src, dst, sym_gain, num_gain, component):
        self.vertices.add(src)
        self.vertices.add(dst)
        self.edges[src].append((dst, sym_gain, num_gain, {component}))

    def merge_parallel(self):
        """Merge parallel edges by summing gains"""
        new_edges = defaultdict(dict)

        for src in self.edges:
            if src[-1] == "v":
                mode = "adm"
            else:
                mode = "imp"

            for dst, sym_gain, num_gain, components in self.edges[src]:
                if dst in new_edges[src]:
                    if mode == "adm":
                        new_edges[src][dst]["sym_gain"] += sym_gain
                        new_edges[src][dst]["num_gain"] += num_gain
                        new_edges[src][dst]["components"].union(components)
                    else:
                        new_imp = 1/(1/new_edges[src][dst]["sym_gain"]+1/sym_gain)
                        new_imp_num = 1/(1/new_edges[src][dst]["num_gain"]+1/num_gain)
                        new_edges[src][dst]["sym_gain"] = new_imp
                        new_edges[src][dst]["num_gain"] = new_imp_num
                        new_edges[src][dst]["components"] = components
                else:
                    new_edges[src][dst] = {"sym_gain": sym_gain, "num_gain": num_gain, "components": components}
        merged = defaultdict(list)

        for src in new_edges:
            for dst, data in new_edges[src].items():
                sym_gain = data["sym_gain"]
                num_gain = data["num_gain"]
                components = data["components"]
                if num_gain != 0:
                    merged[src].append((dst, sym_gain, num_gain, components))

        self.edges = merged

    def collapse_nodes(self, keep, remove):
        """Collapse node 'remove' into node 'keep'"""

        if keep == remove:
            return
        if remove == "0":
            remove = keep
            keep = "0"

        # 1. Redirect outgoing edges from 'remove'
        if remove in self.edges:
            for dst, sym_gain, num_gain, name in self.edges[remove]:
                new_dst = keep if dst == remove else dst
                self.edges[keep].append((new_dst, sym_gain, num_gain, name))

            del self.edges[remove]

        # 2. Redirect incoming edges to 'remove'
        for src in list(self.edges.keys()):
            new_list = []
            for dst, sym_gain, num_gain, name in self.edges[src]:
                if dst == remove:
                    new_list.append((keep, sym_gain, num_gain, name))
                else:
                    new_list.append((dst, sym_gain, num_gain, name))
            self.edges[src] = new_list

        # 3. Remove vertex
        self.remove_vertex(remove)

        # 4. Merge parallel edges (important!)
        self.merge_parallel()

    def remove_vertex(self, vertex="0"):
        """Remove reference node vertices"""
        ref_v = f"{vertex}v"
        ref_i = f"{vertex}i"

        self.vertices.discard(ref_v)
        self.vertices.discard(ref_i)

        if ref_v in self.edges:
            del self.edges[ref_v]
        if ref_i in self.edges:
            del self.edges[ref_i]

        for src in list(self.edges):
            self.edges[src] = [
                (dst, g, ng, n)
                for (dst, g, ng, n) in self.edges[src]
                if dst not in (ref_v, ref_i)
            ]

    def find_all_paths(self, start, end):
        """Return all simple paths with total symbolic and numeric gain"""
        paths = []

        def dfs(current, target, path, elements, sym_gain, num_gain, visited):
            path.append(current)
            visited.add(current)

            if current == target:
                paths.append({
                    "path": path.copy(),
                    "elements": elements.copy(),
                    "sym_gain": sym_gain,
                    "num_gain": num_gain
                })
            else:
                for neighbor, sg, ng, comps in self.edges.get(current, []):
                    if neighbor not in visited:
                        dfs(
                            neighbor,
                            target,
                            path,
                            elements.union(comps),  # log components here
                            sym_gain * sg,
                            num_gain * ng,
                            visited
                        )

            path.pop()
            visited.remove(current)

        dfs(start, end, [], set(), 1, 1, set())
        return paths

    def evaluate_effect(self, points):
        max_gain = 0
        for edges in self.edges.values():
            i = 0
            for edge in edges:
                ng = edge[2]
                num_gain_list = []
                for freq in points:
                    eval_num_gain = utils.evalf(ng, subs={utils.s: 2*utils.pi*freq*utils.j})
                    eval_num_gain_mag = utils.mag(eval_num_gain)
                    num_gain_list.append(eval_num_gain_mag)
                effect = sum(num_gain_list) / len(num_gain_list)
                if effect > max_gain:
                    max_gain = effect
                edges[i] = (edge[0], edge[1], effect, edge[3])
                i+=1
        return max_gain

    def find_most_impactful_path(self, start, end, points):
        """
        Tree search (DFS + branch-and-bound) for the single path between
        start and end with the largest-magnitude overall gain, i.e. the
        path of least impedance / maximal conductance.

        Unlike find_all_paths, this doesn't enumerate every simple path and
        compare afterward - it keeps a running best and prunes any branch
        whose best-case remaining contribution can no longer beat it.

        Returns a dict shaped like the entries from find_all_paths
        ({"path", "elements", "sym_gain", "num_gain"}) for the winning
        path, or None if start and end aren't connected.
        """
        # Upper bound on |gain| any single edge can contribute, used to
        # bound how much a partial path could still grow before pruning.
        max_abs_gain = self.evaluate_effect(points)

        print(max_abs_gain)

        best = {"path": None, "elements": None, "sym_gain": None, "num_gain": 0}

        def upper_bound(current_num_gain, visited):
            # A simple path can extend through at most one edge per
            # remaining unvisited vertex.
            remaining_hops = len(self.vertices) - len(visited)
            return abs(current_num_gain) * (max_abs_gain ** remaining_hops)

        def dfs(current, target, path, elements, sym_gain, num_gain, visited):
            if current == target:
                if abs(num_gain) > abs(best["num_gain"]):
                    best["path"] = path.copy()
                    best["elements"] = elements.copy()
                    best["sym_gain"] = sym_gain
                    best["num_gain"] = num_gain
                return

            # Prune: even at max possible remaining gain, this branch
            # can't beat the best path found so far.
            if upper_bound(num_gain, visited) <= abs(best["num_gain"]):
                return

            for neighbor, sg, ng, comps in self.edges.get(current, []):
                if neighbor not in visited:
                    path.append(neighbor)
                    visited.add(neighbor)

                    dfs(
                        neighbor,
                        target,
                        path,
                        elements.union(comps),
                        sym_gain * sg,
                        num_gain * ng,
                        visited,
                    )

                    path.pop()
                    visited.remove(neighbor)

        dfs(start, end, [start], set(), 1, 1, {start})

        return best if best["path"] is not None else None

    def visualize(self):
        import networkx as nx
        import matplotlib.pyplot as plt

        G = nx.DiGraph()

        for src in self.edges:
            for dst, gain, num_gain, _ in self.edges[src]:
                G.add_edge(src, dst, label=str(gain))

        #pos = nx.spring_layout(G)
        pos = nx.shell_layout(G)

        nx.draw_networkx_nodes(G, pos, node_size=200, node_color="lightblue")
        nx.draw_networkx_labels(G, pos, font_size=9)

        nx.draw_networkx_edges(G, pos, arrows=True)

        edge_labels = nx.get_edge_attributes(G, "label")

        nx.draw_networkx_edge_labels(G, pos, edge_labels=edge_labels, label_pos=0.6, font_size=6)

        plt.title(self.name)
        plt.axis("off")
        plt.show()

