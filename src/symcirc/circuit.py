from typing import List, Union, Any
from collections import defaultdict

from symcirc.component import Component, ParallelAdmittance, SerialAdmittance
from symcirc import parse
from symcirc.graph import Graph
from symcirc.utils import *

from sympy import Symbol


class Circuit:
    """
    When initialized it parses the input netlist, and translates it into a list of components as defined
    in 'symcirc.component'. This class implements methods for further circuit manipulation before analysis.

    :param netlist: string in netlist format which contains the circuit description.
        Note: If you intend to load from a file, use the utils.load_file() function.
    :param operating_point: Use to pass numeric or symbolic operating point values
        for linearized components, if no values are specified the default model with symbolic values is used.
        Example: operating_point = {"Q1": {"gm": 74.5e-3, "gpi": 232e-6, "gmu": 1e-9, "go": 22.5e-6, "gx": 1.66}}
            The Q1 bjt will use these values in the model and expand the basic model accordingly.
    """
    def __init__(self, netlist: str, operating_point: Union[Dict[str, Any], None]=None):
        self.netlist = netlist
        self.nodes = None
        self.adjacency_map = defaultdict(set)
        self.v_graph = None
        self.i_graph = None
        self.components, self.couplings = parse.parse(netlist, operating_point)
        #for component in self.components.values():
        #    print(component.value)
        self._build_graph()

    def __repr__(self):
        body = ",\n    ".join(
            f"{name}={value!r}"
            for name, value in vars(self).items()
        )
        return f"{self.__class__.__name__}(\n    {body}\n)"

    def _build_graph(self):
        self.adjacency_map = defaultdict(set)

        for c in self.components.values():
            nodes = list(c.nodes())

            # Connect all nodes of the component together
            for i in range(len(nodes)):
                for j in range(i + 1, len(nodes)):
                    n1 = nodes[i]
                    n2 = nodes[j]

                    self.adjacency_map[n1].add(n2)
                    self.adjacency_map[n2].add(n1)

        self.nodes = set(self.adjacency_map.keys())

    def build_graph(self, name):
        g = Graph(name)
        for comp in self.components.values():
            if len(comp.nodes()) > 2:
                comp.generate_ports()
                for port in comp.ports():
                    g.add_component(port)
            else:
                g.add_component(comp)
        return g

    def get(self, component_name: str) -> Component:
        return self.components[component_name]

    def add(self, component: Component) -> None:
        if component.name in self.components:
            raise(ValueError('Component already exists'))
        else:
            self.components[component.name] = component
            self._build_graph()

    def delete(self, component_name: str) -> None:
        if component_name in self.components:
            del self.components[component_name]
            self._build_graph()
        else:
            raise(ValueError("Component doesn't exists"))

    def pop(self, component_name: str) -> Component:
        if component_name in self.components:
            return self.components.pop(component_name)
            self._build_graph()
        else:
            raise(ValueError("Component doesn't exists"))

    def change(self, component_name: str, parameter: str, new_value) -> None:
        if component_name not in self.components:
            raise ValueError("Component doesn't exist")
        c = self.components[component_name]
        if hasattr(c, parameter):
            setattr(c, parameter, new_value)
            self._build_graph()
        else:
            raise AttributeError(f"Component '{component_name}' has no attribute '{parameter}'")

    '''def _scan_nodes(self) -> Set[str]:
        node_set = set()

        for c in self.components.values():
            node_set = node_set|c.nodes()
        return node_set'''

    def count_ports(self) -> int:
        """
          Returns the total number of ports in the circuit
          :return int count
        """
        count = 0
        for c in self.components:
            if self.components[c].type in ["a", "e", "g", "f", "h"]:
                count += 2
            elif self.components[c].type in ["k"]:
                pass
            else:
                count += 1
        return count

    def count_components(self) -> int:
        """
          Returns the total number of components in the circuit
          :return int count
        """
        return len(self.components)

    def get_node_dict(self) -> Dict[str, int]:
        nodes = self.nodes
        node_dict = {}
        i = 0
        grounded = False
        for node in nodes:
            if node != "0":
                node_dict[node] = i
                i += 1
            if node == "0":
                grounded = True
        if not grounded:
            print("Circuit not grounded")
        return node_dict

    def get_node_symbols(self) -> List[Symbol]:
        voltage_symbol_list = []
        for node in self.get_node_dict():
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol(f"v({node})"))
        return voltage_symbol_list

    def count_nodes(self) -> int:
        return len(self.nodes) - 1

    def lump_parallel_impedances(self):
        comps = self.components
        lumped_comps = set()
        lumps = set()
        for comp1 in comps.values():
            lump = None
            if comp1 in lumped_comps:
                continue
            elif comp1.type in ["r", "l", "c"]:
                n1 = comp1.node1
                n2 = comp1.node2
                node_set = {n1, n2}
                for comp2 in comps.values():
                    if comp2 is comp1 or comp2 in lumped_comps:
                        continue
                    elif comp2.type in ["r", "l", "c"]:
                        if node_set == {comp2.node1, comp2.node2}:
                            lumped_comps.add(comp2)
                            if lump is None:
                                lumped_comps.add(comp1)
                                lump = ParallelAdmittance(comp1)
                            else:
                                lump.add(comp2)
                if lump is not None:
                    lumps.add(lump)
        for comp in lumped_comps:
            self.pop(comp.name)
        for lump in lumps:
            self.add(lump)

        self._build_graph()

    def lump_serial_impedances(self, fixed_nodes=None):
        g = self.build_graph("graph")
        lumped_comps = []
        lumps = []

        if fixed_nodes is None:
            fixed_nodes = []

        for node in self.nodes:
            if node in fixed_nodes:
                continue
            adjacent = self.adjacency_map[node]
            if len(adjacent) == 2:
                can_be_lumped = True
                comps_to_be_lumped = []
                for n in adjacent:
                    comps_between = g.components_between(node, n)
                    if len(comps_between) == 1:
                        if comps_between[0].type not in ["r", "l", "c"]:
                            can_be_lumped = False
                            break
                        comps_to_be_lumped.append(comps_between[0])
                    else:
                        can_be_lumped = False

                if can_be_lumped:
                    #print(f"LUMPABLE NODE: {node}->{adjacent}")
                    c1 = comps_to_be_lumped[0]
                    c2 = comps_to_be_lumped[1]
                    #print(f"COMPS TO BE LUMPED:\n{c1.node2 if c1.node1 == node else c1.node1}---{c1.name}---{node}---{c2.name}---{c2.node2 if c2.node1 == node else c2.node1}")
                    lumped_comps.append(c1)
                    lumped_comps.append(c2)
                    lump = SerialAdmittance(c1)
                    lump.add(c2)
                    lumps.append(lump)

        for comp in lumped_comps:
            self.pop(comp.name)
        for lump in lumps:
            self.add(lump)
        self._build_graph()

    def lump_impedances(self, fixed_nodes=None):
        self.lump_parallel_impedances()
        self.lump_serial_impedances(fixed_nodes)

    def analyse(self, analysis_type:str = "tf", method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool = False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):
        from symcirc.analysis import DC, AC, TF, TRAN
        analysis_type = analysis_type.lower()
        if analysis_type == "dc":
            analysis = DC(self, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine)
        elif analysis_type == "ac":
            analysis = AC(self, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine)
        elif analysis_type == "tf":
            analysis = TF(self, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine)
        elif analysis_type == "tran":
            analysis = TRAN(self, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine)
        else:
            raise ValueError(f"Nonexistent analysis type: {analysis_type}")
        return analysis