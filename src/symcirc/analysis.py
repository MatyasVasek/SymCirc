import os
import copy, time
from typing import List, Set, Dict, Union

from sympy import Symbol

from symcirc import parse, laplace
from symcirc.component import Component, CurrentSource, VoltageSource
from symcirc.utils import *


class Circuit:
    def __init__(self, netlist):
        self.netlist = netlist
        self.components, self.couplings = parse.parse(netlist)

    def get(self, component_name: str) -> Component:
        return self.components[component_name]

    def add(self, component: Component) -> None:
        if component.name in self.components:
            raise(ValueError('Component already exists'))
        else:
            self.components[component.name] = component

    def remove(self, component_name: str) -> None:
        if component_name in self.components:
            del self.components[component_name]
        else:
            raise(ValueError("Component doesn't exists"))

    def change(self, component_name: str, parameter: str, new_value) -> None:
        if component_name in self.components:
            c = self.components[component_name]
            setattr(c, parameter, new_value)
        else:
            raise(ValueError("Component doesn't exists"))

    def _scan_nodes(self) -> Set[str]:
        node_set = set()

        for c in self.components.values():
            node_set = node_set|c.nodes()
        return node_set

    def count_components(self) -> int:
        """
          Returns the total number of components in the circuit
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

    def get_node_dict(self) -> Dict[str, int]:
        nodes = self._scan_nodes()
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
        return len(self._scan_nodes()) - 1

    def analyse(self, analysis_type:str = "tf", method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):

        analysis_type = analysis_type.lower()
        if analysis_type == "dc":
            analysis = DC(self, method, symbolic, precision, sympy_ilt, use_symengine)
        elif analysis_type == "ac":
            analysis = AC(self, method, symbolic, precision, sympy_ilt, use_symengine)
        elif analysis_type == "tf":
            analysis = TF(self, method, symbolic, precision, sympy_ilt, use_symengine)
        elif analysis_type == "tran":
            analysis = TRAN(self, method, symbolic, precision, sympy_ilt, use_symengine)
        else:
            raise ValueError(f"Nonexistent analysis type: {analysis_type}")
        return analysis

    def dc(self):
        pass

    def ac(self):
        pass

    def tf(self):
        pass

    def tran(self):
        pass

    def pz(self):
        pass


class Analysis:
    """
    Main SymCirc class.
    When initialized it parses input netlist, and conducts the desired analysis which is then stored as a set
    of equations in the equation_matrix variable.

    :param Circuit circuit: Loaded netlist file which contains the circuit description.
        If you intend to load from a file, use the utils.load_file() function.
    :param str analysis_type: Analysis type identifier: "dc", "ac", "tf", "tran".
    :param bool symbolic: False if you want your results evaluated with numerical values from the netlist.

    :raise ValueError: If the analysis_type argument is invalid.
    """
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):

        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        method = method.lower()
        if method not in ["tableau", "two_graph_node", "modified_node"]:
            raise ValueError(f"Nonexistent analysis method: {method}")

        self.is_symbolic: bool = symbolic
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt
        self.circuit: Circuit = circuit
        self.node_voltage_identities: list = []

        self.symbol_dict = {}

        self.node_voltage_symbols: List[sympy.Symbol] = self.circuit.get_node_symbols()
        self.node_dict = self.circuit.get_node_dict()
        self.c_count = self.circuit.count_components()
        self.node_count = self.circuit.count_nodes()

        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}
        self.symbol_dict = self.get_symbols()  # format: {<symbol_name> : <Symbol>}

    def v(self, name: str) -> sympy.Expr:
        """
        Returns the specified voltage
        """
        symbol = self.get_symbols()[f"v({name})"]
        return self.solved_dict[symbol]

    def i(self, name: str) -> sympy.Expr:
        """
        Returns the specified current
        """
        symbol = self.get_symbols()[f"i({name})"]
        return self.solved_dict[symbol]

    def get_symbols(self) -> Dict[str, sympy.Symbol]:
        symbol_dict = {}
        for expr in self.solved_dict.values():
            free_symbols = expr.free_symbols
            for symbol in free_symbols:
                symbol_dict[symbol.name] = symbol
        return symbol_dict

    def get_node_voltage_symbol(self, node: str) -> sympy.Symbol:
        return self.node_voltage_symbols[self.node_dict[node]]

    def component_voltage(self, name: str) -> sympy.Symbol:
        """Has to be implemented in child class"""
        ret = None
        v = f"v({name})"

        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(v)]

        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            n1 = c.node1
            n2 = c.node2
            if n1 == "0":
                vn1 = 0
            else:
                vn1 = self.get_node_voltage(n1)
            if n2 == "0":
                vn2 = 0
            else:
                vn2 = self.get_node_voltage(n2)
            ret = sympy.cancel((vn1 - vn2))
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """Has to be implemented in child class"""
        pass

    def get_node_voltage(self, node: str, force_s_domain: bool=False) -> Union[sympy.Expr, None]:
        """Has to be implemented in child class"""
        func = None
        try:
            func = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return sympy.Expr(0)
                    else:
                        for n in identity:
                            try:
                                func = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        return func

    def _choose_source_val(self, c: Component) -> sympy.Expr:
        """Has to be implemented in child class"""
        pass

    def component_values(self, name: str = "all", default_python_datatypes: bool=False) -> Dict[str, sympy.Expr]:
        """
          Takes a string containing a single component name and returns a dictionary containing the voltage and current
            of the input component.

          :param str name: component id
          :param bool default_python_datatypes:
          :return dict ret: in format {"v(name)" : value, "i(name)" : value}
        """
        ret = {}
        if name == "all":
            ret = self.all_component_values(default_python_datatypes)
        else:
            v = f"v({name})"
            i = f"i({name})"
            ret[v] = self.component_voltage(name)
            ret[i] = self.component_current(name)
        return ret

    def all_component_values(self, default_python_datatypes=False) -> Dict[str, sympy.Expr]:
        """
          Returns a dictionary of all relevant voltages and currents in the circuit.

          :return dict ret: in format {"v(id1)" : value, "i(id2)" : value, ...}
        """
        ret = {}
        for key in self.circuit.components:
            if self.circuit.components[key].type == "k":
                pass
            elif self.circuit.components[key].name[-3:] == "_IC":
                pass
            else:
                if default_python_datatypes:
                    name = self.circuit.components[key].name
                    elem_dict = self.component_values(name)
                    try:
                        for elem_key in elem_dict:
                            ret[elem_key] = float(elem_dict[elem_key])
                    except TypeError:
                        ret.update(elem_dict)
                else:
                    name = self.circuit.components[key].name
                    elem_dict = self.component_values(name)
                    ret.update(elem_dict)
        return ret

    def node_voltages(self):
        """
          Returns a dictionary of all node voltages in the circuit.

          :return dict ret: in format {"v(node1)" : value, ...}
        """
        ret = {}
        for node in self.node_dict:
            if node.find("*ctrl") == -1:
                ret[f"v({node})"] = self.get_node_voltage(node)
            else: # Filter nodes added by virtual current sensors for CCC(V)S
                pass
        return ret

    def transfer_function(self, node1, node2):
        """
          Takes names of two nodes and returns their transfer function

          :param str node1: node id
          :param str node2: node id
          :return sympy_object ret: resulting transfer function
        """
        try:
            if node1 in [0, "0"]:
                voltage1 = 0
            else:
                v1 = sympy.Symbol(f"v({node1})")
                voltage1 = self.solved_dict[v1].simplify()
        except KeyError:
            print("Node {} doesn't exist.".format(node1))
            exit(101)

        try:
            if node2 in [0, "0"]:
                voltage2 = 0
            else:
                v2 = sympy.Symbol(f"v({node2})")
                voltage2 = self.solved_dict[v2].simplify()
        except KeyError:
            print("Node {} doesn't exist.".format(node2))
            exit(101)
        tf = (voltage2/voltage1)
        return tf

    def _analyse(self):
        """
          Implementation of all types of supported analysis.

            eqn_matrix, solved_dict, symbols
          :return sympy.Matrix eqn_matrix: matrix of the system equations
          :return dict solved_dict: dictionary of eqn_matrix solve results
          :return list symbols: list of all used sympy.symbol objects
        """
        eqn_matrix, symbols = self._build_system_eqn()
        t0 = time.time()
        '''solved_dict = sympy.solve_linear_system_LU(eqn_matrix, symbols)
        for i in solved_dict:
            solved_dict[i] = solved_dict[i].cancel()
        print(f"LU solve time: {time.time() - t0}")
        print(solved_dict)'''
        #t0 = time.time()
        solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
        #print(f"Gauss solve time: {time.time()-t0}")
        #print(solved_dict)
        return eqn_matrix, solved_dict, symbols

    def _build_system_eqn(self):
        if self.method == "tableau":
            equation_matrix, symbols = self._build_tableau()
        elif self.method == "two_graph_node":
            equation_matrix, symbols = self._build_tgn()
        return equation_matrix, symbols

    def _build_tableau(self):
        size = self.c_count
        M = sympy.Matrix(sympy.zeros(2 * size + self.node_count))
        R = sympy.Matrix(sympy.zeros(size * 2 + self.node_count, 1))
        for i in range(size):
            M[i, i] = 1
        index = 0
        node_symbols = self.node_voltage_symbols
        voltage_symbols = []
        current_symbols = []
        inductor_index = {}
        couplings = []
        for key in self.circuit.components:
            if self.circuit.components[key].type == "k":
                couplings.append(self.circuit.components[key])
            else:
                c = self.circuit.components[key]
                if c.type in ["r", "l", "c"]:
                    self._add_basic(M, R, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    if c.type == "l":
                        inductor_index[c.name] = index

                if c.type == "v":
                    self._add_voltage_source(M, R, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                if c.type == "i":
                    self._add_current_source(M, R, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                if c.type == "g":
                    self._add_VCT(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "e":
                    self._add_VVT(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "f":
                    self._add_CCT(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "h":
                    self._add_CVT(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "a":
                    self._add_A(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})_control"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})_control"))
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))
                    index += 1
                if c.type == "s":
                    self._add_short(M, c, index)
                    voltage_symbols.append(sympy.Symbol(f"v({c.name})"))
                    current_symbols.append(sympy.Symbol(f"i({c.name})"))

                index += 1
        for coupling in couplings:
            self._add_K(M, R, coupling, inductor_index)

        equation_matrix = M.col_insert(self.c_count * 2 + self.node_count, R)
        symbols = voltage_symbols + current_symbols + node_symbols
        return equation_matrix, symbols

    def _build_tgn(self):
        symbols = []
        v_graph_collapses = []
        i_graph_collapses = []
        v_graph_nodes = []
        i_graph_nodes = []
        matrix_col_expand = 0
        coupled_pairs = []
        coupled_inductors = []

        if len(self.circuit.couplings) > 0:
            for coupling in self.circuit.couplings:
                coupled_inductors.append(coupling.L1)
                coupled_inductors.append(coupling.L2)
                coupled_pairs.append([coupling.L1, coupling.L2, coupling])

        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type in ["r", "l", "c"]:
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)

                if c.type == "l" and c.coupling is not None:
                    matrix_col_expand += 1

            if c.type == "v":
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self._collapse(i_graph_collapses, c.node1, c.node2)

            if c.type == "i":
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)

            if c.type == "g":  # VCT/VCCS
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

            if c.type == "e":  # VVT/VCVS
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

                self._collapse(i_graph_collapses, c.node1, c.node2)

            if c.type == "f":  # CCT/CCCS
                node1 = c.node1
                node2 = c.node2
                if c.node3 is None:
                    c_v = self.circuit.components[c.current_sensor]
                    node3 = c_v.shorted_node
                else:
                    node3 = c.node3
                node4 = c.node4
                self.graph_append(node1, v_graph_nodes)
                self.graph_append(node2, v_graph_nodes)
                self.graph_append(node3, v_graph_nodes)
                self.graph_append(node4, v_graph_nodes)
                self.graph_append(node1, i_graph_nodes)
                self.graph_append(node2, i_graph_nodes)
                self.graph_append(node3, i_graph_nodes)
                self.graph_append(node4, i_graph_nodes)

                self._collapse(v_graph_collapses, node3, node4)
                matrix_col_expand += 1

            if c.type == "h":  # CVT/CCVS
                node1 = c.node1
                node2 = c.node2
                if c.node3 is None:
                    c_v = self.circuit.components[c.current_sensor]
                    node3 = c_v.shorted_node
                else:
                    node3 = c.node3
                node4 = c.node4
                self.graph_append(node1, v_graph_nodes)
                self.graph_append(node2, v_graph_nodes)
                self.graph_append(node3, v_graph_nodes)
                self.graph_append(node4, v_graph_nodes)
                self.graph_append(node1, i_graph_nodes)
                self.graph_append(node2, i_graph_nodes)
                self.graph_append(node3, i_graph_nodes)
                self.graph_append(node4, i_graph_nodes)

                matrix_col_expand += 1
                self._collapse(i_graph_collapses, node1, node2)
                self._collapse(v_graph_collapses, node3, node4)

            if c.type == "a":
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node3, v_graph_nodes)
                self.graph_append(c.node4, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self.graph_append(c.node3, i_graph_nodes)
                self.graph_append(c.node4, i_graph_nodes)

                self._collapse(i_graph_collapses, c.node1, c.node2)
                self._collapse(v_graph_collapses, c.node3, c.node4)

            if c.type == "s":
                pass

        """Collapse nodes based on collapse dictionaries"""
        for collapse_list in i_graph_collapses:
            i = 0
            tmp_i_graph_nodes = copy.copy(i_graph_nodes)
            for n in tmp_i_graph_nodes:
                if n in collapse_list:
                    if "0" in collapse_list:
                        i_graph_nodes.remove(n)
                    else:
                        i_graph_nodes[i] = min(collapse_list)
                i += 1

        for collapse_list in v_graph_collapses:
            i = 0
            tmp_v_graph_nodes = copy.copy(v_graph_nodes)
            for n in tmp_v_graph_nodes:
                if n in collapse_list:
                    if "0" in collapse_list:
                        v_graph_nodes.remove(n)
                    else:
                        v_graph_nodes[i] = min(collapse_list)
                i += 1
        self.node_voltage_identities = v_graph_collapses
        v_graph_nodes = list(set(v_graph_nodes))
        i_graph_nodes = list(set(i_graph_nodes))

        rows = len(i_graph_nodes)
        cols = len(v_graph_nodes)

        m_size = len(v_graph_nodes) + matrix_col_expand

        M = sympy.Matrix(sympy.zeros(m_size))
        S = sympy.Matrix(sympy.zeros(m_size, 1))
        index_row = 0
        index_col = 0
        symbols_to_append = []
        for key in self.circuit.components:
            c = self.circuit.components[key]
            if c.type == "l":
                if c.coupling is not None:
                    raise NotImplementedError(
                        "Coupled inductors not implemented for this method. Use tableau method for coupled inductors.")
                else:
                    self._add_basic_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type in ["r", "c"]:
                self._add_basic_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "v":
                self._add_voltage_source_tgn(M, S, v_graph_nodes, i_graph_nodes, c, index_row, i_graph_collapses,
                                             v_graph_collapses)
                index_row += 1
            if c.type == "i":
                self._add_current_source_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses)
            if c.type == "g":
                self._add_VCT_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "e":
                self._add_VVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_row, i_graph_collapses, v_graph_collapses)
                index_row += 1
            if c.type == "f":
                # raise NotImplementedError("CCCS not supported in 'two_graph_node' analysis")
                self._add_CCT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, i_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                index_col += 1
            if c.type == "h":
                self._add_CVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, index_row, i_graph_collapses,
                                  v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                index_col += 1
                index_row += 1
            if c.type == "a":
                pass
            if c.type == "s":
                raise NotImplementedError("Switch not supported in 'two_graph_node' analysis")

        equation_matrix = M.col_insert(m_size, S)

        for node in v_graph_nodes:
            symbols.append(sympy.Symbol(f"v({node})"))
        for symb in symbols_to_append:
            symbols.append(symb)

        # TODO: experiment with simplification inside matrix - seems like a huge performance upgrade in tgn method!
        for i in range(m_size ** 2):
            expr = equation_matrix[i]
            if expr != 0:
                equation_matrix[i] = sympy.cancel(expr)

        return equation_matrix, symbols

    def _collapse(self, graph_collapses, node1, node2):
        collapsed = False
        node1_in = None
        node2_in = None
        for i in range(len(graph_collapses)):
            if node1 in graph_collapses[i]:
                node1_in = i
            if node2 in graph_collapses[i]:
                node2_in = i
        if (node1_in is None) and (node2_in is None): # collapsed nodes not present in any existing collapse list
            graph_collapses.append([node1, node2])  # set node2 to be collapsed into node1 on the graph
        elif (node1_in is not None) and (node2_in is not None):
            graph_collapses[node2_in] = list(set(graph_collapses[node2_in]) | set(graph_collapses[node1_in]))
            del graph_collapses[node1_in]
        elif node1_in is not None:
            graph_collapses[node1_in].append(node2)  # set node2 to be collapsed into node1 on the graph
        elif node2_in is not None:
            graph_collapses[node2_in].append(node1)  # set node1 to be collapsed into node2 on the graph

    def _add_CVT_tgn(self, M, v_nodes, i_nodes, c, index_col, index_row, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            r = c.sym_value
        else:
            r = c.value

        node1 = c.node1
        node2 = c.node2
        if c.node3 is None:
            c_v = self.circuit.components[c.current_sensor]
            node3 = c_v.shorted_node
            c.node3 = node3
        else:
            node3 = c.node3
        node4 = c.node4

        c_v = self.circuit.components[c.current_sensor]

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)

        col = len(v_nodes) + index_col
        row = len(i_nodes) + index_row

        if n3i is not None:
            M[n3i, col] += 1
        if n4i is not None:
            M[n4i, col] += -1
        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1

        M[row, col] += -r

    def _add_CCT_tgn(self, M, v_nodes, i_nodes, c, index, i_graph_collapses):
        f_gain = None
        if self.is_symbolic:
            f_gain = c.sym_value
        else:
            f_gain = c.value
        node1 = c.node1
        node2 = c.node2
        if c.node3 is None:
            c_v = self.circuit.components[c.current_sensor]
            node3 = c_v.shorted_node
            c.node3 = node3
        else:
            node3 = c.node3
        node4 = c.node4

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)
        col = len(v_nodes) + index
        if n1i is not None:
            M[n1i, col] += f_gain
        if n2i is not None:
            M[n2i, col] += -f_gain
        if n3i is not None:
            M[n3i, col] += 1
        if n4i is not None:
            M[n4i, col] += -1

    def _add_VCT_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            g = c.sym_value
        else:
            g = c.value
        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4
        n1v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node4, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        if n1v is not None:
            if n1i is not None:
                M[n1i, n1v] += +g
            if n2i is not None:
                M[n2i, n1v] += -g
        if n2v is not None:
            if n1i is not None:
                M[n1i, n2v] += -g
            if n2i is not None:
                M[n2i, n2v] += +g

    def _add_VVT_tgn(self, M, v_nodes, i_nodes, c, index, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            e = c.sym_value
        else:
            e = c.value
        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4
        n1v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node4, v_graph_collapses)
        n3v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n4v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        row = len(i_nodes)+index
        if n1v is not None:
            M[row, n1v] += -e
        if n2v is not None:
            M[row, n2v] += e
        if n3v is not None:
            M[row, n3v] += 1
        if n4v is not None:
            M[row, n4v] += -1

    def _add_voltage_source_tgn(self, M, S, v_nodes, i_nodes, c, index, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        val = self._choose_source_val(c)

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        row = len(i_nodes)+index

        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1
        S[row, 0] += val

    def _add_current_source_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        val = self._choose_source_val(c)

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

        if n1i is not None:
            S[n1i, 0] += -val
        if n2i is not None:
            S[n2i, 0] += val

    def _add_basic_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        if c.type == "r":
            y = 1 / val
        if c.type == "l":
            y = 1 / (val*s)
        if c.type == "c":
            y = s * val
        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        if n1v is not None:
            if n1i is not None:
                M[n1i, n1v] += +y
            if n2i is not None:
                M[n2i, n1v] += -y
        if n2v is not None:
            if n1i is not None:
                M[n1i, n2v] += -y
            if n2i is not None:
                M[n2i, n2v] += +y

    @staticmethod
    def index_tgn(nodes, node, collapses):
        try:
            return nodes.index(node)
        except ValueError:
            for collapse_list in collapses:
                if node in collapse_list:
                    if "0" in collapse_list:
                        return None
                    else:
                        return nodes.index(min(collapse_list))
                elif node == "0":
                    return None
    @staticmethod
    def graph_append(node, graph):
        if node == '0':
            pass
        elif node not in graph:
            graph.append(node)
        return graph

    def _incidence_matrix_write(self, N1, N2, matrix, index, y_b=None):
        if N1 == "0":
            pass
        else:
            node_pos = self.node_dict[N1]
            if self.method == "tableau":
                matrix[self.c_count*2 + node_pos, self.c_count + index] += 1
                matrix[index, self.c_count*2 + node_pos] -= 1
        if N2 == "0":
            pass
        else:
            node_pos = self.node_dict[N2]
            if self.method == "tableau":
                matrix[self.c_count*2 + node_pos, self.c_count + index] += -1
                matrix[index, self.c_count*2 + node_pos] -= -1

        return matrix

    def _add_basic(self, matrix, vi_vector, c, index):
        N1 = c.node1
        N2 = c.node2
        y_b = 0
        z_b = 0
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value

        if c.type == "r":
            y_b = 1
            z_b = -val
        elif c.type == "l":
            y_b = 1
            z_b = -s * val
        elif c.type == "c":
            y_b = s * val
            z_b = -1

        matrix[self.c_count+index, index] += y_b
        matrix[self.c_count+index, self.c_count+index] += z_b
        self._incidence_matrix_write(N1, N2, matrix, index)
        return matrix

    def _add_voltage_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        val = self._choose_source_val(c)

        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        result[self.c_count + index, 0] = val
        return matrix

    def _add_current_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        val = self._choose_source_val(c)
        matrix[self.c_count + index, self.c_count + index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        result[self.c_count + index, 0] = val
        return matrix

    def _add_VCT(self, matrix, c, index):  # voltage to current transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[self.c_count + index, self.c_count + index] = 1
        matrix[self.c_count + index + 1, index] = val
        matrix[self.c_count + index + 1, self.c_count + index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index)
        self._incidence_matrix_write(N1, N2, matrix, index + 1)
        return matrix

    def _add_VVT(self, matrix, c, index):  # voltage to voltage transformer
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[self.c_count + index, self.c_count + index] = 1
        matrix[self.c_count + index + 1, index] = val
        matrix[self.c_count + index + 1, index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index)
        self._incidence_matrix_write(N1, N2, matrix, index + 1)
        return matrix

    def _add_CCT(self, matrix, c, index):  # current to current transformer
        N1 = c.node1
        N2 = c.node2
        if c.node3 is None:
            c_v = self.circuit.components[c.current_sensor]
            N3 = c_v.shorted_node
            c.node3 = N3
        else:
            N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[self.c_count + index, index] = 1
        matrix[self.c_count + index + 1, self.c_count + index] = val
        matrix[self.c_count + index + 1, self.c_count + index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index)
        self._incidence_matrix_write(N1, N2, matrix, index + 1)
        return matrix

    def _add_CVT(self, matrix, c, index):  # Current to voltage transformer
        N1 = c.node1
        N2 = c.node2
        if c.node3 is None:
            c_v = self.circuit.components[c.current_sensor]
            N3 = c_v.shorted_node
            c.node3 = N3
        else:
            N3 = c.node3
        N4 = c.node4

        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        matrix[self.c_count + index, index] = 1
        matrix[self.c_count + index + 1, self.c_count + index] = val
        matrix[self.c_count + index + 1, index + 1] = -1
        self._incidence_matrix_write(N3, N4, matrix, index)
        self._incidence_matrix_write(N1, N2, matrix, index + 1)
        return matrix

    def _add_A(self, matrix, c, index):  # Ideal OAMP
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        matrix[self.c_count + index, index] = 1
        matrix[self.c_count + index + 1, self.c_count + index] = 1
        self._incidence_matrix_write(N3, N4, matrix, index)
        self._incidence_matrix_write(N1, N2, matrix, index + 1)
        return matrix

    def _add_K(self, matrix, vi_vector, c, inductor_index):  # coupled inductors
        c_L1 = self.circuit.components[c.L1]
        c_L2 = self.circuit.components[c.L2]
        L1_index = inductor_index[c_L1.name]
        L2_index = inductor_index[c_L2.name]
        N1 = c_L1.node1
        N2 = c_L1.node2
        N3 = c_L2.node1
        N4 = c_L2.node2
        if self.is_symbolic:
            L1 = c_L1.sym_value
            L2 = c_L2.sym_value
            coupling_coeff = c.sym_value
            M = coupling_coeff*sympy.sqrt(L1*L2)
        else:
            L1 = c_L1.value
            L2 = c_L2.value
            coupling_coeff = c.value
            M = coupling_coeff * sympy.sqrt(L1 * L2)

        matrix[self.c_count + L2_index, self.c_count + L1_index] += -s*M
        matrix[self.c_count + L1_index, self.c_count + L2_index] += -s*M
        if c_L1.init_cond != None:
            vi_vector[self.c_count + L2_index, 0] += -M*c_L1.init_cond
        if c_L2.init_cond != None:
            vi_vector[self.c_count + L1_index, 0] += -M*c_L2.init_cond

    def _add_short(self, matrix, c, index):
        N1 = c.node1
        N2 = c.node2
        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        return matrix


class DC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):
        super().__init__(circuit, method, symbolic, precision, sympy_ilt, use_symengine)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        for sym in symbols:
            try:
                solved_dict[sym] = sympy.limit(solved_dict[sym], s, 0)
            except KeyError:
                pass
            except TypeError:
                pass
        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource]) -> sympy.Expr:
        if self.is_symbolic:
            val = c.dc_sym
        else:
            val = c.dc_num
        return val

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = super().component_voltage(name)
        v = f"v({name})"
        if ret is None:
            ret = sympy.Symbol(v)
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = None
        i = f"i({name})"
        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(i)]
        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            if c.type in ["r", "l", "c"]:
                n1 = c.node1
                n2 = c.node2
                if n1 == "0":
                    vn1 = 0
                else:
                    vn1 = self.get_node_voltage(n1)
                if n2 == "0":
                    vn2 = 0
                else:
                    vn2 = self.get_node_voltage(n2)

                if self.is_symbolic:
                    val = c.sym_value
                else:
                    val = c.value

                if c.type == "r":
                    ret = sympy.cancel((vn1 - vn2) / val)
                if c.type == "l":
                    ret = infinity
                if c.type == "c":
                    ret = 0

        if ret is None:
            ret =  sympy.Symbol(i)
        return ret


class TF(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):
        super().__init__(circuit, method, symbolic, precision, sympy_ilt, use_symengine)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource]) -> sympy.Expr:
        if self.is_symbolic:
            if c.ac_num == 0:
                val = 0
            else:
                val = c.ac_sym
        else:
            val = c.ac_num
        return val

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = None
        v = f"v({name})"

        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(v)]

        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            n1 = c.node1
            n2 = c.node2
            if n1 == "0":
                vn1 = 0
            else:
                vn1 = self.get_node_voltage(n1, force_s_domain=True)
            if n2 == "0":
                vn2 = 0
            else:
                vn2 = self.get_node_voltage(n2, force_s_domain=True)
            ret = sympy.cancel((vn1 - vn2))

        if ret is None:
            return sympy.Symbol(v)
        else:
            return ret

    def get_node_voltage(self, node: str, force_s_domain: bool=False) -> Union[sympy.Expr, None]:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        func = None
        try:
            func = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return sympy.Expr(0)
                    else:
                        for n in identity:
                            try:
                                func = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        return func

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = None
        i = f"i({name})"
        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(i)]
        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            if c.type in ["r", "l", "c"]:
                n1 = c.node1
                n2 = c.node2
                if n1 == "0":
                    vn1 = 0
                else:
                    vn1 = self.get_node_voltage(n1)
                if n2 == "0":
                    vn2 = 0
                else:
                    vn2 = self.get_node_voltage(n2)

                if self.is_symbolic:
                    val = c.sym_value
                else:
                    val = c.value

                if c.type == "r":
                    ret = sympy.cancel((vn1 - vn2) / val)
                if c.type == "l":
                    ret = sympy.cancel((vn1 - vn2) / (val * s))
                if c.type == "c":
                    ret = sympy.cancel((vn1 - vn2) * val* s)
        if ret is None:
            return sympy.Symbol(i)
        else:
            return ret

class AC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):
        super().__init__(circuit, method, symbolic, precision, sympy_ilt, use_symengine)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()

        if self.is_symbolic:
            for sym in symbols:
                try:
                    solved_dict[sym] = solved_dict[sym].subs(s, 2 * pi * f * j)
                except KeyError:
                    pass
        else:
            for sym in symbols:
                try:
                    solved_dict[sym] = solved_dict[sym].subs(s, 2 * pi * f * j)
                except KeyError:
                    pass

        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource]) -> sympy.Expr:
        if self.is_symbolic:
            if c.ac_num == 0:
                val = 0
            else:
                val = c.ac_sym
        else:
            val = c.ac_num * c.ac_phase
        return val

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = None
        v = f"v({name})"

        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(v)]

        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            n1 = c.node1
            n2 = c.node2
            if n1 == "0":
                vn1 = 0
            else:
                vn1 = self.get_node_voltage(n1)
            if n2 == "0":
                vn2 = 0
            else:
                vn2 = self.get_node_voltage(n2)
            ret = sympy.cancel((vn1 - vn2))
        if ret is None:
            return sympy.Symbol(v)
        else:
            return ret

    def get_node_voltage(self, node: str, force_s_domain: bool=False) -> Union[sympy.Expr, None]:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        func = None
        try:
            func = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return sympy.Expr(0)
                    else:
                        for n in identity:
                            try:
                                func = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        return func

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = None
        i = f"i({name})"
        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(i)]
        elif self.method == "two_graph_node":
            c = self.circuit.components[name]
            if c.type in ["r", "l", "c"]:
                n1 = c.node1
                n2 = c.node2
                if n1 == "0":
                    vn1 = 0
                else:
                    vn1 = self.get_node_voltage(n1)
                if n2 == "0":
                    vn2 = 0
                else:
                    vn2 = self.get_node_voltage(n2)

                if self.is_symbolic:
                    val = c.sym_value
                else:
                    val = c.value

                if c.type == "r":
                    ret = sympy.cancel((vn1 - vn2) / val)
                if c.type == "l":
                    ret = sympy.cancel((vn1 - vn2) / (val * 2 * pi * f * j))
                if c.type == "c":
                    ret = sympy.cancel((vn1 - vn2) * (val * 2 * pi * f * j))
        if ret is None:
            return sympy.Symbol(i)
        else:
            return ret

class TRAN(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):
        super().__init__(circuit, method, symbolic, precision, sympy_ilt, use_symengine)

    def get_node_voltage(self, node: str, force_s_domain: bool=False) -> Union[sympy.Expr, None]:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        func = super().get_node_voltage(node, force_s_domain)
        if func is None:
            return func
        else:
            res = laplace.iLT(func, self.sympy_ilt)
            res = sympy.factor_terms(res)
            return res

    def _add_basic(self, matrix, vi_vector, c, index):
        matrix = super()._add_basic(matrix, vi_vector, c, index)
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        if c.type == "l":
            vi_vector[self.c_count + index, 0] += -val*c.init_cond
        if c.type == "c":
            vi_vector[self.c_count + index, 0] += val*c.init_cond
        return matrix

    def _choose_source_val(self, c: Component) -> sympy.Expr:
        if self.is_symbolic:
            val = c.tran_sym
        else:
            val = c.tran_num
        return val

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = super().component_voltage(name)
        v = f"v({name})"

        if ret is None:
            return laplace.iLT(sympy.Symbol(v), self.sympy_ilt)
        else:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
            return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        i = f"i({name})"

        if ret is None:
            return laplace.iLT(sympy.Symbol(i), self.sympy_ilt)
        else:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
            return ret


def AnalyseCircuit(netlist: str, analysis_type: str = "DC", method: str = "tableau",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False) -> Analysis:

    circuit = Circuit(netlist)
    analysis_type = analysis_type.lower()
    if analysis_type == "dc":
        analysis = DC(circuit, method, symbolic, precision, sympy_ilt, use_symengine)
    elif analysis_type == "ac":
        analysis = AC(circuit, method, symbolic, precision, sympy_ilt, use_symengine)
    elif analysis_type == "tf":
        analysis = TF(circuit, method, symbolic, precision, sympy_ilt, use_symengine)
    elif analysis_type == "tran":
        analysis = TRAN(circuit, method, symbolic, precision, sympy_ilt, use_symengine)
    else:
        raise ValueError(f"Nonexistent analysis type: {analysis_type}")
    return analysis