import os
import copy
import time
import sympy
from typing import Dict, List
from symcirc import parse, laplace, utils
from symcirc.utils import j,s,t
from symcirc.pole_zero import *
from symcirc.component import Component, Coupling
from sympy.parsing.sympy_parser import T
from sympy.parsing.sympy_parser import (_token_splittable,
                                            standard_transformations, implicit_multiplication,
                                            split_symbols_custom)

class AnalyseCircuit:
    """
    Main SymCirc class.
    When initialized it parses input netlist, and conducts the desired analysis which is then stored as a set
    of equations in the equation_matrix variable.

    :param str netlist: Loaded netlist file which contains the circuit description.
        If you intend to load from a file, use the utils.load_file() function.
    :param str analysis_type: Analysis type identifier: "DC", "AC", "TF", "tran".
    :param bool symbolic: False if you want your results evaluated with numerical values from the netlist.

    :raise ValueError: If the analysis_type argument is invalid.
    """
    def __init__(self, netlist: str, analysis_type: str = "DC", method: str = "tableau", phases: str = "undefined",
                 symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True, use_symengine: bool = False):

        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        if analysis_type not in ["DC", "AC", "TF", "tran"]:
            raise ValueError(f"Nonexistent analysis type: {analysis_type}")

        if method not in ["tableau", "two_graph_node"]:
            raise ValueError(f"Nonexistent analysis method: {method}")

        self.is_symbolic: bool = symbolic
        self.analysis_type: str = analysis_type
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt
        self.netlist: str = netlist
        self.node_voltage_identities: list = []

        if analysis_type == "tran":
            data = parse.parse(netlist, tran=True)
        else:
            data = parse.parse(netlist)

        self.phases = self.parse_phases(phases)

        self.components: Dict[str, Component] = data["components"]   # {<name> : <Component>} (see component.py)
        self.node_dict: Dict[str, int] = data["node_dict"]  # {<node_name>: <index in equation matrix>, ...}
        self.node_count: int = data["node_count"]
        self.couplings: List[Coupling] = data["couplings"]

        self.c_count: int = self.count_components()  # Amount of components
        self.node_voltage_symbols: List[sympy.Symbol] = self._node_voltage_symbols()

        self.eqn_matrix: sympy.Matrix
        self.solved_dict: Dict[sympy.Symbol, sympy.Expr]
        self.symbols: List[sympy.Symbol]
        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}

        self.symbol_dict: Dict[str, sympy.Symbol] = self.generate_symbol_dict()  # format: {<symbol_name> : <Symbol>}

    def v(self, name: str) -> sympy.Expr:
        """
        Returns the specified voltage
        """
        symbol = self.symbol_dict[f"v({name})"]
        return self.solved_dict[symbol]

    def i(self, name: str) -> sympy.Expr:
        """
        Returns the specified current
        """
        symbol = self.symbol_dict[f"i({name})"]
        return self.solved_dict[symbol]

    def generate_symbol_dict(self) -> Dict[str, sympy.Symbol]:
        symbol_dict = {}
        for symbol in self.symbols:
            symbol_dict[symbol.name] = symbol
        return symbol_dict

    def parse_phases(self, phases):
        if phases != "undefined":
            phase_definition = []
            phase_def_syntax_error = ("Invalid phase definition syntax, use 'P=integer' or 'P=[...]' "
                                      "where the list contains lengths of phases written a fraction (the fractions must add up to 1)")
            if phases.startswith("P="):
                phases = phases.replace("P=", "")
                if phases.startswith("[") and phases.endswith("]"):
                    phases = phases.replace("[", "")
                    phases = phases.replace("]", "")
                    phase_definition = sympy.sympify(phases.split(','))
                    phase_sum = sum(phase_definition)
                    if phase_sum != 1:
                        raise ValueError("The sum of phase lengths must be 1")
                    else:
                        phase_definition.insert(0, len(phase_definition))
                elif type(int(phases)) == int:
                    if int(phases) < 2:
                        raise ValueError("The number of phases can't be less than 2")
                    else:
                        phase_definition.append(int(phases[0]))
                        for i in range(phase_definition[0]):
                            phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                else:
                    raise SyntaxError(phase_def_syntax_error)

            else:
                raise SyntaxError(phase_def_syntax_error)
        return phases

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = None
        v = f"v({name})"

        if self.method == "tableau":
            ret = self.solved_dict[sympy.symbols(v)]

        if self.method == "eliminated_tableau":
            if name[0] in ["r", "R", "c", "C", "l", "L"]:
                i = f"i({name})"
                c = self.components[name]
                if self.is_symbolic:
                    if name[0] in ["c", "C"]:
                        impedance = s/c.sym_value
                    else:
                        impedance = c.sym_value
                else:
                    if name[0] in ["c", "C"]:
                        impedance = s/c.value
                    else:
                        impedance = c.value
                ret = self.solved_dict[sympy.symbols(i)]*impedance

        if self.method == "two_graph_node":
            c = self.components[name]
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
        if self.analysis_type == "tran":
            if ret is None:
                return laplace.iLT(sympy.Symbol(v), self.sympy_ilt)
            else:
                return laplace.iLT(ret, self.sympy_ilt)
        elif self.analysis_type == "DC":
            if ret is None:
                return sympy.Symbol(v)
            else:
                try:
                    ret = sympy.limit(ret, s, 0)
                    return ret
                except KeyError:
                    return ret
        else:
            if ret is None:
                return sympy.Symbol(v)
            else:
                return ret

    def get_node_voltage_symbol(self, node: str) -> sympy.Symbol:
        return self.node_voltage_symbols[self.node_dict[node]]

    def get_node_voltage(self, node: str, force_s_domain: bool=False) -> sympy.Expr:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        F = None
        try:
            F = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return 0
                    else:
                        for n in identity:
                            try:
                                F = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        if force_s_domain:
            return F
        elif self.analysis_type == "tran":
            if F is None:
                return F
            else:
                return laplace.iLT(F, self.sympy_ilt)
        else:
            return F

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = None
        i = f"i({name})"
        if self.method in ["tableau", "eliminated_tableau"]:
            ret = self.solved_dict[sympy.symbols(i)]
        if self.method == "two_graph_node":
            c = self.components[name]
            if c.type in ["r", "l", "c"]:
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

                if self.is_symbolic:
                    val = c.sym_value
                else:
                    val = c.value

                if c.type == "r":
                    ret = sympy.cancel((vn1 - vn2) / val)
                if c.type == "l":
                    ret = sympy.cancel((vn1 - vn2) / (val * s))
                if c.type == "c":
                    ret = sympy.cancel((vn1 - vn2) * val*s)

        if self.analysis_type == "tran":
            if ret is None:
                return laplace.iLT(sympy.Symbol(i), self.sympy_ilt)
            else:
                return laplace.iLT(ret, self.sympy_ilt)
        elif self.analysis_type == "DC":
            if ret is None:
                return sympy.Symbol(i)
            else:
                try:
                    ret = sympy.limit(ret, s, 0)
                    return ret
                except KeyError:
                    return ret
        else:
            if ret is None:
                return sympy.Symbol(i)
            else:
                return ret

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
            return ret
        else:
            v = f"v({name})"
            i = f"i({name})"
            ret[v] = self.component_voltage(name)
            ret[i] = self.component_current(name)
        return ret

    def all_component_values(self, default_python_datatypes=False):
        """
          Returns a dictionary of all relevant voltages and currents in the circuit.

          :return dict ret: in format {"v(id1)" : value, "i(id2)" : value, ...}
        """
        ret = {}
        for key in self.components:
            if self.components[key].type == "k":
                pass
            elif self.components[key].name[-3:] == "_IC":
                pass
            else:
                if default_python_datatypes:
                    name = self.components[key].name
                    elem_dict = self.component_values(name)
                    try:
                        for key in elem_dict:
                            ret[key] = float(elem_dict[key])
                    except TypeError:
                        ret.update(elem_dict)
                else:
                    name = self.components[key].name
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
            ret[f"v({node})"] = self.get_node_voltage(node)
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

    def count_components(self):
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

    def _analyse(self):
        """
          Implementation of all types of supported analysis.

            eqn_matrix, solved_dict, symbols
          :return sympy.Matrix eqn_matrix: matrix of the system equations
          :return dict solved_dict: dictionary of eqn_matrix solve results
          :return list symbols: list of all used sympy.symbol objects
        """
        eqn_matrix, symbols = self._build_system_eqn()
        solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)

        if self.analysis_type == "DC":
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = sympy.limit(solved_dict[sym], s, 0)
                    except KeyError:
                        pass
                    except TypeError:
                        pass
            else:
                for sym in symbols:
                    try:
                        solved_dict[sym] = sympy.limit(solved_dict[sym], s, 0)
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["v", "i"]:
                                if c.dc_value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.dc_value)
                            else:
                                if c.value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            #if self.method != "two_graph_node":
                            solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                    except KeyError:
                        pass

        elif self.analysis_type == "AC":
            f = sympy.symbols("f", real=True, positive=True)
            if self.is_symbolic:
                i = 1
                for sym in symbols:
                    try:
                        solved_dict[sym] = solved_dict[sym].subs(s, 2 * sympy.pi * f * j)
                    except KeyError:
                        pass
            else:
                for sym in symbols:
                    try:
                        solved_dict[sym] = solved_dict[sym].subs(s, 2 * sympy.pi * f * j)
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["v", "i"]:
                                if c.ac_value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                            else:
                                if c.value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            #if self.method != "two_graph_node":
                            solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                    except KeyError:
                        pass

        elif self.analysis_type == "TF":
            if not self.is_symbolic:
                for sym in symbols:
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["v", "i"]:
                                if c.ac_value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                            else:
                                if c.value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            #if self.method != "two_graph_node":
                            solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                    except KeyError:
                        pass

        elif self.analysis_type == "tran":
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["i", "v"] and c.name[-3:] != "_IC":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.tran_value)

                    except KeyError:
                        pass
                    try:
                        for name in self.components:
                            c = self.components[name]
                        pass
                    except KeyError:
                        pass
            else:
                for sym in symbols:
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type in ["v", "i"]:
                                if c.tran_value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.tran_value)
                            else:
                                if c.value:
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                    except KeyError:
                        pass

        return eqn_matrix, solved_dict, symbols

    def _node_voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol(f"v({node})"))
        return voltage_symbol_list

    def _build_system_eqn(self):
        if self.method == "tableau":
            size = self.c_count
            M = sympy.Matrix(sympy.zeros(2 * size + self.node_count))
            R = sympy.Matrix(sympy.zeros(size*2 + self.node_count, 1))
            for i in range(size):
                M[i, i] = 1
            index = 0
            node_symbols = self.node_voltage_symbols
            voltage_symbols = []
            current_symbols = []
            inductor_index = {}
            couplings = []
            for key in self.components:
                if self.components[key].type == "k":
                    couplings.append(self.components[key])
                else:
                    c = self.components[key]
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

            equation_matrix = M.col_insert(self.c_count*2 + self.node_count, R)
            symbols = voltage_symbols + current_symbols + node_symbols

        elif self.method == "eliminated_tableau":
            size = self.c_count
            M = sympy.Matrix(sympy.zeros(size + self.node_count))
            R = sympy.Matrix(sympy.zeros(size + self.node_count, 1))
            index = 0
            node_symbols = self.node_voltage_symbols
            voltage_symbols = []
            current_symbols = []
            inductor_index = {}
            couplings = []

            for key in self.components:
                if self.components[key].type == "k":
                    pass
                else:
                    c = self.components[key]
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
                index += 1

            equation_matrix = M.col_insert(self.c_count + self.node_count, R)
            symbols = node_symbols + current_symbols

        elif self.method == "two_graph_node":
            symbols = []
            v_graph_collapses = []
            i_graph_collapses = []
            v_graph_nodes = []
            i_graph_nodes = []
            matrix_col_expand = 0
            coupled_pairs = []
            coupled_inductors = []

            if len(self.couplings) > 0:
                for coupling in self.couplings:
                    coupled_inductors.append(coupling.L1)
                    coupled_inductors.append(coupling.L2)
                    coupled_pairs.append([coupling.L1, coupling.L2, coupling])

            for key in self.components:
                c = self.components[key]
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
                    self.collapse(i_graph_collapses, c.node1, c.node2)

                if c.type == "i":
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)

                if c.type == "g":   # VCT/VCCS
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node3, v_graph_nodes)
                    self.graph_append(c.node4, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)
                    self.graph_append(c.node3, i_graph_nodes)
                    self.graph_append(c.node4, i_graph_nodes)

                if c.type == "e":   # VVT/VCVS
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node3, v_graph_nodes)
                    self.graph_append(c.node4, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)
                    self.graph_append(c.node3, i_graph_nodes)
                    self.graph_append(c.node4, i_graph_nodes)

                    self.collapse(i_graph_collapses, c.node1, c.node2)

                if c.type == "f":   # CCT/CCCS
                    node1 = c.node1
                    node2 = c.node2
                    c_v = self.components[c.control_voltage]
                    node3 = c_v.node2
                    node4 = c_v.shorted_node
                    self.graph_append(node1, v_graph_nodes)
                    self.graph_append(node2, v_graph_nodes)
                    self.graph_append(node3, v_graph_nodes)
                    self.graph_append(node4, v_graph_nodes)
                    self.graph_append(node1, i_graph_nodes)
                    self.graph_append(node2, i_graph_nodes)
                    self.graph_append(node3, i_graph_nodes)
                    self.graph_append(node4, i_graph_nodes)

                    self.collapse(v_graph_collapses, node3, node4)
                    matrix_col_expand += 1

                if c.type == "h":   # CVT/CCVS
                    node1 = c.node1
                    node2 = c.node2
                    c_v = self.components[c.control_voltage]
                    node3 = c_v.node2
                    node4 = c_v.shorted_node
                    self.graph_append(node1, v_graph_nodes)
                    self.graph_append(node2, v_graph_nodes)
                    self.graph_append(node3, v_graph_nodes)
                    self.graph_append(node4, v_graph_nodes)
                    self.graph_append(node1, i_graph_nodes)
                    self.graph_append(node2, i_graph_nodes)
                    self.graph_append(node3, i_graph_nodes)
                    self.graph_append(node4, i_graph_nodes)

                    self.collapse(v_graph_collapses, node1, node2)
                    matrix_col_expand += 1
                    self.collapse(i_graph_collapses, node3, node4)

                if c.type == "a":
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node3, v_graph_nodes)
                    self.graph_append(c.node4, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)
                    self.graph_append(c.node3, i_graph_nodes)
                    self.graph_append(c.node4, i_graph_nodes)

                    self.collapse(i_graph_collapses, c.node1, c.node2)
                    self.collapse(v_graph_collapses, c.node3, c.node4)

            if c.type == "s":
                    pass

            """Collapse nodes based on collapse dictioanaries"""
            for collapse_list in i_graph_collapses:
                i = 0
                tmp_i_graph_nodes = copy.copy(i_graph_nodes)
                for n in tmp_i_graph_nodes:
                    if n in collapse_list:
                        if "0" in collapse_list:
                            i_graph_nodes.remove(n)
                        else:
                            i_graph_nodes[i] = min(collapse_list)
                    i+=1

            for collapse_list in v_graph_collapses:
                i = 0
                tmp_v_graph_nodes = copy.copy(v_graph_nodes)
                for n in tmp_v_graph_nodes:
                    if n in collapse_list:
                        if "0" in collapse_list:
                            v_graph_nodes.remove(n)
                        else:
                            v_graph_nodes[i] = min(collapse_list)
                    i+=1
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
            for key in self.components:
                c = self.components[key]
                if c.type == "l":
                    if c.coupling is not None:
                        raise NotImplementedError("Coupled inductors not implemented for this method. Use tableau method for coupled inductors.")
                    else:
                        self._add_basic_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)

                if c.type in ["r", "c"]:
                    self._add_basic_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)

                if c.type == "v":
                    self._add_V_tgn(M, S, v_graph_nodes, i_graph_nodes, c, index_row, i_graph_collapses, v_graph_collapses)
                    index_row += 1

                if c.type == "i":
                    self._add_I_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses)
                if c.type == "g":
                    self._add_VCT_tgn(M, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
                if c.type == "e":
                    self._add_VVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_row, i_graph_collapses, v_graph_collapses)
                    index_row += 1
                if c.type == "f":
                    #raise TypeError("CCCS not supported in 'two_graph_node' analysis")
                    self._add_CCT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, i_graph_collapses)
                    symbols_to_append.append(sympy.Symbol(f"i({c.control_voltage})"))
                    index_col += 1
                if c.type == "h":
                    self._add_CVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, index_row, i_graph_collapses, v_graph_collapses)
                    symbols_to_append.append(sympy.Symbol(f"i({c.control_voltage})"))
                    index_col += 1
                    index_row += 1
                if c.type == "a":
                    pass
                if c.type == "s":
                    raise TypeError("Switch not supported in 'two_graph_node' analysis")

            equation_matrix = M.col_insert(m_size, S)

            for node in v_graph_nodes:
                symbols.append(sympy.Symbol(f"v({node})"))
            for symb in symbols_to_append:
                symbols.append(symb)

        return equation_matrix, symbols

    def collapse(self, graph_collapses, node1, node2):
        collapsed = False
        for collapse_list in graph_collapses:
            if node1 in collapse_list:
                collapsed = True
                collapse_list.append(node2)  # set node2 to be collapsed into node1 on the current graph
            elif node2 in collapse_list:
                collapsed = True
                collapse_list.append(node1)  # set node1 to be collapsed into node2 on the current graph
        if not collapsed:
            graph_collapses.append([node1, node2])  # set node2 to be collapsed into node1 on the current graph


    def _add_CVT_tgn(self, M, v_nodes, i_nodes, c, index_col, index_row, i_graph_collapses, v_graph_collapses):
        if self.is_symbolic:
            r = c.sym_value
        else:
            r = c.value
        node1 = c.node1
        node2 = c.node2
        c_v = self.components[c.control_voltage]
        node3 = c_v.node2
        node4 = c_v.shorted_node
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        n3v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n4v = self.index_tgn(v_nodes, node4, v_graph_collapses)
        col = len(v_nodes) + index_col
        row = len(i_nodes) + index_row

        if n1i is not None:
            M[n1i, col] += 1
        if n2i is not None:
            M[n2i, col] += -1
        if n3v is not None:
            M[row, n3v] += 1
        if n4v is not None:
            M[row, n4v] += -1

        M[row, col] += -r


    def _add_CCT_tgn(self, M, v_nodes, i_nodes, c, index, i_graph_collapses):
        f = None
        try:
            if self.is_symbolic:
                f = c.sym_value
            else:
                f = c.value
            node1 = c.node1
            node2 = c.node2
            c_v = self.components[c.control_voltage]
            node3 = c_v.node2
            node4 = c_v.shorted_node
            n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
            n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
            n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
            n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)
            col = len(v_nodes) + index
            if n1i is not None:
                M[n1i, col] += f
            if n2i is not None:
                M[n2i, col] += -f
            if n3i is not None:
                M[n3i, col] += 1
            if n4i is not None:
                M[n4i, col] += -1
        except:
            print(f)
            print(type(f))


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



    def _add_V_tgn(self, M, S, v_nodes, i_nodes, c, index, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.analysis_type == "tran":
            val = c.tran_value
        elif self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            else:
                val = c.ac_value

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        row = len(i_nodes)+index

        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1
        S[row, 0] += val

    def _add_I_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.analysis_type == "tran":
            val = c.tran_value
        elif self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            else:
                val = c.ac_value
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

    def index_tgn(self, nodes, node, collapses):
        #i=None
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


    def graph_append(self, node, graph):
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
            elif self.method == "eliminated_tableau":
                matrix[self.c_count + node_pos, self.node_count + index] += 1
                matrix[index, node_pos] += y_b
        if N2 == "0":
            pass
        else:
            node_pos = self.node_dict[N2]
            if self.method == "tableau":
                matrix[self.c_count*2 + node_pos, self.c_count + index] += -1
                matrix[index, self.c_count*2 + node_pos] -= -1
            elif self.method == "eliminated_tableau":
                matrix[self.c_count + node_pos, self.node_count + index] += -1
                matrix[index, node_pos] += -y_b

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

        if self.method == "tableau":
            matrix[self.c_count+index, index] += y_b
            matrix[self.c_count+index, self.c_count+index] += z_b
            self._incidence_matrix_write(N1, N2, matrix, index)
            if c.type == "l" and self.analysis_type == "tran":
                vi_vector[self.c_count + index, 0] += -val*c.init_cond
            if c.type == "c" and self.analysis_type == "tran":
                vi_vector[self.c_count + index, 0] += val*c.init_cond
        elif self.method == "eliminated_tableau":
            matrix[index, self.node_count + index] += z_b
            self._incidence_matrix_write(N1, N2, matrix, index, y_b=y_b)
            #for i in range
            #matrix[index, node_pos] += val
            # matrix[]
            # matrix[]
            #matrix[self.c_count + index, index] += y_b

        return matrix

    def _add_voltage_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value
        if self.method == "tableau":
            matrix[self.c_count + index, index] = 1
            self._incidence_matrix_write(N1, N2, matrix, index)
            result[self.c_count + index, 0] = val
        elif self.method == "eliminated_tableau":
            self._incidence_matrix_write(N1, N2, matrix, index, y_b=1)
            result[index, 0] = val

        return matrix

    def _add_current_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value
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
        c_v = self.components[c.control_voltage]
        N3 = c_v.node2
        N4 = c_v.shorted_node
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
        c_v = self.components[c.control_voltage]
        N3 = c_v.node2
        N4 = c_v.shorted_node
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
        c_L1 = self.components[c.L1]
        c_L2 = self.components[c.L2]
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




