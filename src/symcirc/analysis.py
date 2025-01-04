import os
import copy
import time
import operator
import sympy
from typing import Dict, List
from symcirc import parse, laplace, utils
from symcirc.pole_zero import *
from symcirc.component import Component, Coupling
from symcirc.utils import j, t, z


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
                 scsi: str = "undefined", symbolic: bool = True, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False):

        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        if analysis_type not in ["DC", "AC", "TF", "tran"]:
            raise ValueError(f"Nonexistent analysis type: {analysis_type}")

        if method not in ["tableau", "two_graph_node", "modified_node"]:
            raise ValueError(f"Nonexistent analysis method: {method}")

        self.is_symbolic: bool = symbolic
        self.analysis_type: str = analysis_type
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt
        self.netlist: str = netlist
        self.node_voltage_identities: list = []

        data = parse.parse(netlist, analysis_type)

        self.phases, self.frequency = self.parse_phases(phases)
        self.scsi = scsi
        self.symbol_dict = {}
        self.SCSI_initial_values()

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

    def SCSI_initial_values(self):
        if self.phases != "undefined" and self.analysis_type not in ["AC", "TF", "tran"]:
            raise ValueError("Valid analysis types for SCSI circuits are: 'AC', 'TF' and 'tran'")
        if self.phases != "undefined" and self.method not in ["modified_node", "two_graph_node"]:
            raise ValueError("Valid methods for SCSI circuits are: 'modified_node' and 'two_graph_node'")
        if self.phases != "undefined" and self.scsi not in ["scideal", "siideal"]:
            raise ValueError("Please specify if circuit is in switched capacitor ('scideal') mode or switched current mode ('siideal')")

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

    def parse_phases(self, phases): # used by SCSI analysis
        frequency = 1
        if phases != "undefined":
            phase_definition = []
            phase_def_syntax_error = \
                    (
                    "Invalid phase definition syntax, use 'P=integer', P=(integer,frequency),'P=[...]' or P=([...],frequency), \n"
                    "where the integer is the number of phases "
                    "and the list contains lengths of phases written as a fraction (fractions must add up to 1)")
            if phases.startswith("P="):
                phases = phases.replace("P=", "")
                if phases.startswith("(") and phases.endswith(")"):
                    phases = phases.replace("(", "")
                    phases = phases.replace(")", "")
                    if phases.startswith("[") and "]" in phases:
                        phases = phases.replace("[", "")
                        phases = phases.replace("]", "")
                        phase_definition = sympy.sympify(phases.split(','))
                        frequency = phase_definition[-1]
                        del phase_definition[-1]
                        phase_sum = sum(phase_definition)
                        if phase_sum != 1:
                            raise ValueError("The sum of phase lengths must be 1")
                        else:
                            phase_definition.insert(0, len(phase_definition))
                    elif len(phases.split(',')) == 2:
                        phases = sympy.sympify(phases.split(','))
                        if int(phases[0]) < 2:
                            raise ValueError("The number of phases can't be less than 2")
                        else:
                            phase_definition.append(phases[0])
                            frequency = phases[-1]
                            for i in range(phase_definition[0]):
                                phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                    else:
                        raise SyntaxError(phase_def_syntax_error)
                elif phases.startswith("[") and phases.endswith("]"):
                    phases = phases.replace("[", "")
                    phases = phases.replace("]", "")
                    phase_definition = sympy.sympify(phases.split(','))
                    phase_sum = sum(phase_definition)
                    if phase_sum != 1:
                        raise ValueError("The sum of phase lengths must be 1")
                    else:
                        phase_definition.insert(0, len(phase_definition))
                elif "," not in phases:
                    if int(phases) < 2:
                        raise ValueError("The number of phases can't be less than 2")
                    else:
                        phase_definition.append(int(phases))
                        for i in range(phase_definition[0]):
                            phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                else:
                    raise SyntaxError(phase_def_syntax_error)

            else:
                raise SyntaxError(phase_def_syntax_error)
            phases = phase_definition
        return phases, frequency

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
                        impedance = utils.s/c.sym_value
                    else:
                        impedance = c.sym_value
                else:
                    if name[0] in ["c", "C"]:
                        impedance = utils.s/c.value
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
                ret = laplace.iLT(ret, self.sympy_ilt)
                ret = sympy.factor_terms(ret)
                return ret
        elif self.analysis_type == "DC":
            if ret is None:
                return sympy.Symbol(v)
            else:
                try:
                    ret = sympy.limit(ret, utils.s, 0)
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
                res = laplace.iLT(F, self.sympy_ilt)
                res = sympy.factor_terms(res)
                return res
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
                    if self.analysis_type == "TF":
                        ret = sympy.cancel((vn1 - vn2) / (val * utils.s))
                    elif self.analysis_type == "AC":
                        ret = sympy.cancel((vn1 - vn2) / (val * 2 * sympy.pi * utils.f * j))
                if c.type == "c":
                    if self.analysis_type == "TF":
                        ret = sympy.cancel((vn1 - vn2) * val*utils.s)
                    elif self.analysis_type == "AC":
                        ret = sympy.cancel((vn1 - vn2) * (val * 2 * sympy.pi * utils.f * j))

        if self.analysis_type == "tran":
            if ret is None:
                return laplace.iLT(sympy.Symbol(i), self.sympy_ilt)
            else:
                ret = laplace.iLT(ret, self.sympy_ilt)
                ret = sympy.factor_terms(ret)
                return ret
        elif self.analysis_type == "DC":
            if ret is None:
                return sympy.Symbol(i)
            else:
                try:
                    ret = sympy.limit(ret, utils.s, 0)
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
        if self.scsi == "undefined":
            if name == "all":
                ret = self.all_component_values(default_python_datatypes)
            else:
                v = f"v({name})"
                i = f"i({name})"
                ret[v] = self.component_voltage(name)
                ret[i] = self.component_current(name)
        else:
            if name == "all":
                ret = self.all_component_values(default_python_datatypes)
            else:
                temp = {}
                temp.update(self.SCSI_component_values(name))
                ret = self.SCSI_result_formatter(temp)
        return ret

    def all_component_values(self, default_python_datatypes=False):
        """
          Returns a dictionary of all relevant voltages and currents in the circuit.

          :return dict ret: in format {"v(id1)" : value, "i(id2)" : value, ...}
        """
        ret = {}
        if self.scsi == "undefined":
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
        else:
            temp = {}
            for key in self.components:
                temp.update(self.SCSI_component_values(self.components[key].name))
            ret = self.SCSI_result_formatter(temp)
        return ret

    def node_voltages(self):
        """
          Returns a dictionary of all node voltages in the circuit.

          :return dict ret: in format {"v(node1)" : value, ...}
        """
        ret = {}
        if self.scsi == "undefined":
            for node in self.node_dict:
                if node.find("*ctrl") == -1:
                    ret[f"v({node})"] = self.get_node_voltage(node)
                else: # Filter nodes added by virtual current sensors for CCC(V)S
                    pass
        else:
            num_of_phases = self.phases[0]
            for phase in range(1, num_of_phases + 1):
                for node in self.node_dict:
                    ret[sympy.symbols(f"v({node})_{phase}")] = self.SCSI_get_node_voltage(node, phase)
            ret = self.SCSI_result_formatter(ret)
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
        if self.phases == "undefined":
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)

            if self.analysis_type == "DC":
                if self.is_symbolic:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = sympy.limit(solved_dict[sym], utils.s, 0)
                        except KeyError:
                            pass
                        except TypeError:
                            pass
                else:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = sympy.limit(solved_dict[sym], utils.s, 0)
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["v", "i"]:
                                    if c.dc_value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.dc_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            #solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass

            elif self.analysis_type == "AC":
                if self.is_symbolic:
                    i = 1
                    for sym in symbols:
                        try:
                            solved_dict[sym] = solved_dict[sym].subs(utils.s, 2 * sympy.pi * utils.f * j)
                        except KeyError:
                            pass
                else:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = solved_dict[sym].subs(utils.s, 2 * sympy.pi * utils.f * j)
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["v", "i"]:
                                    if c.ac_value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            #solved_dict[sym] = solved_dict[sym].evalf(self.precision)
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
                            #solved_dict[sym] = solved_dict[sym].evalf(self.precision)
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
                            #solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass

        else: # used by SCSI
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
            self.SCSI_symbol_z_factor(solved_dict)
            self.SCSI_z_pow_inv_sub(solved_dict)
            if self.analysis_type == "TF":
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
                            solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass
            elif self.analysis_type == "AC":
                f = self.frequency
                if self.is_symbolic:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = solved_dict[sym].subs(z, sympy.exp(1) ** (2 * sympy.pi * f * j))
                        except KeyError:
                            pass
                if not self.is_symbolic:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = solved_dict[sym].subs(z, sympy.exp(1) ** (2 * sympy.pi * f * j))
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["v", "i"]:
                                    if c.ac_value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                            solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass
            elif self.analysis_type == "tran":
                raise ValueError("Transient analysis not yet implemented.")

        return eqn_matrix, solved_dict, symbols

    def _node_voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol(f"v({node})"))
        return voltage_symbol_list

    def _build_system_eqn(self):
        if self.method == "tableau" and self.phases == "undefined":
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

        elif self.method == "two_graph_node" and self.phases == "undefined":
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
                    if c.node3 is None:
                        c_v = self.components[c.current_sensor]
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

                    self.collapse(v_graph_collapses, node3, node4)
                    matrix_col_expand += 1

                if c.type == "h":   # CVT/CCVS
                    node1 = c.node1
                    node2 = c.node2
                    if c.node3 is None:
                        c_v = self.components[c.current_sensor]
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
                    self.collapse(i_graph_collapses, node1, node2)
                    self.collapse(v_graph_collapses, node3, node4)

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
                    #raise NotImplementedError("CCCS not supported in 'two_graph_node' analysis")
                    self._add_CCT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, i_graph_collapses)
                    symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                    index_col += 1
                if c.type == "h":
                    self._add_CVT_tgn(M, v_graph_nodes, i_graph_nodes, c, index_col, index_row, i_graph_collapses, v_graph_collapses)
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

        elif self.method == "two_graph_node" and self.phases != "undefined":
            sc_value_error = "Resistors, inductors, current sources, CCVSs and VCCSs aren't permitted in switched capacitor mode"
            num_of_phases = self.phases[0]
            symbols = []
            v_graph_collapses = []
            i_graph_collapses = []
            v_graph_nodes = []
            i_graph_nodes = []
            matrix_col_expand = 0

            for key in self.components:
                c = self.components[key]
                if c.type == "h":
                    raise ValueError("Current controlled voltage sources aren't permitted in SC/SI mode")
                if c.type == "l":
                    raise ValueError("Inductors are not permitted in SC/SI mode")
                if self.scsi == "scideal":
                    if c.type in ["r", "i", "g"]:
                        raise ValueError(sc_value_error)
                # else:
                #     if c.type == "v":
                #         raise ValueError("Inductors, voltage sources and CCVSs aren't permitted in switched current mode")


            for key in self.components:
                c = self.components[key]
                if c.type == "c":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                if c.type == "v":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
                if c.type == "s":
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(c.phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(c.phase), v_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node1 + '_' + str(c.phase), i_graph_nodes)
                    self.SCSI_graph_append_tgn(c.node2 + '_' + str(c.phase), i_graph_nodes)
                    self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(c.phase), c.node2 + '_' + str(c.phase))
                if c.type == "a":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
                        self.SCSI_collapse_tgn(v_graph_collapses, c.node3 + '_' + str(phase), c.node4 + '_' + str(phase))
                if c.type == "e":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_collapse_tgn(i_graph_collapses, c.node1 + '_' + str(phase), c.node2 + '_' + str(phase))
                if c.type == "f": # TODO: fix after bug #8 fix
                    c_v = self.components[c.current_sensor]
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c_v.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c_v.shorted_node + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c_v.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c_v.shorted_node + '_' + str(phase), i_graph_nodes)
                        self.SCSI_collapse_tgn(v_graph_collapses, c_v.node2 + '_' + str(phase),
                                               c_v.shorted_node + '_' + str(phase))
                    matrix_col_expand += num_of_phases
                if c.type == "r":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                if c.type == "i":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                if c.type == "g":
                    for phase in range(1, num_of_phases + 1):
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), v_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node1 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node2 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node3 + '_' + str(phase), i_graph_nodes)
                        self.SCSI_graph_append_tgn(c.node4 + '_' + str(phase), i_graph_nodes)

            i_graph_collapses = self.SCSI_collapse_list_collapser(i_graph_collapses)
            v_graph_collapses = self.SCSI_collapse_list_collapser(v_graph_collapses)

            for collapse_list in i_graph_collapses:
                i = 0
                tmp_i_graph_nodes = copy.copy(i_graph_nodes)
                for n in tmp_i_graph_nodes:
                    if n in collapse_list:
                        if any(node.startswith('0') for node in collapse_list):
                            i_graph_nodes.remove(n)
                        else:
                            i_graph_nodes[i] = min(collapse_list)
                    i += 1
            for collapse_list in v_graph_collapses:
                i = 0
                tmp_v_graph_nodes = copy.copy(v_graph_nodes)
                for n in tmp_v_graph_nodes:
                    if n in collapse_list:
                        if any(node.startswith('0') for node in collapse_list):
                            v_graph_nodes.remove(n)
                        else:
                            v_graph_nodes[i] = min(collapse_list)
                    i+=1
            self.node_voltage_identities = v_graph_collapses

            v_graph_nodes = list(set(v_graph_nodes))
            i_graph_nodes = list(set(i_graph_nodes))

            v_graph_nodes = self.SCSI_node_list_sort(v_graph_nodes)
            i_graph_nodes = self.SCSI_node_list_sort(i_graph_nodes)

            m_size = len(v_graph_nodes) + matrix_col_expand
            M = sympy.Matrix(sympy.zeros(m_size))
            S = sympy.Matrix(sympy.zeros(m_size, 1))
            index_row = 0
            index_col = 0
            symbols_to_append = []

            for key in self.components:
                c = self.components[key]
                if c.type == "c":
                    self.SCSI_add_capacitor_tgn(M, v_graph_nodes, i_graph_nodes,
                                                c, i_graph_collapses, v_graph_collapses,
                                                num_of_phases, matrix_col_expand)
                if c.type == "v":
                    index_row = self.SCSI_add_voltage_source_tgn(M, S, v_graph_nodes, i_graph_nodes,
                                                                 c, index_row, v_graph_collapses, num_of_phases)
                if c.type == "s":
                    index_row = self.SCSI_add_switch_tgn(M, v_graph_nodes, i_graph_nodes,
                                                         c, index_row, v_graph_collapses, num_of_phases)
                if c.type == "a":
                    pass
                if c.type == "e":
                    index_row = self.SCSI_add_VVT_tgn(M, v_graph_nodes, i_graph_nodes,
                                                      c, index_row, v_graph_collapses, num_of_phases)
                if c.type == "f":
                    pass
                if c.type == "r":
                    self.SI_add_resistor_tgn(M, v_graph_nodes, i_graph_nodes,
                                             c, i_graph_collapses, v_graph_collapses, num_of_phases)
                if c.type == "i":
                    self.SI_add_current_source_tgn(M, i_graph_nodes, c, i_graph_collapses, num_of_phases)
                if c.type == "g":
                    self.SI_add_VCT_tgn(M, v_graph_nodes, i_graph_nodes,
                                        c, i_graph_collapses, v_graph_collapses, num_of_phases)

            for key in self.components:
                c = self.components[key]
                if c.type == "f":
                    index_col = self.SCSI_add_QQT_tgn(M, v_graph_nodes, i_graph_nodes,
                                                      c, index_col, i_graph_collapses, num_of_phases)
                    for phase in range(1, num_of_phases + 1):
                        if self.scsi == "scideal":
                            symbols_to_append.append(sympy.Symbol(f"q({c.current_sensor})_{phase}"))
                        else:
                            symbols_to_append.append(sympy.Symbol(f"i({c.current_sensor})_{phase}"))

            equation_matrix = M.col_insert(m_size, S)

            for node in v_graph_nodes:
                symbols.append(sympy.Symbol(f"v({node[0:-2]}){node[-2:]}"))
            for symb in symbols_to_append:
                symbols.append(symb)

        elif self.method == "modified_node" and self.phases != "undefined": # used by SCSI
            if self.scsi == "siideal":
                raise ValueError("Use modified nodal formulation method only for switched capacitor circuits")

            num_of_phases = self.phases[0]
            component_size = self.c_count * num_of_phases
            node_size = self.node_count * num_of_phases
            M = sympy.Matrix(sympy.zeros(component_size + node_size))
            R = sympy.Matrix(sympy.zeros(component_size + node_size, 1))

            component_index = 0
            node_symbols = self.SCSI_node_voltage_symbols()
            charge_symbols = []

            A_matrices = []
            Y_matrices = []
            Z_matrices = []
            for phase_y in range(num_of_phases):
                A_matrices.append([])
                Y_matrices.append([])
                Z_matrices.append([])
                for phase_x in range(num_of_phases):
                    A_matrices[phase_y].append(sympy.Matrix(sympy.zeros(self.node_count, self.c_count)))
                    Y_matrices[phase_y].append(sympy.Matrix(sympy.zeros(self.c_count)))
                    Z_matrices[phase_y].append(sympy.Matrix(sympy.zeros(self.c_count)))

            for key in self.components:
                c = self.components[key]
                if c.type == "v":
                    for phase in range(num_of_phases):
                        self.SCSI_add_voltage_source(A_matrices, Y_matrices, R,
                                                     c, component_index, phase, num_of_phases)
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 1
                elif c.type == "c":
                    for phase_y in range(num_of_phases):
                        for phase_x in range(num_of_phases):
                            self.SCSI_add_capacitor(A_matrices, Y_matrices, Z_matrices,
                                                    c, component_index, phase_y, phase_x, num_of_phases)
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase_y + 1)))
                    component_index += 1
                elif c.type == "s":
                    for phase in range(num_of_phases):
                        self.SCSI_add_switch(A_matrices, Y_matrices, Z_matrices, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 1
                elif c.type == "a":
                    for phase in range(num_of_phases):
                        self.SCSI_add_OpAmp(A_matrices, Y_matrices, Z_matrices, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_in{phase}".format(name=c.name, phase=phase + 1)))
                        charge_symbols.append(sympy.Symbol("q({name})_out{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 2
                elif c.type == "e":
                    for phase in range(num_of_phases):
                        self.SCSI_add_VVT(A_matrices, Y_matrices, Z_matrices, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_in{phase}".format(name=c.name, phase=phase + 1)))
                        charge_symbols.append(sympy.Symbol("q({name})_out{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 2
                elif c.type == "f":
                    for phase in range(num_of_phases):
                        self.SCSI_add_QQT(A_matrices, Y_matrices, Z_matrices, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_in{phase}".format(name=c.name, phase=phase + 1)))
                        charge_symbols.append(sympy.Symbol("q({name})_out{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 2

            for phase_y in range(num_of_phases):
                for phase_x in range(num_of_phases):
                    self.SCSI_submatrix_write(M, A_matrices, 0, self.node_count, phase_y, phase_x)
                    self.SCSI_submatrix_write(M, Z_matrices, self.node_count, self.node_count, phase_y, phase_x)
                    self.SCSI_submatrix_write(M, Y_matrices[phase_y][phase_x] * A_matrices[phase_y][phase_x].T,
                                              self.node_count, 0, phase_y, phase_x)

            self.SCSI_matrix_z_symbol(M)
            equation_matrix = M.col_insert(component_size + node_size, R)
            symbols = self.SCSI_symbol_list_order(node_symbols, charge_symbols)

        return equation_matrix, symbols

    def collapse(self, graph_collapses, node1, node2):
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
            c_v = self.components[c.current_sensor]
            node3 = c_v.shorted_node
            c.node3 = node3
        else:
            node3 = c.node3
        node4 = c.node4

        c_v = self.components[c.current_sensor]

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
        f = None
        if self.is_symbolic:
            f = c.sym_value
        else:
            f = c.value
        node1 = c.node1
        node2 = c.node2
        if c.node3 is None:
            c_v = self.components[c.current_sensor]
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
            M[n1i, col] += f
        if n2i is not None:
            M[n2i, col] += -f
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
            y = 1 / (val*utils.s)
        if c.type == "c":
            y = utils.s * val
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
            z_b = -utils.s * val
        elif c.type == "c":
            y_b = utils.s * val
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
        if c.node3 is None:
            c_v = self.components[c.current_sensor]
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
            c_v = self.components[c.current_sensor]
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

        matrix[self.c_count + L2_index, self.c_count + L1_index] += -utils.s*M
        matrix[self.c_count + L1_index, self.c_count + L2_index] += -utils.s*M
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

    def SCSI_node_voltage_symbols(self):
        num_of_phases = self.phases[0]
        voltage_symbol_list = []
        for phases in range(num_of_phases):
            for node in self.node_dict:
                if node != "0":
                    voltage_symbol_list.append(sympy.Symbol("v({node})_{phases}".format(node=node, phases=phases + 1)))
        return voltage_symbol_list

    def SCSI_symbol_list_order(self, node_symbols, charge_symbols):
        charge_length = len(charge_symbols)
        node_length = len(node_symbols)
        num_of_phases = self.phases[0]
        symbols = []
        for phase in range(1, num_of_phases + 1):
            for symbol in range(node_length):
                if int(str(node_symbols[symbol])[-1]) == phase:
                    symbols.append(node_symbols[symbol])
            for symbol in range(charge_length):
                if int(str(charge_symbols[symbol])[-1]) == phase:
                    symbols.append(charge_symbols[symbol])
        return symbols

    def SCSI_symbol_z_factor(self, symbols):
        num_of_phases = self.phases[0]
        phase_helper_array = []
        for phase in range(num_of_phases):
            phase_helper_array.append(self.phases[phase + 1])
        for symbol in symbols:
            symbol_phase = int(str(symbol)[-1])
            temp = 0
            for phase in range(symbol_phase - 1):
                temp += phase_helper_array[phase]
            if symbol_phase != num_of_phases:
                symbols[symbol] *= z ** (- temp)
            else:
                symbols[symbol] *= z ** phase_helper_array[-1]

    def SCSI_z_pow_inv_sub(self, solved_dict):
        for expression in solved_dict:
            solved_dict[expression] = solved_dict[expression].subs(self.symbol_dict)

    def SCSI_z_power_substitution(self, symbol):
        dict_length = len(self.symbol_dict)
        x = sympy.Symbol("x_{index}".format(index=dict_length))
        self.symbol_dict[x] = symbol
        return x

    def SCSI_submatrix_write(self, matrix, submatrix, start_y, start_x, phase_y, phase_x):
        if isinstance(submatrix, list):
            y_dimension = sympy.shape(submatrix[0][0])[0]
            x_dimension = sympy.shape(submatrix[0][0])[1]
        else:
            y_dimension = sympy.shape(submatrix)[0]
            x_dimension = sympy.shape(submatrix)[1]
        phase_y_offset = (self.c_count + self.node_count) * phase_y
        phase_x_offset = (self.c_count + self.node_count) * phase_x
        for y in range(y_dimension):
            for x in range(x_dimension):
                if isinstance(submatrix, list):
                    matrix[phase_y_offset + start_y + y, phase_x_offset + start_x + x] += submatrix[phase_y][phase_x][y, x]
                else:
                    matrix[phase_y_offset + start_y + y, phase_x_offset + start_x + x] += submatrix[y, x]

    def SCSI_incidence_matrix_generate(self, incidence_matrix, N1, N2, component_index, phase_y, phase_x, capacitor=False):
        if N1 == "0":
            pass
        else:
            node_pos = self.node_dict[N1]
            if capacitor:
                incidence_matrix[phase_y][phase_x][node_pos, component_index] = -1
            else:
                incidence_matrix[phase_y][phase_x][node_pos, component_index] = 1
        if N2 == "0":
            pass
        else:
            node_pos = self.node_dict[N2]
            if capacitor:
                incidence_matrix[phase_y][phase_x][node_pos, component_index] = 1
            else:
                incidence_matrix[phase_y][phase_x][node_pos, component_index] = -1

    def SCSI_matrix_z_symbol(self, matrix):
        num_of_phases = self.phases[0]
        submatrix_dimension = self.node_count + self.c_count
        phase_offset = submatrix_dimension * (num_of_phases - 1)
        for y in range(submatrix_dimension):
            for x in range(submatrix_dimension):
                matrix[phase_offset + y, phase_offset + x] *= z

    def SCSI_add_capacitor(self, A, Y, Z, c, component_index, phase_y, phase_x, num_of_phases):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        y_b = val
        z_b = -1

        Y[phase_y][phase_x][component_index, component_index] = y_b
        if phase_y == phase_x:
            Z[phase_y][phase_x][component_index, component_index] = z_b
            self.SCSI_incidence_matrix_generate(A, N1, N2, component_index, phase_y, phase_x)
        else:
            if (abs(phase_y - phase_x) == 1 and phase_y > phase_x) or (phase_y == 0 and phase_x + 1 == num_of_phases):
                Z[phase_y][phase_x][component_index, component_index] = - z_b
                self.SCSI_incidence_matrix_generate(A, N1, N2, component_index, phase_y, phase_x, True)

    def SCSI_add_voltage_source(self, A, Y, result, c, component_index, phase, num_of_phases):
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
        Y[phase][phase][component_index, component_index] = 1
        self.SCSI_incidence_matrix_generate(A, N1, N2, component_index, phase, phase)

        phase_offset = self.c_count * phase
        node_offset = self.node_count * (phase + 1)
        phase_helper_array = []
        z_symbol_array = []
        temp = 0
        for f in range(num_of_phases):
            phase_helper_array.append(temp)
            temp += self.phases[f + 1]
            z_symbol_array.append(z ** sympy.sympify(phase_helper_array[f]))
        result[phase_offset + node_offset + component_index, 0] = val * self.SCSI_z_power_substitution(
          z_symbol_array[phase])

    def SCSI_add_switch(self, A, Y, Z, c, component_index, phase):
        N1 = c.node1
        N2 = c.node2
        switch_phase = (int(c.phase)) - 1

        self.SCSI_incidence_matrix_generate(A, N1, N2, component_index, phase, phase)
        if phase == switch_phase:
            Y[phase][phase][component_index, component_index] = 1
        else:
            Z[phase][phase][component_index, component_index] = -1

    def SCSI_add_OpAmp(self, A, Y, Z, c, component_index, phase):
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        self.SCSI_incidence_matrix_generate(A, N3, N4, component_index, phase, phase)
        self.SCSI_incidence_matrix_generate(A, N1, N2, component_index + 1, phase, phase)
        Y[phase][phase][component_index, component_index] = 1
        Z[phase][phase][component_index + 1, component_index] = 1

    def SCSI_add_VVT(self, A, Y, Z, c, component_index, phase):
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        self.SCSI_incidence_matrix_generate(A, N3, N4, component_index, phase, phase)
        self.SCSI_incidence_matrix_generate(A, N1, N2, component_index + 1, phase, phase)
        Y[phase][phase][component_index + 1, component_index] = val
        Y[phase][phase][component_index + 1, component_index + 1] = -1
        Z[phase][phase][component_index, component_index] = 1

    def SCSI_add_QQT(self, A, Y, Z, c, component_index, phase):
        # TODO: Fix after bug #8 fix "shorted_node" removed, the node info is now directly in "c"
        N1 = c.node1
        N2 = c.node2
        c_v = self.components[c.current_sensor]
        N3 = c_v.node2
        N4 = c_v.shorted_node
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        self.SCSI_incidence_matrix_generate(A, N3, N4, component_index, phase, phase)
        self.SCSI_incidence_matrix_generate(A, N1, N2, component_index + 1, phase, phase)
        Y[phase][phase][component_index, component_index] = 1
        Z[phase][phase][component_index + 1, component_index] = val
        Z[phase][phase][component_index + 1, component_index + 1] = -1

    def SCSI_graph_append_tgn(self, node, graph):
        if node.startswith('0'):
            pass
        elif node not in graph:
            graph.append(node)

    def SCSI_collapse_tgn(self, graph_collapses, node1, node2):
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

    def SCSI_collapse_list_collapser(self, graph_collapses):
        for collapse_list in graph_collapses:
            for node in collapse_list:
                for cl in graph_collapses:
                    if node in cl and cl!=collapse_list:
                        for n in cl:
                            if n not in collapse_list:
                                collapse_list.append(n)
        for collapse_list in graph_collapses:
            collapse_list.sort()

        temp = []
        for collapse_list in graph_collapses:
            if collapse_list not in temp:
                temp.append(collapse_list)

        return temp

    def SCSI_node_list_sort(self, nodes):
        helper = []
        i = 0
        for element in nodes:
            helper.append([])
            helper[i].append(element[:element.index('_')])
            helper[i].append(element[element.index('_') + 1:])
            i += 1
        helper.sort(key=operator.itemgetter(1, 0))
        nodes_new = []
        for element in helper:
            nodes_new.append(element[0] + '_' + element[1])
        return nodes_new

    def SCSI_index_tgn(self, nodes, node, collapses):
        try:
            return nodes.index(node)
        except ValueError:
            for collapse_list in collapses:
                if node in collapse_list:
                    if any(n.startswith('0') for n in collapse_list):
                        return None
                    else:
                        return nodes.index(min(collapse_list))
                elif node.startswith("0"):
                    return None

    def SCSI_add_capacitor_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases, matrix_col_expand):
        phase_offset = int((M.shape[0] - matrix_col_expand) / num_of_phases)
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1v is not None:
                if n1i is not None:
                    if phase == num_of_phases:
                        M[n1i, n1v] += +val*z
                    else:
                        M[n1i, n1v] += +val
                    if phase == 1:
                        M[n1i, n1v - phase_offset - matrix_col_expand] += -val
                    else:
                        M[n1i, n1v - phase_offset] += -val
                if n2i is not None:
                    if phase == num_of_phases:
                        M[n2i, n1v] += -val*z
                    else:
                        M[n2i, n1v] += -val
                    if phase == 1:
                        M[n2i, n1v - phase_offset - matrix_col_expand] += +val
                    else:
                        M[n2i, n1v - phase_offset] += +val
            if n2v is not None:
                if n1i is not None:
                    if phase == num_of_phases:
                        M[n1i, n2v] += -val*z
                    else:
                        M[n1i, n2v] += -val
                    if phase == 1:
                        M[n1i, n2v - phase_offset - matrix_col_expand] += +val
                    else:
                        M[n1i, n2v - phase_offset] += +val
                if n2i is not None:
                    if phase == num_of_phases:
                        M[n2i, n2v] += +val*z
                    else:
                        M[n2i, n2v] += +val
                    if phase == 1:
                        M[n2i, n2v - phase_offset - matrix_col_expand] += -val
                    else:
                        M[n2i, n2v - phase_offset] += -val

    def SCSI_add_voltage_source_tgn(self, M, S, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value
        phase_helper_array = []
        z_symbol_array = []
        temp = 0
        for f in range(num_of_phases):
            phase_helper_array.append(temp)
            temp += self.phases[f + 1]
            z_symbol_array.append(z ** sympy.sympify(phase_helper_array[f]))
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            row = len(i_nodes) + index_row
            if n1v is not None:
                if phase == num_of_phases:
                    M[row, n1v] += z
                else:
                    M[row, n1v] += 1
            if n2v is not None:
                if phase == num_of_phases:
                    M[row, n2v] += -z
                else:
                    M[row, n2v] += -1
            S[row, 0] += val * self.SCSI_z_power_substitution(z_symbol_array[phase - 1])
            index_row += 1
        return index_row

    def SCSI_add_switch_tgn(self, M, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(c.phase), v_graph_collapses)
        n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(c.phase), v_graph_collapses)
        row = len(i_nodes) + index_row
        if n1v is not None:
            if int(c.phase) == num_of_phases:
                M[row, n1v] += z
            else:
                M[row, n1v] += 1
        if n2v is not None:
            if int(c.phase) == num_of_phases:
                M[row, n2v] += -z
            else:
                M[row, n2v] += -1
        index_row += 1
        return index_row

    def SCSI_add_VVT_tgn(self, M, v_nodes, i_nodes, c, index_row, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            e = c.sym_value
        else:
            e = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node3 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node4 + '_' + str(phase), v_graph_collapses)
            n3v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n4v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            row = len(i_nodes) + index_row
            if n1v is not None:
                if phase == num_of_phases:
                    M[row, n1v] += -e*z
                else:
                    M[row, n1v] += -e
            if n2v is not None:
                if phase == num_of_phases:
                    M[row, n2v] += e*z
                else:
                    M[row, n2v] += e
            if n3v is not None:
                if phase == num_of_phases:
                    M[row, n3v] += z
                else:
                    M[row, n3v] += 1
            if n4v is not None:
                if phase == num_of_phases:
                    M[row, n4v] += -z
                else:
                    M[row, n4v] += -1
            index_row += 1
        return index_row

    def SCSI_add_QQT_tgn(self, M, v_nodes, i_nodes, c, index_col, i_graph_collapses, num_of_phases):
        if self.is_symbolic:
            f = c.sym_value
        else:
            f = c.value
        c_v = self.components[c.current_sensor]
        for phase in range(1, num_of_phases + 1):
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            # TODO: fix after bug #8 fix
            n3i = self.SCSI_index_tgn(i_nodes, c_v.node2 + '_' + str(phase), i_graph_collapses)
            n4i = self.SCSI_index_tgn(i_nodes, c_v.shorted_node + '_' + str(phase), i_graph_collapses)
            col = len(v_nodes) + index_col
            if n1i is not None:
                if phase == num_of_phases:
                    M[n1i, col] += f*z
                else:
                    M[n1i, col] += f
            if n2i is not None:
                if phase == num_of_phases:
                    M[n2i, col] += -f*z
                else:
                    M[n2i, col] += -f
            if n3i is not None:
                if phase == num_of_phases:
                    M[n3i, col] += z
                else:
                    M[n3i, col] += 1
            if n4i is not None:
                if phase == num_of_phases:
                    M[n4i, col] += -z
                else:
                    M[n4i, col] += -1
            index_col += 1
        return index_col

    def SI_add_resistor_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        y = 1 / val
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node1 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node2 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
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

    def SI_add_current_source_tgn(self, S, i_nodes, c, i_graph_collapses, num_of_phases):
        if self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value
        for phase in range(1, num_of_phases + 1):
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
            if n1i is not None:
                S[n1i, 0] += -val
            if n2i is not None:
                S[n2i, 0] += val

    def SI_add_VCT_tgn(self, M, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses, num_of_phases):
        if self.is_symbolic:
            g = c.sym_value
        else:
            g = c.value
        for phase in range(1, num_of_phases + 1):
            n1v = self.SCSI_index_tgn(v_nodes, c.node3 + '_' + str(phase), v_graph_collapses)
            n2v = self.SCSI_index_tgn(v_nodes, c.node4 + '_' + str(phase), v_graph_collapses)
            n1i = self.SCSI_index_tgn(i_nodes, c.node1 + '_' + str(phase), i_graph_collapses)
            n2i = self.SCSI_index_tgn(i_nodes, c.node2 + '_' + str(phase), i_graph_collapses)
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

    def SCSI_get_node_voltage(self, node, phase, force_z = False):
        phase_string = str(phase)
        formatted_node = str(node) + "_" + phase_string
        formatted_voltage = sympy.symbols("v(" + str(node) + ")_" + phase_string)
        try:
            value = self.solved_dict[formatted_voltage]
        except KeyError:
            for identity in self.node_voltage_identities:
                if formatted_node in identity:
                    if any(identity_node.startswith("0") for identity_node in identity):
                        value = 0
                    else:
                        for identity_node in identity:
                            try:
                                if str(identity_node).count("_") == 1:
                                    formatted_voltage = sympy.symbols("v(" + identity_node.split("_")[0]
                                                                      + ")_" + phase_string)
                                else:
                                    formatted_voltage = sympy.symbols("v(" + identity_node.split(")_")[0]
                                                                      + "))_" + phase_string)
                                value = self.solved_dict[formatted_voltage]
                            except KeyError:
                                pass
        if force_z:
            return value
        elif self.analysis_type == "tran":
            pass
            #return z_transform.iZT(value)
        else:
            return value

    def SCSI_component_voltage(self, name, phase):
        value_dict = {}
        c = self.components[name]
        biquad_input_key = sympy.symbols(f"v({name}-in)_{phase}")
        if c.node1 == "0":
            node1_value = 0
        else:
            node1_value = self.SCSI_get_node_voltage(c.node1, phase, force_z=True)
        if c.node2 == "0":
            node2_value = 0
        else:
            node2_value = self.SCSI_get_node_voltage(c.node2, phase, force_z=True)
        if c.type in ["a", "e", "g"]:
            if c.node3 == "0":
                node3_value = 0
            else:
                node3_value = self.SCSI_get_node_voltage(c.node3, phase, force_z=True)
            if c.node4 == "0":
                node4_value = 0
            else:
                node4_value = self.SCSI_get_node_voltage(c.node4, phase, force_z=True)
            value_dict[biquad_input_key] = sympy.cancel(node3_value - node4_value)
        if c.type == "f":
            c_v = self.components[c.current_sensor]
            if c_v.node2 == "0":
                cv_node_value = 0
            else:
                cv_node_value = self.SCSI_get_node_voltage(c_v.node2, phase, force_z=True)
            if c_v.shorted_node == "0": # TODO: fix after bug #8 fix
                cv_short_value = 0
            else:
                cv_short_value = self.SCSI_get_node_voltage(c_v.shorted_node, phase, force_z=True) # TODO: fix after bug #8 fix
            value_dict[biquad_input_key] = sympy.cancel(cv_node_value - cv_short_value)
        if c.type in ["a", "e", "f", "g"]:
            value_dict[sympy.symbols(f"v({name}-out)_{phase}")] = sympy.cancel(node1_value - node2_value)
        else:
            value_dict[sympy.symbols(f"v({name})_{phase}")] = sympy.cancel(node1_value - node2_value)
        # if self.analysis_type == "tran":
        #     for entry in value_dict:
        #         value_dict[entry] = z_transform.IZT(value_dict[entry])
        return value_dict

    def SCSI_component_charge(self, name, phase):
        if self.scsi == "scideal":
            charge = sympy.symbols(f"q({name})_{phase}")
        else:
            charge = sympy.symbols(f"i({name})_{phase}")
        value_dict = {}
        c = self.components[name]
        voltage_key = sympy.symbols(f"v({name})_{phase}")
        if self.is_symbolic:
            value = c.sym_value
        else:
            value = c.value
        if c.type in ["c", "r"]:
            if c.type == "r":
                value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] / value)
            if c.type == "c":
                value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] * value)
        elif c.type == "i":
            if self.is_symbolic:
                value = c.sym_value
            else:
                if self.analysis_type == "DC":
                    value = c.dc_value
                elif self.analysis_type == "tran":
                    value = c.tran_value
                else:
                    value = c.ac_value
            value_dict[charge] = sympy.symbols(value)
            #value_dict[charge] = str(value)
        elif c.type == "g":
            value_dict[charge] = sympy.cancel(self.SCSI_component_voltage(name, phase)[voltage_key] * value)
        elif c.type == "f":
            if self.scsi == "scideal":
                control_charge = sympy.symbols(f"q({c.current_sensor})_{phase}")
            else:
                control_charge = sympy.symbols(f"i({c.current_sensor})_{phase}")
            value_dict[charge] = sympy.cancel(self.solved_dict[control_charge] * value)
        else:
            value_dict[charge] = charge
        # if self.analysis_type == "tran":
        #     for entry in value_dict:
        #         value_dict[entry] = z_transform.IZT(value_dict[entry])
        return value_dict

    def SCSI_component_values(self, name):
        num_of_phases = self.phases[0]
        value_dict = {}
        for phase in range(1, num_of_phases + 1):
            value_dict.update(self.SCSI_component_voltage(name, phase))
            value_dict.update(self.SCSI_component_charge(name, phase))
        return value_dict

    def SCSI_result_formatter(self, dict):
        temp = {}
        for key in dict:
            if list(str(key)).count("_") == 1:
                name = str(key).split("_")[0]
            else:
                name = str(key).split(")_")[0]
                name = "".join([name, ')'])
            temp[name] = []
            for i in range(self.phases[0]):
                temp[name].append(0)
        for key in dict:
            if list(str(key)).count("_") == 1:
                name = str(key).split("_")[0]
                phase = int(str(key).split("_")[1])
            else:
                name = str(key).split(")_")[0]
                name = "".join([name, ')'])
                phase = int(str(key).split(")_")[1])
            temp[name][phase - 1] = dict[key]
        return temp

