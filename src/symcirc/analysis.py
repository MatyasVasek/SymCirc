import os
import copy, time
from typing import List, Union

from symcirc import laplace
from symcirc.circuit import Circuit
from symcirc.component import Component, CurrentSource, VoltageSource, ANALYSIS_DEPENDENT
from symcirc.utils import *


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
    SOLVERS = ['gauss', 'lu', 'ddd']
    METHODS = ["tableau", "two_graph_node"]
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):

        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        method = method.lower()
        if method not in self.METHODS:
            raise ValueError(f"Nonexistent analysis method: {method}\n Pick from: {self.METHODS}")

        if solver not in self.SOLVERS:
            raise ValueError(f"Nonexistent solver: {solver}\nPick from: {self.SOLVERS}")

        self.is_symbolic: bool = symbolic
        self.auto_eval = auto_eval
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt
        self.solver = solver

        self.sbg = True

        self.circuit: Circuit = circuit
        self.node_voltage_identities: list = []

        self.symbol_dict = {}

        self.node_voltage_symbols: List[sympy.Symbol] = self.circuit.get_node_symbols()
        self.node_dict = self.circuit.get_node_dict()
        self.c_count = self.circuit.count_ports()
        self.node_count = self.circuit.count_nodes()

        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}
        self.symbol_dict = self.get_symbols()  # format: {<symbol_name> : <Symbol>}

    def v(self, name: str) -> sympy.Expr:
        """
        Returns the specified voltage
        """
        return self.component_voltage(name)

    def i(self, name: str) -> sympy.Expr:
        """
        Returns the specified current
        """
        return self.component_current(name)

    def get_all_results(self):
        results = self.node_voltages()
        results.update(self.component_values("all"))
        return results

    def get_symbols(self) -> Dict[str, sympy.Symbol]:
        #t0 = time.time()
        symbol_dict = {"s": s, "t": t, "f": f}
        for comp in self.circuit.components.values():
            if comp.type in ANALYSIS_DEPENDENT:
                sym_val = self._choose_source_val(comp, form="sym")
            else:
                sym_val = comp.sym_value

            if sym_val is None:
                continue
            elif hasattr(sym_val, "name"):
                symbol_dict[sym_val.name] = sym_val
            else:
                try:
                    for sym in sym_val.free_symbols:
                        symbol_dict[sym.name] = sym
                except AttributeError:
                    pass
        return symbol_dict
        # Old version for reference, remove later
        '''symbol_dict = {}
        t1 = time.time()
        for expr in self.solved_dict.values():
            try:
                free_symbols = expr.free_symbols
                for symbol in free_symbols:
                    symbol_dict[symbol.name] = symbol
            except AttributeError:
                pass
        print(time.time()-t1)
        return symbol_dict'''

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
                vn1 = self.get_node_voltage(n1, force_s_domain=True)
            if n2 == "0":
                vn2 = 0
            else:
                vn2 = self.get_node_voltage(n2, force_s_domain=True)

            ret = sympy.cancel((vn1 - vn2))
        if ret is None:
            ret = sympy.Symbol(v, )
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """Has to be implemented in child class"""
        ret = None
        i = f"i({name})"
        i_symbol = sympy.symbols(i)

        if self.method == "tableau":
            ret = self.solved_dict[i_symbol]

        if self.method == "two_graph_node":
            if i_symbol in self.solved_dict:
                return self.solved_dict[i_symbol]

            c = self.circuit.components[name]
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
                    ret = sympy.cancel((vn1 - vn2) * val * s)
        return ret

    def get_node_voltage(self, node: str, force_s_domain=False) -> Union[sympy.Expr, None]:
        func = None
        try:
            func = self.solved_dict[self.node_voltage_symbols[self.node_dict[node]]]
        except KeyError:
            for identity in self.node_voltage_identities:
                if node in identity:
                    if "0" in identity:
                        return 0
                    else:
                        for n in identity:
                            try:
                                func = self.solved_dict[self.node_voltage_symbols[self.node_dict[n]]]
                            except:
                                pass
        return func

    def value_map(self, is_float=False):
        value_map = {}
        for comp in self.circuit.components.values():
            if is_float:
                if comp.type in ANALYSIS_DEPENDENT:
                    symbol = self._choose_source_val(comp, form="sym")
                    value_map[symbol] = self._choose_source_val(comp, form="float")
                elif comp.value_float is not None:
                    value_map[comp.sym_value] = comp.value_float
            else:
                if comp.type in ANALYSIS_DEPENDENT:
                    symbol = self._choose_source_val(comp, form="sym")
                    value_map[symbol] = self._choose_source_val(comp, form="rat")
                elif comp.value is not None:
                    value_map[comp.sym_value] = comp.value
        return value_map

    def _choose_source_val(self, c: Component, form="rat") -> sympy.Expr:
        """
        Has to be implemented in child class
        :param form: can be 'rat', 'float', 'sym'
        """
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
            raise ValueError(f"Node {node1} doesn't exist.")

        try:
            if node2 in [0, "0"]:
                voltage2 = 0
            else:
                v2 = sympy.Symbol(f"v({node2})")
                voltage2 = self.solved_dict[v2].simplify()
        except KeyError:
            raise ValueError(f"Node {node2} doesn't exist.")
        tf = (voltage2/voltage1)
        return tf


    def abstract_matrix(self, eqn_matrix, memo):
        i = 0
        for elem in eqn_matrix:
            new_elem = simplify_coeffs(elem, s, f"M{i}", memo)
            eqn_matrix[i] = new_elem
            i += 1

    def unabstract_results(self, solved_dict, memo):
        # Resubstitute
        inverted_memo = {v: k for k, v in memo.items()}
        for key in solved_dict:
            solved_dict[key] = solved_dict[key].subs(inverted_memo)

    def _analyse(self):
        """
          Implementation of all types of supported analysis.

            eqn_matrix, solved_dict, symbols
          :return sympy.Matrix eqn_matrix: matrix of the system equations
          :return dict solved_dict: dictionary of eqn_matrix solve results
          :return list symbols: list of all used sympy.symbol objects
        """

        # TODO: state machine to decide whether to use gauss, LU or DDD

        solved_dict = {}
        if self.solver == "gauss":
            memo = {}
            eqn_matrix, symbols = self._build_system_eqn()
            #self.abstract_matrix(eqn_matrix, memo)
            solved_dict = self._gauss_solve(eqn_matrix, symbols)
            #self.unabstract_results(solved_dict, memo)
        elif self.solver == "lu":
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = self._lu_solve(eqn_matrix, symbols)
            for key in solved_dict:
                solved_dict[key] = general_simplify(solved_dict[key])
        elif self.solver == "ddd":
            # TODO: rework this hack into a more professional solution
            memo = {}
            tmp_symbolic = self.is_symbolic
            self.is_symbolic = True
            eqn_matrix, symbols = self._build_system_eqn()
            #self.abstract_matrix(eqn_matrix, memo)
            self.is_symbolic = tmp_symbolic
            solved_dict = self._ddd_solve(eqn_matrix, symbols)
            for key in solved_dict:
                solved_dict[key] = general_simplify(solved_dict[key])
            #self.unabstract_results(solved_dict, memo)

        t1 = time.time()
        #print(f"Matrix solve time by {self.solver}: {t1-t0}")

        # TODO: test whether this is a speedup or not and neccessity for readability
        #for key in solved_dict:
        #    solved_dict[key] = solved_dict[key].factor()
        #print(f"Gauss solve time: {time.time()-t0}")
        #print(solved_dict)
        # TODO: auto_eval implementation in expression retrieval; not directly here, causes issues with laplace
        #if self.auto_eval:
        #    solved_dict = {key: evalf(func, precision=self.precision) for key, func in solved_dict.items()}
        return eqn_matrix, solved_dict, symbols

    def _gauss_solve(self, eqn_matrix, symbols):
        return sympy.solve_linear_system(eqn_matrix, *symbols)

    def _lu_solve(self, eqn_matrix, symbols):
        return sympy.solve_linear_system_LU(eqn_matrix, symbols)

    def _ddd_solve(self, eqn_matrix, symbols):
        from symcirc.ddd_solve import solve_cramer_ddd, solve_cramer_ddd_semi
        if self.is_symbolic:
            return solve_cramer_ddd(eqn_matrix, symbols, nested=True)
        else:
            return solve_cramer_ddd_semi(eqn_matrix, symbols, self.value_map(is_float=True))

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
        for key in self.circuit.components:
            c = self.circuit.components[key]

            if c.type == "k":
                matrix_col_expand += 2

            if c.type in ["r", "l", "c"]:
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)

            if c.type == "v":
                # Modified topology to also represent current through the source
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                matrix_col_expand += 1

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

            if c.type == "w": # short
                self.graph_append(c.node1, v_graph_nodes)
                self.graph_append(c.node2, v_graph_nodes)
                self.graph_append(c.node1, i_graph_nodes)
                self.graph_append(c.node2, i_graph_nodes)
                self._collapse(i_graph_collapses, c.node1, c.node2)
                self._collapse(v_graph_collapses, c.node1, c.node2)


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
            if c.type in ["r", "c"]:
                self._add_basic_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "l":
                if c.coupling is None:
                    self._add_basic_tgn(M, S, v_graph_nodes, i_graph_nodes, c, i_graph_collapses, v_graph_collapses)
            if c.type == "v":
                self._add_voltage_source_tgn(M, S, v_graph_nodes, i_graph_nodes, c, index_row, index_col,
                                             i_graph_collapses, v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({c.name})"))
                index_row += 1
                index_col += 1
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
                symbols_to_append.append(sympy.Symbol(f"i({c.name}_in)"))
                index_col += 1
                index_row += 1
            if c.type == "a":
                pass
            if c.type == "s":
                raise NotImplementedError("Switch not supported in 'two_graph_node' analysis")

        if len(self.circuit.couplings) > 0:
            for coupling in self.circuit.couplings:
                ind1 = self.circuit.components[coupling.L1]
                ind2 = self.circuit.components[coupling.L2]
                if self.is_symbolic:
                    coef = coupling.sym_value
                else:
                    coef = coupling.value
                self._add_K_tgn(M, S, v_graph_nodes, i_graph_nodes, [ind1, ind2, coef], index_col, index_row,
                                i_graph_collapses, v_graph_collapses)
                symbols_to_append.append(sympy.Symbol(f"i({ind1.name})"))
                symbols_to_append.append(sympy.Symbol(f"i({ind2.name})"))
                index_col += 2
                index_row += 2

        equation_matrix = M.col_insert(m_size, S)

        for node in v_graph_nodes:
            symbols.append(sympy.Symbol(f"v({node})"))
        for symb in symbols_to_append:
            symbols.append(symb)

        # TODO: experiment with simplification inside matrix - seems like a huge performance upgrade in tgn method!
        '''for i in range(m_size ** 2):
            expr = equation_matrix[i]
            if expr != 0:
                equation_matrix[i] = sympy.cancel(expr)'''

        return equation_matrix, symbols

    def _collapse(self, graph_collapses, node1, node2):
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
        node3 = c.node3
        node4 = c.node4

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
        if self.is_symbolic:
            f_gain = c.sym_value
        else:
            f_gain = c.value

        node1 = c.node1
        node2 = c.node2
        node3 = c.node3
        node4 = c.node4


        n1i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node4, i_graph_collapses)
        n3i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        col = len(v_nodes) + index

        if n1i is not None:
            M[n1i, col] += 1/f_gain
        if n2i is not None:
            M[n2i, col] += -1/f_gain
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

    def _add_voltage_source_tgn(self, M, S, v_nodes, i_nodes, c, index_row, index_col, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = self._choose_source_val(c, form)

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

        col = len(v_nodes) + index_col
        row = len(i_nodes) + index_row

        if n1i is not None:
            M[n1i, col] += 1
        if n2i is not None:
            M[n2i, col] += -1
        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1
        S[row, 0] += val

    def _add_current_source_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = self._choose_source_val(c, form)

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

        if n1i is not None:
            S[n1i, 0] += -val
        if n2i is not None:
            S[n2i, 0] += val

    def _add_basic_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        node1 = c.node1
        node2 = c.node2

        y = c.y(self.is_symbolic)
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

    def _add_K_tgn(self, M, S, v_nodes, i_nodes, coupling, index_col, index_row, i_graph_collapses, v_graph_collapses):
        """
        Adds a coupled inductor block to TGN matrix.
        """
        ind1 = coupling[0]
        ind2 = coupling[1]

        if self.is_symbolic:
            L1 = ind1.sym_value
            L2 = ind2.sym_value
        else:
            L1 = ind1.value
            L2 = ind2.value

        coeff = coupling[2] * sympy.sqrt(L1 * L2)

        node1 = ind1.node1
        node2 = ind1.node2
        node3 = ind2.node1
        node4 = ind2.node2

        n1v = self.index_tgn(v_nodes, node1, v_graph_collapses)
        n2v = self.index_tgn(v_nodes, node2, v_graph_collapses)
        n3v = self.index_tgn(v_nodes, node3, v_graph_collapses)
        n4v = self.index_tgn(v_nodes, node4, v_graph_collapses)

        n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
        n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)
        n3i = self.index_tgn(i_nodes, node3, i_graph_collapses)
        n4i = self.index_tgn(i_nodes, node4, i_graph_collapses)

        col1 = len(v_nodes) + index_col
        col2 = len(v_nodes) + index_col + 1
        row1 = len(i_nodes) + index_row
        row2 = len(i_nodes) + index_row + 1

        # L1 KVL row
        if n1v is not None: M[row1, n1v] += 1
        if n2v is not None: M[row1, n2v] += -1
        M[row1, col1] += -s * L1  # self-inductor
        M[row1, col2] += -s * coeff  # mutual term

        # L2 KVL row
        if n3v is not None: M[row2, n3v] += 1
        if n4v is not None: M[row2, n4v] += -1
        M[row2, col2] += -s * L2  # self-inductor
        M[row2, col1] += -s * coeff  # mutual term

        # Node current contributions (KCL)
        if n1i is not None: M[n1i, col1] += 1
        if n2i is not None: M[n2i, col1] += -1
        if n3i is not None: M[n3i, col2] += 1
        if n4i is not None: M[n4i, col2] += -1

        return ind1, ind2, row1, row2, coeff, L1, L2

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
            if val is None:
                val = c.sym_value

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
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = self._choose_source_val(c, form)

        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        result[self.c_count + index, 0] = val
        return matrix

    def _add_current_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            form = "sym"
        else:
            form = "rat"
        val = self._choose_source_val(c, form)
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
        if self.is_symbolic:
            L1 = c_L1.sym_value
            L2 = c_L2.sym_value
            coupling_coeff = c.sym_value
            coeff = coupling_coeff * sympy.sqrt(L1 * L2)
        else:
            L1 = c_L1.value
            L2 = c_L2.value
            coupling_coeff = c.value
            coeff = coupling_coeff * sympy.sqrt(L1 * L2)

        matrix[self.c_count + L2_index, self.c_count + L1_index] += -s * coeff
        matrix[self.c_count + L1_index, self.c_count + L2_index] += -s * coeff

        return coeff, c_L1, c_L2, L1_index, L2_index

    def _add_short(self, matrix, c, index):
        N1 = c.node1
        N2 = c.node2
        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        return matrix


class DC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

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

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.dc_val(form)

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.dc_val(form="sym")
                else:
                    ret = c.dc_val(form="rat")
        else:
            ret = sympy.limit(ret, s, 0)
        return ret


class TF(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        return eqn_matrix, solved_dict, symbols

    def component_current(self, name: str) -> sympy.Symbol:
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.ac_val(form="sym")
                else:
                    ret = c.ac_val(form="rat")
        return ret

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.tf_val(form)


class AC(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def _analyse(self):
        eqn_matrix, solved_dict, symbols = super()._analyse()
        for sym in symbols:
            try:
                solved_dict[sym] = solved_dict[sym].subs(s, 2 * pi * f * j)
            except:
                pass
        return eqn_matrix, solved_dict, symbols

    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.ac_val(form)

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.ac_val(form="sym")
                else:
                    ret = c.ac_val(form="rat")
        else:
            ret = ret.subs(s, 2*pi*f*j)
        return ret


class ACNumeric(AC):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):

        if use_symengine:
            os.environ["USE_SYMENGINE"] = "1"

        method = method.lower()
        if method not in ["tableau", "two_graph_node"]:
            raise ValueError(f"Nonexistent analysis method: {method}")

        self.is_symbolic: bool = symbolic
        self.auto_eval = auto_eval
        self.precision: int = precision
        self.method: str = method
        self.sympy_ilt: bool = sympy_ilt

        self.sbg = True

        self.circuit: Circuit = circuit
        self.node_voltage_identities: list = []

        self.symbol_dict = {}

        self.node_voltage_symbols: List[sympy.Symbol] = self.circuit.get_node_symbols()
        self.node_dict = self.circuit.get_node_dict()
        self.c_count = self.circuit.count_ports()
        self.node_count = self.circuit.count_nodes()

        self.eqn_matrix, self.symbols = self._build_system_eqn()


    def run(self, freqs):
        import numpy as np

        ac_eqn_matrix = self.eqn_matrix.subs(s, 2*j*pi*f)
        A_sym = ac_eqn_matrix[:, :-1]  # all columns except last
        b_sym = ac_eqn_matrix[:, -1]  # last column

        fA = sympy.lambdify(f, A_sym, "numpy")
        fb = sympy.lambdify(f, b_sym, "numpy")

        symbols_str = [str(sym) for sym in self.symbols]
        results = {sym: [] for sym in symbols_str}

        for freq in freqs:
            x = np.linalg.solve(fA(freq), fb(freq)).flatten()

            for sym, val in zip(symbols_str, x):
                results[sym].append(val)

        return results


class TRAN(Analysis):
    def __init__(self, circuit: Circuit, method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, solver: str = "gauss"):
        super().__init__(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)

    def get_node_voltage(self, node: str, force_s_domain=False) -> Union[sympy.Expr, None]:
        """
        Old way to return a node voltage, will be deprecated soon
        """
        func = super().get_node_voltage(node)
        if func is None or force_s_domain:
            return func
        else:
            res = laplace.iLT(func, self.sympy_ilt)
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

    def _add_K(self, matrix, vi_vector, c, inductor_index):
        coeff, c_L1, c_L2, L1_index, L2_index, = super()._add_K(matrix, vi_vector, c, inductor_index)
        if c_L1.init_cond is not None:
            vi_vector[self.c_count + L2_index, 0] += -coeff * c_L1.init_cond
        if c_L2.init_cond is not None:
            vi_vector[self.c_count + L1_index, 0] += -coeff * c_L2.init_cond

    def _add_basic_tgn(self, M, S, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses):
        super()._add_basic_tgn(M, S, v_nodes, i_nodes, c, i_graph_collapses, v_graph_collapses)

        # Only process if there is an actual initial condition
        if c.type in ["c", "l"] and c.init_cond is not None:
            if self.is_symbolic:
                val = c.sym_value
            else:
                val = c.value

            node1 = c.node1
            node2 = c.node2
            n1i = self.index_tgn(i_nodes, node1, i_graph_collapses)
            n2i = self.index_tgn(i_nodes, node2, i_graph_collapses)

            if c.type == "c":
                # Capacitor Norton current: C * V(0)
                if n1i is not None:
                    S[n1i, 0] += val * c.init_cond
                if n2i is not None:
                    S[n2i, 0] += -val * c.init_cond

            elif c.type == "l":
                if c.coupling is None: # If theres no coupling handle IC by Norton current IC model
                    # Inductor Norton current: I(0) / s
                    if n1i is not None:
                        S[n1i, 0] += -(1 / s) * c.init_cond
                    if n2i is not None:
                        S[n2i, 0] += (1 / s) * c.init_cond
                else: # Else let the coupling stamp handle initial conditions
                    pass

    def _add_K_tgn(self, M, S, v_nodes, i_nodes, coupling, index_col, index_row, i_graph_collapses, v_graph_collapses):
        ind1, ind2, row1, row2, m_coeff, L1, L2 = super()._add_K_tgn(M, S, v_nodes, i_nodes, coupling, index_col, index_row, i_graph_collapses, v_graph_collapses)
        # Initial conditions

        ic1 = ind1.init_cond
        ic2 = ind2.init_cond

        if ic1 is not None:
            S[row2, 0] += -m_coeff * ic1 - L2 * ic2
        if ic2 is not None:
            S[row1, 0] += -m_coeff * ic2 - L1 * ic1


    def _choose_source_val(self, c: Union[VoltageSource, CurrentSource], form="rat") -> sympy.Expr:
        return c.tran_val(form)

    def component_voltage(self, name: str) -> sympy.Expr:
        """
        Old way to return a component voltage, will be deprecated soon
        """
        ret = super().component_voltage(name)
        if ret is not None:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
        return ret

    def component_current(self, name: str) -> sympy.Symbol:
        """
        Old way to return a component current, will be deprecated soon
        """
        ret = super().component_current(name)
        if ret is None:
            c = self.circuit.components[name]
            if c.type == "i":
                if self.is_symbolic:
                    ret = c.tran_val(form="sym")
                else:
                    ret = c.tran_val(form="rat")
        if ret is not None:
            ret = laplace.iLT(ret, self.sympy_ilt)
            ret = sympy.factor_terms(ret)
        return ret


def AnalyseCircuit(netlist: str, analysis_type: str = "DC", method: str = "tableau",
                 symbolic: bool = True, auto_eval: bool=False, precision: int = 6, sympy_ilt: bool = True,
                 use_symengine: bool = False, operating_point: Union[Dict[str, float], None] = None,
                   solver: str = "gauss") -> Analysis:

    circuit = Circuit(netlist, operating_point)
    analysis_type = analysis_type.lower()
    if analysis_type == "dc":
        analysis = DC(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "ac":
        analysis = AC(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "tf":
        analysis = TF(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    elif analysis_type == "tran":
        analysis = TRAN(circuit, method, symbolic, auto_eval, precision, sympy_ilt, use_symengine, solver)
    else:
        raise ValueError(f"Nonexistent analysis type: {analysis_type}")
    return analysis