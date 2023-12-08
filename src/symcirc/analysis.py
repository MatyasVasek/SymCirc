import sympy
from symcirc import parse, laplace, utils
from symcirc.utils import j,s,t
from symcirc.pole_zero import *

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
    def __init__(self, netlist, analysis_type="DC", method="tableau", phases="undefined", symbolic=True, precision=6):
        if analysis_type not in ["DC", "AC", "TF", "tran"]:
            raise ValueError("Nonexistent analysis type: {}".format(analysis_type))
        self.is_symbolic = symbolic
        self.analysis_type = analysis_type
        self.precision = precision
        self.method = method
        self.netlist = netlist
        if analysis_type == "tran":
            data = parse.parse(netlist, tran=True)
        else:
            data = parse.parse(netlist)

        self.phases = self.parse_phases(phases)
        if self.phases != "undefined" and analysis_type not in ["AC", "TF", "tran"]:
            raise ValueError("This type of analysis is not allowed in Switched Capacitors/Currents circuits")
        if self.phases != "undefined" and method != "modified_node":
            raise ValueError('The only valid method for Switched Capacitor circuits is "modified_node"')

        self.components = data["components"] # {<name> : <Component>} (see component.py)
        self.node_dict = data["node_dict"]  # {<node_name>: <index in equation matrix>, ...}
        self.matrix_size = data["matrix_size"]
        self.node_count = data["node_count"]
        self.c_count = self.count_components()
        self.node_voltage_symbols = self._node_voltage_symbols()
        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}

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
                        phase_definition.append(int(phases))
                        for i in range(phase_definition[0]):
                            phase_definition.append(sympy.sympify("1/" + str(phase_definition[0])))
                else:
                    raise SyntaxError(phase_def_syntax_error)

            else:
                raise SyntaxError(phase_def_syntax_error)
            phases = phase_definition
        return phases

    def component_voltage(self, name):
        ret = None
        v = "v({})".format(name)
        ret = self.solved_dict[sympy.symbols(v)]
        return ret

    def component_current(self, name):
        ret = None
        i = "i({})".format(name)
        ret = self.solved_dict[sympy.symbols(i)]
        return ret

    def component_values(self, name="all", default_python_datatypes=False):
        """
          Takes a string containing a single component name and returns a dictionary containing the voltage and current
            of the input component.

          :param str name: component id
          :return dict ret: in format {"v(name)" : value, "i(name)" : value}
        """
        ret = {}
        if name == "all":
            ret = self.all_component_values(default_python_datatypes)
            return ret
        else:
            v = "v({})".format(name)
            i = "i({})".format(name)
            ret[v] = self.solved_dict[sympy.symbols(v)]
            ret[i] = self.solved_dict[sympy.symbols(i)]
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
        dic = self.solved_dict
        for key in dic:
            if key in self.node_voltage_symbols and key.name[2] != "*":
                ret[key.name] = dic[key]
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
                v1 = sympy.Symbol("v({})".format(node1))
                voltage1 = self.solved_dict[v1].simplify()
        except KeyError:
            print("Node {} doesn't exist.".format(v1))
            exit(101)

        try:
            if node2 in [0, "0"]:
                voltage2 = 0
            else:
                v2 = sympy.Symbol("v({})".format(node2))
                voltage2 = self.solved_dict[v2].simplify()
        except KeyError:
            print("Node {} doesn't exist.".format(v2))
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
            if self.analysis_type == "DC":
                eqn_matrix, symbols = self._build_system_eqn()
                #print(symbols)
                solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
                #print(solved_dict)
                if self.is_symbolic:
                    for sym in symbols:
                        try:
                            solved_dict[sym] = sympy.limit(solved_dict[sym], s, 0)
                        except KeyError:
                            pass
                        except TypeError:
                            pass
                            #print(solved_dict)
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
                                solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass

            elif self.analysis_type == "AC":
                eqn_matrix, symbols = self._build_system_eqn()
                solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
                f = sympy.symbols("f", real=True, positive=True)
                if self.is_symbolic:
                    i = 1
                    for sym in symbols:
                        try:
                            #print("{}: {}".format(sym, solved_dict[sym]))
                            #solved_dict[sym] = solved_dict[sym].simplify()
                            #print("{}: {}".format(sym, solved_dict[sym]))
                            if i == 1:
                                #sympy.pprint(solved_dict)
                                i = 0
                            solved_dict[sym] = solved_dict[sym].subs(s, 2 * sympy.pi * f * j)
                        except KeyError:
                            pass
                else:
                    for sym in symbols:
                        #print(sym)
                        try:
                            #print(solved_dict[sym])
                            #solved_dict[sym] = solved_dict[sym].simplify()
                            #print(solved_dict[sym])
                            solved_dict[sym] = solved_dict[sym].subs(s, 2 * sympy.pi * f * j)
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["v", "i"]:
                                    if c.ac_value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                    #print(c.ac_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                                solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass

            elif self.analysis_type == "TF":
                eqn_matrix, symbols = self._build_system_eqn()
                solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
                #solved_dict = sympy.solve_linear_system_LU(eqn_matrix, symbols)

                if not self.is_symbolic:
                    for sym in symbols:
                        #print(sym)
                        try:
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["v", "i"]:
                                    #print(c.ac_value)
                                    if c.ac_value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                    #print(c.ac_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                                solved_dict[sym] = solved_dict[sym].evalf(self.precision)
                        except KeyError:
                            pass

            elif self.analysis_type == "tran":
                eqn_matrix, symbols = self._build_system_eqn()
                solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
                #print(solved_dict)
                #print(symbols)
                if self.is_symbolic:
                    #latex_print(solved_dict)
                    for sym in symbols:
                        #print(sym)
                        try:
                            for name in self.components:
                                c = self.components[name]
                                if c.type in ["i", "v"] and c.name[-3:] != "_IC":
                                    #print(c.name)
                                    #print(c.sym_value)
                                    #print("Before: {}".format(solved_dict[sym]))
                                    solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.tran_value)
                                    #print("After: {}".format(solved_dict[sym]))
                                    # print(c.ac_value)

                        except KeyError:
                            pass
                        #solved_dict[sym] = sympy.apart(solved_dict[sym], self.s)
                        #print(solved_dict[sym])
                        #print("{}: {}".format(sym, solved_dict[sym]))

                        inv_l = laplace.iLT(solved_dict[sym])
                        #inv_l = laplace.residue_laplace(solved_dict[sym])

                        #inv_l = sympy.simplify(inv_l)
                        solved_dict[sym] = inv_l
                        #print("{} = {}".format(sym, inv_l))
                        #print("{}: {}".format(sym, solved_dict[sym]))
                        try:
                            for name in self.components:
                                c = self.components[name]

                            #solved_dict[sym] = solved_dict
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
                                    # print(c.ac_value)
                                else:
                                    if c.value:
                                        solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                                #solved_dict[sym] = utils.evaluate(solved_dict[sym], self.precision)
                        except KeyError:
                            pass
                        f = laplace.iLT(solved_dict[sym])
                        solved_dict[sym] = utils.evaluate(f, self.precision)
        else:
            if self.analysis_type == "TF":
                eqn_matrix, symbols = self._build_system_eqn()
                solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
                self.SCSI_symbol_z_factor(solved_dict)
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
                self.SCSI_solved_dict_modifier(solved_dict)
        return eqn_matrix, solved_dict, symbols

    def _node_voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol("v({})".format(node)))
        return voltage_symbol_list

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
        for phase in range(num_of_phases):
            for component in range(int(phase * node_length / num_of_phases), int((phase + 1) * node_length / num_of_phases)):
                symbols.append(node_symbols[component])
            for component in range(phase, charge_length, num_of_phases):
                symbols.append(charge_symbols[component])
        return symbols

    def SCSI_symbol_z_factor(self, symbols):
        num_of_phases = self.phases[0]
        phase_helper_array = []
        z = sympy.Symbol("z")
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
        return symbols

    def SCSI_solved_dict_modifier(self, solved_dict):
        num_of_phases = self.phases[0]
        phase_helper_array = []
        z = sympy.Symbol("z")
        node = "undefined"
        temp = 0
        for phase in range(num_of_phases):
            phase_helper_array.append(self.phases[phase + 1])
        for key in self.components:
            c = self.components[key]
            if c.type == "v":
                node = c.node1
                break
        for key in solved_dict:
            node_index = str(key).find(str(node))
            if node_index != -1 and str(key)[node_index - 1] == "(" and str(key)[node_index + len(str(node))] == ")":
                symbol_phase = int(str(key)[-1])
                for phase in range(symbol_phase - 1):
                    temp += phase_helper_array[phase]
                solved_dict[key] *= z ** temp
                temp = 0
        return solved_dict




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
                        self._add_basic(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        if c.type == "l":
                            inductor_index[c.name] = index

                    if c.type == "v":
                        self._add_voltage_source(M, R, c, index)
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                    if c.type == "i":
                        self._add_current_source(M, R, c, index)
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                    if c.type == "g":
                        self._add_VCT(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({}_control)".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({}_control)".format(c.name)))
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        index += 1
                    if c.type == "e":
                        self._add_VVT(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({}_control)".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({}_control)".format(c.name)))
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        index += 1
                    if c.type == "f":
                        self._add_CCT(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({}_control)".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({}_control)".format(c.name)))
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        index += 1
                    if c.type == "h":
                        self._add_CVT(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({}_control)".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({}_control)".format(c.name)))
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        index += 1
                    if c.type == "a":
                        self._add_A(M, c, index)
                        voltage_symbols.append(sympy.Symbol("v({}_control)".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({}_control)".format(c.name)))
                        voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                        current_symbols.append(sympy.Symbol("i({})".format(c.name)))
                        index += 1

                    index += 1
            for coupling in couplings:
                self._add_K(M, coupling, inductor_index)
            equation_matrix = M.col_insert(self.c_count*2 + self.node_count, R)
            symbols = voltage_symbols + current_symbols + node_symbols

        elif self.method == "modified_node" and self.phases != "undefined":
            num_of_phases = self.phases[0]
            component_size = self.c_count * num_of_phases
            node_size = self.node_count * num_of_phases
            M = sympy.Matrix(sympy.zeros(component_size + node_size))
            R = sympy.Matrix(sympy.zeros(component_size + node_size, 1))

            component_index = 0
            node_symbols = self.SCSI_node_voltage_symbols()
            charge_symbols = []

            for key in self.components:
                c = self.components[key]
                if c.type == "v":
                    for phase in range(num_of_phases):
                        self.SCSI_add_voltage_source(M, R, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 1
                elif c.type == "c":
                    for phase_y in range(num_of_phases):
                        for phase_x in range(num_of_phases):
                            self.SCSI_add_capacitor(M, c, component_index, phase_y, phase_x)
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase_y + 1)))
                    component_index += 1
                elif c.type == "s":
                    self.SCSI_add_switch(M, c, component_index)
                    for phase in range(num_of_phases):
                        charge_symbols.append(sympy.Symbol("q({name})_{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 1
                elif c.type == "a":
                    for phase in range(num_of_phases):
                        self.SCSI_add_OpAmp(M, c, component_index, phase)
                        charge_symbols.append(sympy.Symbol("q({name})_in{phase}".format(name=c.name, phase=phase + 1)))
                        charge_symbols.append(sympy.Symbol("q({name})_out{phase}".format(name=c.name, phase=phase + 1)))
                    component_index += 2

            self.SCSI_matrix_z_symbol(M)
            equation_matrix = M.col_insert(component_size + node_size, R)
            symbols = self.SCSI_symbol_list_order(node_symbols, charge_symbols)
            # print(sympy.shape(equation_matrix))
            # print(len(symbols))
            sympy.pprint(equation_matrix)
            sympy.pprint(symbols)

        elif self.method == "two_graph_node" and self.phases == "undefined":
            symbols = []
            v_graph_collapses = {'0': []}
            i_graph_collapses = {'0': []}
            v_graph_nodes = []
            i_graph_nodes = []

            for key in self.components:
                c = self.components[key]
                if c.type in ["r", "l", "c"]:
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)

                if c.type == "v":
                    self.graph_append(c.node1, v_graph_nodes)
                    self.graph_append(c.node2, v_graph_nodes)
                    self.graph_append(c.node1, i_graph_nodes)
                    self.graph_append(c.node2, i_graph_nodes)
                    if c.node1 in i_graph_collapses:
                        i_graph_collapses[c.node1].append(c.node2)  # set node2 to be collapsed into node1 on the current graph
                    elif c.node2 in i_graph_collapses:
                        i_graph_collapses[c.node2].append(c.node1)    # set node1 to be collapsed into node2 on the current graph
                    else:
                        i_graph_collapses[c.node1] = [c.node2]  # set node2 to be collapsed into node1 on the current graph

                if c.type == "i":
                    pass
                if c.type == "g":
                    pass
                if c.type == "e":
                    pass
                if c.type == "f":
                    pass
                if c.type == "h":
                    pass
                if c.type == "a":
                    pass
                if c.type == "s":
                    pass

            """Collapse nodes based on collapse dictionaries"""
            for main_node in i_graph_collapses:
                collapse_list = i_graph_collapses[main_node]
                i = 0
                for n in i_graph_nodes:
                    if n in collapse_list:
                        if main_node == '0':
                            i_graph_nodes.remove(n)
                        else:
                            i_graph_nodes[i] = main_node
                    i+=1

            v_graph_nodes = list(set(v_graph_nodes))
            i_graph_nodes = list(set(i_graph_nodes))
            # print(f"V-Graph: {v_graph_nodes}")
            # print(f"I-Graph: {i_graph_nodes}")

            M = sympy.Matrix(sympy.zeros(len(v_graph_nodes)))
            S = sympy.Matrix(sympy.zeros(len(v_graph_nodes), 1))

            # sympy.pprint(M)
            # sympy.pprint(S)
            index = 0
            for key in self.components:
                c = self.components[key]
                if c.type in ["r", "l", "c"]:
                    #print(c.value)
                    self._add_basic_tgn(M, v_graph_nodes, i_graph_nodes, c)

                if c.type == "v":
                    self._add_V_tgn(M, S, v_graph_nodes, i_graph_nodes, c, index)
                    index+=1

                if c.type == "i":
                    raise TypeError("Current source not supported in 'two_graph_node' analysis")
                if c.type == "g":
                    raise TypeError("Controlled sources not supported in 'two_graph_node' analysis")
                if c.type == "e":
                    raise TypeError("Controlled sources not supported in 'two_graph_node' analysis")
                if c.type == "f":
                    raise TypeError("Controlled sources not supported in 'two_graph_node' analysis")
                if c.type == "h":
                    raise TypeError("Controlled sources not supported in 'two_graph_node' analysis")
                if c.type == "a":
                    raise TypeError("OpAmp not supported in 'two_graph_node' analysis")
                if c.type == "s":
                    raise TypeError("Switch not supported in 'two_graph_node' analysis")
            # sympy.pprint(M)
            # sympy.pprint(S)
            #print(v_graph_nodes)

            equation_matrix = M.col_insert(len(v_graph_nodes), S)
            #sympy.pprint(equation_matrix)

            for node in v_graph_nodes:
                symbols.append(sympy.Symbol(f"v({node})"))
            # print(symbols)
            sympy.pprint(equation_matrix)

        elif self.method == "two_graph_node" and self.phases != "undefined":  # switched capacitor analysis
            pass

        return equation_matrix, symbols

    def SCSI_submatrix_write(self, matrix, submatrix, start_y, start_x, phase_y, phase_x):
        y_dimension = sympy.shape(submatrix)[0]
        x_dimension = sympy.shape(submatrix)[1]
        phase_y_offset = (self.c_count + self.node_count) * phase_y
        phase_x_offset = (self.c_count + self.node_count) * phase_x
        for y in range(y_dimension):
            for x in range(x_dimension):
                matrix[phase_y_offset + start_y + y, phase_x_offset + start_x + x] += submatrix[y, x]
        return matrix

    def SCSI_incidence_matrix_generate(self, incidence_matrix, N1, N2, component_index):
        if N1 == "0":
            pass
        else:
            node_pos = self.node_dict[N1]
            incidence_matrix[node_pos, component_index] = 1
        if N2 == "0":
            pass
        else:
            node_pos = self.node_dict[N2]
            incidence_matrix[node_pos, component_index] = -1
        return incidence_matrix

    def SCSI_matrix_z_symbol(self, matrix):
        num_of_phases = self.phases[0]
        submatrix_dimension = self.node_count + self.c_count
        phase_offset = submatrix_dimension * (num_of_phases - 1)
        for y in range(submatrix_dimension):
            for x in range(submatrix_dimension):
                matrix[phase_offset + y, phase_offset + x] *= sympy.Symbol("z")
        return matrix

    def SCSI_add_capacitor(self, matrix, c, component_index, phase_y, phase_x):
        num_of_phases = self.phases[0]
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
        y_b = val
        z_b = -1

        Y_2 = sympy.Matrix(sympy.zeros(self.c_count))
        Z_2 = sympy.Matrix(sympy.zeros(self.c_count))
        A_2 = sympy.Matrix(sympy.zeros(self.node_count, self.c_count))
        if phase_y == phase_x:
            Y_2[component_index, component_index] = y_b
            Z_2[component_index, component_index] = z_b
            self.SCSI_incidence_matrix_generate(A_2, N1, N2, component_index)
        else:
            if (abs(phase_y - phase_x) == 1 and phase_y > phase_x) or (phase_y == 0 and phase_x + 1 == num_of_phases):
                Y_2[component_index, component_index] = - y_b
                Z_2[component_index, component_index] = - z_b
                self.SCSI_incidence_matrix_generate(A_2, N1, N2, component_index)
        if phase_y == phase_x:
            self.SCSI_submatrix_write(matrix, A_2, 0, self.node_count, phase_y, phase_x)
        else:
            self.SCSI_submatrix_write(matrix, - A_2, 0, self.node_count, phase_y, phase_x)
        self.SCSI_submatrix_write(matrix, Z_2, self.node_count, self.node_count, phase_y, phase_x)
        self.SCSI_submatrix_write(matrix, Y_2 * A_2.T, self.node_count,0 , phase_y, phase_x)
        return matrix

    def SCSI_add_voltage_source(self, matrix, result, c, component_index, phase):
        phase_offset = self.c_count * phase
        node_offset = self.node_count * (phase + 1)
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
        Y_2 = sympy.Matrix(sympy.zeros(self.c_count))
        A_2 = sympy.Matrix(sympy.zeros(self.node_count, self.c_count))
        Y_2[component_index, component_index] = 1
        self.SCSI_incidence_matrix_generate(A_2, N1, N2, component_index)
        self.SCSI_submatrix_write(matrix, A_2, 0, self.node_count, phase, phase)
        self.SCSI_submatrix_write(matrix, Y_2 * A_2.T, self.node_count, 0, phase, phase)
        result[phase_offset + node_offset + component_index, 0] = val
        return matrix

    def SCSI_add_switch(self, matrix, c, component_index):
        num_of_phases = self.phases[0]
        N1 = c.node1
        N2 = c.node2
        switch_phase = (int(c.phase)) - 1
        for phase in range(num_of_phases):
            A_2 = sympy.Matrix(sympy.zeros(self.node_count, self.c_count))
            self.SCSI_incidence_matrix_generate(A_2, N1, N2, component_index)
            self.SCSI_submatrix_write(matrix, A_2, 0, self.node_count, phase, phase)
            if phase == switch_phase:
                Y_2 = sympy.Matrix(sympy.zeros(self.c_count))
                Y_2[component_index, component_index] = 1
                self.SCSI_submatrix_write(matrix, Y_2 * A_2.T, self.node_count, 0, phase, phase)
            else:
                Z_2 = sympy.Matrix(sympy.zeros(self.c_count))
                Z_2[component_index, component_index] = -1  # ??? why not 1
                self.SCSI_submatrix_write(matrix, Z_2, self.node_count, self.node_count, phase, phase)
        return matrix

    def SCSI_add_OpAmp(self, matrix, c, component_index, phase):
        N1 = c.node1
        N2 = c.node2
        N3 = c.node3
        N4 = c.node4

        A_2 = sympy.Matrix(sympy.zeros(self.node_count, self.c_count))
        self.SCSI_incidence_matrix_generate(A_2, N3, N4, component_index)
        self.SCSI_submatrix_write(matrix, A_2, 0, self.node_count, phase, phase)
        Y_2 = sympy.Matrix(sympy.zeros(self.c_count))
        Y_2[component_index, component_index] = 1
        self.SCSI_submatrix_write(matrix, Y_2 * A_2.T, self.node_count, 0, phase, phase)

        A_2 = sympy.Matrix(sympy.zeros(self.node_count, self.c_count))
        self.SCSI_incidence_matrix_generate(A_2, N1, N2, component_index + 1)
        self.SCSI_submatrix_write(matrix, A_2, 0, self.node_count, phase, phase)
        Z_2 = sympy.Matrix(sympy.zeros(self.c_count))
        Z_2[component_index + 1, component_index] = 1
        self.SCSI_submatrix_write(matrix, Z_2, self.node_count, self.node_count, phase, phase)
        return matrix

    def _add_V_tgn(self, M, S, v_nodes, i_nodes, c, index):
        node1 = c.node1
        node2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value

        n1v = self.index_tgn(v_nodes, node1)
        n2v = self.index_tgn(v_nodes, node2)
        row = len(i_nodes)+index
        #print("ok")
        #print(n1v)
        #print(n2v)
        if n1v is not None:
            M[row, n1v] += 1
        if n2v is not None:
            M[row, n2v] += -1
        S[row, 0] = val

    def _add_basic_tgn(self, M, v_nodes, i_nodes, c):
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
            y = val*s

        n1v = self.index_tgn(v_nodes, node1)
        n2v = self.index_tgn(v_nodes, node2)
        n1i = self.index_tgn(i_nodes, node1)
        n2i = self.index_tgn(i_nodes, node2)
        #print(f"{n1v}, {n2v}, {n1i}, {n2i}")
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

    def index_tgn(self, v_nodes, node):
        try:
            i = v_nodes.index(node)
        except ValueError:
            i = None
        return i

    def graph_append(self, node, graph):
        if node == '0':
            pass
        elif node not in graph:
            graph.append(node)
        return graph

    def _incidence_matrix_write(self, N1, N2, matrix, index):
        if N1 == "0":
            pass
        else:
            val = 1
            node_pos = self.node_dict[N1]
            matrix[self.c_count*2 + node_pos, self.c_count + index] += val
            matrix[index, self.c_count*2 + node_pos] -= val
        if N2 == "0":
            pass
        else:
            val = -1
            node_pos = self.node_dict[N2]
            matrix[self.c_count*2 + node_pos, self.c_count + index] += val
            matrix[index, self.c_count*2 + node_pos] -= val
        return matrix

    def _add_basic(self, matrix, c, index):
        N1 = c.node1
        N2 = c.node2
        y_b = 0
        z_b = 0
        if self.is_symbolic:
            val = c.sym_value
        else:
            val = c.value
            #print(val)
        #print("{}: {}".format(c.name, val))

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
            val = c.sym_value
        else:
            if self.analysis_type == "DC":
                val = c.dc_value
            elif self.analysis_type == "tran":
                val = c.tran_value
            else:
                val = c.ac_value
        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        result[self.c_count + index, 0] = val
        return matrix

    def _add_current_source(self, matrix, result, c, index):
        N1 = c.node1
        N2 = c.node2
        if self.is_symbolic:
            val = c.sym_value
            #print("symbolic: {}".format(val))
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

    def _add_K(self, matrix, c, inductor_index):  # coupled inductors
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
            M = c.sym_value
        else:
            L1 = c_L1.value
            L2 = c_L2.value
            M = c.value

        matrix[self.c_count + L2_index, self.c_count + L1_index] += -s*M*sympy.sqrt(L1*L2)
        matrix[self.c_count + L1_index, self.c_count + L2_index] += -s*M*sympy.sqrt(L1*L2)

    def _add_short(self, matrix, c, index):
        N1 = c.node1
        N2 = c.node2
        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        return matrix