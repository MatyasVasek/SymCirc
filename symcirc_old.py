import sympy
from parse_old import parse
import os
import sys
from laplace import *


def load_file(netlist_addr):
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist


class AnalyseCircuit:
    def __init__(self, netlist, analysis_type="DC", symbolic=True):
        if analysis_type not in ["DC", "AC", "TF", "tran"]:
            #print("ERROR: nonexistent analysis type")
            sys.exit()
        self.is_symbolic = symbolic
        self.analysis_type = analysis_type
        self.s = sympy.symbols("s", real=True, positive=True)
        self.t = sympy.symbols("t", real=True, positive=True)
        self.netlist = netlist
        data = parse(netlist)
        self.components = data["components"]   # {<name> : <Component>} (see component.py)
        self.node_dict = data["node_dict"]  # {<node_name>: <index in equation matrix>, ...}
        self.matrix_size = data["matrix_size"]
        self.node_count = data["node_count"]
        self.eqn_matrix, self.solved_dict = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}

    def component_values(self, name):
        ret = {}
        if name == "all":
            ret = self.all_component_values()
            return ret
        else:
            v = "v({})".format(name)
            i = "i({})".format(name)
            ret[v] = self.component_voltage(name)
            ret[i] = self.component_current(name)
        return ret

    def all_component_values(self):
        ret = {}
        for key in self.components:
            name = self.components[key].name
            val = self.component_values(name)
            ret.update(val)
        return ret

    def node_voltages(self):
        ret = {}
        dic = self.solved_dict
        for key in dic:
            if key.name[0] not in ["i", "I"]:
                ret[key] = dic[key]
        return ret

    def transfer_function(self, node1, node2):
        v1 = sympy.Symbol("v({})".format(node1))
        v2 = sympy.Symbol("v({})".format(node2))
        voltage1 = self.solved_dict[v1].simplify()
        voltage2 = self.solved_dict[v2].simplify()
        tf = (voltage2/voltage1)
        return tf

    def component_voltage(self, name):
        c = self.components[name]
        v_n1 = sympy.Symbol("v({})".format(c.node1))
        v_n2 = sympy.Symbol("v({})".format(c.node2))
        if c.node1 == "0":
            a = 0
        elif v_n1 not in self.solved_dict:
            a = v_n1
        else:
            a = self.solved_dict[v_n1]

        if c.node2 == "0":
            b = 0
        elif v_n2 not in self.solved_dict:
            b = v_n2
        else:
            b = self.solved_dict[v_n2]
        voltage = sympy.simplify(a - b)
        #ret[c.sym_value] = voltage
        return voltage

    def component_current(self, name):
        c = self.components[name]
        current_symbol = sympy.Symbol("i({})".format(name))
        if current_symbol in self.solved_dict:
            current = self.solved_dict[current_symbol]
        else:
            if c.type == "c":
                current = 0
            elif c.type == "l":
                current = "infty"
                #print(current)
            else:
                voltage = self.component_voltage(name)
                if self.is_symbolic:
                    current = sympy.simplify(voltage/c.sym_value)
                else:
                    current = sympy.simplify(voltage/c.value)
        return current

    def freq_to_time(self, frequency_domain_result):
        collected_expr = frequency_domain_result.expand().collect(self.s).apart(self.s)
        #print("F = " + str(collected_expr))
        time_domain_result = sympy.inverse_laplace_transform(collected_expr, self.s, self.t)
        return time_domain_result

    def _analyse(self):
        eqn_matrix, symbols = self._build_system_eqn()
        #latex_print(eqn_matrix)
        #latex_print(symbols)
        solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
        #solved_dict = sympy.solve_linear_system_LU(eqn_matrix, symbols)

        if self.analysis_type == "DC":
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = sympy.limit(solved_dict[sym], self.s, 0).simplify()
                    except KeyError:
                        pass
            else:
                for sym in symbols:
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type == "v":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.dc_value)
                            else:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                    except KeyError:
                        pass

        elif self.analysis_type == "AC":
            f = sympy.symbols("f", real=True, positive=True)
            j = sympy.symbols("j")
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = solved_dict[sym].subs(self.s, 2*sympy.pi*f*j).simplify()
                    except KeyError:
                        pass
            else:
                for sym in symbols:
                    #print(sym)
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type == "v":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                #print(c.ac_value)
                            else:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)

                        solved_dict[sym] = solved_dict[sym].subs(self.s, 2*sympy.pi*f*j).simplify()
                    except KeyError:
                        pass

        elif self.analysis_type == "TF":
            for sym in symbols:
                try:
                    solved_dict[sym] = solved_dict[sym].simplify()
                except KeyError:
                    pass

        elif self.analysis_type == "tran":
            #latex_print(solved_dict)
            for sym in symbols:
                try:
                    #solved_dict[sym] = self.freq_to_time(solved_dict[sym].simplify())
                    pass
                except KeyError:
                    pass

        return eqn_matrix, solved_dict

    def _build_system_eqn(self):
        admittance_matrix = self._admit_matrix_init()
        source_matrix = self._current_matrix_init()
        symbols = self._voltage_symbols()
        i = 0
        priority = []
        non_priority = []
        for key in self.components:
            if key[0] in ["F", "G", "E", "A"]:
                priority.append(key)
            else:
                non_priority.append(key)
        keys = priority + non_priority
        for key in keys:
            c = self.components[key]
            if i == len(self.components):
                pass
            elif c.type == "r":
                self._add_R(admittance_matrix, c)
            elif c.type == "c":
                self._add_C(admittance_matrix, c)
            elif c.type == "l":
                self._add_L(admittance_matrix, c)
            elif c.type == "v":
                self._add_V(admittance_matrix, source_matrix, c)
                symbols.append(sympy.Symbol("i({})".format(c.name)))
            elif c.type == "i":
                self._add_I(source_matrix, c)
            elif c.type == "a":
                self._add_A(admittance_matrix, c)
                symbols.append(sympy.Symbol("i({})".format(c.name)))
            elif c.type == "e":
                self._add_VVT(admittance_matrix, c)
                symbols.append(sympy.Symbol("i({})".format(c.name)))
            elif c.type == "g":
                self._add_VCT(admittance_matrix, c)
            elif c.type == "short":
                self._add_short(admittance_matrix, c)
                symbols.append(sympy.Symbol("i({})".format(c.shorted_component)))
            elif c.type == "f":
                self._add_CCT(admittance_matrix, source_matrix, c)
                symbols.append(sympy.Symbol("i({})".format(c.name)))
            elif c.type == "h":
                self._add_CVT(admittance_matrix, c)
                symbols.append(sympy.Symbol("i({}control)".format(c.name)))
                symbols.append(sympy.Symbol("i({})".format(c.name)))

            i += 1
        equation_matrix = admittance_matrix.col_insert(self.matrix_size + 1, source_matrix)
        n = 0
        for node in self.node_dict:
            if node == "0":
                n = self.node_dict[node]
        #latex_print(symbols)
        #latex_print(equation_matrix)
        equation_matrix.col_del(n)
        equation_matrix.row_del(n)
        source_matrix.row_del(n)
        #latex_print(symbols)
        #latex_print(equation_matrix)
        #sympy.pprint(equation_matrix)
        #print(symbols)

        return equation_matrix, symbols

    def _admit_matrix_init(self):
        return sympy.SparseMatrix(sympy.zeros(self.matrix_size))

    def _voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol("v({})".format(node)))
        return voltage_symbol_list

    def _current_matrix_init(self):
        return sympy.zeros(self.matrix_size, 1)

    # resistor element
    def _add_R(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        admittance = 1/c.sym_value
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # capacitor element
    def _add_C(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        admittance = self.s * c.sym_value
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # inductor element
    def _add_L(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        admittance = 1 / (self.s * c.sym_value)
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # ideal voltage source element
    def _add_V(self, admittance_matrix, source_matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        pos = self.node_count + c.position
        voltage = c.sym_value
        admittance_matrix[pos, N1] += 1
        admittance_matrix[pos, N2] -= 1
        admittance_matrix[N1, pos] += 1
        admittance_matrix[N2, pos] -= 1
        source_matrix[pos, 0] += voltage

    # ideal current source element
    def _add_I(self, source_matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        source_matrix[N1, 0] -= c.sym_value
        source_matrix[N2, 0] += c.sym_value

    # ideal (A = infinity) operational amplifier
    def _add_A(self, admittance_matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        N3 = self.node_dict[c.node3]
        N4 = self.node_dict[c.node4]
        pos = self.node_count + c.position
        admittance_matrix[pos, N3] -= 5
        admittance_matrix[pos, N4] += 5
        admittance_matrix[N1, pos] += 5
        admittance_matrix[N2, pos] -= 5

    def _add_short(self, admittance_matrix, c):
        pos = self.node_count + c.position
        #print("{}, short".format(pos))
        if c.node1 is None:
            s_c = self.components[c.shorted_component]
            n_p = self.node_dict[s_c.node1]
            n_n = self.node_dict[s_c.node2]
        else:
            n_p = self.node_dict[c.node1]
            n_n = self.node_dict[c.node2]

        admittance_matrix[pos, n_p] += 1
        admittance_matrix[pos, n_n] -= 1
        admittance_matrix[n_p, pos] += 1
        admittance_matrix[n_n, pos] -= 1

    def _add_VVT(self, admittance_matrix, c):
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        NC_p = self.node_dict[c.node3]
        NC_n = self.node_dict[c.node4]
        pos = self.node_count + c.position
        gain = c.sym_value
        admittance_matrix[pos, NC_p] -= gain
        admittance_matrix[pos, NC_n] += gain
        admittance_matrix[pos, N_p] += 1
        admittance_matrix[pos, N_n] -= 1
        admittance_matrix[N_p, pos] += 1
        admittance_matrix[N_n, pos] -= 1

    def _add_VCT(self, admittance_matrix, c):
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        NC_p = self.node_dict[c.node3]
        NC_n = self.node_dict[c.node4]
        gain = c.sym_value
        admittance_matrix[N_p, NC_p] -= gain
        admittance_matrix[N_p, NC_n] += gain
        admittance_matrix[N_n, NC_p] -= gain
        admittance_matrix[N_n, NC_n] += gain

    def _add_CCT(self, admittance_matrix, source_matrix, c):
        # TODO finish
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        c_v = self.components[c.control_voltage]
        #print(c_v)
        #print(c_v.node1)
        #print(c_v.node2)
        NC_p = self.node_dict[c_v.node1]
        NC_n = self.node_dict[c_v.node2]
        #print(NC_p)
        #print(NC_n)
        pos = self.node_count + c.position
        #print("{}, CCT".format(pos))
        gain = c.sym_value
        admittance_matrix[pos, NC_p] += 1
        admittance_matrix[pos, NC_n] -= 1
        admittance_matrix[NC_p, pos] += 1
        admittance_matrix[NC_n, pos] -= 1
        admittance_matrix[N_p, pos] += gain
        admittance_matrix[N_n, pos] -= gain

    def _add_CVT(self, admittance_matrix, c):
        # TODO finish
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        NC_p = self.node_dict[c.node3]
        NC_n = self.node_dict[c.node4]
        pos1 = self.node_count + c.position - 1
        pos2 = self.node_count + c.position
        gain = c.sym_value
        admittance_matrix[pos1, NC_p] += 1
        admittance_matrix[pos1, NC_n] -= 1
        admittance_matrix[pos2, N_p] += 1
        admittance_matrix[pos2, N_n] -= 1
        admittance_matrix[NC_p, pos1] += 1
        admittance_matrix[NC_n, pos1] -= 1
        admittance_matrix[N_p, pos2] += 1
        admittance_matrix[N_n, pos2] -= 1
        admittance_matrix[pos2, pos1] -= gain

    # TODO add Coupled Inductors

    # TODO add Transistor


def latex_print(data):
    print("{}".format(sympy.latex(data)))


if __name__ == "__main__":
    dirname = os.path.dirname(__file__)
    netlist = dirname + "/oamp.txt"
    circuit = AnalyseCircuit(netlist)
    latex_print(circuit.eqn_matrix)
    latex_print(circuit.solved_dict)
    latex_print(circuit.transfer_function("1", "2"))



