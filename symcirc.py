import sympy
from parse import parse
import os
import sys
from laplace import *


class AnalyseCircuit:
    def __init__(self, netlist, analysis_type="DC", symbolic=True):
        if analysis_type not in ["DC", "AC", "TF", "tran"]:
            print("ERROR: nonexistent analysis type")
            sys.exit()
        self.symbolic = symbolic
        self.analysis_type = analysis_type
        self.s = sympy.symbols("s", real=True)
        self.t = sympy.symbols("t", real=True)
        self.netlist = netlist
        data = parse(netlist)
        self.components = data["components"]   # {<name> : <Component>} (see component.py)
        self.node_dict = data["node_dict"]  # {<node_name>: <index in equation matrix>, ...}
        self.matrix_size = data["matrix_size"]
        self.node_count = data["node_count"]
        self.eqn_matrix, self.solved_dict = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}

    def analyse_component(self, name):
        ret = {}
        if name == "all":
            for key in self.components:
                name = self.components[key].name
                ret[name] = self.analyse_component(name)
            return ret
        else:
            v = "v{}".format(name)
            i = "i{}".format(name)
            ret[v] = self.component_voltage(name)
            ret[i] = self.component_current(name)
        return ret

    def analyse_all(self):
        ret = {}
        for key in self.components:
            name = self.components[key].name
            ret[name] = self.analyse_component(name)
        return ret

    def transfer_function(self, node1, node2):
        v1 = sympy.Symbol("v{}".format(node1))
        v2 = sympy.Symbol("v{}".format(node2))
        voltage1 = self.solved_dict[v1].simplify()
        voltage2 = self.solved_dict[v2].simplify()
        tf = (voltage2/voltage1)
        return tf

    def component_voltage(self, name):
        c = self.components[name]
        v_n1 = sympy.Symbol("v{}".format(c.node1))
        v_n2 = sympy.Symbol("v{}".format(c.node2))
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
        current_symbol = sympy.Symbol("i{}".format(name))
        if current_symbol in self.solved_dict:
            current = self.solved_dict[current_symbol]
        else:
            voltage = self.component_voltage(name)
            if self.symbolic:
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
        #t(symbols)
        #solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
        solved_dict = sympy.solve_linear_system_LU(eqn_matrix, symbols)

        if self.analysis_type == "DC":
            for sym in symbols:
                try:
                    solved_dict[sym] = sympy.limit(solved_dict[sym], self.s, 0).simplify()
                except KeyError:
                    pass

        elif self.analysis_type == "AC":
            pi, f = sympy.symbols("pi f", real=True, positive=True)
            j = sympy.symbols("j")
            for sym in symbols:
                try:
                    solved_dict[sym] = solved_dict[sym].subs(self.s, 2*pi*f*j).simplify()
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
                    solved_dict[sym] = self.freq_to_time(solved_dict[sym].simplify())
                except KeyError:
                    pass

        return eqn_matrix, solved_dict

    def _build_system_eqn(self):
        admittance_matrix = self._admit_matrix_init()
        source_matrix = self._current_matrix_init()
        symbols = self._voltage_symbols()
        i = 0
        for key in self.components:
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
                symbols.append(sympy.Symbol("i{}".format(c.name)))
            elif c.type == "i":
                self._add_I(source_matrix, c)
            elif c.type == "a":
                self._add_A(admittance_matrix, source_matrix, c)
                symbols.append(sympy.Symbol("i{}".format(c.name)))
            elif c.type == "e":
                self._add_VVT(admittance_matrix, c)
                symbols.append(sympy.Symbol("i{}".format(c.name)))
            elif c.type == "g":
                self._add_VCT(admittance_matrix, c)
            elif c.type == "f":
                self._add_CCT(admittance_matrix, source_matrix, c)
                symbols.append(sympy.Symbol("i{}".format(c.name)))
            elif c.type == "h":
                self._add_CVT(admittance_matrix, c)
                symbols.append(sympy.Symbol("i{}control".format(c.name)))
                symbols.append(sympy.Symbol("i{}".format(c.name)))

            i += 1
        equation_matrix = admittance_matrix.col_insert(self.matrix_size + 1, source_matrix)
        n = 0
        for node in self.node_dict:
            if node == 0:
                n = self.node_dict[node]
        #latex_print(symbols)
        #latex_print(equation_matrix)
        equation_matrix.col_del(n)
        equation_matrix.row_del(n)
        source_matrix.row_del(n)
        #latex_print(symbols)
        #latex_print(equation_matrix)

        return equation_matrix, symbols

    def _admit_matrix_init(self):
        return sympy.SparseMatrix(sympy.zeros(self.matrix_size))

    def _voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != 0:
                voltage_symbol_list.append(sympy.Symbol("v{}".format(node)))
        return voltage_symbol_list

    def _current_matrix_init(self):
        return sympy.zeros(self.matrix_size, 1)

    # resistor element
    def _add_R(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        if self.symbolic:
            admittance = 1/c.sym_value
        else:
            admittance = 1/c.value
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # capacitor element
    def _add_C(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        if self.symbolic:
            admittance = self.s * c.sym_value
        else:
            admittance = self.s * c.value
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # inductor element
    def _add_L(self, matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        if self.symbolic:
            admittance = 1 / (self.s * c.sym_value)
        else:
            admittance = 1/ self.s * c.value
        matrix[N1, N1] += admittance
        matrix[N2, N2] += admittance
        matrix[N1, N2] -= admittance
        matrix[N2, N1] -= admittance

    # ideal voltage source element
    def _add_V(self, admittance_matrix, source_matrix, c):
        N1 = self.node_dict[c.node1]
        N2 = self.node_dict[c.node2]
        pos = self.node_count + c.position
        admittance_matrix[pos, N1] += 1
        admittance_matrix[pos, N2] -= 1
        admittance_matrix[N1, pos] += 1
        admittance_matrix[N2, pos] -= 1
        source_matrix[pos, 0] += c.sym_value

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
        admittance_matrix[pos, N3] -= 1
        admittance_matrix[pos, N4] += 1
        admittance_matrix[N1, pos] += 1
        admittance_matrix[N2, pos] -= 1

    def _add_VVT(self, admittance_matrix, c):
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        NC_p = self.node_dict[c.node3]
        NC_n = self.node_dict[c.node4]
        pos = self.node_count + c.position
        if self.symbolic:
            gain = c.sym_value
        else:
            gain = c.value
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
        if self.symbolic:
            gain = c.sym_value
        else:
            gain = c.value
        admittance_matrix[N_p, NC_p] -= gain
        admittance_matrix[N_p, NC_n] += gain
        admittance_matrix[N_n, NC_p] -= gain
        admittance_matrix[N_n, NC_n] += gain

    def _add_CCT(self, admittance_matrix, source_matrix, c):
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        c_v = self.components[c.control_voltage]
        NC_p = self.node_dict[c_v.node1]
        NC_n = self.node_dict[c_v.node2]
        #print(NC_p)
        #print(NC_n)
        pos = self.node_count + c.position
        if self.symbolic:
            gain = c.sym_value
        else:
            gain = c.value

        admittance_matrix[pos, NC_p] += 1
        admittance_matrix[pos, NC_n] -= 1
        admittance_matrix[NC_p, pos] += 1
        admittance_matrix[NC_n, pos] -= 1
        admittance_matrix[N_p, pos] += gain
        admittance_matrix[N_n, pos] -= gain

    def _add_CVT(self, admittance_matrix, c):
        N_p = self.node_dict[c.node1]
        N_n = self.node_dict[c.node2]
        NC_p = self.node_dict[c.node3]
        NC_n = self.node_dict[c.node4]
        pos1 = self.node_count + c.position - 1
        pos2 = self.node_count + c.position
        if self.symbolic:
            gain = c.sym_value
        else:
            gain = c.value


        admittance_matrix[pos1, NC_p] += 1
        admittance_matrix[pos1, NC_n] -= 1
        admittance_matrix[pos2, N_p] += 1
        admittance_matrix[pos2, N_n] -= 1
        admittance_matrix[NC_p, pos1] += 1
        admittance_matrix[NC_n, pos1] -= 1
        admittance_matrix[N_p, pos2] += 1
        admittance_matrix[N_n, pos2] -= 1
        admittance_matrix[pos2, pos1] -= gain


def latex_print(data):
    print("{}".format(sympy.latex(data)))


if __name__ == "__main__":
    dirname = os.path.dirname(__file__)
    netlist = dirname + "/oamp.txt"
    circuit = AnalyseCircuit(netlist)
    latex_print(circuit.eqn_matrix)
    latex_print(circuit.solved_dict)
    latex_print(circuit.transfer_function("1", "2"))



