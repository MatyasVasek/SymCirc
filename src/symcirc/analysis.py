import sympy
import sys
from symcirc import parse, laplace, utils


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
        data = parse.parse(netlist)
        self.components = data["components"]   # {<name> : <Component>} (see component.py)
        self.node_dict = data["node_dict"]  # {<node_name>: <index in equation matrix>, ...}
        self.matrix_size = data["matrix_size"]
        self.node_count = data["node_count"]
        self.c_count = self.count_components()
        self.node_voltage_symbols = self._node_voltage_symbols()
        self.eqn_matrix, self.solved_dict, self.symbols = self._analyse()  # solved_dict: {sympy.symbols(<vaviable_name>): <value>}

    def simp(self, expr):
        sympy.collect(expr, self.s)
        sympy.expand(expr)
        return expr

    def component_values(self, name):
        ret = {}
        if name == "all":
            ret = self.all_component_values()
            return ret
        else:
            v = "v({})".format(name)
            i = "i({})".format(name)
            ret[v] = self.solved_dict[sympy.symbols(v)]
            ret[i] = self.solved_dict[sympy.symbols(i)]
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
            if key in self.node_voltage_symbols and key.name[2] != "*":
                ret[key] = dic[key]
        return ret

    def transfer_function(self, node1, node2):
        v1 = sympy.Symbol("v({})".format(node1))
        v2 = sympy.Symbol("v({})".format(node2))
        voltage1 = self.solved_dict[v1].simplify()
        voltage2 = self.solved_dict[v2].simplify()
        tf = (voltage2/voltage1)
        return tf

    def count_components(self):
        count = 0
        for c in self.components:
            if self.components[c].type in ["a", "e", "g", "f", "h"]:
                count += 2
            else:
                count += 1
        return count

    def _analyse(self):
        if self.analysis_type == "DC":
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
            if self.is_symbolic:
                for sym in symbols:
                    try:
                        solved_dict[sym] = sympy.limit(solved_dict[sym], self.s, 0)
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
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
            f = sympy.symbols("f", real=True, positive=True)
            j = sympy.symbols("j")
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
                        solved_dict[sym] = solved_dict[sym].subs(self.s, 2 * sympy.pi * f * j)
                    except KeyError:
                        pass
            else:
                for sym in symbols:
                    #print(sym)
                    try:
                        #print(solved_dict[sym])
                        #solved_dict[sym] = solved_dict[sym].simplify()
                        #print(solved_dict[sym])
                        solved_dict[sym] = solved_dict[sym].subs(self.s, 2 * sympy.pi * f * j)
                        for name in self.components:
                            c = self.components[name]
                            if c.type == "v":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                #print(c.ac_value)
                            else:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
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
                            if c.type == "v":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.ac_value)
                                #print(c.ac_value)
                            else:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                    except KeyError:
                        pass

        elif self.analysis_type == "tran":
            eqn_matrix, symbols = self._build_system_eqn()
            solved_dict = sympy.solve_linear_system(eqn_matrix, *symbols)
            if self.is_symbolic:
                #latex_print(solved_dict)
                for sym in symbols:
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type == "v":
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.tran_value)
                                # print(c.ac_value)
                    except KeyError:
                        pass
                    #solved_dict[sym] = sympy.apart(solved_dict[sym], self.s)
                    #print(solved_dict[sym])
                    #print("{}: {}".format(sym, solved_dict[sym]))
                    inv_l = laplace.iLT(solved_dict[sym])
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
                    #print(sym)
                    try:
                        for name in self.components:
                            c = self.components[name]
                            if c.type == "v":

                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.tran_value)
                                # print(c.ac_value)
                            else:
                                solved_dict[sym] = solved_dict[sym].subs(c.sym_value, c.value)
                    except KeyError:
                        pass

                    solved_dict[sym] = laplace.iLT(solved_dict[sym])

        return eqn_matrix, solved_dict, symbols

    def _node_voltage_symbols(self):
        voltage_symbol_list = []
        for node in self.node_dict:
            if node != "0":
                voltage_symbol_list.append(sympy.Symbol("v({})".format(node)))
        return voltage_symbol_list

    def _build_system_eqn(self):
        size = self.c_count
        M = sympy.SparseMatrix(sympy.zeros(2 * size + self.node_count))
        R = sympy.SparseMatrix(sympy.zeros(size*2 + self.node_count, 1))
        for i in range(size):
            M[i, i] = 1
        index = 0
        node_symbols = self.node_voltage_symbols
        voltage_symbols = []
        current_symbols = []
        for key in self.components:
            c = self.components[key]
            if c.type in ["r", "l", "c"]:
                self._add_basic(M, c, index)
                voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                current_symbols.append(sympy.Symbol("i({})".format(c.name)))
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
            if c.type == "s":
                self._add_short(M, c, index)
                voltage_symbols.append(sympy.Symbol("v({})".format(c.name)))
                current_symbols.append(sympy.Symbol("i({})".format(c.name)))

            index += 1
        equation_matrix = M.col_insert(self.c_count*2 + self.node_count, R)
        symbols = voltage_symbols + current_symbols + node_symbols
        return equation_matrix, symbols

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

        if c.type == "r":
            y_b = 1
            z_b = -val
        elif c.type == "l":
            y_b = 1
            z_b = -self.s * val
        elif c.type == "c":
            y_b = self.s * val
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

    def _add_short(self, matrix, c, index):
        N1 = c.node1
        N2 = c.node2
        matrix[self.c_count + index, index] = 1
        self._incidence_matrix_write(N1, N2, matrix, index)
        return matrix


if __name__ == "__main__":
    netlist = "netlists\DC_elem_11.txt"
    c = AnalyseCircuit(utils.load_file(netlist))
    EM, solved_dict, symbols = c._analyse()
    print(symbols)
    sympy.pprint(EM)
    print(solved_dict)
    sympy.pprint(solved_dict)


