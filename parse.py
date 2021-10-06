from component import *
import sympy
from laplace import *

t = sympy.Symbol("t", real=True, positive=True)
s = sympy.Symbol("s", real=True, positive=True)
f = sympy.symbols("f", real=True, positive=True)
j = sympy.symbols("j")
NUMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
UNITS = {"f": 10 ** (-15), "p": 10 ** (-12), "n": 10 ** (-9), "u": 10 ** (-6), "m": 10 ** (-3),
         "k": 10 ** 3, "meg": 10 ** 6, "g": 10 ** 9, "t": 10 ** 12}


def convert_units(val):
    symbolic = False
    if val[0] not in NUMS:
        symbolic = True
        ret = val
    elif val[-3:-1] in UNITS:
        ret = float(val[:-3]) * UNITS["meg"]
    elif val[-1] in UNITS:
        ret = float(val[:-1]) * UNITS[val[-1]]
    else:
        ret = float(val)
    return ret, symbolic


def dc_value(words):
    try:
        if words[3] == "dc":
            dc_value, symbolic = convert_units(words[4])
        else:
            symbolic = True
            dc_value = sympy.Symbol(words[0], real=True)
    except IndexError:
        dc_value = 0
        #print("WARNING: ({}) has no dc value".format(words))
    return dc_value


def ac_value(words):
    try:
        if words[5] == "ac":
            ac_value, symbolic = convert_units(words[6])
            try:
                offset, _ = convert_units(words[7])
            except IndexError:
                offset = 0
        else:
            symbolic = True
            offset = 0
            ac_value = sympy.Symbol(words[0], real=True)
        ac_value = ac_value * sympy.exp(j*offset)
        #print(ac_value)
    except IndexError:
        symbolic = True
        offset = 0
        ac_value = sympy.Symbol(words[0], real=True)
        #print("WARNING: ({}) has no ac value".format(words))
    return ac_value


def value_enum(words, source=False):
    symbolic = False
    if source:
        dc = dc_value(words)
        ac = ac_value(words)
        # TODO add tran value
        return [dc, ac], symbolic
    else:
        # recalculate units
        value, symbolic = convert_units(words[3])
        return value, symbolic


def parse(netlist):
    data = {}
    parsed_netlist = [x.strip() for x in netlist]
    components = {}
    count = 0
    nodes = []
    independent_sources = []
    basic_components = []
    controlled_sources = []
    operational_amplifiers = []
    i_list = []
    i_list = []
    matrix_expansion_coef = 0
    for line in parsed_netlist:
        words = line.split()
        if line == ".end":
            break
        elif count == 0:
            pass
        else:
            # count number of nodes
            name = words[0]
            node1 = words[1]
            node2 = words[2]
            symbolic = False
            if node1 not in nodes:
                nodes.append(node1)
            if node2 not in nodes:
                nodes.append(node2)

            # identify type of component and assign symbolic value

            if name[0] in ["i", "I"]:
                type = "i"
                sym_value = sympy.Symbol(name, real=True)
                value, symbolic = value_enum(words, source=True)
                c = CurrentSource(name, type, node1, node2, sym_value=sym_value, dc_value=value[0], ac_value=value[1])
                independent_sources.append(c)

            elif name[0] in ["v", "V", "u", "U"]:
                type = "v"
                sym_value = sympy.Symbol(name, real=True)
                value, symbolic = value_enum(words, source=True)
                c = VoltageSource(name, type, node1, node2, sym_value=sym_value, position=matrix_expansion_coef,
                                  dc_value=value[0], ac_value=value[1])
                matrix_expansion_coef += 1
                independent_sources.append(c)

            elif name[0] in ["r", "R"]:
                type = "r"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.Symbol(name, real=True)
                c = Resistor(name, type, node1, node2, sym_value=sym_value, value=value)
                basic_components.append(c)

            elif name[0] in ["c", "C"]:
                type = "c"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.Symbol(name, real=True)
                try:
                    init_cond = words[4]
                    c = Capacitor(name, type, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
                except IndexError:
                    #print("No initial condition set for {}".format(name))
                    c = Component(name, type, node1, node2, sym_value=sym_value, value=value)

            elif name[0] in ["l", "L"]:
                type = "l"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.Symbol(name, real=True)
                try:
                    init_cond = words[4]
                    c = Inductor(name, type, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
                except IndexError:
                    #print("No initial condition set for {}".format(name))
                    c = Inductor(name, type, node1, node2, sym_value=sym_value, value=value)

            elif name[0] in ["a", "A"]:
                type = "a"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = OperationalAmplifier(name, type, node1, node2, node3, node4, sym_value,
                                         matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["e", "E"]:  # VVT
                type = "e"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = VoltageControlledSource(name, type, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["g", "G"]:  # VCT
                type = "g"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = VoltageControlledSource(name, type, node1, node2, node3, node4, sym_value=sym_value)

            elif name[0] in ["f", "F"]:  # CCT
                type = "f"
                sym_value = sympy.Symbol(name, real=True)
                v_c = words[3]
                c = CurrentControlledSource(name, type, node1, node2, control_voltage=v_c, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 2

            elif name[0] in ["h", "H"]:  # CVT
                type = "h"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = CurrentControlledSource(name, type, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 3

            components[c.name] = c
        count += 1
    node_dict = {}
    i = 0
    for node in nodes:
        node_dict[node] = i
        i += 1
    #print(i+matrix_expansion_coef)

    data["node_dict"] = node_dict
    data["node_count"] = i
    data["matrix_size"] = i + matrix_expansion_coef
    data["components"] = components

    return data


if __name__ == '__main__':
    circuit = parse("oamp.txt")
    #print(circuit)
    for i in circuit:
        try:
            if i.type == "a":
                print("{}; type:{}; nodes: {} {} {} {}; symbol:{}".format(i.name, i.type, i.node1, i.node2, i.node3,
                                                                          i.node4, i.sym_value))
            else:
                print("{}; type:{}; nodes: {} {}; symbol:{}".format(i.name, i.type, i.node1, i.node2, i.sym_value))
        except:
            pass
