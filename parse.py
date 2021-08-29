from component import Component
import sympy

def source_value(words):
    """
    TODO
    """
    return 10, True

def value_enum(words):
    nums = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
    units = {"f": 10 ** (-15), "p": 10 ** (-12), "n": 10 ** (-9), "u": 10 ** (-6), "m": 10 ** (-3),
             "k": 10 ** 3, "meg": 10 ** 6, "g": 10 ** 9, "t": 10 ** 12}
    symbolic = False
    # recalculate units
    if words[3][0] not in nums:
        symbolic = True
        value = words[3]
    elif words[3][-3:-1] in units:
        value = float(words[3][:-3]) * units["meg"]
    elif words[3][-1] in units:
        value = float(words[3][:-1]) * units[words[3][-1]]
    else:
        value = float(words[3])
    return value, symbolic

def parse(netlist):
    data = {}
    with open(netlist) as f:
        parsed_netlist = f.readlines()
    parsed_netlist = [x.strip() for x in parsed_netlist]
    components = {}
    count = 0
    nodes = []
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
            node1 = int(words[1])
            node2 = int(words[2])
            symbolic = False
            if node1 not in nodes:
                nodes.append(node1)
            if node2 not in nodes:
                nodes.append(node2)

            # identify type of component and assign symbolic value

            if name[0] in ["i", "I"]:
                type = "i"
                sym_value = sympy.Symbol(name, real=True)
                c = Component(name, type, node1, node2, sym_value=sym_value)

            elif name[0] in ["v", "V", "u", "U"]:
                type = "v"
                sym_value = sympy.Symbol(name, real=True)
                c = Component(name, type, node1, node2, sym_value=sym_value, position=matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["r", "R"]:
                type = "r"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.Symbol(name, real=True)
                c = Component(name, type, node1, node2, sym_value=sym_value, value=value)

            elif name[0] in ["c", "C"]:
                type = "c"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.Symbol(name, real=True)
                try:
                    init_cond = words[4]
                    c = Component(name, type, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
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
                    c = Component(name, type, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
                except IndexError:
                    #print("No initial condition set for {}".format(name))
                    c = Component(name, type, node1, node2, sym_value=sym_value, value=value)

            elif name[0] in ["a", "A"]:
                type = "a"
                sym_value = sympy.Symbol(name, real=True)
                node3 = int(words[3])
                node4 = int(words[4])
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = Component(name, type, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["e", "E"]:
                type = "e"
                sym_value = sympy.Symbol(name, real=True)
                node3 = int(words[3])
                node4 = int(words[4])
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = Component(name, type, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["g", "G"]:
                type = "g"
                sym_value = sympy.Symbol(name, real=True)
                node3 = int(words[3])
                node4 = int(words[4])
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = Component(name, type, node1, node2, node3, node4, sym_value=sym_value)

            elif name[0] in ["f", "F"]:
                type = "f"
                sym_value = sympy.Symbol(name, real=True)
                v_c = words[3]
                c = Component(name, type, node1, node2, control_voltage=v_c, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1

            elif name[0] in ["h", "H"]:
                type = "h"
                sym_value = sympy.Symbol(name, real=True)
                node3 = int(words[3])
                node4 = int(words[4])
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = Component(name, type, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 2

            components[c.name] = c
        count += 1
    node_dict = {}
    i = 0
    for node in nodes:
        node_dict[node] = i
        i += 1

    data["node_dict"] = node_dict
    data["node_count"] = i
    data["matrix_size"] = i + matrix_expansion_coef
    data["components"] = components

    return data


if __name__ == '__main__':
    circuit = parse("oamp.txt")
    print(circuit)
    for i in circuit:
        try:
            if i.type == "a":
                print("{}; type:{}; nodes: {} {} {} {}; symbol:{}".format(i.name, i.type, i.node1, i.node2, i.node3,
                                                                          i.node4, i.sym_value))
            else:
                print("{}; type:{}; nodes: {} {}; symbol:{}".format(i.name, i.type, i.node1, i.node2, i.sym_value))
        except:
            pass
