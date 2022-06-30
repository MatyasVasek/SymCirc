from symcirc.component import *

t = sympy.Symbol("t", real=True, positive=True)
s = sympy.Symbol("s", real=True, positive=True)
f = sympy.symbols("f", real=True, positive=True)
j = sympy.symbols("j")
NUMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
UNITS = {"f": sympy.Rational(1, 1000000000000), "p": sympy.Rational(1, 1000000000000),
         "n": sympy.Rational(1, 1000000000), "u": sympy.Rational(1, 1000000), "m": sympy.Rational(1, 1000),
         "k": sympy.Rational(1000, 1), "meg": sympy.Rational(1000000, 1), "g": sympy.Rational(1000000000, 1),
         "t": sympy.Rational(1000000000000, 1)}
OPERATORS = ["+", "-", "*", "/"]
RESERVED = ["sin"]


def check_if_symbolic(val):
    symbolic = False
    for c in val:
        if c in NUMS:
            pass
        elif c in UNITS:
            pass
        elif c in OPERATORS:
            pass
        else:
            symbolic = True
    return symbolic


def convert_units(val, forced_numeric=False):
    if forced_numeric:
        symbolic = False
    else:
        symbolic = check_if_symbolic(val)

    if symbolic:
        symbolic = True
        ret = sympy.parse_expr(val)

    elif val[-3:-1] in UNITS:
        ret = sympy.Rational(sympy.parse_expr(val[:-3])) * UNITS["meg"]
    elif val[-1] in UNITS:
        ret = sympy.Rational(sympy.parse_expr(val[:-1])) * UNITS[val[-1]]
    else:
        ret = sympy.Rational(sympy.parse_expr(val))
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
                if words[7] not in RESERVED:
                    offset, _ = convert_units(words[7])
                else:
                    offset = 0
            except IndexError:
                offset = 0
        else:
            symbolic = True
            offset = 0
            ac_value = sympy.Symbol(words[0], real=True)
        #print(ac_value)
        #print(offset)
        ac_value = ac_value * sympy.exp(j*offset)
        #print(ac_value)
    except IndexError:
        symbolic = True
        offset = 0
        ac_value = sympy.Symbol(words[0], real=True)
        #print("WARNING: ({}) has no ac value".format(words))
    return ac_value

def tran_value(words):
    index = 1
    offset = 0
    amp = 1
    freq = 1
    delay = 0
    damping = 1
    #tran = sympy.Symbol("N/A")

    for word in words:
        if word == "sin":
            break
        else:
            index += 1
    try:
        offset = sympy.Rational(sympy.parse_expr(words[index]))
        amp = sympy.Rational(sympy.parse_expr(words[index+1]))
        freq, _ = convert_units(words[index+2], forced_numeric=True)
        freq = sympy.Rational(freq)
        delay = sympy.Rational(sympy.parse_expr(words[index+3]))
        damping = sympy.Rational(sympy.parse_expr(words[index+4]))
    except IndexError:
        pass
    #offset, amp, freq, delay, damping = sympy.symbols("off, A, f, del, damp", real=True)
    #pi = 3.14159
    pi = sympy.pi
    tran = (offset/s)+amp*sympy.exp(-s*delay)*2*pi*freq/((s+damping)**2+(2*pi*freq)**2)
    #tran = sympy.apart(tran)
    #print(tran)

    return tran



def value_enum(words, source=False):
    symbolic = False
    if source:
        dc = dc_value(words)
        ac = ac_value(words)
        tran = tran_value(words)
        return [dc, ac, tran], symbolic
    else:
        # recalculate units
        value, symbolic = convert_units(words[3])
        return value, symbolic

def nodes_per_element(type):
    if type in ["r", "R", "l", "L", "c", "C", "v", "V", "i", "I", "f", "F", "h", "H"]:
        return 2
    elif type in ["a", "A", "e", "E", "g", "G"]:
        return 4
    elif type in ["k", "K"]:
        return 0

def unpack_subcircuits(parsed_netlist):
    """
    Identifies all subcircuits, unpacks them and returns a unpacked netlist
    Elements inside a subcircuit inherit it's name in the following format: ElementName_SubcircuitName
    :param list parsed_netlist: each cell contains one netlist line in string format
    :return: list unpacked_netlist: each cell contains one netlist line in string format
    """
    unpacked_netlist = [parsed_netlist[0]]
    subckt_instances = []
    subckt_models = {}
    in_subckt = False
    model = ""
    for line in parsed_netlist[1:]:
        words = line.split()
        if line == ".end":
            break
        elif words in [[], "\n", " "]:
            pass
        elif words[0][0] == "*":  # check if line is commentary
            pass

        elif words[0][0] in ["x", "X"]:
            nodes_and_model_id = []
            param_dict = {}
            loading_nodes = True
            for w in words[1:]:
                if w == "PARAMS:":
                    loading_nodes = False
                elif loading_nodes:
                    nodes_and_model_id.append(w)
                else:
                    key, val = w.split("=")
                    param_dict[key] = val

            node_list = nodes_and_model_id[0:-1]
            model_id = nodes_and_model_id[-1]

            c = Subcircuit(words[0], model_id, node_list, param_dict)
            subckt_instances.append(c)

        elif words[0] in [".subckt", ".SUBCKT"]:  # subcircuit model
            in_subckt = True
            loading_nodes = True
            node_list = []
            param_dict = {}
            model_id = words[1]
            for w in words[2:]:
                if w == "PARAMS:":
                    loading_nodes = False
                elif loading_nodes:
                    node_list.append(w)
                else:
                    key, val = w.split("=")
                    param_dict[key] = val
                    model = SubcktModel(model_id, node_list, param_dict)

        elif words[0] in [".ends", ".ENDS"]:
            in_subckt = False

        elif in_subckt:
            line = ""
            for instance in subckt_instances:
                if instance.model_id == model.model_id:
                    node_count = nodes_per_element(words[0][0])
                    node_dict = {}
                    n = 0
                    for node in instance.node_list:
                        node_dict[model.node_list[n]] = node
                        n+=1
                    name = words[0]+"_"+instance.name
                    line += name

                    i = 0
                    for w in words[1:]:
                        if i < node_count:
                            if w in model.node_list:
                                line = line + " " + node_dict[w]
                            else:
                                line = line + " " + w + "_" + instance.name
                        elif w[0] == "{":
                            param = w[1:-1]
                            try:
                                line = line + " " + instance.param_dict[param]
                            except IndexError:
                                try:
                                    line = line + " " + model.param_dict[param]
                                except IndexError:
                                    exit("Subcircuit syntax error")
                        elif w[0] in NUMS:
                            line = line + " " + w
                        elif w in ["dc", "ac", "tran", "sin"]:
                            line = line + " " + w
                        else:
                            line = line + " " + w+"_"+instance.name
                        i+=1


            unpacked_netlist.append(line)
        else:
            unpacked_netlist.append(line)
    return unpacked_netlist


def parse(netlist):
    """
    Translates
    :param str netlist: netlist in a string format
    :return: list data: data contains four items: \n
    * :node_dict: dictionary of circuit nodes
    * :code_count: amount of nodes
    * :matrix_size: matrix size needed to fit all components in
    * :components: list of components
    \n
    Input example: \n
    Circuit AC6
    V1 a 0 dc 0 ac 1 0 sin 0 1 14k 0 0
    R1 a b R1
    L b 0 L1
    R2 b 0 1k
    """
    data = {}
    parsed_netlist = netlist.splitlines() #[x.strip() for x in netlist]
    components = {}
    count = 0
    nodes = []
    independent_sources = []
    basic_components = []
    controlled_sources = []
    operational_amplifiers = []
    add_short = []
    matrix_expansion_coef = 0
    parsed_netlist = unpack_subcircuits(parsed_netlist)

    for line in parsed_netlist:
        words = line.split()
        if line == ".end":
            break
        elif count == 0:  # first line is for title
            pass
        elif words in [[], "\n", " "]:
            pass
        elif words[0][0] == "*":  # check if line is commentary
            pass
        else:


            # count number of nodes
            name = words[0]
            if name[0] not in ['k', 'K']:
                node1 = words[1]
                node2 = words[2]
                symbolic = False
                if node1 not in nodes:
                    nodes.append(node1)
                if node2 not in nodes:
                    nodes.append(node2)

            # identify variant of component and assign symbolic value

            if name[0] in ["i", "I"]:
                variant = "i"
                sym_value = sympy.parse_expr(name)  # sympy.Symbol(name, real=True)
                value, symbolic = value_enum(words, source=True)
                c = CurrentSource(name, variant, node1, node2, sym_value=sym_value, dc_value=value[0], ac_value=value[1],
                                  tran_value=value[2])
                independent_sources.append(c)

            elif name[0] in ["v", "V", "u", "U"]:
                variant = "v"
                sym_value = sympy.parse_expr(name)  # sympy.Symbol(name, real=True)
                value, symbolic = value_enum(words, source=True)
                c = VoltageSource(name, variant, node1, node2, sym_value=sym_value, position=matrix_expansion_coef,
                                  dc_value=value[0], ac_value=value[1], tran_value=value[2])
                matrix_expansion_coef += 1
                independent_sources.append(c)

            elif name[0] in ["r", "R"]:
                variant = "r"
                value, symbolic = value_enum(words)
                if symbolic:
                    #print(type(value))
                    sym_value = value  # sympy.Symbol(value, real=True)

                else:
                    sym_value = sympy.parse_expr(name)  # sympy.Symbol(name, real=True)
                    #print(sym_value.atoms(sympy.Symbol))
                c = Resistor(name, variant, node1, node2, sym_value=sym_value, value=value)
                basic_components.append(c)

            elif name[0] in ["c", "C"]:
                variant = "c"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = value  # sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.parse_expr(name)  # sympy.Symbol(name, real=True)
                try:
                    init_cond = words[4]
                    c = Capacitor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
                except IndexError:
                    #print("No initial condition set for {}".format(name))
                    c = Component(name, variant, node1, node2, sym_value=sym_value, value=value)
                basic_components.append(c)

            elif name[0] in ["l", "L"]:
                variant = "l"
                value, symbolic = value_enum(words)
                if symbolic:
                    sym_value = value  # sympy.Symbol(value, real=True)
                else:
                    sym_value = sympy.parse_expr(name)  # sympy.Symbol(name, real=True)
                try:
                    init_cond = words[4]
                    c = Inductor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
                except IndexError:
                    #print("No initial condition set for {}".format(name))
                    c = Inductor(name, variant, node1, node2, sym_value=sym_value, value=value)
                basic_components.append(c)

            elif name[0] in ["a", "A"]:
                variant = "a"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = OperationalAmplifier(name, variant, node1, node2, node3, node4, sym_value,
                                         matrix_expansion_coef)
                matrix_expansion_coef += 1
                operational_amplifiers.append(c)

            elif name[0] in ["e", "E"]:  # VVT
                variant = "e"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = VoltageControlledSource(name, variant, node1, node2, node3, node4, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1
                controlled_sources.append(c)

            elif name[0] in ["g", "G"]:  # VCT
                variant = "g"
                sym_value = sympy.Symbol(name, real=True)
                node3 = words[3]
                node4 = words[4]
                if node3 not in nodes:
                    nodes.append(node3)
                if node4 not in nodes:
                    nodes.append(node4)
                c = VoltageControlledSource(name, variant, node1, node2, node3, node4, sym_value=sym_value)
                controlled_sources.append(c)

            elif name[0] in ["f", "F"]:  # CCT
                variant = "f"
                sym_value = sympy.Symbol(name, real=True)
                v_c = words[3]
                c = CurrentControlledSource(name, variant, node1, node2, control_voltage=v_c, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1
                add_short.append(v_c)
                controlled_sources.append(c)

            elif name[0] in ["h", "H"]:  # CVT
                variant = "h"
                sym_value = sympy.Symbol(name, real=True)
                v_c = words[3]
                c = CurrentControlledSource(name, variant, node1, node2, control_voltage=v_c, sym_value=sym_value,
                              position=matrix_expansion_coef)
                matrix_expansion_coef += 1
                add_short.append(v_c)
                controlled_sources.append(c)

            elif name[0] in ["k", "K"]:  # coupled inductors
                variant = "k"
                value, symbolic = value_enum(words)
                sym_value = sympy.Symbol(name, real=True)
                L1 = words[1]
                L2 = words[2]
                c = Coupling(name, variant, L1, L2, sym_value, value)
                controlled_sources.append(c)
            components[c.name] = c
        count += 1
    shorts = []
    for key in components:
        if key in add_short:
            shorts.append(key)
    for key in shorts:
        c = components[key]
        new_node = "*short{}".format(c.name)
        nodes.append(new_node)
        c.shorted_node = c.node2
        #short = Short("S{}".format(c.name), "s", new_node, c.node2)
        #components[short.name] = short
        c.node2 = new_node


    node_dict = {}
    i = 0
    for node in nodes:
        if node != "0":
            node_dict[node] = i
            i += 1

    data["node_dict"] = node_dict
    data["node_count"] = i
    data["matrix_size"] = i + matrix_expansion_coef
    data["components"] = components

    return data
