from symcirc.component import *
from symcirc import utils
from symcirc import laplace
from symcirc.utils import t, s, f, j
import sys


NUMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
UNITS = {"f": sympy.Rational(1, 1000000000000000), "p": sympy.Rational(1, 1000000000000),
         "n": sympy.Rational(1, 1000000000), "u": sympy.Rational(1, 1000000), "m": sympy.Rational(1, 1000),
         "k": 1000, "meg": 1000000, "G": 1000000000,
         "T": 1000000000000}
OPERATORS = ["+", "-", "*", "/", "."]
RESERVED = ["sin"]

NETLIST_KEYCHARS = ["R", "r", "C", "c", "L", "l", "V", "v", "U", "u", "I", "i", "A", "a", "F", "f", "H", "h", "G", "g",
                    "E", "e", "K", "k", "S", "s", "X", "x", ".", "*"]

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
        local = {}
        local.update(UNITS)
        local.update(sympy.abc._clash)
        val = val.replace("{", "").replace("}", "")
        ret = sympy.parse_expr(val, local_dict=local)
    elif val[-3:-1] in UNITS:
        ret = sympy.Rational(sympy.parse_expr(val[:-3], local_dict=sympy.abc._clash) * UNITS["meg"])
    elif val[-1] in ["k", "g", "t"]:
        ret = sympy.Rational(sympy.parse_expr(val[:-1], local_dict=sympy.abc._clash) * UNITS[val[-1]])
    elif val[-1] in UNITS:
        ret = sympy.Rational(sympy.parse_expr(val[:-1], local_dict=sympy.abc._clash) * UNITS[val[-1]])
    else:
        #ret = sympy.parse_expr(val)
        try:
            ret = sympy.Rational(sympy.parse_expr(val, local_dict=sympy.abc._clash))
        except TypeError:
            ret = sympy.parse_expr(val, local_dict=sympy.abc._clash)
    return ret, symbolic


def dc_value(words):
    try:
        if words[3] in ["dc", "DC"]:
            dc_value, symbolic = convert_units(words[4])
        elif len(words) == 4:
            dc_value, symbolic = convert_units(words[3])
        else:
            symbolic = True
            dc_value = sympy.Symbol(words[0], real=True)
    except IndexError:
        dc_value = 0
    return dc_value


def ac_value(words):
    try:
        if words[5] in ["ac", "AC"]:
            ac_value, symbolic = convert_units(words[6])
            try:
                if words[7] not in RESERVED:
                    offset, _ = convert_units(words[7])
                else:
                    offset = 0
            except IndexError:
                offset = 0
        else:
            ac_value = 0
    except IndexError:
        ac_value = 0
    return ac_value

def tran_value(words, dc):
    use_DC_val = True
    index = 1
    offset = 0
    amp = 1
    freq = 1
    delay = 0
    damping = 1
    #tran = sympy.Symbol("N/A")

    for word in words:
        if word in ["sin", "SIN"]:
            use_DC_val = False
            break
        else:
            index += 1
    if use_DC_val:
        tran = dc/s
    else:
        try:
            offset, _ = convert_units(words[index])
            amp, _ = convert_units(words[index+1])
            freq, _ = convert_units(words[index+2], forced_numeric=True)
            delay, _ = convert_units(words[index+3])
            damping, _ = convert_units(words[index+4])
        except IndexError:
            pass
        #offset, amp, freq, delay, damping = sympy.symbols("off, A, f, del, damp", real=True)
        #pi = 3.14159
        pi = sympy.pi
        tran = (offset/s)+amp*sympy.exp(-s*delay)*2*pi*freq/((s+damping)**2+(2*pi*freq)**2)
    return tran

def value_enum(words, source=False):
    symbolic = False
    if source:
        dc = dc_value(words)
        ac = ac_value(words)
        tran = tran_value(words, dc)
        return [dc, ac, tran], symbolic
    else:
        # recalculate units
        try:
            value, symbolic = convert_units(words[3])
        except IndexError:
            symbolic = True
            value = sympy.parse_expr(words[0], local_dict=sympy.abc._clash)
        return value, symbolic

def nodes_per_element(type):
    if type in ["r", "R", "l", "L", "c", "C", "v", "V", "i", "I", "f", "F", "h", "H", "s", "S"]:
        return 2
    elif type in ["a", "A", "e", "E", "g", "G"]:
        return 4
    elif type in ["k", "K"]:
        return 0

def parse_subcircuits(netlist):
    subckt_models = {}
    in_model = False
    model_id = ""
    current_model = None
    parsed_netlist = []

    for line in netlist[1:]:
        words = line.split()

        if words in [[], "\n", " "]:
            pass

        elif words[0][0] not in NETLIST_KEYCHARS:
            raise SyntaxError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif words[0] in [".subckt", ".SUBCKT"]:
            in_model = True
            model_id = words[1]
            loading_nodes = True
            node_list = []
            param_dict = {}
            model_id = words[1]
            for w in words[2:]:
                if w in ["PARAMS:", "params:"]:
                    loading_nodes = False
                elif loading_nodes:
                    node_list.append(w)
                else:
                    key, val = w.split("=")
                    param_dict[key] = val
            current_model = SubcktModel(model_id, node_list, param_dict)

        elif words[0][0] == ".":
            if (words[0] in [".ends", ".ENDS"]) and (model_id in words[1]):
                in_model = False
                subckt_models[model_id] = current_model
            elif words[0] in [".end", ".END"]:
                break
            else:
                raise SyntaxError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif in_model:
            current_model.elements.append(line)

        else:
            parsed_netlist.append(line)

    final_netlist = unpack(parsed_netlist, subckt_models)

    return final_netlist


def unpack(parsed_netlist, subckt_models):
    """
    Identifies all subcircuits, unpacks them and returns a unpacked netlist
    Elements inside a subcircuit inherit it's name in the following format: ElementName_SubcircuitName
    :param list parsed_netlist: each cell contains one netlist line in string format
    :return: list unpacked_netlist: each cell contains one netlist line in string format
    """
    final_netlist = []
    for line in parsed_netlist:
        words = line.split()
        if words[0][0] in ["x", "X"]:
            loading_params = False
            nodes = []
            params = {}
            model = None
            for word in words[1:]:
                if word in ["params:", "PARAMS:"]:
                    if model:
                        loading_params = True
                    else:
                        raise SyntaxError(f"Model of element '{words[0]}' is not present in the netlist.")

                elif not loading_params:
                    try:
                        model = subckt_models[word]
                    except KeyError:
                        nodes.append(word)

                else:
                    tmp = word.split("=")
                    if len(tmp) != 2:
                        raise SyntaxError(f"Parameter '{word}' is not in the correct format.")
                    params[tmp[0]], _ = convert_units(tmp[1])
            node_index = 0
            node_dict = {}
            for node in nodes:
                node_dict[model.node_list[node_index]] = node
                node_index += 1

            for elem in model.elements:
                split_elem = elem.split(" ")
                if split_elem[0][0] not in NETLIST_KEYCHARS:
                    raise SyntaxError(
                        f"Keyword/Element '{split_elem[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")
                elif split_elem[0] in ["k", "K"]:
                    split_elem[1] = f"{split_elem[1]}_({words[0]})"
                    split_elem[2] = f"{split_elem[2]}_({words[0]})"
                    if split_elem[3][0] == "{":
                        try:
                            split_elem[3] = str(params[split_elem[3][1:-1]])
                        except KeyError:
                            split_elem[3] = str(model.param_dict[split_elem[3][1:-1]])
                else:
                    node_count = nodes_per_element(split_elem[0][0])
                    for i in range(1, node_count+1): # shift by one to avoid classifying name as node
                        if split_elem[i] in ["0", 0]:
                            pass
                        else:
                            try:
                                split_elem[i] = node_dict[split_elem[i]]
                            except KeyError:
                                split_elem[i] = f"{split_elem[i]}_({words[0]})"
                    index = 0
                    for e in split_elem:
                        if e[0] == "{":
                            local = sympy.abc._clash
                            try:
                                local.update(params)
                                split_elem[index] = str(sympy.parse_expr(e[1:-1], local_dict=local))
                                #split_elem[index] = params[e[1:-1]]
                            except KeyError:
                                local.update(model.param_dict)
                                split_elem[index] = str(sympy.parse_expr(e[1:-1], local_dict=local))
                                #split_elem[index] = model.param_dict[e[1:-1]]

                        index += 1
                split_elem[0] = f"{split_elem[0]}_({words[0]})"
                final_netlist.append(" ".join(split_elem))
        else:
            final_netlist.append(line)

    return final_netlist


def parse(netlist, tran=False):
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
    c = None
    nodes = []
    independent_sources = []
    basic_components = []
    controlled_sources = []
    operational_amplifiers = []
    couplings = []
    SCSI_components = []
    add_short = []
    matrix_expansion_coef = 0
    parsed_netlist = parse_subcircuits(parsed_netlist)

    for line in parsed_netlist:
        words = line.split()
        name = words[0]
        if name[0] == "*":
            continue
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
            sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash)  # sympy.Symbol(name, real=True)
            value, symbolic = value_enum(words, source=True)
            c = CurrentSource(name, variant, node1, node2, sym_value=sym_value, dc_value=value[0], ac_value=value[1],
                              tran_value=value[2])
            independent_sources.append(c)

        elif name[0] in ["v", "V", "u", "U"]:
            variant = "v"
            sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash)  # sympy.Symbol(name, real=True)
            value, symbolic = value_enum(words, source=True)
            c = VoltageSource(name, variant, node1, node2, sym_value=sym_value, position=matrix_expansion_coef,
                              dc_value=value[0], ac_value=value[1], tran_value=value[2])
            matrix_expansion_coef += 1
            independent_sources.append(c)

        elif name[0] in ["r", "R"]:
            variant = "r"
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)

            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash)  # sympy.Symbol(name, real=True)
            c = Resistor(name, variant, node1, node2, sym_value=sym_value, value=value)
            basic_components.append(c)

        elif name[0] in ["c", "C"]:
            variant = "c"
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash)  # sympy.Symbol(name, real=True)
            try:
                init_cond, _ = convert_units(words[4][3:])
                c = Capacitor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            except IndexError:
                init_cond = 0
                c = Capacitor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            basic_components.append(c)

        elif name[0] in ["l", "L"]:
            variant = "l"
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash)  # sympy.Symbol(name, real=True)
            try:
                init_cond, _ = convert_units(words[4][3:])
                c = Inductor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            except IndexError:
                init_cond = 0
                c = Inductor(name, variant, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
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
            value, symbolic = None, None #convert_units(words[5])
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
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
            try:
                value, symbolic = convert_units(words[5])
            except IndexError:
                symbolic = True
                value = sympy.parse_expr(name, local_dict=sympy.abc._clash)

            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)

            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=value, sym_value=sym_value,
                          position=matrix_expansion_coef)
            matrix_expansion_coef += 1
            controlled_sources.append(c)

        elif name[0] in ["g", "G"]:  # VCT
            variant = "g"
            node3 = words[3]
            node4 = words[4]
            if node3 not in nodes:
                nodes.append(node3)
            if node4 not in nodes:
                nodes.append(node4)
            try:
                value, symbolic = convert_units(words[5])
            except IndexError:
                symbolic = True
                value = sympy.parse_expr(name, local_dict=sympy.abc._clash)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=value, sym_value=sym_value)
            controlled_sources.append(c)

        elif name[0] in ["f", "F"]:  # CCT
            variant = "f"
            sym_value = sympy.Symbol(name, real=True)
            v_c = words[3]
            try:
                value, symbolic = convert_units(words[4])
            except IndexError:
                symbolic = True
                value = sympy.parse_expr(name, local_dict=sympy.abc._clash)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = CurrentControlledSource(name, variant, node1, node2, control_voltage=v_c, value=value, sym_value=sym_value,
                          position=matrix_expansion_coef)
            matrix_expansion_coef += 1
            add_short.append(v_c)
            controlled_sources.append(c)

        elif name[0] in ["h", "H"]:  # CVT
            variant = "h"
            sym_value = sympy.Symbol(name, real=True)
            v_c = words[3]
            try:
                value, symbolic = convert_units(words[4])
            except IndexError:
                symbolic = True
                value = sympy.parse_expr(name, local_dict=sympy.abc._clash)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = CurrentControlledSource(name, variant, node1, node2, control_voltage=v_c, value=value, sym_value=sym_value,
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
            couplings.append(c)

        elif name[0] in ["s", "S"]:  # periodic switch used for SC/SI analysis
            variant = "s"
            phase = words[3]
            c = PeriodicSwitch(name, variant, node1, node2, phase=phase)
            SCSI_components.append(c)

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
        c.node2 = new_node

    node_dict = {}
    i = 0
    grounded = False
    for node in nodes:
        if node != "0":
            node_dict[node] = i
            i += 1
        if node == "0":
            grounded = True
    if not grounded:
        exit("Circuit not grounded")

    for couple in couplings:
        L1 = components[couple.L1]
        L2 = components[couple.L2]
        L1.coupling = couple
        L2.coupling = couple

    data["node_dict"] = node_dict
    data["node_count"] = i
    data["components"] = components
    data["basic_components"] = basic_components
    data["independent_sources"] = independent_sources
    data["controlled_sources"] = controlled_sources
    data["oamps"] = operational_amplifiers
    data["couplings"] = couplings
    data["SCSI_components"] = SCSI_components

    return data
