from symcirc.component import *
from symcirc import utils
from symcirc import laplace
from symcirc.utils import *
import sys
from sympy.parsing.sympy_parser import standard_transformations, convert_xor

TRANSFORMS = (standard_transformations + (convert_xor,))

NUMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
UNITS = {"f": sympy.Rational(1, 1000000000000000), "F": sympy.Rational(1, 1000000000000000),
         "p": sympy.Rational(1, 1000000000000), "P": sympy.Rational(1, 1000000000000),
         "n": sympy.Rational(1, 1000000000), "N": sympy.Rational(1, 1000000000),
         "u": sympy.Rational(1, 1000000), "U": sympy.Rational(1, 1000000),
         "m": sympy.Rational(1, 1000), "M": sympy.Rational(1, 1000),
         "k": sympy.Rational(1000), "K": sympy.Rational(1000),
         "meg": sympy.Rational(1000000), "MEG": sympy.Rational(1000000), "Meg": sympy.Rational(1000000),
         "g": sympy.Rational(1000000000), "G": sympy.Rational(1000000000),
         "t": sympy.Rational(1000000000000), "T": sympy.Rational(1000000000000),
         "v": sympy.Rational(1), "V": sympy.Rational(1),
         "a": sympy.Rational(1), "A": sympy.Rational(1)}
OPERATORS = ["+", "-", "*", "/", "."]
RESERVED = ["sin"]

NETLIST_KEYCHARS = ["R", "r", "C", "c", "L", "l", "V", "v", "U", "u", "I", "i", "A", "a", "F", "f", "H", "h", "G", "g",
                    "E", "e", "K", "k", "S", "s", "X", "x", "Q", "q", "M", "m", "D", "d", "J", "j", ".", "*"]

def check_if_symbolic(val):
    return not val.is_number


def convert_units(val, forced_numeric=False, local_dict=None):
    ret = None
    symbolic = True
    if local_dict is None:
        local = {}
    else:
        local = local_dict
    local.update(sympy.abc._clash)
    local.update(global_dict)
    val = val.replace("{", "").replace("}", "")
    if len(val) > 3:
        if (val[-3:] in UNITS) and (val[-4].isnumeric()):
            ret = sympy.parse_expr(val[:-3], local_dict=local, transformations=TRANSFORMS) * UNITS[val[-3:]]
    if len(val) > 1:
        if ret is None:
            if (val[-1] in UNITS) and (val[-2].isnumeric()):
                ret = sympy.parse_expr(val[:-1], local_dict=local, transformations=TRANSFORMS) * UNITS[val[-1]]
    if ret is None:
        ret = sympy.parse_expr(val, local_dict=local, transformations=TRANSFORMS)
    if forced_numeric:
        symbolic = False
    else:
        symbolic = check_if_symbolic(ret)

    try:
        ret = sympy.Rational(ret)
    except:
        pass
    return ret, symbolic


def dc_value(words):
    symbolic = False
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

    if symbolic:
        dc_sym = dc_value
    else:
        dc_sym = sympy.Symbol(words[0], real=True)
    return dc_value, dc_sym


def ac_value(words):
    phase_shift = 0
    symbolic = False
    try:
        if words[5] in ["ac", "AC"]:
            ac_value, symbolic = convert_units(words[6])
            try:
                if words[7] not in RESERVED:
                    phase_shift, _ = convert_units(words[7])
                else:
                    phase_shift = 0
            except IndexError:
                phase_shift = 0
        else:
            ac_value = 0
    except IndexError:
        ac_value = 0

    if symbolic:
        ac_sym = ac_value
    else:
        ac_sym = sympy.Symbol(words[0], real=True)
    return ac_value, phase_shift, ac_sym

def tran_value(words, dc):
    use_DC_val = True
    index = 1
    offset = 0
    amp = 1
    freq = 1
    delay = 0
    damping = 1
    omega = 2*pi*freq
    #tran = sympy.Symbol("N/A")

    for word in words:
        if word in ["sin", "SIN"]:
            use_DC_val = False
            break
        else:
            index += 1
    if use_DC_val:
        tran = dc*(1/s)
    else:
        try:
            offset, _ = convert_units(words[index])
            amp, _ = convert_units(words[index+1])
            freq, _ = convert_units(words[index+2], forced_numeric=True)
            delay, _ = convert_units(words[index+3])
            damping, _ = convert_units(words[index+4])
            omega = 2*pi*freq
        except IndexError:
            pass
        #tran = ((offset/s)+amp*sympy.exp(-s*delay)*2*pi*freq/((s+damping)**2+(2*pi*freq)**2))
        tran = amp*((damping+s)*sympy.sin(delay)+omega*sympy.cos(delay))/(damping**2+2*damping*s+omega**2+s**2)
    return tran

def value_enum(words):
    try:
        value, symbolic = convert_units(words[3])
    except IndexError:
        symbolic = True
        value = sympy.Symbol(words[0], real=True) #value = sympy.parse_expr(words[0], local_dict=sympy.abc._clash, transformations=TRANSFORMS)
    return value, symbolic

def nodes_per_element(type):
    type = type.lower()
    if type in ["r", "l", "c", "v", "i", "f", "h", "d", "s"]:
        return 2
    if type in ["f", "h", "q", "j"]:
        return 3
    elif type in ["a", "e", "g", "m"]:
        return 4
    elif type in ["k"]:
        return 0

def parse_subcircuits(netlist, operating_points):
    subckt_models = {}
    in_model = False
    subckt_model_id = ""
    current_model = None
    parsed_netlist = []

    for line in netlist[1:]:
        words = line.split()

        if words in [[], "\n", " "]:
            pass

        elif words[0][0] not in NETLIST_KEYCHARS:
            raise NotImplementedError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif words[0] in [".subckt", ".SUBCKT"]:
            in_model = True
            loading_nodes = True
            node_list = []
            param_dict = {}
            subckt_model_id = words[1]
            for w in words[2:]:
                if w in ["PARAMS:", "params:"]:
                    loading_nodes = False
                elif loading_nodes:
                    node_list.append(w)
                else:
                    key, val = w.split("=")
                    param_dict[key], _ = convert_units(val)
            current_model = SubcktModel(subckt_model_id, node_list, param_dict)

        elif words[0] in [".model", ".MODEL"]:
            model_id = words[1]
            model_type = words[2].lower()
            # strip optional spice model brackets
            words[3] = words[3].replace("(", "")
            words[-1] = words[-1].replace(")", "")
            # parse arguments

            param_dict = {"nr": 1}
            for w in words[3:]:
                key, val = w.split("=")
                param_dict[key], _ = convert_units(val)
            if model_type in ["npn"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"NPN model not implemented outside of AC or TF analysis.")
                model = NPNModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["pnp"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"PNP model not implemented outside of AC or TF analysis.")
                model = PNPModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["nmos", "njf"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"NFET model not implemented outside of AC or TF analysis.")
                model = NFETModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["pmos", "pjf"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"PFET model not implemented outside of AC or TF analysis.")
                model = PFETModelAC(model_id, param_dict)
                subckt_models[model_id] = model
            elif model_type in ["d"]:
                #if analysis_type not in ["ac", "AC", "tf", "TF"]:
                #    raise NotImplementedError(
                #        f"Diode model not implemented outside of AC or TF analysis.")
                model = DiodeModelAC(model_id, param_dict)
                subckt_models[model_id] = model

        elif words[0][0] == ".":
            if words[0] in [".ends", ".ENDS"]:
                in_model = False
                subckt_models[subckt_model_id] = current_model
            elif words[0] in [".end", ".END"]:
                break
            else:
                raise NotImplementedError(f"Keyword/Element '{words[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")

        elif in_model:
            current_model.elements.append(line)
        else:
            parsed_netlist.append(line)

    final_netlist = unpack(parsed_netlist, subckt_models, operating_points)
    return final_netlist


def unpack(parsed_netlist, subckt_models, operating_points):
    """
    Identifies all subcircuits, unpacks them and returns a unpacked netlist
    Elements inside a subcircuit inherit it's name in the following format: ElementName_SubcircuitName
    :param list parsed_netlist: each cell contains one netlist line in string format
    :return: list unpacked_netlist: each cell contains one netlist line in string format
    """
    final_netlist = []
    for line in parsed_netlist:
        words = line.split()
        name = words[0]
        if name[0].lower() in ["x", "q", "m", "d", "j"]:
            loading_params = False
            nodes = []
            params = {}
            model = None
            for word in words[1:]:
                if word in ["params:", "PARAMS:"]:
                    if model:
                        loading_params = True
                    else:
                        raise NotImplementedError(f"Model of element '{name}' is not present in the netlist.")

                elif not loading_params:
                    try:
                        model = subckt_models[word]
                        params = model.param_dict
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

            if operating_points is not None and name in operating_points:
                op = operating_points[name]
            else:
                op = None

            elements = model.build_instance(op)
            for elem in elements:
                split_elem = elem.split(" ")
                if split_elem[0][0] not in NETLIST_KEYCHARS: # filter unsupported elements and syntax errors
                    raise NotImplementedError(
                        f"Keyword/Element '{split_elem[0]}' not recognized by netlist parser. Check netlist correctness, if your netlist is correct please submit a bug report on GitHub: 'https://github.com/MatyasVasek/SymCirc'.")
                elif split_elem[0] in ["k", "K"]: # parse couplings
                    split_elem[1] = f"{split_elem[1]}_{words[0]}"
                    split_elem[2] = f"{split_elem[2]}_{words[0]}"
                    if split_elem[3][0] == "{":
                        try:
                            split_elem[3] = str(params[split_elem[3][1:-1]])
                        except KeyError:
                            split_elem[3] = str(model.param_dict[split_elem[3][1:-1]])
                else: # parse the rest
                    if split_elem[0][0] in ["f", "F", "h", "H"]: # correct subcircuit CC(C/V)S control voltage source (current sensor) name
                        split_elem[3] = f"{split_elem[3]}_{words[0]}"
                    node_count = nodes_per_element(split_elem[0][0]) # get expected node count
                    for i in range(1, node_count+1): # scan and connect subcircuit nodes with outside circuit nodes
                        if split_elem[i] in ["0", 0]:
                            pass
                        else:
                            try:
                                split_elem[i] = node_dict[split_elem[i]]
                            except KeyError:
                                split_elem[i] = f"{split_elem[i]}_{words[0]}"
                    index = node_count+1
                    for e in split_elem[node_count+1:]: # substitute parameters
                        if e == "":
                            continue

                        name = words[0]
                        if "IC" not in params: # Needed to differentiate between OP of different semiconductor elements
                            params["IC"] = sympy.Symbol(f"IC_{name}")
                        if "ID" not in params: # Needed to differentiate between OP of different semiconductor elements
                            params["ID"] = sympy.Symbol(f"ID_{name}")

                        local = sympy.abc._clash
                        tmp_elem = None

                        tmp_elem = e.replace("{", "")
                        tmp_elem = tmp_elem.replace("}", "")
                        local.update(params)
                        tmp_elem, _ = convert_units(tmp_elem, local_dict=local)
                        string_tmp = str(tmp_elem).replace(" ", "")
                        split_elem[index] = string_tmp

                        index += 1

                device_name = words[0]
                split_elem = append_dev_name(split_elem, device_name)
                final_netlist.append(" ".join(split_elem))
        else:
            final_netlist.append(line)
    return final_netlist

def append_dev_name(split_elem, device_name):
    if len(split_elem) == 3:
        if split_elem[0] == "rpi":
            split_elem.append(f"1/gpi_{device_name}")
        elif split_elem[0] == "rd":
            split_elem.append(f"1/gd_{device_name}")
        elif split_elem[0] == "ro":
            split_elem.append(f"1/go_{device_name}")

    split_elem[0] = f"{split_elem[0]}_{device_name}"
    return split_elem

def preparse(netist_lines):
    preparsed_netlist_lines = [netist_lines[0]]
    for line in netist_lines[1:]:
        split_line = line.split()
        if split_line:
            first_char = split_line[0][0]
            if first_char in ["m", "M"]:  # add 'PARAMS:' between mos model name and 'W' and 'L' params
                if split_line[6] != "PARAMS:":
                    split_line.insert(6, "PARAMS:")
                preparsed_line = " ".join(split_line)
                preparsed_netlist_lines.append(preparsed_line)
            elif first_char == "*":  # Strip commentaries
                pass
            else:
                preparsed_line = " ".join(split_line)
                preparsed_netlist_lines.append(preparsed_line)
    return preparsed_netlist_lines

def parse(netlist, operating_points=None):
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
    parsed_netlist = preparse(parsed_netlist)
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
    add_short = {}
    parsed_netlist = parse_subcircuits(parsed_netlist, operating_points)

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
            #sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)  # sympy.Symbol(name, real=True)
            dc_num, dc_sym = dc_value(words)
            ac_num, ac_phase, ac_sym = ac_value(words)
            tran_num = tran_value(words, dc_num)
            tran_sym = dc_sym / s
            c = CurrentSource(name, node1, node2,
                              dc_num=dc_num, dc_sym=dc_sym,
                              ac_num=ac_num, ac_sym=ac_sym, ac_phase=ac_phase,
                              tran_num=tran_num, tran_sym=tran_sym)
            independent_sources.append(c)

        elif name[0] in ["v", "V", "u", "U"]:
            #sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)  # sympy.Symbol(name, real=True)

            dc_num, dc_sym = dc_value(words)
            ac_num, ac_phase, ac_sym = ac_value(words)
            tran_num = tran_value(words, dc_num)
            tran_sym = dc_sym/s

            c = VoltageSource(name, node1, node2,
                              dc_num=dc_num, dc_sym=dc_sym,
                              ac_num=ac_num, ac_sym=ac_sym, ac_phase=ac_phase,
                              tran_num=tran_num, tran_sym=tran_sym)
            independent_sources.append(c)

        elif name[0] in ["r", "R"]:
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)

            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)  # sympy.Symbol(name, real=True)
            c = Resistor(name, node1, node2, sym_value=sym_value, value=value)
            basic_components.append(c)

        elif name[0] in ["c", "C"]:
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)  # sympy.Symbol(name, real=True)
            try:
                init_cond, _ = convert_units(words[4][3:])
                c = Capacitor(name, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            except IndexError:
                init_cond = 0
                c = Capacitor(name, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            basic_components.append(c)

        elif name[0] in ["l", "L"]:
            value, symbolic = value_enum(words)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)  # sympy.Symbol(name, real=True)
            try:
                init_cond, _ = convert_units(words[4][3:])
                c = Inductor(name, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            except IndexError:
                init_cond = 0
                c = Inductor(name, node1, node2, sym_value=sym_value, init_cond=init_cond, value=value)
            basic_components.append(c)

        elif name[0] in ["a", "A"]:
            sym_value = sympy.Symbol(name, real=True)
            node3 = words[3]
            node4 = words[4]
            if node3 not in nodes:
                nodes.append(node3)
            if node4 not in nodes:
                nodes.append(node4)
            c = IdealOperationalAmplifier(name, node1, node2, node3, node4)
            operational_amplifiers.append(c)

        elif name[0] in ["e", "E"]:  # VVT (VCVS)
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
                value = sympy.Symbol(name, real=True) #value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)

            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)

            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=value, sym_value=sym_value)
            controlled_sources.append(c)

        elif name[0] in ["g", "G"]:  # VCT (VCCS)
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
                value = sympy.Symbol(name, real=True) #sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = VoltageControlledSource(name, variant, node1, node2, node3, node4, value=value, sym_value=sym_value)
            controlled_sources.append(c)

        elif name[0] in ["f", "F"]:  # CCT (CCCS)
            variant = "f"
            sym_value = sympy.Symbol(name, real=True)
            v_c = words[3]
            try:
                value, symbolic = convert_units(words[4])
            except IndexError:
                symbolic = True
                value = sympy.Symbol(name, real=True) #value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = CurrentControlledSource(name, variant, node1, node2, current_sensor=v_c, value=value, sym_value=sym_value,)
            try:
                add_short[v_c].append(name)
                c.node3 = f"*ctrl{v_c}{add_short[v_c][-2]}"
                c.node4 = f"*ctrl{v_c}{add_short[v_c][-1]}"
            except:
                add_short[v_c] = [name]
                c.node3 = None
                c.node4 = f"*ctrl{v_c}{add_short[v_c][0]}"
            controlled_sources.append(c)

        elif name[0] in ["h", "H"]:  # CVT (CCVS)
            variant = "h"
            sym_value = sympy.Symbol(name, real=True)
            v_c = words[3]
            try:
                value, symbolic = convert_units(words[4])
            except IndexError:
                symbolic = True
                value = sympy.Symbol(name, real=True) #value = sympy.parse_expr(name, local_dict=sympy.abc._clash, transformations=TRANSFORMS)
            if symbolic:
                sym_value = value  # sympy.Symbol(value, real=True)
            else:
                sym_value = sympy.Symbol(name, real=True)
            c = CurrentControlledSource(name, variant, node1, node2, current_sensor=v_c, value=value, sym_value=sym_value)
            try:
                add_short[v_c].append(name)
                c.node3 = f"*ctrl{v_c}{add_short[v_c][-2]}"
                c.node4 = f"*ctrl{v_c}{add_short[v_c][-1]}"
            except:
                add_short[v_c] = [name]
                c.node3 = None
                c.node4 = f"*ctrl{v_c}{add_short[v_c][0]}"
            controlled_sources.append(c)

        elif name[0] in ["k", "K"]:  # coupled inductors
            value, symbolic = value_enum(words)
            sym_value = sympy.Symbol(name, real=True)
            L1 = words[1]
            L2 = words[2]
            c = Coupling(name, L1, L2, sym_value, value)
            couplings.append(c)

        elif name[0] in ["s", "S"]:  # periodic switch used for SC/SI analysis
            variant = "s"
            phase = words[3]
            c = PeriodicSwitch(name, variant, node1, node2, phase=phase)
            SCSI_components.append(c)

        components[c.name] = c
        count += 1

    for key in add_short:
        c_s = components[key]
        new_node = None

        shorted_node = c_s.node1
        for c_name in add_short[key]:
            c = components[c_name]
            if c.node3 is None:
                c.node3 = shorted_node
            new_node = f"*ctrl{c_s.name}{c_name}"
            nodes.append(new_node)

        c_s.node1 = new_node

    for couple in couplings:
        L1 = components[couple.L1]
        L2 = components[couple.L2]
        L1.coupling = couple
        L2.coupling = couple

    data["components"] = components
    data["basic_components"] = basic_components
    data["independent_sources"] = independent_sources
    data["controlled_sources"] = controlled_sources
    data["oamps"] = operational_amplifiers
    data["couplings"] = couplings
    data["SCSI_components"] = SCSI_components
    return components, couplings
