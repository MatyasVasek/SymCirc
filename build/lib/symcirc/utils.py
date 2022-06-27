import sympy


def load_file(netlist_addr):
    """
    File to string netlist converter, provided for user convenience.

    :param str netlist_addr: string containing netlist file address
    :return str netlist: netlist formatted for AnalyseCircuit() input
    """
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist


def to_latex(dictionary):
    """
    Python dict to latex converter.

    :param dict dictionary: arbitrary dictionary
    :return str ret: string containing resulting latex code
    """
    ret = {}
    for key in dictionary:
        ret.update({str(key): sympy.latex(dictionary[key])})
    return ret


def latex_print(data):
    """
    Print data in latex code.

    :param str data: arbitrary string
    """
    print("{}".format(sympy.latex(data)))