import sympy

def load_file(netlist_addr):
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist
    
def to_latex(dictionary):
    ret = {}
    for key in dictionary:
        ret.update({str(key): sympy.latex(dictionary[key])})
    return ret


def latex_print(data):
    print("{}".format(sympy.latex(data)))