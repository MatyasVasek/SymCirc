import sympy

t = sympy.Symbol("t", real=True, positive=True)
s = sympy.Symbol("s", real=True, positive=True)
f = sympy.symbols("f", real=True, positive=True)
j = sympy.symbols("j", real=False)


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


def xpoints(start, stop, points):
    x = []
    step = (stop-start)/points
    for i in range(points):
        x.append(start+i*step)
    return x


def ypoints(func, xpoints, var):
    y = []
    lambda_func = sympy.lambdify(var, func)
    for i in xpoints:
        y.append(lambda_func(i))
    return y


def evaluate(func, precision=6):
    func = sympy.N(func, precision)
    arg_list = []
    for arg in func.args:
        #print(arg)
        arg_list.append(evaluate(arg))
    #print(arg_list)
    if len(arg_list) > 0:
        func = func.func(*arg_list)
    return func

