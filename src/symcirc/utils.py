import sympy
from sympy import expand, factor, simplify, limit, diff, solve, parse_expr
from sympy import oo as infinity
from sympy import log, exp, sin, cos, tan, cot
from sympy import I as j
from sympy import pi

from typing import Dict

f = sympy.symbols("f", real=True, positive=True)
s = sympy.symbols("s", real=False)
t = sympy.symbols("t", real=True, positive=True)
z = sympy.symbols("z", real=True, positive=True)

global_dict = {"pi": pi,
               "Pi": pi,
               "PI": pi,
               "pI": pi,
               "f": f,
               "s": s,
               "t": t,
               "z": z,
               "oo": infinity,
               "inf": infinity,
               "Inf": infinity,
               "i": j,
               "j": j}

def numer(H):
    N, _ = sympy.fraction(H)
    return N

def denom(H):
    _, D = sympy.fraction(H)
    return D

def evalf(H, subs:dict={}, precision:int=6):
    subs_dict = {**subs, **global_dict}
    return H.evalf(subs=subs_dict, n=precision)

def load_file(netlist_addr):
    """
    File to string netlist converter, provided for user convenience.

    :param str netlist_addr: string containing netlist file address
    :return str netlist: netlist formatted for AnalyseCircuit() input
    """
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist


def to_latex(dictionary: Dict[str, sympy.Expr], symbol_names: Dict[sympy.Symbol, str]=None) -> Dict[str, str]:
    """
    Python dict to latex converter.

    :param dict dictionary: arbitrary dictionary
    :param symbol_names dictionary: contains desired representations of symbols to override default sympy interpretation
    :return str ret: string containing resulting latex code
    """
    ret = {}
    for key in dictionary:
        if symbol_names:
            ret.update({key: sympy.latex(dictionary[key], imaginary_unit="j", symbol_names=symbol_names)})
        else:
            ret.update({key: sympy.latex(dictionary[key], imaginary_unit="j")})
    return ret


def latex_print(data: sympy.Expr, symbol_names: Dict[sympy.Symbol, str]=None):
    """
    Print data in latex code.

    :param str data: arbitrary string
    """
    print("{}".format(sympy.latex(data), imaginary_unit="j", symbol_names=symbol_names))


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
    """
    Deep numeric evaluation via recursion
    """
    try:
        func = sympy.N(func, precision)
        arg_list = []
        for arg in func.args:
            #print(arg)
            arg_list.append(evaluate(arg))
        #print(arg_list)
        if len(arg_list) > 0:
            func = func.func(*arg_list)
    except RecursionError:
        pass
    return func

