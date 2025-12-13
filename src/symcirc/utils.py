import sympy
import random
from sympy import expand, factor, simplify, limit, diff, solve, parse_expr
from sympy import oo as infinity
from sympy import log, exp, sin, cos, tan, cot
from sympy import I as j
from sympy import pi

from typing import Dict
from logging import warning

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
    return factor(H.evalf(subs=subs_dict, n=precision))

def load_file(netlist_addr):
    """
    File to string netlist converter, provided for user convenience.

    :param str netlist_addr: string containing netlist file address
    :return str netlist: netlist formatted for AnalyseCircuit() input
    """
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist


def to_latex(dictionary: Dict[str, sympy.Expr], symbol_names: Dict[str, sympy.Symbol]=None) -> Dict[str, str]:
    """
    Python dict to latex converter.

    :param dict dictionary: arbitrary dictionary
    :param symbol_names dictionary: contains desired representations of symbols to override default sympy interpretation
    :return str ret: string containing resulting latex code
    """
    ret = {}
    for key in dictionary:
        if symbol_names:
            symbol_names_inversed = dict((v, k) for k, v in symbol_names.items())
            ret.update({key: sympy.latex(dictionary[key], imaginary_unit="j", symbol_names=symbol_names_inversed)})
        else:
            ret.update({key: sympy.latex(dictionary[key], imaginary_unit="j")})
    return ret


def latex_print(data: sympy.Expr, symbol_names: Dict[str, sympy.Symbol]=None):
    """
    Print data in latex code.

    :param str data: arbitrary string
    """
    latex_dict = to_latex({'expr': data}, symbol_names=symbol_names)
    print(f"{latex_dict['expr']}")

def xpoints(start, stop, points, log=False):
    """
    Generate a list of points between start and stop.

    Parameters:
        start (float): start value (must be >0 for log scale)
        stop (float): stop value
        points (int): number of points
        x_log (bool): if True, return logarithmically spaced points

    Returns:
        list of floats
    """
    import numpy as np
    if log:
        if start <= 0:
            raise ValueError("Logarithmic spacing requires start > 0")
        x = np.logspace(np.log10(start), np.log10(stop), points)
    else:
        x = np.linspace(start, stop, points)
    return list(x)


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
            arg_list.append(evaluate(arg))
        if len(arg_list) > 0:
            func = func.func(*arg_list)
    except RecursionError:
        pass
    return func

def gauss_points(mean, sd, points):
    res = []
    for i in range(points):
        res.append(random.gauss(mean, sd))
    return res

def plot(func, var, start, stop, points, title="", x_log=False):
    import matplotlib.pyplot as plt
    from numpy import array
    x = xpoints(start, stop, points, log=x_log)
    y = ypoints(func, x, var)
    arrx = array(x)
    arry = array(y)
    plt.plot(arrx, arry)
    plt.title(title)
    plt.show()

def mag(func, decibel=False):
    if func != 0:
        func = abs(func)
        if decibel:
            func = 20*sympy.log(func, 10)
        return func
    else:
        warning("Cannot plot 20*log(mag(func)), because func == 0!")
        return -infinity

def arg(func, deg=False):
    func = sympy.arg(func)
    if deg:
        func = sympy.deg(func)
        func = (func + 180) % 360 - 180
    return func

def plot_mag(func, var, start, stop, points, title="Mag Plot"):
    import matplotlib.pyplot as plt
    func = mag(func, decibel=True)
    plt.xscale("log")
    plt.xlabel("f (Hz)")
    plt.ylabel(r'$20 \log_{10} |H(j\omega)| \;(\mathrm{dB})$')
    plt.grid(True)
    plt.grid(which='both', color='gray', alpha=0.7)
    plot(func, var, start, stop, points, title, x_log=True)

def plot_phase(func, var, start, stop, points, title="Phase Plot"):
    import matplotlib.pyplot as plt
    func = arg(evalf(func))
    plt.xscale("log")
    plt.xlabel("f (Hz)")
    plt.ylabel(r'$\arg(H(j\omega))\;(\mathrm{deg})$')
    plt.grid(True)
    plt.grid(which='both', color='gray', alpha=0.7)
    plot(func, var, start, stop, points, title, x_log=True)

def plot_bode(func, var, start, stop, points=500, title="Bode Plot"):
    import matplotlib.pyplot as plt
    from numpy import array
    # Evaluate magnitude and phase
    func_arg = arg(func, deg=True)  # phase in radians
    func_mag = mag(func, decibel=True)  # magnitude in dB

    # Frequency points (log scale)
    x = xpoints(start, stop, points, log=True)
    y_mag = ypoints(func_mag, x, var)
    y_arg = ypoints(func_arg, x, var)

    arrx = array(x)
    arry_mag = array(y_mag)
    arry_arg = array(y_arg)

    # Create figure
    fig, ax1 = plt.subplots()
    ax1.set_xscale("log")
    ax1.set_xlabel("f (Hz)")
    ax1.set_ylabel(r'$20 \log_{10} |H(j\omega)| \;(\mathrm{dB})$', color='tab:blue')
    ax1.plot(arrx, arry_mag, color='tab:blue', label='Magnitude')
    ax1.tick_params(axis='y', labelcolor='tab:blue')
    ax1.grid(which='both', color='gray', alpha=0.7)

    # Create second y-axis for phase
    ax2 = ax1.twinx()
    ax2.set_ylabel(r'$\arg(H(j\omega)) \;(\mathrm{deg})$', color='tab:orange')
    ax2.plot(arrx, arry_arg, color='tab:orange', label='Phase')  # convert rad -> deg
    ax2.tick_params(axis='y', labelcolor='tab:orange')

    # Title
    plt.title(title)
    fig.tight_layout()
    plt.show()