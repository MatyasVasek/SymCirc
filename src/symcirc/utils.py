import sympy
import random
from sympy import expand, factor, simplify, limit, diff, solve, parse_expr, factor_terms
from sympy import oo as infinity
from sympy import log, exp, sin, cos, tan, cot
from sympy import I as j
from sympy import pi
from sympy import pretty #pprint

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

def pprint(A):
    sympy.pprint(A, wrap_line=False)

def numer(H):
    N, _ = sympy.fraction(H)
    return N

def denom(H):
    _, D = sympy.fraction(H)
    return D

'''def evalf(H, subs:Dict={}, precision:int=6):
    subs_dict = {**subs, **global_dict}
    num, den = sympy.fraction(H)
    H = sympy.simplify(num / den.as_leading_term(f))
    print(H)
    print(sympy.N(H))
    return factor_terms(H.evalf(subs=subs_dict, n=precision))
'''
def normalize(H):
    num, den = sympy.fraction(H)
    # leading coefficient
    try:
        leading_coeff = sympy.LC(den)
    except:
        return H
    scale = leading_coeff.as_coeff_Mul()[0]
    H = (num / scale) / (den / scale)
    return H

def evalf(H, subs: Dict = {}, precision: int = 6, norm=True):
    subs_dict = {**subs, **global_dict}
    H = H.subs(subs_dict)
    H = sympy.together(H)

    H = H.evalf(precision)
    if norm:
        H = normalize(H)
    return H

def load_file(netlist_addr):
    """
    File to string netlist converter, provided for user convenience.

    :param str netlist_addr: string containing netlist file address
    :return str netlist: netlist formatted for AnalyseCircuit() input
    """
    with open(netlist_addr) as f:
        netlist = f.read()
    return netlist

def latex(expr: sympy.Expr, symbol_names: Dict[str, sympy.Symbol]=None):
    if symbol_names:
        symbol_names_inversed = dict((v, k) for k, v in symbol_names.items())
        ret = sympy.latex(expr, imaginary_unit="j", symbol_names=symbol_names_inversed)
    else:
        ret = sympy.latex(expr, imaginary_unit="j")
    return ret

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

def simplify_poly_coeffs(poly, var, prefix, memo, start_index=0):
    """
    Extracts coefficients, creates dummy symbols,
    and returns the new expression and the mapping.
    """
    coeffs = poly.all_coeffs()  # Gets all coeffs including zeros
    new_terms = []
    degree = poly.degree()

    for i, coeff in enumerate(coeffs):
        if len(coeff.free_symbols) < 2:
            symbol = coeff
        elif coeff in memo:
            symbol = memo[coeff]
            #print("MEMO HIT!", coeff)
        else:
            symbol = sympy.Symbol(f"{prefix}{start_index + i}", real=True)
            memo[coeff] = symbol

        # Reconstruct the polynomial term with the new symbol
        power = degree - i
        new_terms.append(symbol * var ** power)
    return sum(new_terms)

def simplify_coeffs(expr, var, prefix, memo, start_index=0):
    """
    Wraps simplify_poly_coeffs() for more general use on expressions.
    """
    if expr == 1 or expr == 0 or len(expr.free_symbols) < 2:
        return expr
    num, den = sympy.fraction(sympy.cancel(expr))
    num_poly = sympy.Poly(num, var)
    den_poly = sympy.Poly(den, var)
    new_num = simplify_poly_coeffs(num_poly, var, f"{prefix}_n", memo, start_index)
    new_den = simplify_poly_coeffs(den_poly, var, f"{prefix}_d", memo, start_index)
    new_expr = new_num / new_den
    return new_expr

def gauss_points(mean, sd, points):
    res = []
    for i in range(points):
        res.append(random.gauss(mean, sd))
    return res


def plot(func, var, start, stop, points, title="", x_label="x_axis", y_label="y_axis",
         x_log=False, param_list=None, param_symbol=None) -> "matplotlib.figure":
    import matplotlib.pyplot as plt
    from numpy import array

    fig, ax = plt.subplots()
    if x_log:
        ax.set_xscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.tick_params(axis='y')
    ax.grid(which='both', color='gray', alpha=0.7)

    if param_list is not None:
        i = 0
        for param_val in param_list:
            func_temp = func.subs(param_symbol, param_val)
            x = xpoints(start, stop, points, log=x_log)
            y = ypoints(func_temp, x, var)
            arrx = array(x)
            arry = array(y)
            ax.plot(arrx, arry, label=f'{param_symbol}={param_val}')
            i+=1
    else:
        x = xpoints(start, stop, points, log=x_log)
        y = ypoints(func, x, var)
        arrx = array(x)
        arry = array(y)
        ax.plot(arrx, arry)

    # Title
    ax.legend(loc="upper left")
    plt.title(title)
    fig.tight_layout()
    plt.show()

    return fig


def mag(func, decibel=False, do_evalf=True):
    if do_evalf:
        func = evalf(func, precision=20)
    if func != 0:
        # some error is accrued during evalf, sympy for some reason returns a negligible Im from abs
        func = sympy.re(abs(func))
        if decibel:
            func = 20*sympy.log(func, 10)
            if do_evalf:
                func = evalf(func, precision=20)
        return func
    else:
        return -infinity

def arg(func, deg=False, do_evalf=True):
    if do_evalf:
        func = evalf(func, precision=20)
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

def bode(func, start, stop, points=500):
    # Evaluate magnitude and phase
    func_arg = arg(func, deg=True, do_evalf=True)  # phase in radians
    func_mag = mag(func, decibel=True, do_evalf=True)  # magnitude in dB

    # Frequency points (log scale)
    x = xpoints(start, stop, points, log=True)
    y_mag = ypoints(func_mag, x, f)
    y_arg = ypoints(func_arg, x, f)

    return x, y_mag, y_arg


def plot_bode(func, start, stop, points=500, title="Bode Plot", param_list=None, param_symbol=None, label_list=None) -> "matplotlib.figure":
    import matplotlib.pyplot as plt
    from numpy import array

    # Create figure
    fig, ax_mag = plt.subplots()
    ax_mag.set_xscale("log")
    ax_mag.set_xlabel("f (Hz)")
    ax_mag.set_ylabel(r'$20 \log_{10} |H(j\omega)| \;(\mathrm{dB})$')
    ax_mag.tick_params(axis='y')
    ax_mag.grid(which='both', color='gray', alpha=0.7)

    # Create second y-axis for phase
    ax_phase = ax_mag.twinx()
    ax_phase.set_ylabel(r'$\arg(H(j\omega)) \;(\mathrm{deg})$')
    ax_phase.tick_params(axis='y')

    if param_list is not None:
        i = 0
        for param_val in param_list:
            func_temp = func.subs(param_symbol, param_val)
            x, y_mag, y_arg = bode(func_temp, start, stop, points)
            arrx = array(x)
            arry_mag = array(y_mag)
            arry_arg = array(y_arg)
            ax_mag.plot(arrx, arry_mag, label=f'mag_{i}')
            ax_phase.plot(arrx, arry_arg, label=f'arg_{i}', linestyle='--')
            i+=1
    elif isinstance(func, list):
        i = 0
        if label_list is None:
            labels = [x for x in range(len(func))]
        else:
            labels = label_list
        for func_temp in func:
            x, y_mag, y_arg = bode(func_temp, start, stop, points)
            arrx = array(x)
            arry_mag = array(y_mag)
            arry_arg = array(y_arg)
            ax_mag.plot(arrx, arry_mag, label=f'mag_{label_list[i]}')
            ax_phase.plot(arrx, arry_arg, label=f'arg_{label_list[i]}', linestyle='--')
            i += 1
    else:
        x, y_mag, y_arg = bode(func, start, stop, points)
        arrx = array(x)
        arry_mag = array(y_mag)
        arry_arg = array(y_arg)
        ax_mag.plot(arrx, arry_mag, label=f'mag')
        ax_phase.plot(arrx, arry_arg, label=f'arg', linestyle='--')

    # Title
    ax_mag.legend(loc="lower left")
    ax_phase.legend(loc="lower right")
    plt.title(title)
    fig.tight_layout()
    plt.show()

    return fig


def plot_bode_num(values, control_freqs):
    import matplotlib.pyplot as plt
    from numpy import array
    from sympy import sympify

    # Create figure
    fig, ax_mag = plt.subplots()
    ax_mag.set_xscale("log")
    ax_mag.set_xlabel("f (Hz)")
    ax_mag.set_ylabel(r'$20 \log_{10} |H(j\omega)| \;(\mathrm{dB})$')
    ax_mag.tick_params(axis='y')
    ax_mag.grid(which='both', color='gray', alpha=0.7)

    # Create second y-axis for phase
    ax_phase = ax_mag.twinx()
    ax_phase.set_ylabel(r'$\arg(H(j\omega)) \;(\mathrm{deg})$')
    ax_phase.tick_params(axis='y')

    x = control_freqs
    y_mag = [mag(sympify(v), decibel=True) for v in values]
    y_arg = [arg(sympify(v), deg=True) for v in values]
    arrx = array(x)
    arry_mag = array(y_mag)
    arry_arg = array(y_arg)
    ax_mag.plot(arrx, arry_mag, label=f'mag')
    ax_phase.plot(arrx, arry_arg, label=f'arg', linestyle='--')

    # Title
    ax_mag.legend(loc="lower left")
    ax_phase.legend(loc="lower right")
    plt.title("Bode of bjt lab")
    fig.tight_layout()
    plt.show()


def plot_latex(latex_str, title=None):
    import matplotlib.pyplot as plt
    plt.figure()
    plt.axis('off')
    if title:
        plt.title(title)
    plt.text(0.5, 0.5, f"${latex_str}$", ha='center', va='center', fontsize=9)
    plt.show()
