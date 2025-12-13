from logging import warning

import sympy
import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
from symcirc import utils

def plot(func, var, start, stop, points, title="", x_log=False):
    import numpy
    import matplotlib.pyplot as plt
    x = utils.xpoints(start, stop, points, log=x_log)
    y = utils.ypoints(func, x, var)
    arrx = numpy.array(x)
    arry = numpy.array(y)
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
        return -utils.infinity

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
    func = arg(utils.evalf(func))
    plt.xscale("log")
    plt.xlabel("f (Hz)")
    plt.ylabel(r'$\arg(H(j\omega))\;(\mathrm{deg})$')
    plt.grid(True)
    plt.grid(which='both', color='gray', alpha=0.7)
    plot(func, var, start, stop, points, title, x_log=True)

def plot_bode(func, var, start, stop, points=500, title="Bode Plot"):
    import numpy
    import matplotlib.pyplot as plt
    # Evaluate magnitude and phase
    func_arg = arg(func, deg=True)  # phase in radians
    func_mag = mag(func, decibel=True)  # magnitude in dB

    # Frequency points (log scale)
    x = utils.xpoints(start, stop, points, log=True)
    y_mag = utils.ypoints(func_mag, x, var)
    y_arg = utils.ypoints(func_arg, x, var)

    arrx = numpy.array(x)
    arry_mag = numpy.array(y_mag)
    arry_arg = numpy.array(y_arg)

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


