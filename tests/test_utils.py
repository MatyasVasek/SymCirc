import matplotlib.pyplot as plt
import numpy
import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
from symcirc.utils import xpoints, ypoints

def plot(func, var, start, stop, points, title=""):
    x = xpoints(start, stop, points)
    y = ypoints(func, x, var)
    arrx = numpy.array(x)
    arry = numpy.array(y)
    plt.plot(arrx, arry)
    plt.title(title)
    plt.show()
