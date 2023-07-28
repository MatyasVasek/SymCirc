import matplotlib.pyplot as plt
import numpy
import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
from symcirc.utils import xpoints, ypoints

def plot(func, var, start, stop, points, title=""):
    x = xpoints(start, stop, points)
    y = ypoints(func, x, var)
    plt.plot(numpy.array(x), numpy.array(y))
    plt.title(title)
    plt.show()
