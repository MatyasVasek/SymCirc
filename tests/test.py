import sys, os
import time

# project root path
sys.path.append(os.path.dirname(__file__)+"/../")
from src.symcirc import *

if __name__ == '__main__':
    t0 = time.time()
    netlist = "netlists\AC6.txt"
    s = sympy.symbols("s", real=True)
    circuit = AnalyseCircuit(load_file(netlist), "tran", symbolic=False)
    #print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
    #latex_print(circuit.solved_dict)
    t1 = time.time()
    print(t1 - t0)
    all = circuit.component_values("all")
    print("---------------------------------------------------------")
    print("All components: {}".format(all))
    #latex_print(all)
    nodes = circuit.node_voltages()
    #print(nodes)
    #print(all["i(V1)"].subs(t, 1))
    #print("RunTime: {}".format(t1 - t0))



