import sys, os
import time
import sympy

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")

from symcirc import *

if __name__ == '__main__':
    netlist = "netlists\\coupled.txt"
    """n = utils.load_file(netlist)
    circuit = parse.unpack_subcircuit(n)"""

    t0 = time.time()
    s = sympy.symbols("s", real=True)
    circuit = AnalyseCircuit(load_file(netlist), "tran", symbolic=False)
    print("TEST -- print matrix: start")
    sympy.pprint(circuit.eqn_matrix)
    print("TEST -- print matrix: end")
    #print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
    #latex_print(circuit.solved_dict)
    t1 = time.time()
    print("run time: {}".format(t1 - t0))
    all = circuit.component_values("all")
    print("---------------------------------------------------------")
    print("All components: {}".format(all))
    #latex_print(all)
    nodes = circuit.node_voltages()
    print(nodes)
    #print(all["i(V1)"].subs(t, 1))
    #print("RunTime: {}".format(t1 - t0))




