import sys
from symcirc import *
from laplace import *

if __name__ == '__main__':
    netlist = "netlists\AC5.txt"
    s = sympy.symbols("s", real=True)
    circuit = AnalyseCircuit(load_file(netlist), "AC", symbolic=True)
    #print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
    #latex_print(circuit.solved_dict)
    all = circuit.component_values("all")
    print("---------------------------------------------------------")
    print("All components: {}".format(all))
    latex_print(all)
    nodes = circuit.node_voltages()
    print(nodes)



