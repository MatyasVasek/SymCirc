import sys, os
import json
sys.path.append(os.path.dirname(__file__)+"/../src/")
from symcirc import *

if __name__ == '__main__':
    netlist = """
    C 1 0 10m IC=5
    I1 0 1 dc 5m ac 0 """
    circuit = AnalyseCircuit(netlist, "tran", symbolic=False, use_symengine=True)
    result = circuit.node_voltages()
    result.update(circuit.component_values("all"))
    print(json.dumps(to_latex(result)))