import sys, os
import json

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")

from symcirc import *


netlist = """
*circ
R1 1 0 2k
R2 3 0 (1/G)
V1 2 1 dc 1 ac 1
R3 3 4 6k

C1 3 4 1n
R4 4 0 10k

V2 4 0 dc 5
I1 3 2 dc 1m ac 0"""

circuit = AnalyseCircuit(netlist, "DC", symbolic=True)
print(json.dumps(to_latex(circuit.component_values("all"))))
