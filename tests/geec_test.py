import sys, os
import json

import sympy

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")

from symcirc import *

'''
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
'''
netlist_1799 = """1799
RC C 0 1k
V1 1 0 dc 1 ac 1
R1 1 2 1k
RE E 0 1k
R2 1 C 1k
F1 C E Vaux beta
Vaux 2 E 0
"""
netlist_1837 = """1837
V1 1 0 dc 0 ac 1 0 sin 0 1 1591 0 0
R2 2 0 {R/b}
C1 1 2 C
R1 1 3 R
C2 2 o {C*a}
R3 3 o {R/a}
C3 3 0 {C*b}
"""

circuit = AnalyseCircuit(netlist_1799, "TF", symbolic=True)
result = circuit.component_values()
result.update(circuit.node_voltages())
# Rovnice 1
lhs1 = circuit.parse_expr("(UE-V1)/R1+UE/RE-beta*(V1-UE)/R1")
rov1 = sympy.Eq(lhs1, 0)
print("ROV1")
sympy.pprint(rov1)
# Rovnice 2
lhs2 = circuit.parse_expr("(UC-V1)/R2+beta*(V1-UE)/R1+UC/RC")
rov2 = sympy.Eq(lhs2, 0)
print("\nROV2")
sympy.pprint(rov2)
# Reseni 12
UE = sympy.symbols("UE")
UC = sympy.symbols("UC")
reseni12 = sympy.solve([rov1, rov2], [UE, UC])
print("\nRESENI12")
sympy.pprint(reseni12)

# v(C)-UC
tmp = circuit.parse_expr("v(C)-UC")
for i in result:
    print(f"{i}: {result[i]}")

tmp = tmp.subs(result)
tmp = tmp.subs(reseni12)
tmp = tmp.subs("Vaux", 0)
print("\nv(C)-UC")
sympy.pprint(sympy.simplify(tmp))

tf1 = circuit.parse_expr("v(C)/v(1)")
tf1 = tf1.subs(circuit.solved_dict)
tf1 = tf1.subs("Vaux", 0)
tf1 = limit(tf1, sympy.Symbol("beta"), sympy.oo)
tf1 = limit(tf1, sympy.Symbol("R2"), sympy.oo)
print("\nv(C)/v(1), R2=infinity, beta=infinity")
sympy.pprint(tf1)

# ------------------------------------------------------------------
# Rovnice 3
lhs3 = circuit.parse_expr("-V1+R1*IB+RE*(IB+beta*IB)")
rov3 = sympy.Eq(lhs3, 0)
print("\nROV3")
sympy.pprint(rov3)
# Rovnice 4
lhs4 = circuit.parse_expr("-V1+R2*IA+RC*(IA-beta*IB)")
rov4 = sympy.Eq(lhs4, 0)
print("\nROV4")
sympy.pprint(rov4)
# Reseni 34
IA = sympy.symbols("IA")
IB = sympy.symbols("IB")
reseni34 = sympy.solve([rov3, rov4], [IA, IB])
print("\nRESENI34")
sympy.pprint(reseni34)
# ------------------------------------------------------------------

