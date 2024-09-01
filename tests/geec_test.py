import sys, os
import json

import sympy
from sympy import real_root

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")

from symcirc import *

def parse_expression(expr, result_dict):
    tmp = {
        "v": lambda x: result_dict[f"v({x})"],
        "i": lambda x: result_dict[f"i({x})"]
    }
    tmp.update(result_dict)
    tmp.update(sympy.abc._clash)
    return sympy.parse_expr(expr, local_dict=tmp)

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


def netlist_1799_custom_results():
    netlist_1799 = """1799
    RC C 0 1k
    V1 1 0 dc 1 ac 1
    R1 1 2 1k
    RE E 0 1k
    R2 1 C 1k
    F1 C E Vaux beta
    Vaux 2 E 0
    """
    print("CUSTOM RESULTS FOR CIRCUIT 1799:\n\nANALYSIS RESULTS:")
    circuit = AnalyseCircuit(netlist_1799, "TF", symbolic=True)
    result = circuit.component_values()
    result.update(circuit.node_voltages())
    for key in result:
        print(f"{key}: {result[key]}")

    print("\n")

    # Rovnice 1
    lhs1 = parse_expression("(UE-V1)/R1+UE/RE-beta*(V1-UE)/R1", result)
    rov1 = sympy.Eq(lhs1, 0)
    print("ROV1: ")
    sympy.pprint(rov1)
    print("-----------------------")
    # Rovnice 2
    lhs2 = parse_expression("(UC-V1)/R2+beta*(V1-UE)/R1+UC/RC", result)
    rov2 = sympy.Eq(lhs2, 0)
    print("ROV2")
    sympy.pprint(rov2)
    print("-----------------------")
    # Reseni 12
    UE = sympy.symbols("UE")
    UC = sympy.symbols("UC")
    reseni12 = sympy.solve([rov1, rov2], [UE, UC])
    print("RESENI12")
    sympy.pprint(reseni12)
    print("-----------------------")

    # v(C)-UC
    r4 = parse_expression("simplify(v(C)-UC)", result)
    r4 = r4.subs(reseni12)
    r4 = r4.subs("Vaux", 0)  # Není mi jasné jak to je obecně s "Vaux", má být vždy nulované i v symbolice?
    print("v(C)-UC")
    sympy.pprint(r4)
    print("-----------------------")

    # Rovnice 3
    lhs3 = parse_expression("-V1+R1*IB+RE*(IB+beta*IB)", result)
    rov3 = sympy.Eq(lhs3, 0)
    print("ROV3")
    sympy.pprint(rov3)
    print("-----------------------")

    # Rovnice 4
    lhs4 = parse_expression("-V1+R2*IA+RC*(IA-beta*IB)", result)
    rov4 = sympy.Eq(lhs4, 0)
    print("ROV4")
    sympy.pprint(rov4)
    print("-----------------------")

    # Reseni 34
    IA = sympy.symbols("IA")
    IB = sympy.symbols("IB")
    reseni34 = sympy.solve([rov3, rov4], [IA, IB])
    print("RESENI34")
    sympy.pprint(reseni34)
    print("-----------------------")

    # v(E)-(IB+beta*IB)*RE
    r8 = parse_expression("v(E)-(IB+beta*IB)*RE", result)
    r8 = r8.subs("Vaux", 0)
    r8 = r8.subs(reseni34)
    print("v(E)-(IB+beta*IB)*RE")
    sympy.pprint(sympy.simplify(r8))
    print("-----------------------")

    # ------------------------------------------------------------------
    # simplify(v(E))
    r9 = parse_expression("simplify(v(E))", result)
    r9 = r9.subs('Vaux', 0)
    print("simplify(v(E))")
    sympy.pprint(r9)
    print("-----------------------")

    # v(1)/i(R1)
    r10 = parse_expression("v(1)/i(R1)", result)
    r10 = r10.subs("Vaux", 0)
    print("v(1)/i(R1)")
    sympy.pprint(r10)
    print("-----------------------")

    # v(C), R2=infinity
    r11 = parse_expression("v(C)", result)
    r11 = r11.subs("Vaux", 0)
    r11 = limit(r11, sympy.Symbol("R2"), sympy.oo)
    print("v(C), R2=infinity")
    sympy.pprint(r11)
    print("-----------------------")

    # v(C)/v(1), R2=infinity, beta=infinity
    r12 = parse_expression("v(C)/v(1)", result)
    r12 = r12.subs(circuit.solved_dict)
    r12 = r12.subs("Vaux", 0)
    r12 = limit(r12, sympy.Symbol("beta"), sympy.oo)
    r12 = limit(r12, sympy.Symbol("R2"), sympy.oo)
    print("v(C)/v(1), R2=infinity, beta=infinity")
    sympy.pprint(r12)
    print("-----------------------")


def netlist_1837_custom_results():
    netlist_1837 = """1837
    V1 1 0 dc 0 ac 1 0 sin 0 1 1591 0 0
    R2 2 0 {R/b}
    C1 1 2 C
    R1 1 3 R
    C2 2 o {C*a}
    R3 3 o {R/a}
    C3 3 0 {C*b}
    """
    print("CUSTOM RESULTS FOR CIRCUIT 1837:\n\nANALYSIS RESULTS:")
    circuit = AnalyseCircuit(netlist_1837, "TF", symbolic=True)
    result = circuit.component_values()
    result.update(circuit.node_voltages())
    for key in result:
        print(f"{key}: {result[key]}")
    print("\n")

    # H
    print("H")
    H = parse_expression("simplify(v(o)/v(1))", result)
    sympy.pprint(H)
    print("-----------------------")
    result["H"] = H

    # DH
    #result["C"] = sympy.Symbol("C", positive=True)  # místo assuming C>0
    #result["R"] = sympy.Symbol("R", positive=True)  # místo assuming R>0
    # asumpce je potreba dat do "local_dict" pri pouziti "sympy.parse_expr(expr, local_dict)"
    print("DH")
    DH = parse_expression("collect(expand(denom(H)/b/C**2/R**2), s)", result)
    sympy.pprint(DH)
    print("-----------------------")
    result["DH"] = DH

    # NH
    print("NH")
    NH = parse_expression("collect(expand(numer(H)/b/C**2/R**2), s)", result)
    sympy.pprint(NH)
    print("-----------------------")
    result["NH"] = NH

    # o0
    print("o0")
    o0 = parse_expression("sqrt(DH.coeff(s, 0))", result)
    o0_pos, changed_symbol_dict_o0 = sympy.posify(o0)
    sympy.pprint(o0_pos)
    result["o0"] = o0
    print("-----------------------")

    # Q
    print("Q")
    Q = parse_expression("o0/DH.coeff(s, 1)", result)
    Q_pos, changed_symbol_dict_Q = sympy.posify(Q)
    sympy.pprint(sympy.cancel(Q_pos))
    result["Q"] = Q
    print("-----------------------")

if __name__ == "__main__":
    netlist_1799_custom_results()
    netlist_1837_custom_results()