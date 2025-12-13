import sys, os
import json
sys.path.append(os.path.dirname(__file__)+"/../src/")
# DC/rz2
import sys
import json
import re
import sympy

from symcirc import *
from symcirc import utils

from importlib.metadata import version
sys.stderr.write(f"Current version of Symcirc: {version('symcirc')}\n")

def get_voltage(result_dict, x, y = None):
	if(y):
		return sympy.simplify(result_dict[f"v({x})"] - result_dict[f"v({y})"])
	else:
		return result_dict[f"v({x})"]

def parse_expression(expr, result_dict):
    tmp = {
        "v": lambda x, y=None: get_voltage(result_dict, x, y),
        "i": lambda x: result_dict[f"i({x})"],
        "q": lambda x: result_dict[f"q({x})"],
        "ph": lambda x: sympy.arg(x)*180/sympy.pi,
		"degtorad": lambda x: x*sympy.pi/180,
		"radtodeg": lambda x: x*180/sympy.pi
    }
    tmp.update(result_dict)
    tmp.update(sympy.abc._clash)
    return sympy.parse_expr(expr, local_dict=tmp)

netlist = """

v3101 node1 0 DC 1 AC 1 0
e1401 node2 0 node1 0 A

r2701 node1 node2 1k 


"""

circuit = AnalyseCircuit(netlist, "TF", symbolic=True, method="tableau")

results = circuit.node_voltages()
results.update(circuit.component_values("all"))



# design var (assign numeric value if checkbox is checked && is not used below in subs!!)
# TODO: temporary implementation



# Apply global substitutions




output_dict = {
	"result_0": parse_expression("v(node1)", results),
	"result_1": parse_expression("v(node2)", results),
	"result_2": parse_expression("v(v3101)", results),
	"result_3": parse_expression("i(v3101)", results),
	"result_4": parse_expression("v(v3101)*i(v3101)", results),
	"result_5": parse_expression("(v(node2)-0)", results),
	"result_6": parse_expression("v(r2701)", results),
	"result_7": parse_expression("i(r2701)", results)
}


symbol_names = circuit.get_symbols()

final_results = to_latex(output_dict, symbol_names)




f = open("{{OUT}}", "w")
f.write(json.dumps(final_results))
f.close()
