import json
from typing import Dict
from symcirc import *

if __name__ == "__main__":

    netlist = """

    V1 1 0 DC 10 AC 0 0
    R1 1 2 1k 
    R2 2 0 1k 

    """

    circuit = AnalyseCircuit(netlist, "DC", symbolic=True, method="tableau")
    results = circuit.node_voltages()
    results.update(circuit.component_values("all"))

    print(results)

    symbol_names: Dict[sympy.Symbol, str] = {}
    for key in circuit.components:
        component = circuit.components[key]
        symbol_names[component.sym_value] = component.name

    final_results = to_latex(results, symbol_names)

    output_dict = {
        "result_0": final_results["v(V1)"],
        "result_1": final_results["v(R2)"],
        # ...
    }

    print(json.dumps(output_dict))
