from symcirc import *

netlist = """
I1 0 a dc 6m
R1 a 0 10k
R2 a 0 5k
"""

circuit = Circuit(netlist)

# You can pass a symbolic value as the 'value' of the component.
# Semisymbolic analysis always uses the 'value' parameter, symbolic analysis uses the 'sym_value' parameter.
R1 = sympy.Symbol("R1")
circuit.change("R1", "value", R1)
R2 = sympy.Symbol("R2")
circuit.change("R2", "value", R2)

dc_analysis = DC(circuit, symbolic=False)
node_a = dc_analysis.get_node_voltage("a")
print(node_a)

param_list = [1, 10, 100, 1000, 3000, 50000, 10000]


plot(node_a, R1, 1, 10e3, 1000, title=f"R1 sweep", x_label="Resistance (Ω)", y_label="Voltage (V)",
     param_list=param_list, param_symbol=R2)
