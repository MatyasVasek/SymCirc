from symcirc import *

netlist = """Basic HP
V1 1 0 DC 1 AC 1
R1 2 0
C1 1 2 10n
"""

circuit = Circuit(netlist)

# Semisymbolic analysis with all components evaluated to their numeric values
ac_analysis = AC(circuit, symbolic=False)
vnode2 = ac_analysis.get_node_voltage("2")
print(vnode2)

# You can pass a symbolic value as the 'value' of the component.
# Semisymbolic analysis always uses the 'value' parameter, symbolic analysis uses the 'sym_value' parameter.
R1 = sympy.Symbol("R1")
circuit.change("R1", "value", R1)

ac_analysis = AC(circuit, symbolic=False)
vnode2 = ac_analysis.get_node_voltage("2")
print(vnode2)

# Parametric bode plot:
param_list = [1, 17, 100, 1342, 10000, 100000, 10**6]
plot_bode(vnode2, 1, 10e6, 1000, param_list=param_list, param_symbol=R1, title=f"R1 sweep")




