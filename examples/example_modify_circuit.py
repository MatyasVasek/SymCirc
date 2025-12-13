from symcirc import *

netlist = """Basic HP
V1 1 0 DC 1 AC 1
R1 2 0
C1 1 2 10n
"""

circuit = Circuit(netlist)

# Any component parameter can be changed
print(f'Old R1 semisymbolic value: {circuit.get("R1").value}')
circuit.change("R1", "value", 10**6)
print(f'New R1 semisymbolic value: {circuit.get("R1").value}')
circuit.change("V1", "ac_num", 6)
print(f'New V1 ac_val: {circuit.get("V1").ac_num}')

# AC analysis without phase shift
ac_analysis = AC(circuit, symbolic=False)
vnode2 = ac_analysis.get_node_voltage("2")
print("Plotting bode...")
plot_bode(vnode2, f, 1, 10**6, 10000, f"Bode plot of v(2)")

# Set the phase of V1 to 3*pi/4
circuit.change("V1", "ac_phase", parse_expr("3*pi/4"))
print(f'New V1 ac_phase: {circuit.get("V1").ac_phase}')

# AC analysis with phase shift
ac_analysis = AC(circuit, symbolic=False)
vnode2_shift = ac_analysis.get_node_voltage("2")
print("Plotting bode...")
plot_bode(vnode2_shift, f, 1, 10**6, 10000, f"Bode plot of v(2) with 3*pi/4 shift")




