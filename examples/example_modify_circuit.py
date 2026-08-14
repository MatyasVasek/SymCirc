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
print(f'New V1 ac_val: {circuit.get("V1").ac_val("num")}')
circuit.change("V1", "ac_value", ["6", "0"])
print(f'New V1 ac_val: {circuit.get("V1").ac_val("num")}')

# AC analysis without phase shift
ac_analysis = AC(circuit, symbolic=False)
vnode2 = ac_analysis.get_node_voltage("2")
print("Plotting bode...")
plot_bode(vnode2, 1, 10**6, 10000, f"Bode plot of v(2)")

# Set the phase of V1 to 3*pi/4
circuit.change("V1", "ac_value", ["6", "3*pi/4"])
print(f'New V1 ac_phase: {circuit.get("V1").ac_val("num")}')

# AC analysis with phase shift
ac_analysis = AC(circuit, symbolic=False)
vnode2_shift = ac_analysis.get_node_voltage("2")
print("Plotting bode...")
plot_bode(vnode2_shift, 1, 10**6, 10000, f"Bode plot of v(2) with 3*pi/4 shift")

circuit.change("V1", "tran_value", ["SIN", "0", "1", "20", "0", "0"])
tran_analysis = TRAN(circuit, symbolic=False)
vnode2_tran = tran_analysis.get_node_voltage("2")
plot(vnode2_tran, t, 0, 1, 1000, f"Tran plot of v(2) with sin source", "t (s)", "V")

