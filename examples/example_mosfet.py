from symcirc import *

netlist = """Simple PMOS circuit
C1 i 2 10u
C2 3 o 10n
Rz o 0 1Meg
Vi i 0 dc 0 ac 1 0 sin 0 10m 10k 0 0
VN 1 0 dc -10
RD 3 1 5k
I1 4 2 dc 1m ac 0
V1 4 0 dc 10
M1 3 0 2 2 PMOSmodel L=10u W=160u
.model PMOSmodel PMOS VTO=-1.9 KP=0.0016
"""
# Dictionary should contain an operating point for each transistor
# To improve model precision at the cost of performance you can add 'cbd', 'cbs'...
op_dict = {"M1": {"gm": 7.15e-3}}

circuit = Circuit(netlist, operating_points=op_dict)

ac_analysis = AC(circuit, symbolic=False)

results = ac_analysis.get_all_results()

for i in results:
    tmp = utils.evalf(results[i])
    print(f"{i}: {tmp}")

plot_bode(results["v(o)"], 1, 10 ** 6, 10000, f"Bode plot of v(o)")
