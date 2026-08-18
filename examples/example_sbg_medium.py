from symcirc.sbg import *
from symcirc import AC
import time
from symcirc.utils import *

# netlist = utils.load_file("..\\tests\\netlists\\geec_1652.txt")

netlist = """output stage AB class - GEEC #1652
Rz nodeo 0 10
D1 node3 node5 Dmod
D2 node6 node4 Dmod
D3 nodeo node5 modZD
D4 node6 nodeo modZD
RG node3 node4 1k
M1 nodeVc node3 nodeo nodeo NMOSmodelP L=10u W=500u
M2 nodeVe node4 nodeo nodeo PMOSmodelP L=10u W=500u
M3 node3 node1 nodeVc2 nodeVc2 PMOSmodel L=10u W=160u
M4 node4 node2 nodeVe2 nodeVe2 NMOSmodel L=10u W=160u
M5 node1 node1 nodeVc2 nodeVc2 PMOSmodel L=10u W=160u
M6 node2 node2 nodeVe2 nodeVe2 NMOSmodel L=10u W=160u
IB node1 node2 DC 4.5m AC 0 {0*Pi/180}
C1 nodei node1 100u
Vi nodei 0 DC 0 AC 1 {0*Pi/180} SIN 0 50m 10k 0 0
V1 nodeVc2 nodeVc DC 5
V2 nodeVc 0 DC 15
V3 nodeVe nodeVe2 DC 5
V4 0 nodeVe DC 15
.model Dmod d IS=1e-14
.model modZD d bv=7.5
.model NMOSmodelP NMOS VTO=2.2 KP=0.003 CGSO=2e-7 CGDO=2e-8
.model PMOSmodelP PMOS VTO=-2.2 KP=0.003 CGSO=2e-7 CGDO=2e-8
.model PMOSmodel PMOS VTO=-1.9 KP=0.0016 lambda=0.0125 CGSO=6e-9 CGDO=6e-10
.model NMOSmodel NMOS VTO=1.9 KP=0.0016 lambda=0.0125 CGSO=6e-9 CGDO=6e-10
"""

# Dictionary should contain an operating point for each transistor
# To improve model precision at the cost of performance you can add 'cbd', 'cbs'...
op_dict = {"gd:D1": 94.8e-12,
           "gd:D2": 94.8e-12,
           "gd:D3": 1e-12,
           "gd:D4": 1e-12,
           "gm:M1": 68.3e-3, "cgs:M1": 100e-12, "cgd:M1": 10e-12, "gds:M1": 54.5e-6,
           "gm:M2": 68.3e-3, "cgs:M2": 100e-12, "cgd:M2": 10e-12, "gds:M2": 54.5e-6,
           "gm:M3": 18.1e-3, "cgs:M3": 0.96e-12, "cgd:M3": 0.096e-12, "gds:M3": 54.5e-6,
           "gm:M4": 18.1e-3, "cgs:M4": 0.96e-12, "cgd:M4": 0.096e-12, "gds:M4": 54.5e-6,
           "gm:M5": 15.4e-3, "cgs:M5": 0.96e-12, "cgd:M5": 0.096e-12, "gds:M5": 54.5e-6,
           "gm:M6": 15.4e-3, "cgs:M6": 0.96e-12, "cgd:M6": 0.096e-12, "gds:M6": 54.5e-6,
           }

out_node = "nodeo"
noi = [out_node, "nodei"]

metric = "mag"
#metric = "s_plane_dist"

solver = "gauss"

res_list = []
label_list = []
circuit = Circuit(netlist, operating_point=op_dict)

t0 = time.time()

# Semisymbolic reference plot
ac_analysis = AC(circuit, method="two_graph_node", symbolic=False)
out = ac_analysis.get_node_voltage(out_node)
res_list.append(out)
label_list.append("Semisymbol ref")

# Optimize model with relative error in 100-10K Hz range accuracy
# Vo = -Rz*Vi*gm_M1*gm_M3*rds_M3/(Rz*gm_M1 + Rz*gm_M2 + 1)
max_err = 0.05 # 5%
rel_err = True
circ1 = Circuit(netlist, operating_point=op_dict)
circ1.approximate(noi, max_err, 1e2, 1e4, 10, rel_err, metric, solver=solver)
ac_circ1 = AC(circ1, method="two_graph_node", symbolic=False, solver=solver)
res1 = evalf(ac_circ1.get_node_voltage(out_node), precision=20)
res_list.append(res1)
label_list.append("5%; 100-10k")

# Optimize model with relative error in 100-10K Hz range accuracy
# Vo = -Vi*gm_M1*gm_M3*rds_M3/(gm_M1 + gm_M2)
max_err = 0.80 # 80%
rel_err = True
circ2 = Circuit(netlist, operating_point=op_dict)
circ2.approximate(noi, max_err, 1e2, 1e4, 10, rel_err, metric)
ac_circ2 = AC(circ2, method="two_graph_node", symbolic=False, solver=solver)
res2 = evalf(ac_circ2.get_node_voltage(out_node), precision=20)
res_list.append(res2)
label_list.append("80%; 100-10k")

# Optimize model with relative error in 1-1M Hz range accuracy
# Vo = -C1*Rz*Vi*gm_M3*rds_M3*s*(cgs_M1*s + gm_M1)/((C1*s + gm_M5)*(Rz*cgs_M1*s + Rz*gm_M1 + cgs_M1*rds_M3*s + 1))
max_err = 0.05 # 5%
rel_err = True
circ3 = Circuit(netlist, operating_point=op_dict)
circ3.approximate(noi, max_err, 1e0, 1e6, 10, rel_err, metric)
ac_circ3 = AC(circ3, method="two_graph_node", symbolic=False, solver=solver)
res3 = evalf(ac_circ3.get_node_voltage(out_node), precision=20)
res_list.append(res3)
label_list.append("5%; 1-1M")

# Optimize model for 1-1M Hz range accuracy
# Vo = -C1*Rz*Vi*gm_M3*rds_M3*s*(cgs_M1*s + gm_M1)/((C1*s + gm_M5)*(Rz*cgs_M1*gm_M2*rds_M3*s + Rz*cgs_M1*s + Rz*gm_M1 + Rz*gm_M2 + cgs_M1*rds_M3*s + 1))
max_err = 0.3
rel_err = True
circ4 = Circuit(netlist, operating_point=op_dict)
circ4.approximate(noi, max_err, 1e0, 1e6, 10, rel_err, metric)
ac_circ4 = AC(circ4, method="two_graph_node", symbolic=False, solver=solver)
res4 = evalf(ac_circ4.get_node_voltage(out_node), precision=20)
res_list.append(res4)
label_list.append("30%; 1-1M")

# Optimize model with relative error in 1-1M Hz range accuracy
# Vo = -C1*Rz*Vi*gm_M3*rds_M3*s*(cgs_M1*s + gm_M1)/((C1*s + gm_M5)*(Rz*cgs_M1*s + Rz*gm_M1 + cgs_M1*rds_M3*s + 1))
max_err = 0.8 # 80%
rel_err = True
circ5 = Circuit(netlist, operating_point=op_dict)
circ5.approximate(noi, max_err, 1e0, 1e6, 10, rel_err, metric)
ac_circ5 = AC(circ5, method="two_graph_node", symbolic=False, solver=solver)
res5 = evalf(ac_circ5.get_node_voltage(out_node), precision=20)
res_list.append(res5)
label_list.append("80%; 1-1M")


plot_bode(res_list, 1e-2, 1e8, 1000, "Bode plot of all approximations", label_list=label_list)
