from symcirc.analysis import AC, ACNumeric, Circuit
from symcirc.utils import *
from symcirc.sbg import *
import time

netlist = load_file("..\\tests\\netlists\\nonlinear\\geec_1684.txt")

# Dictionary should contain an operating point for each transistor
# To improve model precision at the cost of performance you can add 'cbd', 'cbs'...
op_dict = {"gd:D1": 94.8e-12,
           "gd:D2": 94.8e-12,
           "gd:D3": 1e-12,
           "gd:D4": 1e-12,
           "gm:M1": 49.2e-3, "cgs:M1": 100e-12, "cgd:M1": 10e-12, "gds:M1": 54.5e-6,
           "gm:M2": 85.1e-3, "cgs:M2": 100e-12, "cgd:M2": 10e-12, "gds:M2": 54.5e-6,
           "gm:M3": 18.1e-3, "cgs:M3": 0.96e-12, "cgd:M3": 0.096e-12, "gds:M3": 54.5e-6,
           "gm:M4": 18.1e-3, "cgs:M4": 0.96e-12, "cgd:M4": 0.096e-12, "gds:M4": 54.5e-6,
           "gm:M5": 11.8e-3, "cgs:M5": 0.96e-12, "cgd:M5": 0.096e-12, "gds:M5": 54.5e-6,
           "gm:M6": 15.4e-3, "cgs:M6": 0.96e-12, "cgd:M6": 0.096e-12, "gds:M6": 54.5e-6,
           "gm:M7": 13e-3, "cgs:M7": 0.96e-12, "cgd:M7": 0.096e-12, "gds:M7": 54.5e-6,
           "gm:M8": 13e-3, "cgs:M8": 0.96e-12, "cgd:M8": 0.096e-12, "gds:M8": 54.5e-6,
           "gm:M9": 18.2e-3, "cgs:M9": 0.96e-12, "cgd:M9": 0.096e-12, "gds:M9": 54.5e-6,
           "gm:M10": 11.8e-3, "cgs:M10": 0.96e-12, "cgd:M10": 0.096e-12, "gds:M10": 54.5e-6,
           }

fixed_comps = ["Vi"]

out_node = "nodeo"
noi = [out_node, "nodei"]

metric = "mag"
#metric = "s_plane_dist"

res_list = []
label_list = []
circuit = Circuit(netlist, operating_point=op_dict)

# Semisymbolic reference plot
ac_analysis = AC(circuit, method="two_graph_node", symbolic=False, solver="gauss")
out = ac_analysis.get_node_voltage(out_node)

ac_analysis_ddd = AC(circuit, method="two_graph_node", symbolic=False, solver="ddd")
out_ddd = evalf(ac_analysis_ddd.get_node_voltage(out_node), precision=20)

res_list.append(out)
label_list.append("Semisymbol ref gauss")
res_list.append(out_ddd)
label_list.append("Semisymbol ref ddd")
'''t0 = time.time()
plot_bode(res_list, 1, 1e8, 1000, "Bode plot of all approximations", label_list=label_list)
print(f"plot time ddd: {time.time()-t0}")'''


solver = 'ddd'

# Optimize model with relative error in 100-10K Hz range accuracy
# Vo = -Rz*Vi*gm_M1*gm_M3*rds_M3/(Rz*gm_M1 + Rz*gm_M2 + 1)
max_err = 0.05 # 5%
rel_err = True
circ1 = Circuit(netlist, operating_point=op_dict)
#res0, log0 = run_optimization(circ1, noi, max_err, 1e2, 1e4, 10, rel_err, metric, solver=solver)
circ1.approximate(noi, max_err, 1e2, 1e4, 10, rel_err, metric, solver=solver)
ac_circ1 = AC(circ1, method="two_graph_node", symbolic=False, solver="ddd")
res1 = evalf(ac_circ1.get_node_voltage(out_node), precision=20)
res_list.append(res1)
label_list.append("5%; 100-100k")

# Optimize model with relative error in 100-10K Hz range accuracy
# Vo = -Rz*Vi*gm_M1*gm_M3*rds_M3/(Rz*gm_M1 + Rz*gm_M2 + 1)
max_err = 0.5 # 50%
rel_err = True
circ2 = Circuit(netlist, operating_point=op_dict)
circ2.approximate(noi, max_err, 1e2, 1e4, 10, rel_err, metric, solver=solver)
ac_circ2 = AC(circ2, method="two_graph_node", symbolic=False, solver="ddd")
res2 = evalf(ac_circ2.get_node_voltage(out_node), precision=20)
res_list.append(res2)
label_list.append("50%; 100-10k")

# Optimize model with relative error in 100-100k Hz range accuracy
# Vo = -Vi*gm_M1*gm_M3*rds_M3/(gm_M1 + gm_M2)
max_err = 0.8 # 80%
rel_err = True
circ3 = Circuit(netlist, operating_point=op_dict)
circ3.approximate(noi, max_err, 1e2, 1e4, 10, rel_err, metric, solver=solver)
ac_circ3 = AC(circ3, method="two_graph_node", symbolic=False, solver="ddd")
res3 = evalf(ac_circ3.get_node_voltage(out_node), precision=20)
res_list.append(res3)
label_list.append("100%; 100-10k")

# Optimize model for 1-100M Hz range accuracy
# Vo = -C1*Rz*Vi*gm_M3*rds_M3*s*(cgs_M1*s + gm_M1)/((C1*s + gm_M5)*(Rz*cgs_M1*gm_M2*rds_M3*s + Rz*cgs_M1*s + Rz*gm_M1 + Rz*gm_M2 + cgs_M1*rds_M3*s + 1))
max_err = 0.5 # 50%
rel_err = True
circ4 = Circuit(netlist, operating_point=op_dict)
circ4.approximate(noi, max_err, 1, 1e8, 10, rel_err, metric, solver=solver)
ac_circ4 = AC(circ4, method="two_graph_node", symbolic=False, solver="ddd")
res4 = evalf(ac_circ4.get_node_voltage(out_node), precision=20)
res_list.append(res4)
label_list.append("50%; 1-100M")

# Optimize model with relative error in 1-100M Hz range accuracy
# Vo = -Vi*gm_M1*gm_M3*rds_M3/(gm_M1 + gm_M2)
max_err = 0.8 # 80%
rel_err = True
circ5 = Circuit(netlist, operating_point=op_dict)
circ5.approximate(noi, max_err, 1, 1e8, 10, rel_err, metric, solver=solver)
ac_circ5 = AC(circ5, method="two_graph_node", symbolic=False, solver="ddd")
res5 = evalf(ac_circ5.get_node_voltage(out_node), precision=20)
res_list.append(res5)
label_list.append("100%; 1-100M")

plot_bode(res_list, 1, 1e8, 1000, "Bode plot of all approximations", label_list=label_list)
