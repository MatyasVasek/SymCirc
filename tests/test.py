import sys, os
import time
import sympy
# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
import symcirc
from symcirc import utils
import test_utils

if __name__ == '__main__':
    plots = True
    test_prints = True
    parser_test = False
    analysis_test = True
    netlist = symcirc.utils.load_file("netlists\\RLLC.txt")

    #method = "two_graph_node"
    method = "tableau"
    #method = "eliminated_tableau"


    if parser_test:
        data = symcirc.parse.parse(netlist)

    if analysis_test:

        """n = utils.load_file(netlist)
        circuit = parse.unpack_subcircuit(n)"""

        analysis = "TF"
        t0 = time.time()
        #s = sympy.symbols("s", real=True)
        t0 = time.time()
        circuit = symcirc.analysis.AnalyseCircuit(netlist, analysis, symbolic=False, precision=6, method=method)
        t1 = time.time()
        print(t1-t0)
        """
        c = circuit.components["I1"]
        F = circuit.component_voltage("C").subs(c.sym_value, c.tran_value)
        print("F = {}".format(F))
        f = residue_laplace(F)
        print("f = {}".format(f))
        """

    if test_prints:
        #print(circuit.components)
        t1 = time.time()

        print("TEST -- print matrix: start")
        sympy.pprint(circuit.eqn_matrix)
        print("TEST -- print matrix: end")

        #print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
        #latex_print(circuit.solved_dict)
        print("run time: {}".format(t1 - t0))
        #all = circuit.solved_dict
        all = circuit.component_values()
        print("---------------------------------------------------------")
        print("All components: {}".format(all))
        #print(f"Node voltages: {circuit.node_voltages()['v(2)']}")
        #latex_print(all)
        #vR = circuit.component_voltage("R")
        #print(vR)
        #vRt = sympy.integrals.transforms.inverse_laplace_transform(vR, symcirc.s, symcirc.t, plane=None)
        #print(vRt)


    if plots:
        xpoints = []
        ypoints = []
        all_voltages = circuit.component_values()
        node_voltages = circuit.node_voltages()
        n = 0
        for symbol_eqn in node_voltages:
            #try:
            n += 1
            voltage = node_voltages[symbol_eqn]

            for symbol in circuit.components:
                value = circuit.components[symbol].value
                voltage = voltage.subs(symbol, value)
            t_symbol = voltage.free_symbols
            print(voltage)

            #print(str(node)+": "+str(voltage))
            if analysis == "tran":
                test_utils.plot(voltage, utils.t, 0, 0.001, 10000, title=symbol_eqn)
            else:
                test_utils.plot(voltage, utils.s, 0, 0.01, 10000, title=symbol_eqn)
            """
            except Exception as e:
            print(e)
            """
    #print(all["i(V1)"].subs(t, 1))
    #print("RunTime: {}".format(t1 - t0))





