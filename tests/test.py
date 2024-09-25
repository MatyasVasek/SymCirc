import sys, os
import time
import sympy
# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
import symcirc
from symcirc import utils
import test_utils

if __name__ == '__main__':
    plots = False
    test_prints = True
    parser_test = False
    analysis_test = True
    #netlist = symcirc.utils.load_file("netlists\\symbulator\\NR11_13_7_tran.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_1685.txt")
    netlist = symcirc.utils.load_file("netlists\\ic_test_1.txt")
    #netlist = symcirc.utils.load_file("netlists\\DC_elem_10.txt")

    method = "two_graph_node"
    #method = "tableau"
    #method = "eliminated_tableau"

    if parser_test:
        data = symcirc.parse.parse(netlist)

    if analysis_test:

        """n = utils.load_file(netlist)
        circuit = parse.unpack_subcircuit(n)"""

        analysis = "tran"
        t0 = time.time()
        circuit = symcirc.analysis.AnalyseCircuit(netlist, analysis, symbolic=False, precision=6, method=method, sympy_ilt=True)
        print(time.time() - t0)
        all = circuit.component_values()
        t1 = time.time()
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

        print("run time: {}".format(t1 - t0))
        print(circuit.node_voltage_symbols)
        print(f"Node voltages: {circuit.node_voltages()}")
        all = circuit.component_values()
        print("---------------------------------------------------------")
        print("All components: {}".format(all))
        print(f"Node voltages: {circuit.node_voltages()}")
        print(circuit.symbols)
        #utils.latex_print(all)
        #utils.latex_print(circuit.node_voltages())


    if plots:
        xpoints = []
        ypoints = []
        all_values = circuit.component_values()
        all_voltages = circuit.component_values()
        node_voltages = circuit.node_voltages()
        n = 0
        for symbol_eqn in all_values:
            #try:
            n += 1
            func = all_values[symbol_eqn]

            for symbol in circuit.components:
                value = circuit.components[symbol].value
                func = func.subs(symbol, value)
            t_symbol = func.free_symbols
            print(func)

            #print(str(node)+": "+str(func))
            if analysis == "tran":
                test_utils.plot(func, utils.t, 0, 2, 10000, title=symbol_eqn)
            else:
                test_utils.plot(func, utils.s, 0, 0.01, 10000, title=symbol_eqn)
            """
            except Exception as e:
            print(e)
            """
    #print(all["i(V1)"].subs(t, 1))
    #print("RunTime: {}".format(t1 - t0))





