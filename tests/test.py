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
    netlist = "netlists\\3OpAmp.txt"


    if parser_test:
        data = symcirc.parse.parse(symcirc.utils.load_file(netlist))

    if analysis_test:
        """n = utils.load_file(netlist)
        circuit = parse.unpack_subcircuit(n)"""
        analysis = "DC"
        t0 = time.time()
        s = sympy.symbols("s", real=True)
        circuit = symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), analysis, symbolic=False, precision=6)

        """
        c = circuit.components["I1"]
        F = circuit.component_voltage("C").subs(c.sym_value, c.tran_value)
        print("F = {}".format(F))
        f = residue_laplace(F)
        print("f = {}".format(f))
    
        """
    if test_prints:
        print(circuit.components)
        t1 = time.time()

        print("TEST -- print matrix: start")
        sympy.pprint(circuit.eqn_matrix)
        print("TEST -- print matrix: end")

        #print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
        #latex_print(circuit.solved_dict)
        print("run time: {}".format(t1 - t0))
        all = circuit.component_values()
        print("---------------------------------------------------------")
        print("All components: {}".format(all))
        print(f"Node voltages: {circuit.node_voltages()['v(2)']}")
        #latex_print(all)
    if plots:
        xpoints = []
        ypoints = []
        all_voltages = circuit.component_values()
        n = 0
        for symbol_eqn in all_voltages:
            n += 1
            voltage = all_voltages[symbol_eqn]

            for symbol in circuit.components:
                value = circuit.components[symbol].value
                voltage = voltage.subs(symbol, value)
            t_symbol = voltage.free_symbols
            print(voltage)

            #print(str(node)+": "+str(voltage))
            if analysis == "tran":
                test_utils.plot(voltage, utils.t, 0, 0.001, 10000)
            else:
                test_utils.plot(voltage, utils.s, 0, 0.001, 10000)
    #print(all["i(V1)"].subs(t, 1))
    #print("RunTime: {}".format(t1 - t0))





