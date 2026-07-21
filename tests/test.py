import sys, os
import time
import sympy
from sympy import simplify

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
    #netlist = symcirc.utils.load_file("netlists\\ic_test_1.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_355.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_1536.txt")
    #netlist = symcirc.utils.load_file("netlists\\coupled2.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_1800.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_1529.txt")
    #netlist = symcirc.utils.load_file("netlists\\geec_371.txt")
    netlist = symcirc.utils.load_file("netlists\\AC1.txt")

    analysis_type = "AC"
    symbolic = False

    method = "two_graph_node"
    #method = "tableau"
    solver = "gauss"
    #solver = "ddd"

    if parser_test:
        data = symcirc.parse.parse(netlist)

    if analysis_test:

        """n = utils.load_file(netlist)
        circuit = parse.unpack_subcircuit(n)"""

        t0 = time.time()
        #circuit = symcirc.analysis.Circuit(netlist)
        #op_dict = {"gm:Q1": 74.5e-3, "gpi:Q1": 232e-6, "gmu:Q1": 1e-9, "go:Q1": 22.5e-6, "gx:Q1": 1.66}}
        #op_dict = {"gm:M1": 7.15e-3}
        op_dict = {"gpi:Q1": 345e-6, "gpi:Q2": 373e-6, "gpi:Q3": 386e-6, "gpi:Q4": 386e-6}
        #op_dict = None
        analysis = symcirc.analysis.AnalyseCircuit(netlist, analysis_type, symbolic=symbolic, auto_eval=True, precision=6, method=method, sympy_ilt=True, operating_point=op_dict, solver=solver)
        #analysis = circuit.analyse("dc")
        print(time.time() - t0)

        M = analysis.eqn_matrix
        t2 = time.time()
        #sol = symcirc.solver_DDD_symengine.cramer_ddd_solve(M)
        #print(f"Cramer+DDD solve time: {time.time() - t2}")
        #print(sol)

    if test_prints:
        #print(circuit.components)
        t1 = time.time()

        print("TEST -- print matrix: start")
        print(analysis.eqn_matrix)
        sympy.pprint(analysis.eqn_matrix)
        print("TEST -- print matrix: end")

        print("run time: {}".format(t1 - t0))
        #print(analysis.node_voltage_symbols)
        #print(analysis.solved_dict)
        print(f"Node voltages: {analysis.solved_dict}")
        analysis_tgn = symcirc.analysis.AnalyseCircuit(netlist, analysis_type, symbolic=symbolic, auto_eval=True,
                                                       precision=6, method="two_graph_node", sympy_ilt=True,
                                                       operating_point=op_dict, solver=solver)
        utils.pprint(analysis_tgn.eqn_matrix)
        voltages = analysis_tgn.node_voltages()
        print(f"Node voltages tgn: {voltages}")
        tableau_volt = analysis.node_voltages()
        #for volt in voltages:
        #    print(voltages[volt] == tableau_volt[volt])
        all = analysis.component_values()
        print("---------------------------------------------------------")
        #print(f'v(node5): {sympy.simplify(analysis.get_node_voltage("node5"))}')
        #print("All components: {}".format(all))
        #print(f"Node voltages: {analysis.node_voltages()}")
        #print(analysis.get_all_results())
        #print(analysis.symbols)
        #sympy.pprint(all)
        #utils.latex_print(all)
        #utils.latex_print(circuit.node_voltages())


    if plots:
        if symbolic is False:
            xpoints = []
            ypoints = []
            all_values = analysis.component_values()
            all_voltages = analysis.component_values()
            node_voltages = analysis.node_voltages()
            n = 0
            for symbol_eqn in all_values:
                #try:
                n += 1
                func = all_values[symbol_eqn]

                if func is None:
                    continue

                for symbol in analysis.circuit.components:
                    value = analysis.circuit.components[symbol].value
                    func = func.subs(symbol, value)
                t_symbol = func.free_symbols
                print(func)

                #print(str(node)+": "+str(func))
                if analysis_type.lower() == "tran":
                    test_utils.plot(func, utils.t, 0, 2, 10000, title=symbol_eqn)
                elif analysis_type.lower() == "ac":
                    #test_utils.plot_mag(func, utils.f, 1, 1000000, 100000, title=f"Magnitude plot of {symbol_eqn}")
                    #test_utils.plot_phase(func, utils.f, 1, 1000000, 10000, title=f"Phase plot of {symbol_eqn})")
                    test_utils.plot_bode(func, utils.f, 1, 1000000, 10000, title=f"Bode plot of {symbol_eqn}")
                else:
                    test_utils.plot(func, utils.s, 0, 0.01, 10000, title=symbol_eqn)
                """
                except Exception as e:
                print(e)
                """
        #print(all["i(V1)"].subs(t, 1))
        #print("RunTime: {}".format(t1 - t0))





