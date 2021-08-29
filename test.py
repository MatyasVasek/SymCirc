import sys
from symcirc import *
from laplace import *

if __name__ == '__main__':
    #netlist = "DC_elem_11.txt"
    netlist = "DC_sup_1.txt"
    s = sympy.symbols("s", real=True)
    circuit = AnalyseCircuit(netlist, "DC", symbolic=True)
    #latex_print(circuit.eqn_matrix)
    #latex_print(circuit.solved_dict)
    '''
    latex_print(circuit.eqn_matrix)
    latex_print(circuit.solved_dict)
    tf = sympy.simplify(circuit.transfer_function("1", "o"))
    print("F = " + str(tf))
    func = circuit.freq_to_time(tf)
    latex_print(tf)
    print(sympy.simplify(func))
    '''
    """res = circuit.solved_dict
    keys = []
    for key in res:
        keys.append(key)
    v1 = sympy.limit(res[keys[0]], s, 0)
    print("{}: {}".format(keys[0], v1))
    try:
        v = circuit.component_voltage("R2")
        i = circuit.component_current("R2")
        print(v)
        print(i)
    except:
        pass"""
    print("Dictionary of solved V/C: {}".format(circuit.solved_dict))
    latex_print(circuit.solved_dict)
    #comp = circuit.analyse_component("R2")
    #print("Single component: {}".format(comp))
    all = circuit.analyse_component("all")
    print("All components: {}".format(all))
    #transfer = circuit.transfer_function("1", "o")
    #print("Transfer: {}".format(transfer))



