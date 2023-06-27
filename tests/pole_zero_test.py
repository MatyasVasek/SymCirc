import sys, os
# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
import symcirc
import sympy

if __name__ == '__main__':
    netlist = "netlists\\AC6.txt"
    """n = utils.load_file(netlist)
    circuit = parse.unpack_subcircuit(n)"""
    analysis = "TF"
    s = sympy.symbols("s", real=True)
    circuit = symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), analysis, symbolic=True, precision=6)

    tf = circuit.transfer_function(1, 2)
    print("Transfer function: {}".format(tf))
    p, z = symcirc.analysis.pole_zero(tf)
    print("Poles: {}".format(p))
    print("Zeros: {}".format(z))