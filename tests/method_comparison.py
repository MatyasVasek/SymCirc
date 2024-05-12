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

    netlists = ["netlists\\RC.txt"]

    #method = "two_graph_node"
    #method = "tableau"
    method = "eliminated_tableau"

    for netlist in netlists:


        t_t = 0
        t_tgn = 0
        t_et = 0

        tmp = time.time()
        circ_t = (symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), "TF", symbolic=True, precision=6, method="tableau"))
        t_t = (time.time()-tmp)

        tmp1 = time.time()
        circ_tgn = (symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), "TF", symbolic=True, precision=6, method="two_graph_node"))
        t_tgn = (time.time() - tmp1)

        tmp2 = time.time()
        circ_et = (symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), "TF", symbolic=True, precision=6, method="eliminated_tableau"))
        t_et = (time.time() - tmp2)

        print(t_t)
        print(t_tgn)
        print(t_et)

        print(circ_t.solved_dict)
        print(circ_tgn.solved_dict)
        print(circ_et.solved_dict)

