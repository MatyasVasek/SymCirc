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

    netlists = ["netlists\\L2C3R.txt"]
    #netlists = ["netlists\\3OpAmp.txt"]

    for netlist in netlists:


        t_t = 0
        t_tgn = 0
        t_et = 0

        #tmp = time.time()
        #circ_t = (symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), "TF", symbolic=True, precision=6, method="tableau"))
        #t_t = (time.time()-tmp)

        #tmp1 = time.time()
        circ_tgn = (symcirc.analysis.AnalyseCircuit(symcirc.utils.load_file(netlist), "TF", symbolic=True, precision=6, method="two_graph_node"))
        circ_tgn.all_component_values()
        #t_tgn = (time.time() - tmp1)


        #print(t_t)
        #print(t_tgn)

        #print(circ_t.solved_dict)
        #print(circ_tgn.solved_dict)
