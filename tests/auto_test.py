import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src")
from symcirc import *
from symcirc import utils

def test_analysis(analysis_type, analysis_method="tableau", is_symbolic=True):
    netlists = []
    for filename in os.listdir("{}/netlists".format(os.getcwd())):
        netlists.append(filename)

    for filename in netlists:
        try:
            circuit = AnalyseCircuit(utils.load_file("netlists/{}".format(filename)), analysis_type=analysis_type, symbolic=is_symbolic, method=analysis_method, use_symengine=True)
            print("{}[OK]: {} analysis of {} successful".format('\033[96m', analysis_type, filename))
        except:
            print("{}[FAIL]: {} analysis of {} ended with an error".format('\033[91m', analysis_type, filename))
    return 1


if __name__ == "__main__":
    #method = "two_graph_node"
    method = "tableau"
    #numeric_test("tran")
    print("{}DC:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=False, analysis_method=method)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("DC", is_symbolic=True, analysis_method="two_graph_node")

    print("{}TF:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=False, analysis_method=method)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("TF", is_symbolic=True, analysis_method=method)

    print("{}tran:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=False, analysis_method=method)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("tran", is_symbolic=True, analysis_method=method)




