import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src")
from symcirc import *
from symcirc import utils

def test_analysis(analysis_type, is_symbolic=True):
    netlists = []
    for filename in os.listdir("{}/netlists".format(os.getcwd())):
        netlists.append(filename)

    for filename in netlists:
        try:
            circuit = AnalyseCircuit(utils.load_file("netlists/{}".format(filename)), analysis_type=analysis_type, symbolic=is_symbolic)
            print("{}[OK]: {} analysis of {} successful".format('\033[96m', analysis_type, filename))
        except:
            print("{}[FAIL]: {} analysis of {} ended with an error".format('\033[91m', analysis_type, filename))
    return 1


def numeric_test(analysis_type, time=1):
    netlists = []
    for filename in os.listdir("{}/netlists".format(os.getcwd())):
        netlists.append(filename)

    for filename in netlists:
        circuit = AnalyseCircuit(utils.load_file("netlists/{}".format(filename)), type=analysis_type, symbolic=False)
        node_voltages = circuit.node_voltages()
        #print(filename)
        for i in node_voltages:
            try:
                value = node_voltages[i].subs(t, time)
            except:
                pass
            value = value.evalf( )
            #print("{}: {}".format(i, value))
    return 1


if __name__ == "__main__":
    #numeric_test("tran")
    print("{}DC:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=False)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("DC", is_symbolic=True)

    print("{}TF:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=False)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("TF", is_symbolic=True)

    print("{}tran:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=False)
    #print("{} Symbolic test:".format('\033[95m'))
    #test_analysis("tran", is_symbolic=True)




