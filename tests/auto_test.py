import sys, os

# project root path
sys.path.append(os.path.dirname(__file__)+"/../")
from symcirc import *

def test_analysis(analysis_type, is_symbolic=True):
    netlists = []
    for filename in os.listdir("{}/netlists".format(os.getcwd())):
        netlists.append(filename)

    for filename in netlists:
        try:
            circuit = AnalyseCircuit(load_file("netlists/{}".format(filename)), analysis_type, is_symbolic)
            print("{}[OK]: {} analysis of {} successful".format('\033[96m', analysis_type, filename))
        except:
            print("{}[FAIL]: {} analysis of {} ended with an error".format('\033[91m', analysis_type, filename))
    return 1


if __name__ == "__main__":
    print("{}DC:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=False)
    print("{} Symbolic test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=True)

    print("{}TF:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=False)
    print("{} Symbolic test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=True)

    print("{}AC:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("AC", is_symbolic=False)
    print("{} Symbolic test:".format('\033[95m'))
    test_analysis("AC", is_symbolic=True)

    print("{}tran:".format('\033[93m'))
    print("{} Numeric test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=False)
    print("{} Symbolic test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=True)




