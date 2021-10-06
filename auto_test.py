from symcirc import *
import os


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
    print("{}SYMBOLIC:".format('\033[93m'))
    print("{} DC analysis test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=True)
    print("{} AC analysis test:".format('\033[95m'))
    test_analysis("AC", is_symbolic=True)
    print("{} TF analysis test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=True)
    print("{} tran analysis test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=True)

    print("{}SEMISYMBOLIC:".format('\033[93m'))
    print("{} DC analysis test:".format('\033[95m'))
    test_analysis("DC", is_symbolic=False)
    print("{} AC analysis test:".format('\033[95m'))
    test_analysis("AC", is_symbolic=False)
    print("{} TF analysis test:".format('\033[95m'))
    test_analysis("TF", is_symbolic=False)
    print("{} tran analysis test:".format('\033[95m'))
    test_analysis("tran", is_symbolic=False)




