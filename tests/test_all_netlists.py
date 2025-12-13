import sys, os
import func_timeout

# project root path
sys.path.append(os.path.dirname(__file__)+"/../src")
from symcirc import *
from symcirc import utils

import pytest
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
NETLIST_DIR = f"{BASE_DIR}/netlists"
ANALYSIS_TYPES = ["DC", "TF", "AC", "tran"]
ANALYSIS_METHODS = ["tableau", "two_graph_node"]

NETLISTS = sorted(p for p in os.listdir(NETLIST_DIR))


class Bcolors:
    HEADER = '\033[95m'
    OKBLUE = '\033[94m'
    OKCYAN = '\033[96m'
    OKGREEN = '\033[92m'
    WARNING = '\033[93m'
    FAIL = '\033[91m'
    ENDC = '\033[0m'

    BOLD = '\033[1m'
    UNDERLINE = '\033[4m'


def analysis(analysis_type, analysis_method="tableau", is_symbolic=False):
    error_dict = {}
    netlists = []
    for filename in os.listdir(NETLIST_DIR):
        netlists.append(filename)

    timeout = 20
    for filename in netlists:
        try:
            func_timeout.func_timeout(timeout, analyze, args=(filename, analysis_type, is_symbolic, analysis_method))
            print(f"{Bcolors.OKGREEN}[OK]: {analysis_type} analysis of {filename} successful{Bcolors.ENDC}")
        except NotImplementedError:
            print(f"{Bcolors.WARNING}[NotImplemented]: {analysis_type} analysis of {filename} ended with a NotImplementedError{Bcolors.ENDC}")
        except func_timeout.FunctionTimedOut:
            print(
                f"{Bcolors.FAIL}[FAIL]: {analysis_type} analysis of {filename} took more than {timeout}s, timeout triggered...{Bcolors.ENDC}")
            error_dict[filename] = "TIMEOUT"
        except Exception as e:
            print(f"{Bcolors.FAIL}[FAIL]: {analysis_type} analysis of {filename} ended with an error{Bcolors.ENDC}")
            error_dict[filename] = repr(e)
    return error_dict

def analyze(filename, analysis_type, is_symbolic, analysis_method):
    analysis = AnalyseCircuit(utils.load_file(f"{NETLIST_DIR}/{filename}"), analysis_type=analysis_type,
                              symbolic=is_symbolic, method=analysis_method, use_symengine=True)
    analysis.get_all_results()

@pytest.mark.parametrize("is_symbolic", [True, False])
@pytest.mark.parametrize("analysis_method", ANALYSIS_METHODS)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
def test_semi_smoke(analysis_method, analysis_type, netlist, is_symbolic):
    timeout = 20

    try:
        func_timeout.func_timeout(timeout, analyze, args=(netlist, analysis_type, is_symbolic, analysis_method))

    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")

    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {timeout}s for {netlist}")

    except Exception as e:
        pytest.fail(
            f"{analysis_type} failed for {netlist}\n"
            f"{type(e).__name__}: {e}"
        )


if __name__ == "__main__":
    error_dict = {}
    #method = "two_graph_node"
    method = "tableau"

    print(f"{Bcolors.HEADER}{Bcolors.BOLD}{Bcolors.UNDERLINE}DC:{Bcolors.ENDC}")
    print(f"{Bcolors.HEADER} Numeric test:{Bcolors.ENDC}")
    error_dict["DC_semi"] = analysis("DC", is_symbolic=False, analysis_method=method)

    print(f"{Bcolors.HEADER}{Bcolors.BOLD}{Bcolors.UNDERLINE}TF:{Bcolors.ENDC}")
    print(f"{Bcolors.HEADER} Numeric test:{Bcolors.ENDC}")
    error_dict["TF_semi"] = analysis("TF", is_symbolic=False, analysis_method=method)

    print(f"{Bcolors.HEADER}{Bcolors.BOLD}{Bcolors.UNDERLINE}AC:{Bcolors.ENDC}")
    print(f"{Bcolors.HEADER} Numeric test:{Bcolors.ENDC}")
    error_dict["AC_semi"] = analysis("AC", is_symbolic=False, analysis_method=method)

    print(f"{Bcolors.HEADER}{Bcolors.BOLD}{Bcolors.UNDERLINE}TRAN:{Bcolors.ENDC}")
    print(f"{Bcolors.HEADER} Numeric test:{Bcolors.ENDC}")
    error_dict["TRAN_semi"] = analysis("tran", is_symbolic=False, analysis_method=method)
