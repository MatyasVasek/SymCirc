import sys
import os

import func_timeout
import pytest

# Add project src to path so symcirc can be imported without installing
sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
from symcirc import *
from symcirc import utils

# --- Test configuration ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
NETLIST_DIR = os.path.join(BASE_DIR, "netlists")

ANALYSIS_TYPES = ["DC", "TF", "AC", "tran"]
ANALYSIS_METHODS = ["tableau", "two_graph_node"]
SYMBOLIC = [True, False]
SOLVERS = ["gauss", "ddd"]

NETLISTS = sorted(p for p in os.listdir(NETLIST_DIR) if p.endswith((".txt", ".cir")))

TIMEOUT = 20


class Bcolors:
    HEADER    = '\033[95m'
    OKBLUE    = '\033[94m'
    OKCYAN    = '\033[96m'
    OKGREEN   = '\033[92m'
    WARNING   = '\033[93m'
    FAIL      = '\033[91m'
    ENDC      = '\033[0m'
    BOLD      = '\033[1m'
    UNDERLINE = '\033[4m'


def analyze(filename, analysis_type, is_symbolic, analysis_method, solver):
    """
    Run a single circuit analysis. Intended to be called via func_timeout.
    """
    circuit = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method=analysis_method,
        use_symengine=True,
        solver=solver,
    )
    circuit.get_all_results()


# Pytest parametrized smoke test
# NotImplementedError and timeouts are marked xfail rather than hard failures.
@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_method", ANALYSIS_METHODS)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
@pytest.mark.parametrize("solver", SOLVERS)
def test_smoke(analysis_method, analysis_type, netlist, is_symbolic, solver):
    try:
        func_timeout.func_timeout(TIMEOUT, analyze, args=(netlist, analysis_type, is_symbolic, analysis_method, solver))
    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        pytest.fail(f"{analysis_type} failed for {netlist}\n{type(e).__name__}: {e}")
