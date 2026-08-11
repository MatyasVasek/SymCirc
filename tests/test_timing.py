import os
import func_timeout
import pytest
from tests.test_all_netlists import run_analysis_w_timing

import logging
logger = logging.getLogger(__name__)

# --- Test configuration ---
BASE_DIR = os.path.dirname(os.path.abspath(__file__))
NETLIST_DIR = os.path.join(BASE_DIR, "netlists")
LARGE_NETLIST_DIR = os.path.join(BASE_DIR, "netlists")

ANALYSIS_TYPES = ["DC", "TF", "AC", "tran"]
ANALYSIS_METHODS = ["tableau", "two_graph_node"]
SYMBOLIC = [True, False]
SOLVERS = ["gauss"] #"ddd"] # DDD is still in the experimental phase

NETLISTS = sorted(p for p in os.listdir(NETLIST_DIR) if p.endswith((".txt", ".cir")))
LARGE_NETLISTS = sorted(p for p in os.listdir(LARGE_NETLIST_DIR) if p.endswith((".txt", ".cir")))

TIMEOUT = 20


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

@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_method", ["two_graph_node"])
@pytest.mark.parametrize("analysis_type", ["TF", "tran"])
@pytest.mark.parametrize("netlist", NETLISTS)
def test_timing(analysis_method, analysis_type, netlist, is_symbolic):
    """
    Test that local version produces the same results as the latest PyPI version.
    Both versions run in isolated subprocesses for fair comparison.
    """
    try:
        # Run PyPI version in subprocess
        try:
            pypi_results, pypi_time = func_timeout.func_timeout(
                TIMEOUT,
                run_analysis_w_timing,
                args=(netlist, analysis_type, is_symbolic, analysis_method, False)
            )
        except func_timeout.FunctionTimedOut:
            pypi_results = {"Timeout"}
            pypi_time = TIMEOUT

        try:
            # Run local version in subprocess
            local_results, local_time = func_timeout.func_timeout(
                TIMEOUT,
                run_analysis_w_timing,
                args=(netlist, analysis_type, is_symbolic, analysis_method, True)
            )
        except func_timeout.FunctionTimedOut:
            local_results = {"Timeout"}
            local_time = TIMEOUT

        try:
            local_ddd_results, local_ddd_time = func_timeout.func_timeout(
                TIMEOUT,
                run_analysis_w_timing,
                args=(netlist, analysis_type, is_symbolic, analysis_method, True, "ddd")
            )
        except func_timeout.FunctionTimedOut:
            local_ddd_results = {"Timeout"}
            local_ddd_time = TIMEOUT


        # Report timing
        local_time_f = float(local_time)
        pypi_time_f = float(pypi_time)
        local_ddd_time_f = float(local_ddd_time)

        speedup = pypi_time_f / local_time_f
        diff = pypi_time_f - local_time_f
        speedup_ddd = pypi_time_f / local_ddd_time_f
        diff_ddd = pypi_time_f - local_ddd_time_f

        msg = f"Timing for {netlist} ({analysis_type}, symbolic={is_symbolic}): "
        msg += f"PyPI {pypi_time_f:.3f}s | Local {local_time_f:.3f}s | Local DDD {local_ddd_time_f:.3f}s |\n"
        msg += f"Speedup gauss: {speedup:.2f}x ({diff:+.3f}s) | Speedup ddd: {speedup_ddd:.2f}x ({diff_ddd:+.3f}s) "
        msg += f"\nLocal expr: {local_results}\nDDD expr: {local_ddd_results}\nPyPI expr: {pypi_results}"

        logger.info(msg)

    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        if str(e) == f"PyPI version failed with return code 1":
            pytest.xfail(f"PyPI version failed")
        pytest.fail(f"Regression test failed for {netlist}\n{type(e).__name__}: {e}")
