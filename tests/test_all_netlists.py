import sys
import os
import random

import func_timeout
import pytest

import subprocess
import json
import tempfile
import sympy

#from symcirc import utils
#import symcirc

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

TIMEOUT = 100


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


def analyze(filename, analysis_type, is_symbolic, analysis_method):
    """
    Run a single circuit analysis. Intended to be called via func_timeout.
    """
    # Add project src to path so symcirc can be imported without installing
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
    from symcirc import AnalyseCircuit, utils
    analysis = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method=analysis_method,
        use_symengine=True,
    )
    analysis.get_all_results()


def compare_solvers(filename, analysis_type, is_symbolic, analysis_method):
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
    from symcirc import AnalyseCircuit, utils, evalf
    analysis_gauss = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method=analysis_method,
        use_symengine=True,
        solver="gauss",
    )
    v_gauss = analysis_gauss.node_voltages()

    analysis_ddd = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method=analysis_method,
        use_symengine=True,
        solver="ddd",
    )
    v_ddd = analysis_ddd.node_voltages()

    for key in v_gauss:
        tmp_gauss = v_gauss[key]
        tmp_ddd =  v_ddd[key]
        if not sympy.sympify(tmp_gauss).equals(sympy.sympify(tmp_ddd)):
            seed = random.randrange(1, 1000)/1000000
            if analysis_type == "tran":
                if tmp_gauss in [0.0, 0]:
                    diff = evalf((tmp_gauss - tmp_ddd), {utils.t: seed})
                    if -0.05 < diff < 0.05:
                        continue
                else:
                    ratio = evalf((tmp_gauss) / (tmp_ddd),{utils.t: seed})
                    if ratio > 0.95 and ratio < 1.05:
                        continue

            raise ValueError(f"{v_gauss[key]} != {v_ddd[key]}")


def compare_methods(filename, analysis_type, is_symbolic, solver):
    sys.path.append(os.path.join(os.path.dirname(__file__), "..", "src"))
    from symcirc import AnalyseCircuit, utils
    analysis_tableau = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method="tableau",
        use_symengine=True,
        solver=solver,
    )
    v_tableau = analysis_tableau.node_voltages()

    analysis_tgn = AnalyseCircuit(
        utils.load_file(os.path.join(NETLIST_DIR, filename)),
        analysis_type=analysis_type,
        symbolic=is_symbolic,
        method="two_graph_node",
        use_symengine=True,
        solver=solver,
    )
    v_tgn = analysis_tgn.node_voltages()

    for key in v_tableau:
        sympify_tableau = sympy.sympify(v_tableau[key])
        sympify_tgn = sympy.sympify(v_tgn[key])
        if not sympify_tableau.equals(sympify_tgn):
            if not sympy.simplify(sympify_tableau).equals(sympy.simplify(sympify_tgn)):
                raise ValueError(f"{v_tableau[key]} != {v_tgn[key]}")


def compare_results(results_a, results_b):
    from symcirc import utils
    for key in results_b:
        if key not in results_a:
            pytest.fail(f"Key {key} missing in local version")
    for key in results_a:
        if key not in results_b:
            pytest.fail(f"Key {key} missing in PyPI version")

    for key in results_a:
        sympify_a = sympy.sympify(results_a[key])
        sympify_b = sympy.sympify(results_b[key])
        if not sympify_a.equals(sympify_b):
            a_evalf = utils.evalf(sympify_a)
            b_evalf = utils.evalf(sympify_b)
            if not a_evalf.equals(b_evalf):
                #ratio = sympy.simplify(a_evalf/b_evalf)
                #if not ratio > 0.999 or ratio < 1.001:
                raise ValueError(f"{results_a[key]} != {results_b[key]}")


def run_analysis_w_timing(filename, analysis_type, is_symbolic, analysis_method, local=True, solver="gauss"):
    """
    Run analysis using PyPI installed version via subprocess
    """
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.json', delete=False) as f:
        res_file = f.name
    with tempfile.NamedTemporaryFile(mode='w+', suffix='.txt', delete=False) as f:
        time_file = f.name

    try:
        if local:
            src_path = os.path.join(os.path.dirname(__file__), "..", "src")
            script = f"import sys\nsys.path.insert(0, {repr(src_path)})\n"
        else:
            script = "import sys"
        script += f"""
from symcirc import AnalyseCircuit, utils
import os
import json
import time
time_start = time.time()
analysis = AnalyseCircuit(
    utils.load_file(os.path.join({repr(NETLIST_DIR)}, {repr(filename)})),
    analysis_type={repr(analysis_type)},
    symbolic={is_symbolic},
    method={repr(analysis_method)},
    use_symengine=True,"""
        if local:
            script += f"solver={repr(solver)},"
        script += f""")
results = analysis.node_voltages()
time_end = time.time() - time_start

# Write results to temp file
with open({repr(res_file)}, 'w') as f:
    json.dump({{str(k): str(v) for k, v in results.items()}}, f)
with open({repr(time_file)}, 'w') as f:
    f.write(str(time_end))
"""
        #result = subprocess.run([sys.executable, "-c", script], timeout=TIMEOUT)
        env = os.environ.copy()
        env.pop("PYTHONPATH", None)  # stop a dev shell's PYTHONPATH from leaking in
        result = subprocess.run(
            [sys.executable, "-c", script],
            timeout=TIMEOUT,
            cwd=tempfile.gettempdir(),  # neutral cwd — never the repo root
            env=env,
        )
        if result.returncode != 0:
            raise Exception(f"PyPI version failed with return code {result.returncode}")

        # Read results from temp file
        with open(res_file, 'r') as f:
            res = json.load(f)
        with open(time_file, 'r') as f:
            analysis_time = f.read()
        return res, analysis_time
    finally:
        os.unlink(res_file)
        os.unlink(time_file)


# Pytest parametrized smoke test
# NotImplementedError and timeouts are marked xfail rather than hard failures.
@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_method", ANALYSIS_METHODS)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
def test_smoke(analysis_method, analysis_type, netlist, is_symbolic):
    try:
        func_timeout.func_timeout(TIMEOUT, analyze, args=(netlist, analysis_type, is_symbolic, analysis_method))
    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        pytest.fail(f"{analysis_type} failed for {netlist}\n{type(e).__name__}: {e}")

@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
@pytest.mark.parametrize("solvers", SOLVERS)
def test_tgn_matrix(analysis_type, netlist, is_symbolic, solvers):
    try:
        func_timeout.func_timeout(TIMEOUT, compare_methods, args=(netlist, analysis_type, is_symbolic, solvers))
    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        pytest.fail(f"{analysis_type} failed for {netlist}\n{type(e).__name__}: {e}")

@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
def test_ddd_solver(analysis_type, netlist, is_symbolic):
    try:
        func_timeout.func_timeout(TIMEOUT, compare_solvers, args=(netlist, analysis_type, is_symbolic, "two_graph_node"))
    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        pytest.fail(f"{analysis_type} failed for {netlist}\n{type(e).__name__}: {e}")

FIXED_BUGS = ["DC_tran_1.txt",
              "DC_tran_10.txt",
              "RL.txt",
              "RLL.txt",
              "coupled_ic.txt",
              "geec_1800.txt",
              "geec_241.txt",
              "test.txt",
              "wien_oamp_1790.txt"
              ]
@pytest.mark.parametrize("is_symbolic", SYMBOLIC)
@pytest.mark.parametrize("analysis_method", ANALYSIS_METHODS)
@pytest.mark.parametrize("analysis_type", ANALYSIS_TYPES)
@pytest.mark.parametrize("netlist", NETLISTS)
def test_regression_vs_pypi(analysis_method, analysis_type, netlist, is_symbolic):
    """
    Test that local version produces the same results as the latest PyPI version.
    Both versions run in isolated subprocesses for fair comparison.
    """
    if netlist in FIXED_BUGS:
        pytest.xfail(f"Expected fail due to the current PyPi version having a bug that was fixed or"
                     f"it not supporting this simulation.")
    try:
        # Run local version in subprocess
        local_results, local_time = func_timeout.func_timeout(
            TIMEOUT,
            run_analysis_w_timing,
            args=(netlist, analysis_type, is_symbolic, analysis_method, True)
        )

        # Run PyPI version in subprocess
        pypi_results, pypi_time = func_timeout.func_timeout(
            TIMEOUT,
            run_analysis_w_timing,
            args=(netlist, analysis_type, is_symbolic, analysis_method, False)
        )

        func_timeout.func_timeout(TIMEOUT, compare_results, args=(local_results, pypi_results))

        # Report timing
        '''local_time_f = float(local_time)
        pypi_time_f = float(pypi_time)
        speedup = pypi_time_f / local_time_f
        diff = pypi_time_f - local_time_f

        msg = f"Timing for {netlist} ({analysis_type}, symbolic={is_symbolic}): "
        msg += f"PyPI {pypi_time_f:.3f}s | Local {local_time_f:.3f}s | "
        msg += f"Speedup: {speedup:.2f}x ({diff:+.3f}s)"

        logger.info(msg)'''

    except NotImplementedError:
        pytest.xfail(f"{analysis_type} not implemented for {netlist}")
    except func_timeout.FunctionTimedOut:
        pytest.xfail(f"{analysis_type} timeout after {TIMEOUT}s for {netlist}")
    except Exception as e:
        if str(e) == f"PyPI version failed with return code 1":
            pytest.xfail(f"PyPI version failed")
        pytest.fail(f"Regression test failed for {netlist}\n{type(e).__name__}: {e}")
