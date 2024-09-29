import sys, os
import sympy
import random
import json
# project root path
sys.path.append(os.path.dirname(__file__)+"/../src/")
import symcirc
from symcirc import utils


class ErrData:
    def __init__(self):
        self.count = 0
        self.data = []


class NumericTestErrors:
    def __init__(self):
        self.total = 0
        self.not_in_ref = ErrData()
        self.index_err = ErrData()
        self.type_err = ErrData()
        self.key_err = ErrData()
        self.incorrect = ErrData()

    def update(self, err_type, msg=""):
        self.total += 1
        if err_type == IndexError:
            self.index_err.count += 1
        if err_type == TypeError:
            self.type_err.count += 1
        if err_type == KeyError:
            self.key_err.count += 1
        if err_type == "incorrect":
            self.incorrect.count += 1
            self.incorrect.data.append(msg)
        if err_type == "not_in_ref":
            self.not_in_ref.count += 1
            self.not_in_ref.data.append(msg)


def numeric_test(filename, analysis_type, precision=6, test_data_creation=False, method="tableau", runs=10):
    err = NumericTestErrors()
    divisor = 10**precision
    circuit = symcirc.AnalyseCircuit(utils.load_file("netlists/{}".format(filename)), analysis_type, symbolic=False, method=method)
    result_dict = circuit.component_values()
    reference_string = utils.load_file("reference_data/{}/{}".format(analysis_type, filename))
    symbol_dict = {"t": utils.t, "f": utils.f, "j": utils.j, "s": utils.s}
    for item in circuit.symbols:
        symbol_dict.update({str(item): item})

    reference_dict = sympy.parse_expr(reference_string, symbol_dict)

    if test_data_creation:
        print("{}--------Analysis result for {}:--------".format('\33[37m', filename))
        print(result_dict)
        print("----------------Reference:---------------")
        print(reference_dict)
        print("-----------------------------------------")
    #print(result_dict)
    for res_key in result_dict:
        for i in range(runs):
            try:
                if result_dict[res_key] != None:
                    ratio = None
                    seed = random.random()
                    if reference_dict[res_key] in [0.0, 0]:
                        ratio = sympy.simplify(result_dict[res_key] - reference_dict[res_key])
                    elif analysis_type == "tran":
                        ratio = (sympy.simplify((result_dict[res_key])/(reference_dict[res_key])).subs(utils.t, seed))
                    elif analysis_type == "TF":
                        ratio = (sympy.simplify((result_dict[res_key])/(reference_dict[res_key])).subs(utils.s, seed))
                    elif analysis_type == "AC":
                        ratio = sympy.simplify((result_dict[res_key])/(reference_dict[res_key]))
                    elif analysis_type == "DC":
                        ratio = sympy.simplify((result_dict[res_key])/(reference_dict[res_key]))

                    if res_key not in reference_dict:
                        msg = "{} is not in reference file".format(res_key)
                        err.update("not_in_ref", msg)
                    elif result_dict[res_key] == reference_dict[res_key]:
                        pass
                    elif ratio > 0.98 or ratio < 1.02:
                        pass
                    else:
                        msg = "{} is {} in result and {} in reference".format(res_key, result_dict[res_key], reference_dict[res_key])
                        err.update("incorrect", msg)

            except IndexError:
                err.update(IndexError)
            except TypeError:
                msg = "{} is {} in result and {} in reference".format(res_key, result_dict[res_key],
                                                                      reference_dict[res_key])
                err.update("incorrect", msg)
            except KeyError:
                err.update(KeyError)

    return err


def print_test_results(filename, ref, analysis_type, err, full=False):
    if err.total != 0:
        print("{}[FAIL]: {} analysis of {} ended with {} errors"
              .format('\033[91m', analysis_type, filename, err.total))
        if full:
            print("     ----------ERROR ANALYSIS FOR {}----------".format(filename))
            print("     Total errors found: {}".format(err.total))
            print("     -----------------------------------------")
            print("     Numerical errors: {}".format(err.incorrect.count))
            for e in err.incorrect.data:
                print("     " + e)
            print("     -----------------------------------------")
            print("     Value not in reference file: {}".format(err.not_in_ref.count))
            for e in err.not_in_ref.data:
                print("     "+e)
            print("     -----------------------------------------")
            print("     KeyError: {}".format(err.key_err.count))
            print("     TypeError: {}".format(err.type_err.count))
            print("     IndexError: {}".format(err.index_err.count))
    else:
        print("{}[OK]: {} analysis of {} successful!"
              .format('\033[96m', analysis_type, filename))


def analysis_test(analysis="all", w=False, test_data_creation=False, full=True, method="tableau", runs=10):
    if analysis == "all":
        a_types = ["DC", "TF", "tran"]
    elif type(analysis) == list:
        a_types = analysis
    else:
        a_types = [analysis]

    for analysis_type in a_types:
        netlists = []
        ref = []
        for filename in os.listdir("{}/netlists".format(os.getcwd())):
            netlists.append(filename)
        for f in os.listdir("{}/reference_data/{}".format(os.getcwd(), analysis_type)):
            ref.append(f)

        for filename in netlists:
            try:
                if filename in ref:
                    err = numeric_test(filename, analysis_type, test_data_creation=test_data_creation, precision=6, method=method, runs=10)
                    print_test_results(filename, ref, analysis_type, err, full=full)
                elif w:
                    print("{}[Warning]: Reference for {} analysis of {} not found".format('\033[93m', analysis_type, filename))
            except TypeError:
                print("{}[FAIL]: {} analysis of {} ended with an error".format('\033[91m', analysis_type, filename))
    return True


if __name__ == "__main__":
    state = analysis_test(analysis="all", w=False, test_data_creation=False, full=True, method="tableau", runs=10)
    if state:
        print("")
        print("Test successful!")

