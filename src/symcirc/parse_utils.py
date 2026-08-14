import traceback
from symcirc.utils import *
from sympy.parsing.sympy_parser import standard_transformations, convert_xor

TRANSFORMS = (standard_transformations + (convert_xor,))

NUMS = ["1", "2", "3", "4", "5", "6", "7", "8", "9", "0"]
UNITS = {"f": sympy.Rational(1, 1000000000000000), "F": sympy.Rational(1, 1000000000000000),
         "p": sympy.Rational(1, 1000000000000), "P": sympy.Rational(1, 1000000000000),
         "n": sympy.Rational(1, 1000000000), "N": sympy.Rational(1, 1000000000),
         "u": sympy.Rational(1, 1000000), "U": sympy.Rational(1, 1000000),
         "m": sympy.Rational(1, 1000), "M": sympy.Rational(1, 1000),
         "k": sympy.Rational(1000), "K": sympy.Rational(1000),
         "meg": sympy.Rational(1000000), "MEG": sympy.Rational(1000000), "Meg": sympy.Rational(1000000),
         "g": sympy.Rational(1000000000), "G": sympy.Rational(1000000000),
         "t": sympy.Rational(1000000000000), "T": sympy.Rational(1000000000000),
         "v": sympy.Rational(1), "V": sympy.Rational(1),
         "a": sympy.Rational(1), "A": sympy.Rational(1)}
OPERATORS = ["+", "-", "*", "/", "."]
RESERVED = ["sin"]

NETLIST_KEYCHARS = ["R", "r", "C", "c", "L", "l", "V", "v", "U", "u", "I", "i", "A", "a", "F", "f", "H", "h", "G", "g",
                    "E", "e", "K", "k", "S", "s", "X", "x", "Q", "q", "M", "m", "D", "d", "J", "j", ".", "*"]


def dc_value(name, dc_sig):
    if dc_sig is None:
        symbolic = True
        dc_val = sympy.Symbol(name, real=True)
    else:
        dc_val, symbolic = convert_units(dc_sig)

    dc_rat, dc_flt, dc_sym = process_value(name, dc_val, symbolic)
    return dc_rat, dc_flt, dc_sym

def ac_value(name, ac_sig):
    if ac_sig is not None:
        if len(ac_sig) not in [1, 2]:
            raise ValueError(f"{ac_sig} is not a correct AC signature. Use the following format: ['<mag>', '<phase>']")
        mag, phase = ac_sig
        ac_val, val_symbolic = convert_units(mag)
        if phase is not None:
            phase_rad, _ = convert_units(phase)
            #phase_rad = sympy.rad(phase_deg)
        else:
            phase_rad = Integer(0)
    else:
        val_symbolic = False
        ac_val = Integer(0)
        phase_rad = Integer(0)

    ac_rat, ac_flt, ac_sym = process_value(name, ac_val, val_symbolic)
    if ac_rat == 0:  # GEEC requires this
        ac_sym = ac_rat
    return ac_rat, ac_flt, ac_sym, phase_rad

def tran_value(name, tran_sig):
    if tran_sig is not None:
        if tran_sig[0].lower() == "sin":
            return parse_sin(tran_sig[1:])

    tran_rat, tran_flt = [None], [None]
    return tran_rat, tran_flt


def parse_sin(tran_sig):
    offset = Integer(0)
    amp = Integer(1)
    freq = Integer(1)
    delay = Integer(0)
    damping = Integer(1)
    #omega = Integer(2) * pi * freq
    try:
        offset, _ = convert_units(tran_sig[0])
        amp, _ = convert_units(tran_sig[1])
        freq, _ = convert_units(tran_sig[2], forced_numeric=True)
        delay, _ = convert_units(tran_sig[3])
        damping, _ = convert_units(tran_sig[4])
    except IndexError:
        pass
    tran_params = [offset, amp, freq, delay, damping]
    tran_rat = ["sin"] + [nsimplify(i, rational=True) for i in tran_params]
    tran_flt = ["sin"] + [evalf(i) for i in tran_params]
    return tran_rat, tran_flt
    #tran = amp*((damping+s)*sympy.sin(delay)+omega*sympy.cos(delay))/(damping**2+2*damping*s+omega**2+s**2)

def check_if_symbolic(val):
    return not val.is_number

def convert_units(val: str, forced_numeric: bool=False, local_dict=None):
    ret = None
    if local_dict is None:
        local = {}
    else:
        local = local_dict
    local.update(sympy.abc._clash)
    local.update(global_dict)
    val = val.replace("{", "").replace("}", "")
    if len(val) > 3:
        if (val[-3:] in UNITS) and (val[-4].isnumeric()):
            ret = sympy.parse_expr(val[:-3], local_dict=local, transformations=TRANSFORMS) * UNITS[val[-3:]]
    if len(val) > 1:
        if ret is None:
            if (val[-1] in UNITS) and (val[-2].isnumeric()):
                ret = sympy.parse_expr(val[:-1], local_dict=local, transformations=TRANSFORMS) * UNITS[val[-1]]
    if ret is None:
        ret = sympy.parse_expr(val, local_dict=local, transformations=TRANSFORMS)
    if forced_numeric:
        symbolic = False
    else:
        symbolic = check_if_symbolic(ret)
    return ret, symbolic

def process_value(name, val, symbolic):
    if symbolic:
        sym = val
        rat = val
        flt = val
    else:
        sym = sympy.Symbol(name, real=True)
        try:
            rat = nsimplify(val, rational=True)
        except:
            traceback.print_exc()
            rat = val
        flt = evalf(val)
    return rat, flt, sym

def value_enum(name: str, value: str, local_dict=None):
    try:
        value, is_symbolic = convert_units(value)
    except IndexError:
        is_symbolic = True
        if local_dict is None:
            local_dict = {}
        value = sympy.parse_expr(name, local_dict=sympy.abc._clash|local_dict, transformations=TRANSFORMS)
    rat, flt, sym = process_value(name, value, is_symbolic)
    return rat, flt, sym, is_symbolic

def nodes_per_element(type):
    type = type.lower()
    if type in ["r", "l", "c", "v", "i", "f", "h", "d", "s"]:
        return 2
    if type in ["f", "h", "q", "j"]:
        return 3
    elif type in ["a", "e", "g", "m"]:
        return 4
    elif type in ["k"]:
        return 0