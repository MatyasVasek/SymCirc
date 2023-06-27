import sympy
from symcirc.utils import s

def poles(TF):
    _, d = TF.as_numer_denom()
    p = sympy.roots(d, s)
    return p

def zeros(TF):
    n, _ = TF.as_numer_denom()
    z = sympy.roots(n, s)
    return z

def pole_zero(TF):
    n, d = TF.as_numer_denom()
    p = sympy.roots(d, s)
    z = sympy.roots(n, s)
    return p, z
