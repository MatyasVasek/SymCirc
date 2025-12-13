import sympy
import time
from sympy import exp, sin, cos, sqrt, factorial, DiracDelta, abc
from symcirc.utils import s, t, pi


def poles(F):
    #print("poles F: {}".format(F))
    N, D = F.as_numer_denom()
    roots = sympy.roots(D, s)
    return roots


def residue(F, pole, order):
    #print(F)
    if order == 1:
        tmp = (s - pole) * F * sympy.exp(s*t)
        #print(tmp)
        func = sympy.limit(tmp, s, pole, "+")
        #print(func)
        #print("hit 1-order")
    elif order == 2:
        func = sympy.limit(sympy.simplify(sympy.diff(sympy.simplify((s - pole)**2 * F * sympy.exp(s*t)), s)), s, pole)
        #print("hit 2-order")
    else:
        #print("hit n-order")
        func = (1 / sympy.factorial(order - 1)) * sympy.limit(sympy.diff(((s - pole)**order) * F * sympy.exp(s*t), s, order-1), s, pole)
    return func


def residue_laplace(F):
    f = 0
    pole_dict = poles(F)
    #print(pole_dict)
    for p in pole_dict:
        order = pole_dict[p]
        f += residue(F, p, order)
    return f


def latex_print(data):
    print("{}".format(sympy.latex(data)))

def laplace(func):
    return sympy.laplace_transform(func, t, s, noconds=True)

def visualise(func):
    func_laplace = laplace(func)
    func_inv_laplace = iLT(func_laplace)
    print("Laplace: " + str(func) + " => " + str(func_laplace))
    print("Inverse Laplace: " + str(func_laplace) + " => " + str(func_inv_laplace))


def _split_by_numerator(f, list):
    N, D = f.as_numer_denom()
    N = N.expand()
    N = N.collect(s)
    #print("PARTFRAC: P, Q: {}, {}".format(N, D))
    if N.func == sympy.Add:
        for arg in N.args:
            #print(arg/D)
            list.append(arg / D)
    else:
        list.append(f)
    return list


def split_parts(f):
    part_list = []
    t0 = time.time()
    f = sympy.apart(f, s)
    t1 = time.time()
    #print(f)
    if f.func == sympy.Add:
        for arg in f.args:
            #arg = sympy.factor(arg)
            part_list = _split_by_numerator(arg, part_list)
    else:
        part_list = _split_by_numerator(f, part_list)

    return part_list

def separate_s(f):
    f = sympy.Poly(f, s)
    coeff = f.all_coeffs()
    return coeff

def iLT(F, sympy_ilt=True):
    if sympy_ilt:
        try:
            return sympy.integrals.transforms.inverse_laplace_transform(F, s, t, plane=None)
        except:
            return symcirc_iLT(F)
    else:
        return symcirc_iLT(F)


def symcirc_iLT(F):
    f = 0
    F = sympy.apart(F, s)
    part_list = split_parts(F)
    for Func in part_list:
        func = table_inverse_laplace_transform(Func)
        if func is None:
            func = residue_laplace(F)
        f += func
    return f


def table_inverse_laplace_transform(F):
    """
    This function performs the inverse laplace transform via comparison of factors and known table values, if possible.
    If this method fails, the ILT is computed by definition using sympy library.
    """
    f = None
    #print("F = {}".format(F))
    N, D = F.as_numer_denom()
    N = N.expand()
    N = N.collect(s)

    D = D.expand()
    D = D.collect(s)
    #print("PARTFRAC: N, D: {}, {}".format(N, D))
    N_coeff = separate_s(N)
    D_coeff = separate_s(D)
    #print("N_coeff: {}".format(N_coeff))
    #print("D_coeff: {}".format(D_coeff))
    d_size = len(D_coeff)
    n_size = len(N_coeff)

    #print("F:{}".format(F))
    if F == 0:
        f = 0
    elif s not in F.atoms():
        f = DiracDelta(t, 0)*F
    elif d_size == 1 and n_size == 1 and N_coeff[1] == 0:
        f = DiracDelta(t, 1)*N_coeff[0]
        #print(f)
    elif d_size == 3 and D_coeff[1] == 0 and D_coeff[2] != 0:  # second order polynomial of type: (s**2 + a)
        if n_size == 1:  # sin form:  c*a/(s**2+a**2) --> c*sin(a*t)
            a = sqrt(sympy.cancel(D_coeff[2]/D_coeff[0]))
            c = sympy.cancel(N_coeff[0]/(a*D_coeff[0]))
            f = c * sin(a * t)
        elif n_size == 2 and N_coeff[1] == 0:  # cos form:  c*s/(s**2+a**2) --> c*cos(a*t)
            #print("cos form")
            a = sqrt(sympy.cancel(D_coeff[2]/D_coeff[0]))
            c = sympy.cancel(N_coeff[0]/D_coeff[0])
            f = c * cos(a * t)
    else:
        non_zero_c = 0
        for i in D_coeff:
            if i != 0:
                non_zero_c += 1
        if non_zero_c == 1:  # t form:  n!/s**(n+1) --> t**n
            #print("t form")
            n = d_size-2
            c = N_coeff[0]/(factorial(n)*D_coeff[0])
            f = c*t**n
        if non_zero_c == d_size:  # exp form:  c*n!/(s+a)**(n+1) --> c*t**n*exp(-a*t)
            #print("EXP")
            try:
                check = D_coeff[1]/(2*sqrt(D_coeff[2]))
            except IndexError:
                check = 1
            if check in [1, -1]:
                #print(D_coeff)
                n = d_size-2
                c = N_coeff[0]/(factorial(n)*D_coeff[0])
                a = sympy.cancel(D_coeff[-1]/D_coeff[0])
                #print("a: {}".format(a))
                f = c*t**n*exp(-a*t)
            else:
                if n_size == 1: # exp*sin form: c*omega/((s+a)**2 + omega**2)  --> c*exp(-a*t)*sin(omega*t)
                    #print("exp*sin form")
                    a = D_coeff[1]/(D_coeff[0]*2)
                    omega = sqrt(D_coeff[2]//D_coeff[0] - a**2)
                    c = N_coeff[0]/(omega*D_coeff[0])
                    f = c*exp(-a*t)*sin(omega*t)
                elif n_size == 2 and N_coeff[1] == 0:
                    #print("exp*cos*sin form")
                    a = D_coeff[1] / (D_coeff[0] * 2)
                    omega = sqrt(D_coeff[2]/D_coeff[0] - a ** 2)
                    c_cos = N_coeff[0]/D_coeff[0]
                    c_sin = a*N_coeff[0]/(omega*D_coeff[0])
                    f = exp(-a*t)*(c_cos*cos(omega*t)+c_sin*sin(omega*t))
        else:
            pass
    #print("result: {}".format(f))
    return f


if __name__ == "__main__":
    F = sympy.parse_expr('100000000000/(s**3 + 10000*s**2 + 100000000000*s)')
    s, t = sympy.symbols('s, t')
    sympy.integrals.transforms.inverse_laplace_transform(F, s, t)
