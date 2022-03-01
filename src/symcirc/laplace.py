import sympy
import time
from sympy import exp, sin, cos, sqrt, factorial

def latex_print(data):
    print("{}".format(sympy.latex(data)))

t, s, a = sympy.symbols('t, s, a', positive=True)
R, R2, C, C2 = sympy.symbols('R, R2, C, C2', positive=True)

def laplace(func):
    return sympy.laplace_transform(func, t, s, noconds=True)


def inv_laplace(func):
    return sympy.inverse_laplace_transform(func, s, t)


def visualise(func):
    func_laplace = laplace(func)
    func_inv_laplace = inv_laplace(func_laplace)
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


def iLT(F):
    f = 0
    #print("F = {}".format(F))
    #F = sympy.apart(F, s)
    part_list = split_parts(F)
    #print(part_list)
    for Func in part_list:
        func = table_inverse_laplace_transform(Func)
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


    if F == 0:
        f = 0
    elif s not in F.atoms():
        f = F
    elif d_size == 3 and D_coeff[1] == 0 and D_coeff[2] != 0:  # second order polynomial of type: (s**2 + a)
        if n_size == 1:  # sin form:  c*a/(s**2+a**2) --> c*sin(a*t)
            #print("sin form")
            a = sqrt(sympy.cancel(D_coeff[2]/D_coeff[0]))
            #print("a: {}".format(a))
            c = sympy.cancel(N_coeff[0]/(a*D_coeff[0]))
            f = c * sin(a * t)
            #print("f : {}".format(f))
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
            c = N_coeff[0]/factorial(n)
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




    #print(f)
    return f


if __name__ == "__main__":
    L, R1, R2, C1= sympy.symbols('L, R1, R2, C1', real=True, positive=True)
    omega = sympy.Symbol('omega', real=True)
    print("--------------------")
    pi = sympy.symbols("pi", real=True, positive=True)
    a = sympy.Rational(1, 5)
    #F = (s+a)/(s*(s+a) + 1)
    #F = -2*pi*C1*R1*s/((s - 1)**2 + 4*pi**2)
    #F = (-20*pi*(s-7))/((s - 7)**2 + 4*pi**2)
    #F = 6 / (s - 32000*pi)**2
    #F = (36000*pi*C*R1*R2*s/(s**2 + 1296000000*pi**2) + 36000*pi*R1/(s**2 + 1296000000*pi**2))/(C*R1*R2*s + R1 + R2)
    F = (3*s)/(7*s**2+11*s+13)
    time0 = time.time()
    f = sympy.inverse_laplace_transform(F, s, t)
    time1 = time.time()
    print(time1 - time0)
    print("f = " + str(f))
    f = iLT(F)
    time2 = time.time()
    print(time2 - time1)
    simp_f = sympy.simplify(f)
    print("f = "+str(simp_f))
    #print(sympy.simplify(inv_laplace(3*s/(5*s**2 + s + 5))))
