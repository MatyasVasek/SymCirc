import sympy

t, s = sympy.symbols('t, s', real=True, positive=True)
R1, R2, C1, C2 = sympy.symbols('R1, R2, C1, C2', real=True, positive=True)

def laplace(func):
    return sympy.laplace_transform(func, t, s, noconds=True)


def inv_laplace(func):
    return sympy.inverse_laplace_transform(func, s, t)


def visualise(func):
    func_laplace = laplace(func)
    func_inv_laplace = inv_laplace(func_laplace)
    print("Laplace: " + str(func) + " => " + str(func_laplace))
    print("Inverse Laplace: " + str(func_laplace) + " => " + str(func_inv_laplace))


if __name__ == "__main__":
    a = sympy.symbols('a', real=True, positive=True)
    omega = sympy.Symbol('omega', real=True)

    print("--------------------")
    freq_func = (C1*R2*s)/(C1*C2*R1*R2*s**2 + C1*R1*s + C1*R2*s + C2*R2*s + 1)
    freq_func = 1 / (R1*s**2)
    print("F = "+str(freq_func))
    partial_fractions = sympy.apart(freq_func, s)
    print("F = "+str(partial_fractions))
    func = inv_laplace(partial_fractions)
    print("f = "+str(func))
