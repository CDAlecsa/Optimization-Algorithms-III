'''
                                                    Packages
'''

from numpy import array
from numpy.linalg import norm
from scipy.optimize import fsolve


'''
                                                    Newton Solver for the Implicit Method
'''


def F_Newton(x, u_n, v_n, h, c, coeff_1, coeff_2):
    sol1 = x[0] - u_n[0] - v_n[0] * (h/c) + (h ** 2 / 2) * coeff_1 * (x[0] + u_n[0])
    sol2 = x[1] - u_n[1] - v_n[1] * (h/c) + (h ** 2 / 2) * coeff_2 * (x[1] + u_n[1])
    return [sol1, sol2]


'''
                                                    Strang Splitting
'''


def Strang_2D(h, u0, v0, E, gradE, coeff_1, coeff_2, E_min, eps=1e-8):

    c = (1+h/2)/(1-h/2)

    u, v = [], []

    u.append(u0)
    v.append(v0)

    n, ok = 0, 0

    while ok == 0:

        guess = u[-1]
        args = (u[-1], v[-1], h, c, coeff_1, coeff_2)
        u.append(fsolve(F_Newton, guess, args, xtol=1e-10))

        C1 = [(1/c**2) * i for i in v[-1]]
        C2 = [i + j for i, j in zip(u[-1], u[-2])]
        C3 = [(h/c) * i for i in C2]
        C4 = [i * j for i, j in zip(C3, [coeff_1, coeff_2])]

        res = [i - j for i, j in zip(C1, C4)]
        v.append(res)

        if norm(gradE(u[-1], coeff_1, coeff_2)) <= eps:
            ok = 1
        n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = [E(i, coeff_1, coeff_2) for i in u]
    H = [0.5 * norm(i) ** 2 + j - E_min for i, j in zip(v, E)]

    return it, u, v, E, H

