'''
                                                    Packages
'''

from numpy import array
from numpy.linalg import norm
from scipy.optimize import fsolve


'''
                                                    Newton Solver for the Implicit Method
'''

'''
                                                    Newton Solver for Crank-Nicolson Method
'''


def F_NewtonCN(x, u_n, v_n, h, c, gradE):

    x1 = x[0]
    x2 = x[1]

    sol1 = x1 - u_n[0] - v_n[0] * ((h / 2) * (1 + c)) + (h ** 2 / (4 * (1 + h))) * (
                gradE(u_n)[0] + gradE(x)[0])

    sol2 = x2 - u_n[1] - v_n[1] * ((h / 2) * (1 + c)) + (h ** 2 / (4 * (1 + h))) * (
                gradE(u_n)[1] + gradE(x)[1])

    return array([sol1, sol2])



'''
                                                    Strang Splitting
'''


def CN_2D(h, u0, v0, E, gradE, E_min, eps=1e-8):

    c = (1-h)/(1+h)

    u, v = [], []

    u.append(u0)
    v.append(v0)

    n, ok = 0, 0

    while ok == 0:

        guess = u[-1]
        args = (u[-1], v[-1], h, c, gradE)
        u.append(fsolve(F_NewtonCN, guess, args, xtol=1e-10))

        v.append([c * v[-1][0] - (h/(2*(1+h))) * (gradE(u[-1])[0] + gradE(u[-2])[0]),
                  c * v[-1][1] - (h/(2*(1+h))) * (gradE(u[-1])[1] + gradE(u[-2])[1])])

        if norm(gradE(u[-1])) <= eps:
            ok = 1
        n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = [E(i) for i in u]
    H = [0.5 * norm(i) ** 2 + j - E_min for i, j in zip(v, E)]

    Change = [(i - j)/h for i, j in zip(E[1:], E[0:-1])]
    Dissip = [(i - j)/h for i, j in zip(H[1:], H[0:-1])]


    return it, u, v, E, H, Change, Dissip
