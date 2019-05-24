'''
                                                    Packages
'''

from numpy import array
from numpy.linalg import norm
from scipy.optimize import fsolve


'''
                                                    Newton Solver for the Implicit Method
'''


def F_Newton(x, u_n, v_n, h, c):

    A = (1 / 3) * (300 * u_n[0] ** 2 - 200 * u_n[1] - 400 * x[1] + 3) * x[0]
    B = (1 / 3) * (-400 * u_n[1] - 200 * x[1] + 3) * u_n[0]

    auxNewton1 = 100 * x[0] ** 3 + 100 * u_n[0] * x[0] ** 2 - 2 + 100 * u_n[0] ** 3 + A + B
    auxNewton2 = (-200/3) * u_n[0] ** 2 + (-200 / 3) * u_n[0] * x[0] + (-200/3) * x[0] ** 2 + 100 * u_n[1] + 100 * x[1]

    sol1 = x[0] - u_n[0] - v_n[0] * (h / c) + (h ** 2 / 2) * auxNewton1
    sol2 = x[1] - u_n[1] - v_n[1] * (h / c) + (h ** 2 / 2) * auxNewton2

    return array([sol1, sol2])



'''
                                                    Strang Splitting
'''


def Strang_2D(h, u0, v0, E, gradE, E_min, eps=1e-8):

    c = (1+h/2)/(1-h/2)

    u, v = [], []

    u.append(u0)
    v.append(v0)

    n, ok = 0, 0

    while ok == 0:

        guess = u[-1]
        args = (u[-1], v[-1], h, c)
        u.append(fsolve(F_Newton, guess, args, xtol=1e-10))

        A = (1 / 3) * (300 * u[-2][0] ** 2 - 200 * u[-2][1] - 400 * u[-1][1] + 3) * u[-1][0]
        B = (1 / 3) * (-400 * u[-2][1] - 200 * u[-1][1] + 3) * u[-2][0]

        aux1 = 100 * u[-1][0] ** 3 + 100 * u[-2][0] * u[-1][0] ** 2 - 2 + 100 * u[-2][0] ** 3 + A + B
        aux2 = (-200/3) * u[-2][0] ** 2 + (-200/3) * u[-2][0] * u[-1][0] + \
               (-200/3) * u[-1][0] ** 2 + 100 * u[-2][1] + 100 * u[-1][1]

        v.append([(1 / c ** 2) * v[-1][0] - (h / c) * aux1, (1 / c ** 2) * v[-1][1] - (h / c) * aux2])

        if norm(gradE(u[-1])) <= eps:
            ok = 1
        n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = [E(i) for i in u]
    H = [0.5 * norm(i) ** 2 + j - E_min for i, j in zip(v, E)]

    Change = [(i-j)/h for i, j in zip(E[1:], E[0:-1])]
    Dissip = [(i-j)/h for i, j in zip(H[1:], H[0:-1])]

    return it, u, v, E, H, Change, Dissip
