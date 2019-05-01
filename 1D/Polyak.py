'''
                                                    Packages
'''

from numpy import array
from numpy.linalg import norm
from scipy.optimize import fsolve




'''
                                                    Polyak Method
'''


def Polyak(h, u0, v0, condStop, E, gradE, E_min, eps=1e-8):

    u, v = [], []
    u.append(u0)
    v.append(v0)
    n, ok = 0, 0

    if condStop == 'gradFct':
        while ok == 0:

            u.append(u[-1] + h * v[-1])
            v.append((1 - 2 * h) * v[-1] - h * gradE(u[-1]))

            if norm(gradE(u[-1])) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'objFct':
        while ok == 0:

            u.append(u[-1] + h * v[-1])
            v.append((1 - 2 * h) * v[-1] - h * gradE(u[-1]))

            if norm(E(u[-1]-E(u[-2]))) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'itFct':
        while ok == 0:

            u.append(u[-1] + h * v[-1])
            v.append((1 - 2 * h) * v[-1] - h * gradE(u[-1]))

            if norm(u[-1]-u[-2]) <= eps:
                ok = 1
            n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = E(u)
    H = 0.5 * abs(v) ** 2 + E - E_min

    Change = E[1:] - E[0:-1]
    Dissip = H[1:] - H[0:-1]

    return it, u, v, E, H, Change, Dissip

