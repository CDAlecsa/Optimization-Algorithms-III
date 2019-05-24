'''
                                                    Packages
'''

from numpy import array
from numpy.linalg import norm



'''
                                                    Strang Splitting - Explicit Method with Predictor-Corrector
'''


def Strang_PC(h, u0, v0, condStop, E, gradE, E_min, eps=1e-8):

    c = (1+h/2)/(1-h/2)
    u, v = [], []
    u.append(u0)
    v.append(v0)
    n, ok = 0, 0

    if condStop == 'gradFct':
        while ok == 0:

            u.append(u[-1] + (h / c) * v[-1] - (h ** 2 / 2) * gradE(u[-1] + (h / c) * v[-1]))

            A1 = (h ** 3 * (1 - h / 2) ** 3 * v[-1] ** 3) / (4 * (1 + h / 2) ** 3)
            A2 = (u[-2] * h ** 2 * (1 - h / 2) ** 2 * v[-1] ** 2) / ((1 + h / 2) ** 2)
            A3 = (3 * u[-2] ** 2 * h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A4 = (-h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A5 = u[-2] ** 3 - u[-2]

            v.append((1 / c ** 2) * v[-1] - (h / c) * (A1 + A2 + A3 + A4 + A5))

            if norm(gradE(u[-1])) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'objFct':
        while ok == 0:

            u.append(u[-1] + (h / c) * v[-1] - (h ** 2 / 2) * gradE(u[-1] + (h / c) * v[-1]))

            A1 = (h ** 3 * (1 - h / 2) ** 3 * v[-1] ** 3) / (4 * (1 + h / 2) ** 3)
            A2 = (u[-2] * h ** 2 * (1 - h / 2) ** 2 * v[-1] ** 2) / ((1 + h / 2) ** 2)
            A3 = (3 * u[-2] ** 2 * h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A4 = (-h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A5 = u[-2] ** 3 - u[-2]

            v.append((1 / c ** 2) * v[-1] - (h / c) * (A1 + A2 + A3 + A4 + A5))

            if norm(E(u[-1]-E(u[-2]))) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'itFct':
        while ok == 0:

            u.append(u[-1] + (h / c) * v[-1] - (h ** 2 / 2) * gradE(u[-1] + (h / c) * v[-1]))

            A1 = (h ** 3 * (1 - h / 2) ** 3 * v[-1] ** 3) / (4 * (1 + h / 2) ** 3)
            A2 = (u[-2] * h ** 2 * (1 - h / 2) ** 2 * v[-1] ** 2) / ((1 + h / 2) ** 2)
            A3 = (3 * u[-2] ** 2 * h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A4 = (-h * (1 - h / 2) * v[-1]) / (2 * (1 + h / 2))
            A5 = u[-2] ** 3 - u[-2]

            v.append((1 / c ** 2) * v[-1] - (h / c) * (A1 + A2 + A3 + A4 + A5))

            if norm(u[-1]-u[-2]) <= eps:
                ok = 1
            n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = E(u)
    H = 0.5 * abs(v) ** 2 + E - E_min

    Change = E[1:] - E[0:-1]
    Change = [i/h for i in Change]
    
    Dissip = H[1:] - H[0:-1]
    Dissip = [i/h for i in Dissip]

    return it, u, v, E, H, Change, Dissip


