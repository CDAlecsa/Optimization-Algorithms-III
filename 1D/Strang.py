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
    return x - u_n - v_n * (h/c) + (h ** 2 / 8) * (x + u_n) * (x ** 2 + u_n ** 2 - 2)


'''
                                                    Strang Splitting
'''


def Strang(h, u0, v0, condStop, E, gradE, E_min, eps=1e-8):

    c = (1+h/2)/(1-h/2)
    u, v = [], []
    u.append(u0)
    v.append(v0)
    n, ok = 0, 0

    if condStop == 'gradFct':
        while ok == 0:

            guess = u[-1]
            args = (u[-1], v[-1], h, c)
            u.append(fsolve(F_Newton, guess, args, xtol=1e-10))
            v.append((1/c**2) * v[-1] - (h/(4*c)) * (u[-1]+u[-2]) * (u[-1]**2 + u[-2]**2 - 2))

            if norm(gradE(u[-1])) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'objFct':
        while ok == 0:

            guess = u[-1]
            args = (u[-1], v[-1], h, c)
            u.append(fsolve(F_Newton, guess, args, xtol=1e-10))
            v.append((1/c**2) * v[-1] - (h/(4*c)) * (u[-1]+u[-2]) * (u[-1]**2 + u[-2]**2 - 2))

            if norm(E(u[-1]-E(u[-2]))) <= eps:
                ok = 1
            n = n + 1

    elif condStop == 'itFct':
        while ok == 0:

            guess = u[-1]
            args = (u[-1], v[-1], h, c)
            u.append(fsolve(F_Newton, guess, args, xtol=1e-10))
            v.append((1/c**2) * v[-1] - (h/(4*c)) * (u[-1]+u[-2]) * (u[-1]**2 + u[-2]**2 - 2))

            if norm(u[-1]-u[-2]) <= eps:
                ok = 1
            n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = E(u)
    H = 0.5 * abs(v)**2 + E - E_min

    Change = E[1:] - E[0:-1]
    Change = [i/h for i in Change]
    
    Dissip = H[1:] - H[0:-1]
    Dissip = [i/h for i in Dissip]

    return it, u, v, E, H, Change, Dissip

