'''
                                                    Packages
'''

from numpy import array, cos, sin, sqrt
from numpy.linalg import norm
from scipy.optimize import fsolve


'''
                                                    Newton Solver for the Implicit Method
'''


def F_Newton(x, u_n, v_n, h, c):

#    x[0] = nan_to_num(x[0])
#    x[1] = nan_to_num(x[1])
    
    A = 1 / ( 6 * ( u_n[0] + u_n[1] - x[0] - x[1] ) )

    
    B = 9 * cos(sqrt(2)*x[0]) * cos(sqrt(2)*x[1]) - \
    9 * sin(sqrt(2)*x[0]) * sin(sqrt(2)*x[1]) - \
    9 * cos(sqrt(2)*u_n[0]) * cos(sqrt(2)*u_n[1]) + \
    9 * sin(sqrt(2)*u_n[0]) * sin(sqrt(2)*u_n[1])
    
    C1 = 4 * ( u_n[0] + x[0] ) * ( u_n[0] + u_n[1] - x[0] - x[1] )
    C2 = 4 * ( u_n[1] + x[1] ) * ( u_n[0] + u_n[1] - x[0] - x[1] )

    auxNewton1 = A * ( B + C1 )
    auxNewton2 = A * ( B + C2 )

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

        
        guess = [0, 0]    
        args = (u[-1], v[-1], h, c)    
        u.append(fsolve(F_Newton, guess, args, xtol=1e-10))
        

        A = 1 / ( 6 * ( u[-2][0] + u[-2][1] - u[-1][0] - u[-1][1] ) )
    
        B = 9 * cos(sqrt(2)*u[-1][0]) * cos(sqrt(2)*u[-1][1]) - \
        9 * sin(sqrt(2)*u[-1][0]) * sin(sqrt(2)*u[-1][1]) - \
        9 * cos(sqrt(2)*u[-2][0]) * cos(sqrt(2)*u[-2][1]) + \
        9 * sin(sqrt(2)*u[-2][0]) * sin(sqrt(2)*u[-2][1])
    
        C1 = 4 * ( u[-2][0] + u[-1][0] ) * ( u[-2][0] + u[-2][1] - u[-1][0] - u[-1][1] )
        C2 = 4 * ( u[-2][1] + u[-1][1] ) * ( u[-2][0] + u[-2][1] - u[-1][0] - u[-1][1] )

        aux1 = A * ( B + C1 )
        aux2 = A * ( B + C2 )

        v.append([(1 / c ** 2) * v[-1][0] - (h / c) * aux1, (1 / c ** 2) * v[-1][1] - (h / c) * aux2])

        if norm(gradE(u[-1])) <= eps:
            ok = 1
        n = n + 1

    it = n + 1

    u = array(u)
    v = array(v)

    E = [E(i) for i in u]
    H = [0.5 * norm(i) ** 2 + j - E_min for i, j in zip(v, E)]

    Change = [i-j for i, j in zip(E[1:], E[0:-1])]
    Dissip = [i-j for i, j in zip(H[1:], H[0:-1])]

    return it, u, v, E, H, Change, Dissip
