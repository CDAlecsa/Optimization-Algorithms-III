'''
                                                    Packages
'''

from numpy import array, real
from numpy.linalg import norm
from scipy.integrate import ode


'''
                                                    Runge Kutta integrator
'''


def RK_2D(h, t_final, u0, v0, E, gradE, E_min):

    u, v = [], []
    u.append(u0)
    v.append(v0)

    y0, t0 = [u0[0], v0[0], u0[1], v0[1]], 0

    def f(t, y):
        return [y[1], -2 * y[1] - gradE([y[0], y[2]])[0], y[3], -2 * y[3] - gradE([y[0], y[2]])[1]]

    r = ode(f).set_integrator('dopri5')
    r.set_initial_value(y0, t0)

    lst_t = [t0]  # time steps
    lst_u1 = [y0[0]]
    lst_v1 = [y0[1]]
    lst_u2 = [y0[2]]
    lst_v2 = [y0[3]]

    while r.successful() and r.t <= t_final:

        r.integrate(r.t + h)
        lst_t.append(r.t)
        lst_u1.append(r.y[0])
        lst_v1.append(r.y[1])
        lst_u2.append(r.y[2])
        lst_v2.append(r.y[3])

    # Description
    # u -> [u1, u2]
    # v -> [v1, v2]

    lst_t = array(lst_t)
    lst_u1 = array(lst_u1)
    lst_v1 = array(lst_v1)
    lst_u2 = array(lst_u2)
    lst_v2 = array(lst_v2)

    u = [array([real(i), real(j)]) for i, j in zip(lst_u1, lst_u2)]
    v = [array([real(i), real(j)]) for i, j in zip(lst_v1, lst_v2)]

    u = array(u)
    v = array(v)

    E = [E(i) for i in u]
    H = [0.5 * norm(i) ** 2 + j - E_min for i, j in zip(v, E)]

    Change = [(i - j)/h for i, j in zip(E[1:], E[0:-1])]
    Dissip = [(i - j)/h for i, j in zip(H[1:], H[0:-1])]

    return lst_t, u, v, E, H, Change, Dissip
