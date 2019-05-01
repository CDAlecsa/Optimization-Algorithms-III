'''
                                        Import Packages
'''

import numpy as np
from numpy.linalg import norm


import seaborn as sns
sns.set(palette='deep')



'''
                                        Import Modules
'''

from Strang_2D_Rosenbrock import Strang_2D
from CN_2D_Rosenbrock import CN_2D
from RK_2D_Rosenbrock import RK_2D


'''
                                                    Objective Function
'''


def E(u):
    return 100 * (u[1] - u[0] ** 2) ** 2 + 1 * (1 - u[0]) ** 2


'''
                                                    Gradient of the Objective Function
'''


def gradE(u):
    return np.array([(-4) * 100 * u[0] * (u[1] - u[0]**2) - 2 * (1 - u[0]), 200 * 1 * (u[1] - u[0] ** 2)])



'''
                                        Reading the Parameters 
'''

file_path = 'Parameters.txt'
dict_values = {}

with open(file_path) as fp:
    line = fp.readline()
    index_line = 1
    while line:
        print('Line {} : Value {}'.format(index_line, line.strip()))

        key, value = line.split('=')

        value = value.rstrip()
        if len(value.split()) != 1:
            dict_values[key] = [float(i) for i in value.split()]
        else:
            dict_values[key] = float(value)

        line = fp.readline()
        index_line += 1


u0, v0, tol, h_Strang, h_CN, h_RK = dict_values.values()
E_min = 0


'''
                                            Strang and CN Optimizers
'''

itStrang, u_Strang, v_Strang, E_Strang, H_Strang, Change_Strang, Dissip_Strang = Strang_2D(h_Strang, u0, v0, E, gradE, E_min, tol)
itCN, u_CN, v_CN, E_CN, H_CN, Change_CN, Dissip_CN = CN_2D(h_CN, u0, v0, E, gradE, E_min, tol)



'''
                                            Time Grid for the Optimizers
'''

tStrang = (np.linspace(1, itStrang, itStrang)-1) * h_Strang
tCN = (np.linspace(1, itCN, itCN)-1) * h_CN


'''
                                            Runge Kutta Method
'''

t_final = tStrang[-1]
tRK, u_RK, v_RK, E_RK, H_RK, Change_RK, Dissip_RK = RK_2D(h_RK, t_final, u0, v0, E, gradE, E_min)

