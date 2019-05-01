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

from Strang_2D import Strang_2D




'''
                                        Objective Function & Gradient of the Objective Function
'''


def E(u, coeff_1, coeff_2):
    return coeff_1 * u[0] ** 2 + coeff_2 * u[1] ** 2


def gradE(u, coeff_1, coeff_2):
    return [2 * coeff_1 * u[0], 2 * coeff_2 * u[1]]




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


u0, v0, tol, h_Strang, coeff_1, coeff_2 = dict_values.values()
E_min = 0


'''
                                            Strang Optimizer
'''


itStrang, u_Strang, v_Strang, E_Strang, H_Strang = Strang_2D(h_Strang, u0, v0, E, gradE, coeff_1, coeff_2, E_min, tol)


'''
                                            Time Grid for the Strang Optimizer
'''

tStrang = (np.linspace(1, itStrang, itStrang)-1) * h_Strang


'''
                                            Exact Solution
'''

Coeff1_List = [coeff_1] * len(tStrang)
Coeff2_List = [coeff_2] * len(tStrang)

sol1 = [(1/2)*(1/(2*j-1))*(((2*j-1-np.sqrt(1-2*j))*np.exp(i*(-1+np.sqrt(1-2*j)))) + ((2*j-1+np.sqrt(1-2*j))*np.exp(i*(-1-np.sqrt(1-2*j))))) for i, j in zip(tStrang, Coeff1_List)]
sol2 = [(1/2)*(1/(2*j-1))*(((2*j-1-np.sqrt(1-2*j))*np.exp(i*(-1+np.sqrt(1-2*j)))) + ((2*j-1+np.sqrt(1-2*j))*np.exp(i*(-1-np.sqrt(1-2*j))))) for i, j in zip(tStrang, Coeff2_List)]

v1 = [(1/2)*(1/(2*j-1))*(((-1+2*j-np.sqrt(1-2*j))*(-1+np.sqrt(1-2*j))*np.exp(i*(-1+np.sqrt(1-2*j)))) + ((2*j-1+np.sqrt(1-2*j))*(-1-np.sqrt(1-2*j))*np.exp(i*(-1-np.sqrt(1-2*j))))) for i, j in zip(tStrang, Coeff1_List)]
v2 = [(1/2)*(1/(2*j-1))*(((-1+2*j-np.sqrt(1-2*j))*(-1+np.sqrt(1-2*j))*np.exp(i*(-1+np.sqrt(1-2*j)))) + ((2*j-1+np.sqrt(1-2*j))*(-1-np.sqrt(1-2*j))*np.exp(i*(-1-np.sqrt(1-2*j))))) for i, j in zip(tStrang, Coeff2_List)]

sol = np.column_stack((sol1, sol2))
v = np.column_stack((v1, v2))

E = [E(i, coeff_1, coeff_2) for i in sol]
H = [0.5 * i ** 2 + j - E_min for i, j in zip(abs(v), E)]


'''
                                                Compute Absolute Error
'''

err = [norm(i-j) for i, j in zip(u_Strang, sol)]

