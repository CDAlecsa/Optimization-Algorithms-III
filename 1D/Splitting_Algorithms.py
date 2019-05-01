'''
                                        Import Packages
'''

import numpy as np


import seaborn as sns
sns.set(palette='deep')



'''
                                        Import Modules
'''

from Heun import Heun
from Polyak import Polyak
from CrankNicolson import CN
from Strang import Strang
from Strang_PC import Strang_PC



'''
                                        Objective Function & Gradient of the Objective Function
'''


def E(u):
    return u ** 4 / 4 - u ** 2 / 2 + 1 / 4


def gradE(u):
    return u ** 3 - u


'''
                                        Critical Points of the Objective Function
'''

optPoints = np.array([-1, 0, 1])
optValues = np.array([E(-1), E(0), E(1)])

E_min1 = optValues[0]
E_max = optValues[1]
E_min2 = optValues[2]



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
        if index_line != 1:
            dict_values[key] = float(value)
        else:
            dict_values[key] = value
            dict_values[key] = value.strip()

        line = fp.readline()
        index_line += 1


condStop, u0, v0, tol, h_Heun, h_CN, h_Strang, h_Strang_PC, h_Polyak = dict_values.values()

if u0 < 0:
    E_min = E_min1
elif u0 > 0:
    E_min = E_min2
else:
    E_min = E_max


'''
                                            Optimizers
'''


itStrang, u_Strang, v_Strang, E_Strang, H_Strang, Change_Strang, Dissip_Strang = Strang(h_Strang, u0, v0, condStop, E, gradE, E_min, tol)
itStrang_PC, u_Strang_PC, v_Strang_PC, E_Strang_PC, H_Strang_PC, Change_Strang_PC, Dissip_Strang_PC = Strang_PC(h_Strang_PC, u0, v0, condStop, E, gradE, E_min, eps=1e-8)

itCN, u_CN, v_CN, E_CN, H_CN, Change_CN, Dissip_CN = CN(h_CN, u0, v0, condStop, E, gradE, E_min, eps=1e-8)
itP, u_P, v_P, E_P, H_P, Change_P, Dissip_P = Polyak(h_Polyak, u0, v0, condStop, E, gradE, E_min, tol)
itHeun, u_Heun, v_Heun, E_Heun, H_Heun, Change_Heun, Dissip_Heun = Heun(h_Heun, u0, v0, condStop, E, gradE, E_min, tol)


'''
                                            Time Grid for the Optimizers
'''

tStrang = (np.linspace(1, itStrang, itStrang)-1) * h_Strang
tStrang_PC = (np.linspace(1, itStrang_PC, itStrang_PC)-1) * h_Strang_PC
tPolyak = (np.linspace(1, itP, itP)-1) * h_Polyak
tCN = (np.linspace(1, itCN, itCN)-1) * h_CN
tHeun = (np.linspace(1, itHeun, itHeun)-1) * h_Heun


