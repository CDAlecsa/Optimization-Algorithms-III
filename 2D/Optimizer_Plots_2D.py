import numpy as np
import Splitting_Algorithms as SA
import matplotlib.pyplot as plt


'''
                                            Components - 2D Plots 
'''

minPoint = np.array([0, 0])

plt.figure()
plt.plot(SA.u_Strang[:, 0], SA.u_Strang[:, 1], 'kv-')
plt.plot(SA.sol[:, 0], SA.sol[:, 1], 'g-')
plt.plot(minPoint[0], minPoint[1], 'bo')
plt.xlabel(r'$u_{1}$', fontsize=18)
plt.ylabel(r'$u_{2}$', fontsize=18)
plt.title('Contours of the Objective Function and Iterations for the Strang Splitting', fontsize=14)



'''
                                            Time Decay in the Regularization Function
'''


plt.figure()
plt.plot(SA.tStrang, SA.H_Strang, 'kv-')
plt.plot(SA.tStrang, SA.H, 'b--')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$H \left( u , v \right)$', fontsize=18)
plt.title('Energy decay', fontsize=14)


'''
                                            Absolute Error between the Exact and the Numerical Solution
'''


plt.figure()
plt.plot(SA.tStrang, SA.err, 'r*-')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$Error$', fontsize=18)
plt.title('Absolute error', fontsize=14)
plt.ticklabel_format(style='sci', axis='y', scilimits=(0, 0))
plt.show()

