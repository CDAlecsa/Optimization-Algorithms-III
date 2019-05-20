import numpy as np
import Splitting_Algorithms_NonConvex as SA
import matplotlib.pyplot as plt


'''
                                            Components - 2D Plots 
'''

minPoint = np.array([0, 0])

plt.figure()

plt.plot(SA.u_Strang[:, 0], SA.u_Strang[:, 1], 'ro-', markersize=10)
plt.plot(SA.u_Strang[0, 0], SA.u_Strang[0, 1], 'ks-', markersize=10)

plt.plot(SA.u_CN[:, 0], SA.u_CN[:, 1], 'bo-', markersize=10)
plt.plot(SA.u_RK[:, 0], SA.u_RK[:, 1], 'y*-', markersize=10)

plt.plot(minPoint[0], minPoint[1], 'go', markersize=10)
plt.plot(SA.u_Strang[-1, 0], SA.u_Strang[-1, 1], 'ys-', markersize=4)

plt.xlabel(r'$u_{1}$', fontsize=18)
plt.ylabel(r'$u_{2}$', fontsize=18)
plt.title('Contours of the Objective Function and Iterations for the Strang Splitting', fontsize=14)



'''
                                            Time Decay in the Regularization Function
'''


plt.figure()
plt.plot(SA.tStrang, SA.H_Strang, 'ro-')
plt.plot(SA.tRK, SA.H_RK, 'y*-')
plt.plot(SA.tCN, SA.H_CN, 'bv-')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$H \left( u , v \right)$', fontsize=18)
plt.title('Energy decay', fontsize=14)



'''
                                            Time Decay in the Energy Function
'''


plt.figure()
plt.plot(SA.tStrang, SA.E_Strang, 'ro-')
plt.plot(SA.tRK, SA.E_RK, 'y*-')
plt.plot(SA.tCN, SA.E_CN, 'bv-')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$E \left( u \right)$', fontsize=18)
plt.title('Energy function', fontsize=14)



'''
                                            Energy Dissipation
'''

plt.figure()
plt.plot(SA.tStrang[1:], SA.Dissip_Strang, color='red', marker='o', label='Strang Splitting')
plt.plot(SA.tCN[1:], SA.Dissip_CN, color='blue', marker='v', label='CN')
plt.plot(SA.tRK[1:], SA.Dissip_RK, color='yellow', marker='*', label='RK')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r' $ \frac{H_{n+1} - H_{n}}{h} $ ', fontsize=18)
plt.title('Energy dissipation', fontsize=14)
plt.show()

