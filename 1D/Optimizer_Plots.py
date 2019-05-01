import Splitting_Algorithms as SA
import matplotlib.pyplot as plt


'''
                                            Decay in the Objective Function
'''

plt.figure()
plt.plot(SA.tStrang, SA.E_Strang, color='black', marker='v', label='Strang')
plt.plot(SA.tPolyak, SA.E_P, color='yellow', marker='v', label='Polyak')
plt.plot(SA.tCN, SA.E_CN,  color='blue', marker='o', label='CN')
plt.plot(SA.tHeun, SA.E_Heun, color='red', marker='o', label='Heun')
plt.plot(SA.tStrang_PC, SA.E_Strang_PC, color='green', marker='*', label='Strang Predictor Corrector')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$E \left( u \right)$', fontsize=18)
plt.legend()
plt.title('Energy Decay')



'''
                                            Decay in the Regularization Function
'''

plt.figure()
plt.plot(SA.tStrang, SA.H_Strang, color='black', marker='v', label='Strang Splitting')
plt.plot(SA.tPolyak, SA.H_P, color='yellow', marker='v', label='Polyak')
plt.plot(SA.tCN, SA.H_CN, color='blue', marker='o', label='CN')
plt.plot(SA.tHeun, SA.H_Heun, color='red', marker='o', label='Heun')
plt.plot(SA.tStrang_PC, SA.H_Strang_PC, color='green', marker='*', label='Strang Predictor Corrector')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r'$H \left( u , v \right)$', fontsize=18)
plt.legend()
plt.title('Regularization Decay')



'''
                                            Phase Portrait
'''

plt.figure()
plt.plot(SA.u_Strang, SA.v_Strang, color='black', marker='v', label='Strang Splitting')
plt.plot(SA.u_CN, SA.v_CN, color='blue', marker='o', label='CN')
plt.plot(SA.u_P, SA.v_P, color='yellow', marker='v', label='Polyak')
plt.plot(SA.u_Heun, SA.v_Heun, color='red', marker='o', label='Heun')
plt.plot(SA.u_Strang_PC, SA.v_Strang_PC, color='green', marker='*', label='Strang Predictor Corrector')
plt.xlabel(r'$ u $', fontsize=18)
plt.ylabel(r'$ v $', fontsize=18)
plt.legend()
plt.title('Phase Portrait')



'''
                                            Energy Dissipation
'''

plt.figure()
plt.plot(SA.tStrang[1:], SA.Dissip_Strang, color='black', marker='v', label='Strang Splitting')
plt.plot(SA.tCN[1:], SA.Dissip_CN, color='blue', marker='o', label='CN')
plt.plot(SA.tPolyak[1:], SA.Dissip_P, color='yellow', marker='v', label='Polyak')
plt.plot(SA.tHeun[1:], SA.Dissip_Heun, color='red', marker='o', label='Heun')
plt.plot(SA.tStrang_PC[1:], SA.Dissip_Strang_PC, color='green', marker='*', label='Strang Predictor Corrector')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r' $ \frac{H_{n+1} - H_{n}}{h} $ ', fontsize=18)
plt.legend()



'''
                                            Change in the Objective Function    
'''

plt.figure()
plt.plot(SA.tStrang[1:], SA.Change_Strang, color='black', marker='v', label='Strang Splitting')
plt.plot(SA.tCN[1:], SA.Change_CN, color='blue', marker='o', label='CN')
plt.plot(SA.tPolyak[1:], SA.Change_P, color='yellow', marker='v', label='Polyak')
plt.plot(SA.tHeun[1:], SA.Change_Heun, color='red', marker='o', label='Heun')
plt.plot(SA.tStrang_PC[1:], SA.Change_Strang_PC, color='green', marker='*', label='Strang Predictor Corrector')
plt.xlabel(r'$t = nh$', fontsize=18)
plt.ylabel(r' $ \frac{E_{n+1} - E_{n}}{h} $ ', fontsize=18)
plt.legend()
plt.show()

