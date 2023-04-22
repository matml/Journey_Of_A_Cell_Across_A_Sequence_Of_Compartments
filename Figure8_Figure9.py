

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from pylab import step,plot,show,legend,xlabel,xticks,ylabel,xlim,ylim,text,savefig,title


''' Function definition '''

def sys_ODE(Ec, t=0):
    Ec = np.array([-(mu_[0]+nu_[0]+s_[0]-l_[0])*Ec[0],
                    (2*s_[0]+a_[0]+nu_[0])*Ec[0]-(mu_[1]+nu_[1]+s_[1]-l_[1])*Ec[1],
                    (2*s_[1]+a_[1]+nu_[1])*Ec[1]-(mu_[2]+nu_[2]+s_[2]-l_[2])*Ec[2],
                    (2*s_[2]+a_[2]+nu_[2])*Ec[2]-(mu_[3]+nu_[3]+s_[3]-l_[3])*Ec[3],
                    (2*s_[3]+a_[3]+nu_[3])*Ec[3]-(mu_[4]+nu_[4]+s_[4]-l_[4])*Ec[4],])
    return Ec

''' Initial conditions '''

EcIC = [1, 0, 0, 0, 0]
timeODE = np.linspace(0, 3000, 500)
labels = (r'$E_1(t)$',r'$E_2(t)$',r'$E_3(t)$',r'$E_4(t)$',r'$E_5(t)$')


''' Scenario Symm1 '''

fig_Symm1, axs_Symm1 = plt.subplots(1, 2)
fig_Symm1.suptitle('Cell dynamics, Scenario Symm1')

ind = 0 
for i in [0.01, 0.1]:
    
    ''' Parameter definition '''
    s_ = np.zeros(5)
    s_[0] = i
    a_ = np.zeros(5)
    phi_ = np.array([0.019063,0.04334,0.057295,0.14557,0.3385]) + a_ + s_
    lambda_ = np.array([0.016497,0.007847,0.032834,0.16113,0]) + a_ + 2*s_
    delta_ = np.array([(-0.0075899-0.0016495)/2,0.0017357,0.0044844,0.01556,0.0293]) + s_
    
    mu_ = delta_- lambda_ + phi_
    nu_ = lambda_ - a_ - 2*s_
    l_ = phi_ - a_ - s_
    
    ''' Solve ODEs system '''
    pop = integrate.odeint(sys_ODE, EcIC, timeODE)
    E1,E2,E3,E4,E5 = pop.T
    
    ''' Plot section '''
    axs_Symm1[ind].plot(timeODE,E1)
    axs_Symm1[ind].plot(timeODE,E2)
    axs_Symm1[ind].plot(timeODE,E3)
    axs_Symm1[ind].plot(timeODE,E4)
    axs_Symm1[ind].plot(timeODE,E5)
    
    axs_Symm1[0].set_ylabel('cell number')
    axs_Symm1[0].set_xlabel('time (days)')
    axs_Symm1[1].set_xlabel('time (days)')
    axs_Symm1[0].set_xlim(0,3000)
    axs_Symm1[1].set_xlim(0,3000)
    axs_Symm1[0].set_ylim(0,100)
    axs_Symm1[1].set_ylim(0,100)
    fig_Symm1.legend(labels)
    
    ind+=1
    
    
''' Scenario SymmAll '''

fig_SymmAll, axs_SymmAll = plt.subplots(1, 2)
fig_SymmAll.suptitle('Cell dynamics, Scenario SymmAll')

ind = 0   
for i in [0.01, 0.1]:

    ''' Parameter definition '''
    s_ = np.ones(5)*i
    a_ = np.zeros(5)
    phi_ = np.array([0.019063,0.04334,0.057295,0.14557,0.3385]) + a_ + s_
    lambda_ = np.array([0.016497,0.007847,0.032834,0.16113,0]) + a_ + 2*s_
    delta_ = np.array([(-0.0075899-0.0016495)/2,0.0017357,0.0044844,0.01556,0.0293]) + s_
    
    mu_ = delta_- lambda_ + phi_
    nu_ = lambda_ - a_ - 2*s_
    l_ = phi_ - a_ - s_
    
    ''' Solve ODEs system '''
    pop = integrate.odeint(sys_ODE, EcIC, timeODE)
    E1,E2,E3,E4,E5 = pop.T
    
    ''' Plot section '''
    axs_SymmAll[ind].plot(timeODE,E1)
    axs_SymmAll[ind].plot(timeODE,E2)
    axs_SymmAll[ind].plot(timeODE,E3)
    axs_SymmAll[ind].plot(timeODE,E4)
    axs_SymmAll[ind].plot(timeODE,E5)
    
    axs_SymmAll[0].set_ylabel('cell number')
    axs_SymmAll[0].set_xlabel('time (days)')
    axs_SymmAll[1].set_xlabel('time (days)')
    axs_SymmAll[0].set_xlim(0,1000)
    axs_SymmAll[1].set_xlim(0,1000)
    axs_SymmAll[0].set_ylim(0,100)
    axs_SymmAll[1].set_ylim(0,100)
    fig_SymmAll.legend(labels)
    
    ind+=1
        


''' Scenario AsymmAll '''

fig_AsymmAll, axs_AsymmAll = plt.subplots(1, 4)
fig_AsymmAll.suptitle('Cell dynamics, Scenario AsymmAll')

ind = 0
for i in [0.0001, 0.001, 0.01, 0.1]:
    
    ''' Parameter definition '''
    a_ = np.ones(5)*i
    s_ = np.zeros(5)
    s_[0] = 0.005
    phi_ = np.array([0.019063,0.04334,0.057295,0.14557,0.3385]) + a_ + s_
    lambda_ = np.array([0.016497,0.007847,0.032834,0.16113,0]) + a_ + 2*s_
    delta_ = np.array([(-0.0075899-0.0016495)/2,0.0017357,0.0044844,0.01556,0.0293]) + s_
    
    mu_ = delta_- lambda_ + phi_
    nu_ = lambda_ - a_ - 2*s_
    l_ = phi_ - a_ - s_
    
    ''' Solve ODEs system '''
    pop = integrate.odeint(sys_ODE, EcIC, timeODE)
    E1,E2,E3,E4,E5 = pop.T
    
    ''' Plot section '''
    axs_AsymmAll[ind].plot(timeODE,E1)
    axs_AsymmAll[ind].plot(timeODE,E2)
    axs_AsymmAll[ind].plot(timeODE,E3)
    axs_AsymmAll[ind].plot(timeODE,E4)
    axs_AsymmAll[ind].plot(timeODE,E5)
    
    axs_AsymmAll[0].set_ylabel('cell number')
    axs_AsymmAll[0].set_xlabel('time (days)')
    axs_AsymmAll[1].set_xlabel('time (days)')
    axs_AsymmAll[2].set_xlabel('time (days)')
    axs_AsymmAll[3].set_xlabel('time (days)')
    axs_AsymmAll[0].set_yscale('log')
    axs_AsymmAll[1].set_yscale('log')
    axs_AsymmAll[2].set_yscale('log')
    axs_AsymmAll[3].set_yscale('log')
    axs_AsymmAll[0].set_ylim(1e-3,1e5)
    axs_AsymmAll[1].set_ylim(1e-3,1e5)
    axs_AsymmAll[2].set_ylim(1e-3,1e5)
    axs_AsymmAll[3].set_ylim(1e-3,1e5)
    fig_AsymmAll.legend(labels)
    
    ind+=1





plt.show()
