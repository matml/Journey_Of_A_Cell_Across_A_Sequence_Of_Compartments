
from pylab import step,plot,show,legend,xlabel,xticks,ylabel,xlim,ylim
import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate

''' Function definition '''

def sys_ODE(Ec, t=0):
    Ec = np.array([-delta_[0]*Ec[0],
                    lambda_[0]*Ec[0]-delta_[1]*Ec[1],
                    lambda_[1]*Ec[1]-delta_[2]*Ec[2],
                    lambda_[2]*Ec[2]-delta_[3]*Ec[3],
                    lambda_[3]*Ec[3]-delta_[4]*Ec[4],])
    return Ec

''' Parameter values and initial conditions from Barile et al. 2020, for population of HSC1, HSC2, MPP1+2, MPP3, HPC1 '''

phi_ = np.array([0.019063,0.04334,0.057295,0.14557,0.3385])
lambda_ = np.array([0.016497,0.007847,0.032834,0.16113,0])
delta_ = np.array([(-0.0075899-0.0016495)/2,0.0017357,0.0044844,0.01556,0.0293])

EcIC = [890, 1370, 1540, 2020, 1.5*10**4]

timeDots = [50,75,100,170,200,250]


''' ODE system - definition and solution '''

timeODE = np.linspace(50, 1000, 100)

pop = integrate.odeint(sys_ODE, EcIC, timeODE)
E1,E2,E3,E4,E5 = pop.T


''' Plot section '''

fig = plt.figure(1)    
ax = fig.add_subplot(1,1,1)

plt.plot(timeODE,E1, label=r'$C_1(t)$')
plt.plot(timeODE,E2, label=r'$C_2(t)$')
plt.plot(timeODE,E3, label=r'$C_3(t)$')
plt.plot(timeODE,E4, label=r'$C_4(t)$')
plt.plot(timeODE,E5, label=r'$C_5(t)$')

xlabel('time (days)')
ylabel('cell number')
plt.yscale('log')
legend()

plt.show()



