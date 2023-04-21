import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# System of ODEs
def model(C,t,pSR,pSD,pAD):
    dc1dt = -(mu1+nu1+omega*pSD-omega*pSR)*C[0]
    dc2dt = (nu1+omega*pAD+2*omega*pSD)*C[0]-(mu2+nu2+omega*pSD-omega*pSR)*C[1]
    dc3dt = (nu2+omega*pAD+2*omega*pSD)*C[1]-(mu3+nu3+omega*pSD-omega*pSR)*C[2]
    dc4dt = (nu3+omega*pAD+2*omega*pSD)*C[2]
    return [dc1dt,dc2dt,dc3dt,dc4dt]


# Parameter values
mu1=1
mu2=1
mu3=1

nu1=0.5
nu2=0.5
nu3=0.5

omega=0.9

# initial condition
C0 = [100,0,0,0]

# time points
t = np.linspace(0,10)

# solve ODEs
C_Scenario1 = odeint(model,C0,t,args=(1.0,0,0))
C_Scenario2 = odeint(model,C0,t,args=(0.8,0.1,0.1))
C_Scenario3 = odeint(model,C0,t,args=(0.1,0.8,0.1))
C_Scenario4 = odeint(model,C0,t,args=(0.1,0.1,0.8))


# Dynamics over time
fig, axs = plt.subplots(2, 2)
axs[0, 0].plot(t,C_Scenario1[:,0],'r-',label=r'$C1$')
axs[0, 0].plot(t,C_Scenario1[:,1],'g-',label=r'$C2$')
axs[0, 0].plot(t,C_Scenario1[:,2],'y-',label=r'$C3$')
axs[0, 0].plot(t,C_Scenario1[:,3],'b-',label=r'$C4$')
axs[0,0].set_title('Only-SR')
axs[0,0].set_ylabel('Number of cells')

axs[0, 1].plot(t,C_Scenario2[:,0],'r-',label=r'$C1$')
axs[0, 1].plot(t,C_Scenario2[:,1],'g-',label=r'$C2$')
axs[0, 1].plot(t,C_Scenario2[:,2],'y-',label=r'$C3$')
axs[0, 1].plot(t,C_Scenario2[:,3],'b-',label=r'$C4$')
axs[0,1].set_title('Dominant-SR')

axs[1, 0].plot(t,C_Scenario3[:,0],'r-',label=r'$C1$')
axs[1, 0].plot(t,C_Scenario3[:,1],'g-',label=r'$C2$')
axs[1, 0].plot(t,C_Scenario3[:,2],'y-',label=r'$C3$')
axs[1, 0].plot(t,C_Scenario3[:,3],'b-',label=r'$C4$')
axs[1,0].set_title('Dominant-SD')
axs[1,0].set_ylabel('Number of cells')
axs[1,0].set_xlabel('time')

axs[1, 1].plot(t,C_Scenario4[:,0],'r-',label=r'$C1$')
axs[1, 1].plot(t,C_Scenario4[:,1],'g-',label=r'$C2$')
axs[1, 1].plot(t,C_Scenario4[:,2],'y-',label=r'$C3$')
axs[1, 1].plot(t,C_Scenario4[:,3],'b-',label=r'$C4$')
axs[1,1].set_title('Dominant-AD')
axs[1,1].set_xlabel('time')

lgd=fig.legend(loc='right')
fig.tight_layout(pad=2.0)

plt.show()


