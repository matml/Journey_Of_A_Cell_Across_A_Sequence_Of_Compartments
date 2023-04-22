
import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

''' Function definition '''

def cells_in_genealogy(pSR,pSD,pAD):
    Delta1 = mu1+nu1+omega*pSD-omega*pSR
    Delta2 = mu2+nu2+omega*pSD-omega*pSR
    Delta3 = mu3+nu3+omega*pSD-omega*pSR
    
    Lambda1 = nu1+omega*pAD+2*omega*pSD
    Lambda2 = nu2+omega*pAD+2*omega*pSD
    Lambda3 = nu3+omega*pAD+2*omega*pSD
    
    m11= Delta1**(-1)*(2*omega*pSR+omega*pAD)
    m12= Delta1**(-1)*(Delta2**(-1)*Lambda1*(2*omega*pSR+omega*pAD)+omega*pAD+2*omega*pSD)
    m13= Delta1**(-1)*Delta2**(-1)*Lambda1*(2*omega*pSD+omega*pAD+Lambda2*Delta3**(-1)*(2*omega*pSR+omega*pAD))
    m14= Delta1**(-1)*Delta2**(-1)*Delta3**(-1)*Lambda1*Lambda2*(2*omega*pSD+omega*pAD)

    return [m11,m12,m13,m14]

''' Parameter values '''

mu1=1
mu2=1
mu3=1

nu1=0.5
nu2=0.5
nu3=0.5

omega=0.9

''' Calculation of the mean number of cells in the genealogy '''

Scenario_OnlySR_CellsGenealogy=cells_in_genealogy(1.0,0.0,0.0)
Scenario_DominantSR_CellsGenealogy=cells_in_genealogy(0.8,0.1,0.1)
Scenario_DominantSD_CellsGenealogy=cells_in_genealogy(0.1,0.8,0.1)
Scenario_DominantAD_CellsGenealogy=cells_in_genealogy(0.1,0.1,0.8)


''' Plot section '''


labels = (r'$m_1(1)$',r'$m_1(2)$',r'$m_1(3)$',r'$m_1(4)$')
y_pos = 2*np.arange(len(labels))

fig, axs = plt.subplots(2, 2)

axs[0,0].bar(y_pos,Scenario_OnlySR_CellsGenealogy,align='center',alpha=0.5,width=1.0,color=['red','green','yellow','blue'])
axs[0,0].set_xticks(y_pos)
axs[0,0].set_xticklabels(labels)
axs[0,0].set_ylabel('Number of cells')
axs[0,0].set_title('Only-SR')
axs[0,0].set_ylim(0,4)

axs[0,1].bar(y_pos,Scenario_DominantSR_CellsGenealogy,align='center',alpha=0.5,width=1.0,color=['red','green','yellow','blue'])
axs[0,1].set_xticks(y_pos)
axs[0,1].set_xticklabels(labels)

axs[0,1].set_title('Dominant-SR')
axs[0,1].set_ylim(0,4)

axs[1,0].bar(y_pos,Scenario_DominantSD_CellsGenealogy,align='center',alpha=0.5,width=1.0,color=['red','green','yellow','blue'])
axs[1,0].set_xticks(y_pos)
axs[1,0].set_xticklabels(labels)
axs[1,0].set_ylabel('Number of cells')
axs[1,0].set_title('Dominant-SD')
axs[1,0].set_ylim(0,4)

axs[1,1].bar(y_pos,Scenario_DominantAD_CellsGenealogy,align='center',alpha=0.5,width=1.0,color=['red','green','yellow','blue'])
axs[1,1].set_xticks(y_pos)
axs[1,1].set_xticklabels(labels)
axs[1,1].set_title('Dominant-AD')
axs[1,1].set_ylim(0,4)

fig.tight_layout(pad=2.0)

plt.show()
