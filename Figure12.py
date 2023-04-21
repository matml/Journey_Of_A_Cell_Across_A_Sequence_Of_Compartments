import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt

# System of ODEs
def model(C,t,pSR,pSD,pAD):
    dc1dt = -(mu1+nu1+l*pSD-l*pSR)*C[0]
    dc2dt = (nu1+l*pAD+2*l*pSD)*C[0]-(mu2+nu2+l*pSD-l*pSR)*C[1]
    dc3dt = (nu2+l*pAD+2*l*pSD)*C[1]-(mu3+nu3+l*pSD-l*pSR)*C[2]
    dc4dt = (nu3+l*pAD+2*l*pSD)*C[2]
    return [dc1dt,dc2dt,dc3dt,dc4dt]

# Parameter values
mu1 = 0.263
nu1 = 0.137
mu2 = 1.369

nu24 = 0.07
nu28 = 0.054

l4 = 0.216
l8 = 0.093

mu4 = 0.04
mu8 = 0.11

nu4 = 0.21
nu8 = 0.14

# Cellular fate probabilities
b1 = mu1/(mu1+nu1)
b2 = nu1/(mu1+nu1) * mu2/(mu2+nu24+nu28)
b4 = nu1/(mu1+nu1) * nu24/(mu2+nu24+nu28) * mu4/(mu4+nu4)
b8 = nu1/(mu1+nu1) * nu28/(mu2+nu24+nu28) * mu8/(mu8+nu8)
b4p = nu1/(mu1+nu1) * nu24/(mu2+nu24+nu28) * nu4/(mu4+nu4)
b8p = nu1/(mu1+nu1) * nu28/(mu2+nu24+nu28) * nu8/(mu8+nu8)

# Bar chart with cellular fate probabilities

labels = ('pre-DP\n death','post-DP\n death','CD4 SP\n death', 'CD8 SP\n death','CD4\n Periphery', 'CD8\n Periphery')
y_pos = 2*np.arange(len(labels))
fate_probabilities=[b1,b2,b4,b8,b4p,b8p]

plt.bar(y_pos,fate_probabilities,align='center',alpha=0.5,log=True,width=1.0,color=['orange','green','blue','red','cyan','pink'])
plt.xticks(y_pos,labels)
plt.ylabel('Probability')
plt.title('Cellular fate probabilities')
plt.show()

lifespan = 1/(mu1+nu1) * (nu1/(mu2+nu24+nu28)*(nu24/(mu4+nu4)+nu28/(mu8+nu8)+1)+1)
print("\n Average lifespan: ",lifespan)

eta4 = nu1*nu24*l4/((mu1+nu1)*(mu2+nu24+nu28)*(mu4+nu4))
eta8 = nu1*nu28*l8/((mu1+nu1)*(mu2+nu24+nu28)*(mu8+nu8))

print("\n Average number of division events: ",eta4,eta8,eta4+eta8)

print("\n Cellular fate probabilities: ",fate_probabilities)

print("\n Sensitivities:")

print("dBeta_1(1)/dmu1: ",(nu1)/(mu1+nu1)**2)
print("dBeta_1(1)/dnu1: ",-mu1/(mu1+nu1)**2)
print("dBeta_1(2)/dmu1: ",-mu2*nu1/((mu1+nu1)**2*(mu2+nu24+nu28)))
print("dBeta_1(2)/dnu1: ",(2*nu1+mu1)*mu2/((mu1+nu1)**2*(mu2+nu24+nu28)))
print("dBeta_1(2)/dmu2: ",(2*mu2+nu24+nu28)*nu1/((mu1+nu1)*(mu2+nu24+nu28)**2))
print("dBeta_1(2)/dnu24: ",-mu2*nu1/((mu1+nu1)*(mu2+nu24+nu28)**2))
print("dBeta_1(2)/dnu28: ",-mu2*nu1/((mu1+nu1)*(mu2+nu24+nu28)**2))
print("dBeta_1(4)/dmu1: ",-nu1*nu24*mu4/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu4+nu4)))
print("dBeta_1(4)/dnu1: ",(2*nu1+mu1)*nu24+mu4/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu4+nu4)))
print("dBeta_1(4)/dmu2: ",-nu1*nu24*mu4/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4)/dnu24: ",nu1*mu4*(2*nu24+mu2+nu28)/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4)/dnu28: ",-nu1*nu24*mu4/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4)/dmu4: ",nu1*nu24*(2*mu4*nu4)/((mu1+nu1)*(mu2+nu24+nu28)*(mu4+nu4)**2))
print("dBeta_1(4)/dnu4: ",-nu1*nu24*mu4/((mu1+nu1)*(mu2+nu24+nu28)*(mu4+nu4)**2))
print("dBeta_1(8)/dmu1: ",-nu1*nu28*mu8/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu8+nu8)))
print("dBeta_1(8)/dnu1: ",(2*nu1+mu1)*nu28*mu8/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu8+nu8)))
print("dBeta_1(8)/dmu2: ",-nu1*nu28*mu8/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8)/dnu24: ",-nu1*nu28*mu8/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8)/dnu28: ",nu1*(2*nu28+mu2+nu24)*mu8/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8)/dmu8: ",(2*mu8+nu8)*nu1*nu28/((mu1+nu1)*(mu2+nu24+nu28)*(mu8+nu8)**2))
print("dBeta_1(8)/dnu8: ",-nu1*nu28*mu8/((mu1+nu1)*(mu2+nu24+nu28)*(mu8+nu8)**2))
print("dBeta_1(4P)/dmu1: ",-nu1*nu24*nu4/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu4+nu4)))
print("dBeta_1(4P)/dnu1: ",(2*nu1+mu1)*nu24*nu4/((mu1+nu1)**2*(mu2+nu24+nu28)*(mu4+nu4)))
print("dBeta_1(4P)/dmu2: ",-nu1*nu24*nu4/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4P)/dnu24: ",(2*nu24+mu2+nu28)*nu1*nu4/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4P)/dnu28: ",-nu1*nu24*nu4/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu4+nu4)))
print("dBeta_1(4P)/dmu4: ",-nu1*nu24*nu4/((mu1+nu1)*(mu2+nu24+nu28)*(mu4+nu4)**2))
print("dBeta_1(4P)/dnu4: ",(2*nu4+mu4)*nu1*nu24/((mu1+nu1)*(mu2+nu24+nu28)*(mu4+nu4)**2))
print("dBeta_1(8P)/dmu1: ",-nu1*nu28*nu8/((mu1+nu1)**2*(mu4+nu24+nu28)*(mu8+nu8)))
print("dBeta_1(8P)/dnu1: ",(2*nu1+mu1)*nu28*nu8/((mu1+nu1)**2*(mu4+nu24+nu28)*(mu8+nu8)))
print("dBeta_1(8P)/dmu2: ",-nu1*nu28*nu8/((mu1+nu1)*(mu4+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8P)/dnu24: ",-nu1*nu28*nu8/((mu1+nu1)*(mu4+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8P)/dnu28: ",(2*nu28+mu2+nu24)*nu1*nu8/((mu1+nu1)*(mu2+nu24+nu28)**2*(mu8+nu8)))
print("dBeta_1(8P)/dmu8: ",-nu1*nu28*nu8/((mu1+nu1)*(mu4+nu24+nu28)*(mu8+nu8)**2))
print("dBeta_1(8P)/dnu8: ",(2*nu8+mu8)*nu1*nu28/((mu1+nu1)*(mu4+nu24+nu28)*(mu8+nu8)**2))
