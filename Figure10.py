

import numpy as np
import matplotlib.pyplot as plt
from scipy import integrate
from pylab import step,plot,show,legend,xlabel,xticks,ylabel,xlim,ylim,text,savefig,title


''' Function definition '''

def cells_in_genealogy_HSC1(a_,s_):
    
    m11= delta_[0]**(-1)*(2*l_[0]+a_[0])
    m12= delta_[0]**(-1)*(delta_[1]**(-1)*lambda_[0]*(2*l_[1]+a_[1])+2*s_[0]+a_[0])
    m13= delta_[0]**(-1)*delta_[1]**(-1)*lambda_[0]*(2*s_[1]+a_[1]+lambda_[1]*delta_[2]**(-1)*(2*l_[2]+a_[2]))
    m14= delta_[0]**(-1)*delta_[1]**(-1)*delta_[2]**(-1)*lambda_[0]*lambda_[1]*(2*s_[2]+a_[2]+lambda_[2]*delta_[3]**(-1)*(2*l_[3]+a_[3]))
    m15= delta_[0]**(-1)*delta_[1]**(-1)*delta_[2]**(-1)*delta_[3]**(-1)*lambda_[0]*lambda_[1]*lambda_[2]*(2*s_[3]+a_[3]+lambda_[3]*delta_[4]**(-1)*(2*l_[4]+a_[4]))

    return [m11,m12,m13,m14,m15]

def cells_in_genealogy_HSC2(a_,s_):
    
    m22= delta_[1]**(-1)*(2*l_[1]+a_[1])
    m23= delta_[1]**(-1)*(delta_[2]**(-1)*lambda_[1]*(2*l_[2]+a_[2])+2*s_[1]+a_[1])
    m24= delta_[1]**(-1)*delta_[2]**(-1)*lambda_[1]*(2*s_[2]+a_[2]+lambda_[2]*delta_[3]**(-1)*(2*l_[3]+a_[3]))
    m25= delta_[1]**(-1)*delta_[2]**(-1)*delta_[3]**(-1)*lambda_[1]*lambda_[2]*(2*s_[3]+a_[3]+lambda_[3]*delta_[4]**(-1)*(2*l_[4]+a_[4]))

    return [m22,m23,m24,m25]




''' Scenario Symm1 '''

gen_Symm1_HSC1 = np.zeros([5,2])

ind1, ind2 = 0, 0
for i in [0.0001, 0.001, 0.01, 0.1]:
    
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
    
    ''' Genealogy computation '''
    HSC1_CellsGenealogy=cells_in_genealogy_HSC1(a_,s_)
    if HSC1_CellsGenealogy[0]>=0:
        gen_Symm1_HSC1[:,ind1] = HSC1_CellsGenealogy
        ind1+=1

''' Plot section '''

s_vec = [0.01, 0.1]

fig_Symm1,ax_Symm1 = plt.subplots(1,1)

plt.plot(s_vec,gen_Symm1_HSC1[0,:],'^:', color='midnightblue', markersize=7)
plt.plot(s_vec,gen_Symm1_HSC1[1,:],'^:', color='firebrick', markersize=7)
plt.plot(s_vec,gen_Symm1_HSC1[2,:],'^:', color='seagreen', markersize=7)
plt.plot(s_vec,gen_Symm1_HSC1[3,:],'^:', color='gold', markersize=7)
plt.plot(s_vec,gen_Symm1_HSC1[4,:],'^:', color='royalblue', markersize=7, label='Scenario Symm1')

plt.xlabel('$s_1$',fontsize=15)
plt.ylabel('$m_1(j)$',fontsize=15, rotation=0)
ax_Symm1.set_yscale('log')
ax_Symm1.set_xlim([0.000071, 0.14])
ax_Symm1.set_ylim([0.1, 1e9])
ax_Symm1.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_Symm1.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_Symm1.set_xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)




''' Scenario SymmAll '''

gen_SymmAll_HSC1 = np.zeros([5,2])
gen_SymmAll_HSC2 = np.zeros([4,4])

ind1, ind2 = 0, 0

for i in [0.0001, 0.001, 0.01, 0.1]:
    
    ''' Parameter definition '''
    s_ = np.ones(5)*i
    a_ = np.zeros(5)
    phi_ = np.array([0.019063,0.04334,0.057295,0.14557,0.3385]) + a_ + s_
    lambda_ = np.array([0.016497,0.007847,0.032834,0.16113,0]) + a_ + 2*s_
    delta_ = np.array([(-0.0075899-0.0016495)/2,0.0017357,0.0044844,0.01556,0.0293]) + s_
    
    mu_ = delta_- lambda_ + phi_
    nu_ = lambda_ - a_ - 2*s_
    l_ = phi_ - a_ - s_
    
    ''' Genealogy computation '''
    
    HSC1_CellsGenealogy=cells_in_genealogy_HSC1(a_,s_)
    if HSC1_CellsGenealogy[0]>=0:
        gen_SymmAll_HSC1[:,ind1] = HSC1_CellsGenealogy
        ind1+=1
    
    HSC2_CellsGenealogy=cells_in_genealogy_HSC2(a_,s_)
    if HSC2_CellsGenealogy[0]>=0:
        gen_SymmAll_HSC2[:,ind2] = HSC2_CellsGenealogy
        ind2+=1
 
''' Plot section '''

s_vec = [0.01, 0.1]
s_vec2 = [0.0001, 0.001, 0.01, 0.1]

# HSC1 cells, m_1(j)
fig_SymmAll_HSC1,ax_SymmAll_HSC1 = plt.subplots(1,1)

plt.plot(s_vec,gen_SymmAll_HSC1[0,:],'o:', color='midnightblue', markersize=7)
plt.plot(s_vec,gen_SymmAll_HSC1[1,:],'o:', color='firebrick', markersize=7)
plt.plot(s_vec,gen_SymmAll_HSC1[2,:],'o:', color='seagreen', markersize=7)
plt.plot(s_vec,gen_SymmAll_HSC1[3,:],'o:', color='gold', markersize=7)
plt.plot(s_vec,gen_SymmAll_HSC1[4,:],'o:', color='royalblue', markersize=7, label='Scenario SymmAll')

plt.xlabel('$s_k$',fontsize=15)
plt.ylabel('$m_1(j)$',fontsize=15, rotation=0)
ax_SymmAll_HSC1.set_yscale('log')
ax_SymmAll_HSC1.set_xlim([0.000071, 0.14])
ax_SymmAll_HSC1.set_ylim([0.1, 1e9])
ax_SymmAll_HSC1.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_SymmAll_HSC1.set_xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# HSC2 cells - m_2(j)
fig_SymmAll_HSC2,ax_SymmAll_HSC2 = plt.subplots(1,1)

plt.plot(s_vec2,gen_SymmAll_HSC2[0,:],'o:', color='firebrick', markersize=7)
plt.plot(s_vec2,gen_SymmAll_HSC2[1,:],'o:', color='seagreen', markersize=7)
plt.plot(s_vec2,gen_SymmAll_HSC2[2,:],'o:', color='gold', markersize=7)
plt.plot(s_vec2,gen_SymmAll_HSC2[3,:],'o:', color='royalblue', markersize=7, label='Scenario SymmAll')

plt.xlabel('$s_k$',fontsize=15)
plt.ylabel('$m_2(j)$',fontsize=15, rotation=0)
ax_SymmAll_HSC2.set_yscale('log')
ax_SymmAll_HSC2.set_ylim([0.1, 1e9])
ax_SymmAll_HSC2.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_SymmAll_HSC2.set_xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)




''' Scenario AsymmAll '''

gen_AsymmAll_HSC1 = np.zeros([5,4])
gen_AsymmAll_HSC2 = np.zeros([4,4])

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
    
    ''' Genealogy computation '''
    
    HSC1_CellsGenealogy_a=cells_in_genealogy_HSC1(a_,s_)
    if HSC1_CellsGenealogy[0]>=0:
        gen_AsymmAll_HSC1[:,ind] = HSC1_CellsGenealogy_a
    HSC2_CellsGenealogy_a=cells_in_genealogy_HSC2(a_,s_)
    if HSC2_CellsGenealogy[0]>=0:
        gen_AsymmAll_HSC2[:,ind] = HSC2_CellsGenealogy_a
    ind+=1 



''' Plot section '''
s_vec2 = [0.0001, 0.001, 0.01, 0.1]

# HSC1 cells, m_1(j)
fig_AsymmAll_HSC1,ax_AsymmAll_HSC1 = plt.subplots(1,1)

plt.plot(s_vec2,gen_AsymmAll_HSC1[0,:],'s:', color='midnightblue', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC1[1,:],'s:', color='firebrick', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC1[2,:],'s:', color='seagreen', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC1[3,:],'s:', color='gold', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC1[4,:],'s:', color='royalblue', markersize=7, label='Scenario AsymmAll')

plt.xlabel('$a_k$',fontsize=15)
plt.ylabel('$m_1(j)$',fontsize=15, rotation=0)
ax_AsymmAll_HSC1.set_yscale('log')
ax_AsymmAll_HSC1.set_ylim([0.1, 1e9])
ax_AsymmAll_HSC1.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_AsymmAll_HSC1.set_xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

# HSC2 cells, m_2(j)
fig_AsymmAll_HSC2,ax_AsymmAll_HSC2 = plt.subplots(1,1)

plt.plot(s_vec2,gen_AsymmAll_HSC2[0,:],'s:', color='firebrick', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC2[1,:],'s:', color='seagreen', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC2[2,:],'s:', color='gold', markersize=7)
plt.plot(s_vec2,gen_AsymmAll_HSC2[3,:],'s:', color='royalblue', markersize=7, label='Scenario AsymmAll')

plt.xlabel('$a_k$',fontsize=15)
plt.ylabel('$m_2(j)$',fontsize=15, rotation=0)
ax_AsymmAll_HSC2.set_yscale('log')
ax_AsymmAll_HSC2.set_ylim([0.1, 1e9])
ax_AsymmAll_HSC2.set_yticks([10**1,10**3,10**5,10**7,10**9])
ax_AsymmAll_HSC2.set_xscale('log')
plt.xticks(fontsize=15)
plt.yticks(fontsize=15)

plt.show()





