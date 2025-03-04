#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 21 12:25:39 2024

@author: s5257291
"""

import numpy as np
from matplotlib import pyplot as plt
# from qutip import *
# from matplotlib import cm
# import matplotlib as mpl
import matplotlib.ticker as mticker

plt.rcParams.update({"mathtext.fontset": 'cm',
                      'font.family':'times new roman'})
plt.rc('font', size=11) #controls default text size
plt.rc('axes', titlesize=11) #fontsize of the title
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rc('xtick', labelsize=11) #fontsize of the x tick labels
plt.rc('ytick', labelsize=11) #fontsize of the y tick labels
plt.rc('legend', fontsize=11) #fontsize of the legend

# # Return Analytic and Numerical ss distributions for R.P. laser
# def distributions(D,p,q,lam):

#     mu = (D-1)/2    
#     k = np.pi**2*p/(4*mu**2)    

#     pn_num = np.zeros(D, dtype='f16') # cav photon distribution (numerical)
#     pn_an = np.zeros(D, dtype='f16') # cav photon distribution (analytical)
#     pn_lin = np.zeros(D, dtype='f16') # cav photon distribution (linearized)
#     G = np.zeros((D,D), dtype='f16') # gain operator
#     L = np.zeros((D,D), dtype='f16') # loss operator

#     if q==0:    
        
#         # ss cav photon distribution (analytical)
#         for i in range(D):
#             pn_an[i] = (np.sin(np.pi*(i+1)/(D+1)))**p
#             pn_lin[i] = np.exp(-k*(i-mu)**2/2)
#         pn_an = pn_an/np.sum(pn_an) # normalise distribution
#         pn_lin = pn_lin/np.sum(pn_lin) # normalise distribution
        
#         for i in range(D):
#             if i > 0:
#                 G[i,i-1] = (pn_an[i]/pn_an[i-1])**(lam/2)
#                 L[i-1,i] = (pn_an[i-1]/pn_an[i])**((1-lam)/2)
#         Gdag = np.conjugate(np.transpose(G)) # adjoint of G
#         Ldag = np.conjugate(np.transpose(L)) # adjoint of L
        
#         pn_num = pn_an

#     else:    
        
#         # ss cav photon distribution (analytical)
#         for i in range(D):
#             pn_an[i] = (np.sin(np.pi*(i+1)/(D+1)))**p
#             pn_lin[i] = np.exp(-k*(i-mu)**2/2)
#         pn_an = pn_an/np.sum(pn_an) # normalise distribution
#         pn_lin = pn_lin/np.sum(pn_lin) # normalise distribution
        
#         for i in range(D):
#             if i > 0:
#                 G[i,i-1] = 1
#                 L[i-1,i] = (pn_an[i-1]/pn_an[i])**((1+q/2)/2)
#         Gdag = np.conjugate(np.transpose(G)) # ctranspose of G
#         Ldag = np.conjugate(np.transpose(L)) # ctranspose of L
        
#         # ss cav photon distribution (numerical)
#         for i in range(D):
#             if i == 0:
#                 pn_num[i] = 1.
#             elif i == 1:
#                 pn_num[i] = (1 - q/2)*pn_num[i-1]/L[i-1,i]**2
#             elif i == 2:
#                 pn_num[i] = ((q-1)*pn_num[i-2] + (1 + L[i-2,i-1]**2 - q/2)*pn_num[i-1])/L[i-1,i]**2
#             else:
#                 pn_num[i] = (-(q/2)*pn_num[i-3] + (q-1)*pn_num[i-2] + (1 + L[i-2,i-1]**2 - q/2)*pn_num[i-1])/L[i-1,i]**2
#         for i in range(D):
#             if pn_num[i] < 0:
#                 pn_num[i] = 0 # Get rid of negative values
#         pn_num = pn_num/np.sum(pn_num) # normalise distribution
    
#     return pn_an, pn_lin, pn_num, G, Gdag, L, Ldag

# # Action of Liouvillain
# def Action(pss_mat, G, Gdag, L, Ldag, q, lam):
    
#     if q==0:
#         change_gain1 = G@pss_mat@Gdag - 0.5*(Gdag@G@pss_mat + pss_mat@Gdag@G)
#         change_loss = L@pss_mat@Ldag - 0.5*(Ldag@L@pss_mat + pss_mat@Ldag@L)
#         b = np.trace(Ldag@L@pss_mat)
#         return (change_gain1 + change_loss)/b

#     else:
    
#         change_gain1 = G@pss_mat@Gdag - 0.5*(Gdag@G@pss_mat + pss_mat@Gdag@G)
#         change_gain2 = G@change_gain1@Gdag - 0.5*(Gdag@G@change_gain1 + change_gain1@Gdag@G)
#         change_loss = L@pss_mat@Ldag - 0.5*(Ldag@L@pss_mat + pss_mat@Ldag@L)
#         b = np.trace(Ldag@L@pss_mat)
#         return (change_gain1 + (q/2)*change_gain2 + change_loss)/b

# # Define pure cavity state
# def Pure_cav(pss,D):
#     pcav = np.zeros([D,D], dtype='f16')
#     for i in range(D):
#         for j in range(D):
#             pcav[i,j] = np.sqrt(pss[i]*pss[j])
#     return pcav

# # Define phase diffusion state change
# def Pure_diff(pcav,lw,D):
#     pchange = np.zeros([D,D], dtype='f16')
#     for i in range(D):
#         for j in range(D):
#             pchange[i,j] = -(lw/2)*(i-j)**2*pcav[i,j]
#     return pchange

# ###############################################################################
# ##### HS Calculations #########################################################
# ###############################################################################

# D_arr = np.arange(500,1001,25)
# p_arr = np.arange(10,51,2.5)

# HS_num_lam = np.zeros([len(D_arr),len(p_arr)], dtype = 'f16')
# HS_num_q = np.zeros([len(D_arr),len(p_arr)], dtype = 'f16')
# HS_num_0 = np.zeros([len(D_arr),len(p_arr)], dtype = 'f16')

# for i in range(len(D_arr)):
#     for j in range(len(p_arr)):
        
#         lam = 0.5
#         q = 0
        
#         pss_an, pss_lin, pss_num, G, Gdag, L, Ldag = distributions(int(D_arr[i]),p_arr[j],q,lam)
        
#         # Pure cavity states based on analytic and numeric ss
#         pcav_num = Pure_cav(pss_num,int(D_arr[i]))
    
#         # Action of Liouvillian on above pure cavity states
#         Acav_num = Action(pcav_num, G, Gdag, L, Ldag, q, lam)
        
#         # Pure diffusion change to states
#         lw = (2*lam**2 - 2*lam + 1) * np.pi**4 * p_arr[j]**2 / ( 64 * (D_arr[i]/2)**4 )
#         pchange_num = Pure_diff(pcav_num,lw,int(D_arr[i]))
        
#         print('Action lam: Done')
        
#         # Compute HS norms
        
#         HS_num_lam[i,j] = np.linalg.norm(Acav_num - pchange_num)/np.linalg.norm(pchange_num)
        
#         print('HS lam: Done')
        
#         lam = 0
#         q = -1
        
#         pss_an, pss_lin, pss_num, G, Gdag, L, Ldag = distributions(int(D_arr[i]),p_arr[j],q,lam)
        
#         # Pure cavity states based on analytic and numeric ss
#         pcav_num = Pure_cav(pss_num,int(D_arr[i]))
    
#         # Action of Liouvillian on above pure cavity states
#         Acav_num = Action(pcav_num, G, Gdag, L, Ldag, q, lam)
        
#         # Pure diffusion change to states
#         lw = (1 + q/2)**2 * np.pi**4 * p_arr[j]**2 / ( 64 * (D_arr[i]/2)**4 )
#         pchange_num = Pure_diff(pcav_num,lw,int(D_arr[i]))
        
#         print('Action lam: Done')
        
#         # Compute HS norms
        
#         HS_num_q[i,j] = np.linalg.norm(Acav_num - pchange_num)/np.linalg.norm(pchange_num)
        
#         print('HS lam: Done')
        
#         lam = 0
#         q = 0
        
#         pss_an, pss_lin, pss_num, G, Gdag, L, Ldag = distributions(int(D_arr[i]),p_arr[j],q,lam)
        
#         # Pure cavity states based on analytic and numeric ss
#         pcav_num = Pure_cav(pss_num,int(D_arr[i]))
    
#         # Action of Liouvillian on above pure cavity states
#         Acav_num = Action(pcav_num, G, Gdag, L, Ldag, q, lam)
        
#         # Pure diffusion change to states
#         lw = (1 + q/2)**2 * np.pi**4 * p_arr[j]**2 / ( 64 * (D_arr[i]/2)**4 )
#         pchange_num = Pure_diff(pcav_num,lw,int(D_arr[i]))
        
#         print('Action lam: Done')
        
#         # Compute HS norms
        
#         HS_num_0[i,j] = np.linalg.norm(Acav_num - pchange_num)/np.linalg.norm(pchange_num)
        
#         print('HS 0: Done')
        
#         print(i,j)
        
# np.save('p_arr',p_arr)
# np.save('D_arr',D_arr)
# np.save('HS_num_0',HS_num_0)
# np.save('HS_num_lam',HS_num_lam)
# np.save('HS_num_q',HS_num_q)

p_arr = np.load('p_arr.npy')
D_arr = np.load('D_arr.npy')
HS_num_0 = np.load('HS_num_0.npy')
HS_num_lam = np.load('HS_num_lam.npy')
HS_num_q = np.load('HS_num_q.npy')

fit_0 = np.polyfit(np.log(p_arr[:]), np.array(np.log(HS_num_0),dtype = 'float64')[-1,:], 1)
fit_lam = np.polyfit(np.log(p_arr[:]), np.array(np.log(HS_num_lam),dtype = 'float64')[-1,:], 1)
fit_q = np.polyfit(np.log(p_arr[:5]), np.array(np.log(HS_num_q),dtype = 'float64')[-1,:5], 1)

fig = plt.figure(figsize = (1*(3+3/8),1.75*(3+3/8)))

m = 0.1 # amount to take off axes width, in units of figure width
axs1 = fig.add_axes([0 + m, 2/3 + m, 1 - m, 1/3 - m]) # [left, bottom, width, height]

axs1.loglog(p_arr,np.exp(fit_0[1])*p_arr**fit_0[0], 'k-', label = r'$4.96p^{-1.11}$')
axs1.loglog(p_arr,HS_num_0[0,:],'x', color = 'C0', label = r'$\mu = 250$')
axs1.loglog(p_arr,HS_num_0[10,:],'x', color = 'C2', label = r'$\mu = 375$')
axs1.loglog(p_arr,HS_num_0[20,:],'x', color = 'C3', label = r'$\mu = 500$')
axs1.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs1.xaxis.set_major_formatter(mticker.ScalarFormatter())
axs1.set_xticks([10,20,30,40,50])
axs1.legend(loc = 'upper right', bbox_to_anchor=(1.01, 1.05), labelspacing = 0.2, handlelength = 1.0)
axs1.set_ylabel(r'Relative Distance')
axs1.set_xlabel(r'$p$')
axs1.text(0.025, 0.7, r'(a)', transform=axs1.transAxes)
# axs1.grid('on')
    
axs2 = fig.add_axes([0 + m, 1/3 + m, 1 - m, 1/3 - m]) # [left, bottom, width, height]

axs2.loglog(p_arr,np.exp(fit_lam[1])*p_arr**fit_lam[0], 'k-', label = r'$5.24p^{-1.13}$')
axs2.loglog(p_arr,HS_num_lam[0,:],'x', color = 'C0')
axs2.loglog(p_arr,HS_num_lam[10,:],'x', color = 'C2')
axs2.loglog(p_arr,HS_num_lam[20,:],'x', color = 'C3')
axs2.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs2.xaxis.set_major_formatter(mticker.ScalarFormatter())
axs2.set_xticks([10,20,30,40,50])
axs2.set_ylabel(r'Relative Distance')
axs2.set_xlabel(r'$p$')
axs2.legend(loc = 'upper right', labelspacing = 0.2, handlelength = 1.0)
axs2.text(0.025, 0.7, r'(b)', transform=axs2.transAxes)
    
axs3 = fig.add_axes([0 + m, 0/3 + 0*m, 1 - m, 1/3 - 0*m]) # [left, bottom, width, height]

axs3.loglog(p_arr,np.exp(fit_q[1])*p_arr**fit_q[0], 'k-', label = r'$R_0 = 5.90p^{-1.17}$')
axs3.loglog(p_arr,HS_num_q[0,:],'x', color = 'C0')
axs3.loglog(p_arr,HS_num_q[10,:],'x', color = 'C2')
axs3.loglog(p_arr,HS_num_q[20,:],'x', color = 'C3')
axs3.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs3.xaxis.set_major_formatter(mticker.ScalarFormatter())
axs3.set_xticks([10,20,30,40,50])
axs3.set_xlabel(r'$p$')
axs3.legend(loc = 'lower left', labelspacing = 0.2, handlelength = 1.0)
axs3.set_ylabel(r'Relative Distance')
axs3.text(0.025, 0.8, r'(c)', transform=axs3.transAxes)
axs3.set_xticks([10,20,30,40,50])

axsin = fig.add_axes([1/2 + m, 1/3 - 1/4.7 + m, 1/2.15 - m, 1/5 - m]) # [left, bottom, width, height]
x = D_arr/2
y = HS_num_q - np.exp(fit_q[1])*p_arr**fit_q[0]
fit_q2 = np.polyfit(np.log(x)[12:21], np.array(np.log(y[:,-1])[12:21],dtype = 'float64'), 1)
axsin.loglog(x[12:21],np.exp(fit_q2[1])*x[12:21]**fit_q2[0], '-', color = 'C1')
axsin.plot(x[12:21],y[:,-1][12:21],'.', color = 'C0')
axsin.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axsin.xaxis.set_major_formatter(mticker.ScalarFormatter())
axsin.set_xticks([400,450,500])
axsin.yaxis.set_minor_formatter(mticker.ScalarFormatter())
axsin.yaxis.set_major_formatter(mticker.ScalarFormatter())
axsin.set_yticks([0.016,0.018,0.020,0.022,0.024])
axsin.minorticks_off()
axsin.set_xlabel(r'$\mu$',labelpad = -1)
axsin.set_ylabel(r'$R-R_0$',labelpad = -1)
axsin.text(0.55, 0.75, r'$p=50$', transform=axsin.transAxes)


# fig.savefig('diff_verify.pdf',bbox_inches = 'tight')
















