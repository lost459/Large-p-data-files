#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  8 14:16:57 2025

@author: s5257291
"""

import numpy as np
from matplotlib import pyplot as plt
# from qutip import *
# from matplotlib import cm
# import matplotlib as mpl
import matplotlib.ticker as mticker
from scipy.stats import linregress

plt.rcParams.update({"mathtext.fontset": 'cm',
                      'font.family':'times new roman'})
plt.rc('font', size=11) #controls default text size
plt.rc('axes', titlesize=11) #fontsize of the title
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rc('xtick', labelsize=11) #fontsize of the x tick labels
plt.rc('ytick', labelsize=11) #fontsize of the y tick labels
plt.rc('legend', fontsize=11) #fontsize of the legend

# Return Analytic and Numerical ss distributions for R.P. laser
def distributions(D,p,q,lam):

    mu = (D-1)/2    
    k = np.pi**2*p/(4*mu**2)    

    pn_num = np.zeros(D, dtype='float128') # cav photon distribution (numerical)
    pn_an = np.zeros(D, dtype='float128') # cav photon distribution (analytical)
    pn_lin = np.zeros(D, dtype='float128') # cav photon distribution (linearized)
    G = np.zeros((D,D), dtype='float128') # gain operator
    L = np.zeros((D,D), dtype='float128') # loss operator

    if q==0:    
        
        # ss cav photon distribution (analytical)
        for i in range(D):
            pn_an[i] = (np.sin(np.pi*(i+1)/(D+1)))**p
            pn_lin[i] = np.exp(-k*(i-mu)**2/2)
        pn_an = pn_an/np.sum(pn_an) # normalise distribution
        pn_lin = pn_lin/np.sum(pn_lin) # normalise distribution
        
        for i in range(D):
            if i > 0:
                G[i,i-1] = (pn_an[i]/pn_an[i-1])**(lam/2)
                L[i-1,i] = (pn_an[i-1]/pn_an[i])**((1-lam)/2)
        Gdag = np.conjugate(np.transpose(G)) # adjoint of G
        Ldag = np.conjugate(np.transpose(L)) # adjoint of L
        
        pn_num = pn_an

    else:    
        
        # ss cav photon distribution (analytical)
        for i in range(D):
            pn_an[i] = (np.sin(np.pi*(i+1)/(D+1)))**p
            pn_lin[i] = np.exp(-k*(i-mu)**2/2)
        pn_an = pn_an/np.sum(pn_an) # normalise distribution
        pn_lin = pn_lin/np.sum(pn_lin) # normalise distribution
        
        for i in range(D):
            if i > 0:
                G[i,i-1] = 1
                L[i-1,i] = (pn_an[i-1]/pn_an[i])**((1+q/2)/2)
        Gdag = np.conjugate(np.transpose(G)) # ctranspose of G
        Ldag = np.conjugate(np.transpose(L)) # ctranspose of L
        
        # ss cav photon distribution (numerical)
        for i in range(D):
            if i == 0:
                pn_num[i] = 1.
            elif i == 1:
                pn_num[i] = (1 - q/2)*pn_num[i-1]/L[i-1,i]**2
            elif i == 2:
                pn_num[i] = ((q-1)*pn_num[i-2] + (1 + L[i-2,i-1]**2 - q/2)*pn_num[i-1])/L[i-1,i]**2
            else:
                pn_num[i] = (-(q/2)*pn_num[i-3] + (q-1)*pn_num[i-2] + (1 + L[i-2,i-1]**2 - q/2)*pn_num[i-1])/L[i-1,i]**2
        pn_num = pn_num/np.sum(pn_num) # normalise distribution
    
    return pn_an, pn_lin, pn_num, G, Gdag, L, Ldag

p_arr = np.array([4,6,8,10,20,30])
p_lst = [r'$p=4$',r'$p=6$',r'$p=8$',r'$p=10$',r'$p=20$',r'$p=30$']
# D_arr = np.arange(21,5001,50)
D_arr = np.arange(21,2001,50)
q = -1
lam = 0

# neg_sums = np.zeros((len(D_arr),len(p_arr)), dtype='float128')
neg_sums = np.zeros((len(D_arr),len(p_arr)))

for i in range(len(D_arr)):
    for j in range(len(p_arr)):
        
        _, _, pnum, _, _, _, _ = distributions(D_arr[i],p_arr[j],q,lam)
        neg_sums[i,j] = np.sum(pnum[pnum < 0])
        if abs(neg_sums[i,j]) < 1e-16:
            neg_sums[i,j] = 1
        
        print(i,j)

# D = 50
# p = 4
# q = -1
# lam = 0

# n = np.arange(D)

# _, _, pnum, _, _, _, _ = distributions(D,p,q,lam)

# negs = np.sum(pnum[pnum < 0])
# print(negs)
# print(np.sum(pnum))

fig = plt.figure(figsize = (1*(3+3/8),0.9*(3+3/8)))

m = 0.1 # amount to take off axes width, in units of figure width
axs1 = fig.add_axes([0 + m, 0 + m, 1 - m, 1 - m]) # [left, bottom, width, height]

# a = [[] for i in range(len(p_arr))]
a = [plt.Line2D([], [], alpha=0) for i in range(len(p_arr))]
a_labels = [[] for i in range(len(p_arr))]
# b = [[] for i in range(4)]
b = [plt.Line2D([], [], alpha=0) for i in range(len(p_arr))]
b_labels = ['' for i in range(len(b))]
clrs = ['C0', 'C1',  'C2',  'C3',  'C4', 'k'] 

for i in range(len(p_arr)):
    
    a[i], = axs1.loglog((D_arr-1)/2, np.abs(neg_sums[:,i]), '.', color = clrs[i], markersize = 4.0)
    a_labels[i] = p_lst[i]
    
    if i <= 3:
    
        idx = np.where((np.abs(neg_sums[:,i]) <= 10**-10) & (np.abs(neg_sums[:,i]) >= 10**-16))[0]
        # axs1.loglog((np.array([D_arr[idx[0]],D_arr[idx[-1]]])-1)/2, np.array([np.abs(neg_sums[:,i])[idx[0]],np.abs(neg_sums[:,i])[idx[-1]]]), '-', label = p_lst[i])
        f = np.polyfit(np.log((D_arr[idx[0]:idx[-1]]-1)/2),  np.log(np.abs(neg_sums[:,i])[idx[0]:idx[-1]]), 1)
        mu_short = (D_arr[idx[0]:idx[-1]]-1)/2
        b[i], = axs1.loglog(mu_short, np.exp(f[1])*mu_short**f[0],'-', label = 'test', color = clrs[i+1])
        b_labels[i] = f'$O(\mu^{{{f[0]:.1f}}})$'
        print(np.exp(f[1]),f[0])

axs1.set_ylim([1e-16,1e-2])
axs1.set_xlim([1e1,1e3])
axs1.legend(np.append(a,b), 
            np.append(a_labels,b_labels), 
            ncol=2, 
            loc='lower left', 
            bbox_to_anchor=(0, 0),
            handletextpad=0.3,     # tighter symbol-to-text spacing
            columnspacing=0.5,     # tighter column spacing
            labelspacing=0.2,      # tighter vertical spacing
            borderpad=0.3,         # smaller border around content
            borderaxespad=0.3,     # closer to the axes
            handlelength=1.5,
            framealpha=1)
axs1.grid('on')
axs1.set_xlabel(r'$\mu$')
axs1.set_ylabel(r'Summed Elements')


# fig.savefig('Apx1.pdf', bbox_inches = 'tight')




