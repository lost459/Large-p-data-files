#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 13 11:32:52 2024

@author: s5257291
"""

import numpy as np
import matplotlib.pyplot as plt
import math

plt.rcParams.update({"mathtext.fontset": 'cm',
                     'font.family':'times new roman'})

plt.rc('font', size=9) #controls default text size
plt.rc('axes', titlesize=9) #fontsize of the title
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rc('xtick', labelsize=9) #fontsize of the x tick labels
plt.rc('ytick', labelsize=9) #fontsize of the y tick labels
plt.rc('legend', fontsize=9) #fontsize of the legend

# Return Analytic and Numerical ss distributions for R.P. laser
def distributions(D,p,q):

    mu = (D-1)/2    
    k = np.pi**2*p/(4*mu**2)

    pn_num = np.zeros(D) # cav photon distribution (numerical)
    pn_an = np.zeros(D) # cav photon distribution (analytical)
    pn_lin = np.zeros(D) # cav photon distribution (linearized)
    G = np.zeros((D,D)) # gain operator
    L = np.zeros((D,D)) # loss operator
    
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
    Gdag = np.conjugate(np.transpose(G)) # adjoint of G
    Ldag = np.conjugate(np.transpose(L)) # adjoint of L
    
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
    for i in range(D):
        if pn_num[i] < 0:
            pn_num[i] = 0 # Get rid of negative values
    pn_num = pn_num/np.sum(pn_num) # normalise distribution

    
    return pn_an, pn_lin, pn_num, G, Gdag, L, Ldag

# Function for exact ss.
def rho(p,D):
    
    arr = np.zeros(D)
    
    for i in range(D):
        arr[i] = (np.sin(np.pi*(i+1)/(D+1)))**p
        
    arr = arr/np.sum(arr)
        
    return arr

# Function for approximate (Gaussian) ss.
def rho_apx(p,D):
    
    arr = np.zeros(D)
    mu = (D - 1)/2
    sigma = np.sqrt(4*mu**2/(p*np.pi**2))
    
    for i in range(D):
        arr[i] = np.exp(-(i - mu)**2/(2*sigma**2))/(sigma*np.sqrt(2*np.pi))
        
    arr = arr/np.sum(arr)
        
    return arr

# Function for nonzero coefficients of loss operator.
def Loss(p,lam,D):
    
    arr = np.zeros(D-1)
    
    for i in range(1,D):
        # arr[i-1] = (np.sin(np.pi*(i)/(D + 1))/np.sin(np.pi*(i+1)/(D + 1)))**(p*(1 - lam)/2)
        arr[i-1] = (np.sin(np.pi*(i)/(D + 1))/np.sin(np.pi*(i+1)/(D + 1)))**(p*(1 - lam))
        
    return arr

# Function for nonzero approximate (linearized) coefficients of loss operator.
def Loss_apx(p,lam,D):
    
    arr = np.zeros(D-1)
    mu = (D - 1)/2
    
    for i in range(1,D):
        # arr[i-1] = 1 - np.pi**2*p*(1 - lam)*(mu - i - 1/2)/(8*mu**2)
        arr[i-1] = 1 - np.pi**2*p*(1 - lam)*(mu - i)/(4*mu**2)
        
    return arr

# Function for nonzero coefficients of gain operator.
def Gain(p,lam,D):
    
    arr = np.zeros(D-1)
    
    for i in range(1,D):
        # arr[i-1] = (np.sin(np.pi*(i + 1)/(D + 1))/np.sin(np.pi*(i)/(D + 1)))**(p*lam/2)
        arr[i-1] = (np.sin(np.pi*(i + 1)/(D + 1))/np.sin(np.pi*(i)/(D + 1)))**(p*lam)
        
    return arr

# Function for nonzero approximate (linearized) coefficients of Gain operator.
def Gain_apx(p,lam,D):
    
    arr = np.zeros(D-1)
    mu = (D - 1)/2
    
    for i in range(1,D):
        # arr[i-1] = 1 + np.pi**2*p*lam*(mu - i + 1/2)/(8*mu**2)
        arr[i-1] = 1 + np.pi**2*p*lam*(mu - i)/(4*mu**2)
        
    return arr

##################
#### Small p #####
##################

D = 101
p = 4
lam = 0.5

lw = 1.0
lw2 = 1.5

n = np.arange(0,D)
n_minus_one = np.arange(1,D)
rho_arr = rho(p,D)
rho_apx_arr = rho_apx(p,D)
Loss_arr = Loss(p,lam,D)
Loss_apx_arr = Loss_apx(p,lam,D)
Gain_arr = Gain(p,lam,D)
Gain_apx_arr = Gain_apx(p,lam,D)
_, _, rho_reg, _, _, _, _ = distributions(D,p,-1)

fig = plt.figure(figsize = (3+3/8,3+3/8))

axs1 = fig.add_subplot(2,2,1)
axs1.text(0.05, 0.9, r'(a)', transform=axs1.transAxes)
axs1.text(0.4, 0.7, r'$\mu = 50$' '\n' r'$p = 4$', transform=axs1.transAxes)
axs1.plot(n_minus_one,Loss_arr,'C3-',linewidth = lw2)
axs1.plot(n_minus_one,Loss_apx_arr,'k--',linewidth = lw)
axs1.plot(n_minus_one,Gain_arr,'C2-',linewidth = lw2)
axs1.plot(n_minus_one,Gain_apx_arr,'k--',linewidth = lw)
# axs1.set_xlabel(r'$n$')
axs1.set_ylabel(r'Matrix Coefficient')
axs1.grid('on')
axs1.set_xlim([0,100])
axs1.set_ylim([min([min(Loss_arr),min(Gain_arr)]),max([max(Loss_arr),max(Gain_arr)])])
axs1.tick_params(axis='y', which='major', pad=1)

axs2 = fig.add_subplot(2,2,2)
axs2.text(0.05, 0.9, r'(b)', transform=axs2.transAxes)
axs2.text(0.35, 0.15, r'$\mu = 50$' '\n' r'$p = 4$', transform=axs2.transAxes)
axs2.plot(n,rho_arr,'C1-',linewidth = 2.5)
axs2.plot(n,rho_reg,'C0-',linewidth = lw)
axs2.plot(n,rho_apx_arr,'k--',linewidth = lw)
axs2.ticklabel_format(axis = 'y', style='sci', scilimits=(0,0))
axs2.yaxis.get_offset_text().set_visible(False)
axs2.text(0.0, 1.015, r'$10^{-2}$', transform=axs2.transAxes)
# axs2.set_xlabel(r'$n$')
axs2.grid('on')
axs2.set_xlim([min(n),max(n)+1])
axs2.set_ylim([0,1.05*max([max(rho_arr),max(rho_apx_arr)])])
axs2.tick_params(axis='y', which='major', pad=1)

##################
#### Large p #####
##################

D = 501
p = 50
# lam = 0.5

lw = 1.0

n = np.arange(0,D)
n_minus_one = np.arange(1,D)
rho_arr = rho(p,D)
rho_apx_arr = rho_apx(p,D)
Loss_arr = Loss(p,lam,D)
Loss_apx_arr = Loss_apx(p,lam,D)
Gain_arr = Gain(p,lam,D)
Gain_apx_arr = Gain_apx(p,lam,D)
_, _, rho_reg, _, _, _, _ = distributions(D,p,-1)

axs3 = fig.add_subplot(2,2,3)
axs3.text(0.4, 0.7, r'$\mu = 250$' '\n' r'$p = 50$', transform=axs3.transAxes)
axs3.text(0.1, 0.9, r'(c)', transform=axs3.transAxes)
axs3.plot(n_minus_one,Loss_arr,'C3-',linewidth = lw2)
axs3.plot(n_minus_one,Loss_apx_arr,'k--',linewidth = lw)
axs3.plot(n_minus_one,Gain_arr,'C2-',linewidth = lw2)
axs3.plot(n_minus_one,Gain_apx_arr,'k--',linewidth = lw)
axs3.set_xlabel(r'$n$')
axs3.set_ylabel(r'Matrix Coefficient')
axs3.grid('on')
axs3.set_xlim([175,325])
axs3.set_ylim([1-.1,1+.1])
axs3.tick_params(axis='y', which='major', pad=1)

axs4 = fig.add_subplot(2,2,4)
axs4.text(0.05, 0.9, r'(d)', transform=axs4.transAxes)
axs4.text(0.325, 0.15, r'$\mu = 250$' '\n' r'  $p = 50$', transform=axs4.transAxes)
axs4.plot(n,rho_arr,'C1-',linewidth = 2.5)
axs4.plot(n,rho_reg,'C0-',linewidth = lw)
axs4.plot(n,rho_apx_arr,'k--',linewidth = lw)
axs4.ticklabel_format(axis = 'y', style='sci', scilimits=(0,0))
axs4.yaxis.get_offset_text().set_visible(False)
axs4.text(0.0, 1.015, r'$10^{-2}$', transform=axs4.transAxes)
axs4.set_xlabel(r'$n$')
axs4.grid('on')
axs4.set_xlim([175,325])
axs4.set_ylim([0,1.05*max([max(rho_arr),max(rho_apx_arr)])])
axs4.tick_params(axis='y', which='major', pad=1)

fig.tight_layout(pad=0.25)

plt.subplots_adjust(wspace = 0.3)

# fig.savefig('Linearization_Visualisation.pdf', bbox_inches = 'tight')