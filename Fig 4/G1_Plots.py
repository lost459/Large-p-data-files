#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 20 12:24:12 2024

@author: s5257291
"""

import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({"mathtext.fontset": 'cm',
                     'font.family':'times new roman'})

plt.rc('font', size=9) #controls default text size
plt.rc('axes', titlesize=9) #fontsize of the title
plt.rc('axes', labelsize=11) #fontsize of the x and y labels
plt.rc('xtick', labelsize=9) #fontsize of the x tick labels
plt.rc('ytick', labelsize=9) #fontsize of the y tick labels
plt.rc('legend', fontsize=9) #fontsize of the legend
plt.rc('lines', linewidth=1) #linewidth of plot lines
plt.rc('lines', markersize=3.75) #size of datapoint
plt.rc('axes', labelpad=0.75)

###############################################################################

# load data, rows sorted according to coherence, D, p, lam, etc...
data = np.loadtxt(r'/Users/s5257291/Documents/Large_p_Paper/2024/Data Files/Fig 4/G1_final.txt', delimiter='\t', skiprows=1)
last_rows_0 = data[4*121:]
last_rows_0 = last_rows_0[last_rows_0[:, 3].argsort()]
last_rows_1 = last_rows_0[:121]
last_rows_1 = last_rows_1[last_rows_1[:, 5].argsort()]
last_rows_2 = last_rows_0[121:]
last_rows_2 = last_rows_2[last_rows_2[:, 5].argsort()]
data = np.vstack((data[:4*121], last_rows_1,last_rows_2))

taus_short = data[0:1*121:6,5]
taus_long = data[1*121:2*121:6,5]


G1_short = data[0:1*121:6,0]
G1_long = data[1*121:2*121:6,0]

G1lam_short = data[5*121:6*121:6,0]
G1lam_long = data[2*121:3*121:6,0]

G1q_short = data[4*121:5*121:6,0]
G1q_long = data[3*121:4*121:6,0]


###############################################################################

ts_long = np.linspace(0,taus_long[-1],1000)
ts_short = np.linspace(0,taus_short[-1],1000)

def G1_lam_an(mu,p,lam,ts):
    k = np.pi**2 * p / (4 * mu**2)
    g2 = 1 + lam * (lam - 1) * k * np.exp(- k * ts)
    ell = (2*lam**2 - 2*lam + 1) * np.pi**4 * p**2 / (64 * mu**4)
    G1 = np.exp(- ell * ts / 2) * (1 + (1/4) * (g2 - g2[0]))
    G1_pure = np.exp(- ell * ts / 2) * (1 + (1/4) * (1 - g2[0]))
    return G1, G1_pure

def G1_q_an(mu,p,q,ts):
    k = np.pi**2 * p / (4 * mu**2)
    g2 = 1 + q * (1/2 + q/4) * k * np.exp(- k * (1 + q/2) * ts)
    ell = (1 + q/2)**2 * np.pi**4 * p**2 / (64 * mu**4)
    G1 = np.exp(- ell * ts / 2) * (1 + (1/4) * (g2 - g2[0]))
    G1_pure = np.exp(- ell * ts / 2) * (1 + (1/4) * (1 - g2[0]))
    return G1, G1_pure

G1lam_longan, _ = G1_lam_an(250,50,0.5,ts_long)
G1q_longan, _ = G1_q_an(250,50,-1,ts_long)
G1_longan, _ = G1_lam_an(250,50,0,ts_long)
G1lam_shortan, G1lam_pure = G1_lam_an(250,50,0.5,ts_short)
G1q_shortan, G1q_pure = G1_q_an(250,50,-1,ts_short)
G1_shortan, _ = G1_lam_an(250,50,0,ts_short)

k = 4*250**2 / (np.pi**2 * 50)

###############################################################################

fig = plt.figure(figsize = (1*(3+3/8),1.5*(3+3/8)))

m = 0.1 # amount to take off axes width, in units of figure width
axs1 = fig.add_axes([0 + m, 1/2 + m, 1 - m, 1/2 - m]) # [left, bottom, width, height]
axs2 = fig.add_axes([0 + m, 0 + m, 1 - m, 1/2 - m]) # [left, bottom, width, height]

axs1.plot(ts_short/k,G1lam_shortan-1, 'k')
axs1.plot(ts_short/k,G1q_shortan-1, 'k')
axs1.plot(ts_short/k,G1_shortan-1, 'k')
axs1.plot(taus_short/k,G1q_short-1, 'C3x', label = r'$p, q$-fam.; $q=-1$')
axs1.plot(taus_short/k,G1lam_short-1, 'C2x', label = r'$p,\lambda$-fam.; $\lambda=0.5$')
axs1.plot(taus_short/k,G1_short-1, 'C0x', label = r'Poissonian Limit')
axs1.legend()
axs1.set_xlabel(r'$\mathcal{N}t/k$')
axs1.set_ylabel(r'$g^{(1)}(s+t,s)-1$')
axs1.grid('on')
axs1.set_xlim([0,taus_short[-1]/k])
axs1.set_ylim([-10.1e-4,0.1e-4])
axs1.ticklabel_format(axis = 'y', style='sci', scilimits=(0,0))
axs1.yaxis.get_offset_text().set_visible(False)
axs1.text(0.0, 1.025, r'$\times10^{-4}$', transform=axs1.transAxes)

axs2.plot(ts_long/k**2,G1lam_longan, 'k')
axs2.plot(ts_long/k**2,G1q_longan, 'k')
axs2.plot(ts_long/k**2,G1_longan, 'k')
axs2.plot(taus_long/k**2,G1lam_long, 'C2x')
axs2.plot(taus_long/k**2,G1q_long, 'C3x')
axs2.plot(taus_long/k**2,G1_long, 'C0x')
axs2.set_xlabel(r'$\mathcal{N}t/k^2$')
axs2.set_ylabel(r'$g^{(1)}(s+t,s)$')
axs2.grid('on')
axs2.set_xlim([0,taus_long[-1]/k**2])
axs2.set_ylim([0,1.025])

# fig.savefig('G1_plot.pdf', bbox_inches = 'tight')



