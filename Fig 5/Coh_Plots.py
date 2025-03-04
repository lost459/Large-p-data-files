#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 22 15:09:48 2024

@author: s5257291
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import gamma as g
import matplotlib.ticker as mticker

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
data = np.loadtxt(r'/Users/s5257291/Documents/Large_p_Paper/2024/Data Files/Fig 5/coh.txt', delimiter='\t', skiprows=1)

ps = data[0:47,2]

Cohlam = data[0:47,0]
Cohq = data[47:94,0]
Cohp = data[94:141,0]

###############################################################################

pdense = np.linspace(4,50,1000)

def coh_laman(mu,p,lam):
    ell = (2*lam**2 - 2*lam + 1) * np.pi**4 * p**2 / (64 * mu**4)
    ell_h = ell * g((p+2)/2)*g((p-3)/2)/( g((p+1)/2)*g((p-2)/2) )
    return 4/ell, 4/ell_h

def coh_qan(mu,p,q):
    ell = (1 + q/2)**2 * np.pi**4 * p**2 / (64 * mu**4)
    ell_h = ell * g((p+2)/2)*g((p-3)/2)/( g((p+1)/2)*g((p-2)/2) )
    return 4/ell, 4/ell_h

mu = 250

cohlaman, cohlaman_h = coh_laman(mu,pdense,0.5)
cohqan, cohqan_h = coh_qan(mu,pdense,-1)
cohpan, cohpan_h = coh_laman(mu,pdense,0)

cohlaman_sp, _ = coh_laman(mu,ps,0.5)
cohqan_sp, _ = coh_qan(mu,ps,-1)
cohpan_sp, _ = coh_laman(mu,ps,0)

###############################################################################

fig = plt.figure(figsize = (1*(3+3/8),2*(3+3/8)))

m = 0.1 # amount to take off axes width, in units of figure width
axs1 = fig.add_axes([0 + m, 1/2 + m, 1 - m, 1/2 - m]) # [left, bottom, width, height]

axs1.plot(pdense[100:],cohlaman[100:], 'k')
axs1.plot(pdense[100:],cohqan[100:], 'k')
axs1.plot(pdense[100:],cohpan[100:], 'k')
axs1.plot(pdense,cohlaman_h, '-', color = [0.85,0.85,0.85])
axs1.plot(pdense,cohqan_h, '-', color = [0.85,0.85,0.85])
axs1.plot(pdense,cohpan_h, '-', color = [0.85,0.85,0.85])
axs1.plot(ps,Cohlam, 'C2.')
axs1.plot(ps,Cohq, 'C3.')
axs1.plot(ps,Cohp, 'C0.')
axs1.grid('on')
axs1.set_xlim([4,50])
axs1.set_xticks([4,10,20,30,40,50])
axs1.set_ylim([0,1.0e9])
axs1.ticklabel_format(axis = 'y', style='sci', scilimits=(0,0))
axs1.yaxis.get_offset_text().set_visible(False)
axs1.text(0.0, 1.025, r'$\times10^{9}$', transform=axs1.transAxes)
axs1.set_xlabel(r'$p$')
axs1.set_ylabel(r'$\mathfrak{C}$')


axsin = fig.add_axes([1/4 + m, 4.9/8 + m, 2.9/4 - m, 2.8/8 - m]) # [left, bottom, width, height]
axsin.plot(pdense[100:],cohlaman[100:], 'k', label = r'Analytical')
axsin.plot(pdense[100:],cohqan[100:], 'k')
axsin.plot(pdense[100:],cohpan[100:], 'k')
axsin.plot(pdense,cohlaman_h, '-', color = [0.85,0.85,0.85])
axsin.plot(pdense,cohqan_h, '-', color = [0.85,0.85,0.85])
axsin.plot(pdense,cohpan_h, '-', color = [0.85,0.85,0.85])
axsin.plot(ps,Cohq, 'C3.', label = r'$p,q$-fam.; $q = -1$')
axsin.plot(ps,Cohlam, 'C2.', label = r'$p,\lambda$-fam.; $\lambda = 0.5$')
axsin.plot(ps,Cohp, 'C0.', label = r'Poissonian Limit')
axsin.grid('on')
axsin.set_xlim([30,50])
axsin.set_ylim([0,5e7])
axsin.ticklabel_format(axis = 'y', style='sci', scilimits=(0,0))
axsin.yaxis.get_offset_text().set_visible(False)
axsin.text(0.0, 1.015, r'$\times10^{7}$', transform=axsin.transAxes)
axsin.legend(loc='lower left', bbox_to_anchor=(0.35,0.56), framealpha = 1.0)


###############################################################################
# Fit to error in p > 20
fit2 = np.polyfit(np.log(ps[16:]), np.log(np.abs(Cohq - cohqan_sp)/Cohq)[16:], 1)

###############################################################################

axs2 = fig.add_axes([0 + m, 1.08/6 + m, 1 - m, 1/3 - m]) # [left, bottom, width, height]
axs2.loglog(ps,np.exp(fit2[1])*ps**fit2[0], 'k-', label = r'$7.74p^{-1.45}$')
axs2.loglog(ps,np.abs(Cohq - cohqan_sp)/Cohq, 'C3.')
axs2.xaxis.set_minor_formatter(mticker.ScalarFormatter())
axs2.xaxis.set_major_formatter(mticker.ScalarFormatter())
axs2.set_ylim([.01,1])
axs2.set_xlim([5,50])
axs2.set_xticks([5,10,20,30,40,50])
axs2.set_xlabel(r'$p$')
axs2.set_ylabel(r'$\epsilon$')
axs2.legend()


# fig.savefig('coh_plot_1.pdf', bbox_inches = 'tight')


